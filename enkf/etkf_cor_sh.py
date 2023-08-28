from pylab import *
from numpy import *
import numpy.linalg as nlg
import sys
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse.linalg import spsolve, factorized
import numpy as np

class etkf_c:
    """ class for Ensemble Transform Kalman Filter 
    """
    def __init__(self, x_mean, y_mean, x, y, r, xref=1, xnorm=1.0, do_debug=False):
        """ initialize the class 
            Arguments:  
                x   :    array of dimension [nx, ne]  with nx = the number of the  variables,  ne=the number of ensemble 
                y   :   array  of dimension [nobs, ne]   ensemble of model forecasts 
                r   :   observation error covariances given in 1 dimension
            Notes:  
                    self.x=(x-mean(x))/sqrt(2) 
                   self.y=(y-mean(y))/sqrt(2) 
        """
        self.nx, self.ne=shape(x)
        self.x=array(x)
        print('etkf: nx, ne', self.nx, self.ne)
        self.x_mean=array(x_mean)
        print('etkf: x_mean', self.x_mean[0:6])
        diag_r=diag(r)
        print('diag--r', shape(r), diag_r[0:5])
        
        max_y, min_y=max(y.flat), min(y.flat)
        max_y=max([abs(max_y), abs(min_y)])
        ulimt=1.0e-8*max_y
        new_y=array(y) # where(abs(y)<ulimt, 0.0, y)
        
        self.xnorm=xnorm
        # self.x=self.x[:,:]
        self.x=self.x/sqrt(self.xnorm)
        
        self.y=new_y/sqrt(self.xnorm)
        self.sp_y=csr_matrix(self.y)
        
        self.ny=size(y[:,0])  # ---> nobs 
        self.y_mean=array(y_mean)
        self.yy=None
        # observation error 
        self.r=r
        self.sp_r=csr_matrix(r)

        np.save('test_r', self.r)

        # s=r+dy*dy.t
#        self.sp_b=self.sp_y*self.sp_y.T
        self.sp_b = dot(self.sp_y, transpose(self.sp_y))
        self.sp_s=self.sp_b + self.sp_r
        
    
        
    def reset(self, x_mean, y_mean, x, y, r, xref=1.0, xnorm=1.0, do_debug=False):
        """ initialize the class 
        Arguments:  
        x   :    array of dimension [nx, ne]  with nx = the number of the  variables,  ne=the number of ensemble 
                y   :   array  of dimension [nobs, ne]   ensemble of model forecasts 
                r   :   observation error covariances given in 1 dimension
            Notes:  
                    self.x=(x-mean(x))/sqrt(2) 
                   self.y=(y-mean(y))/sqrt(2) 
        """
        
        


        self.nx, self.ne=shape(x)
        self.x=array(x)
        print('etkf: nx, ne', self.nx, self.ne)
        self.x_mean=array(x_mean)
        print('etkf: x_mean', self.x_mean[0:6])

        max_y, min_y=max(y.flat), min(y.flat)
        max_y=max([abs(max_y), abs(min_y)])
        ulimt=1.0e-8*max_y
        new_y=where(abs(y)<ulimt, 0, y)
        
        self.xnorm=xnorm
        self.x=self.x[:,:]
        self.x=self.x/sqrt(self.xnorm)
        
        self.y=new_y/sqrt(self.xnorm)
        self.sp_y=csr_matrix(self.y)
        
        self.ny=size(y[:,0])  # ---> nobs 
        self.y_mean=array(y_mean)
        
        # observation error 
        self.r=r
        self.sp_r=csr_matrix(r)
        
        # s=r+dy*dy.t
        self.sp_b=self.sp_y*self.sp_y.T
        self.sp_s=self.sp_b+self.sp_r
        
        
        
        return 0
    
    
    def get_postprior(self,yobs):
        """ calculate postprior from assimilation of observations 
            Arguments:
                yobs    :   array of observations
            Returns:
                postprior
            Notes:
                xa=xf+k*(yobs-f(x))=xf+xf'*inc_m
        
        """
        inc_m=self.get_inc_m(yobs)
        x=self.x_mean+dot(self.x, transpose(inc_m))
        return x
    
    def get_tm(self):
        """ the analysis increment on the variables  
                Arguments:
                    None
                Returns:
                    tm   :   matrix of [ne, ne] for transform from X'(f)--> X'(a)      
                Notes: 
                    tm=U*(1+W*WT)^-0.5)*UT
        """
        tm=list()
#        asf=linsolve.factorized(self.sp_r)
        self.sp_r = csc_matrix(self.sp_r)
        asf=factorized(self.sp_r)
        print(shape(self.r), shape(self.y), self.ne)
        
        
        for ie in range(self.ne):
            #  t1=linsolve.spsolve(self.sp_r, self.y[:,ie])
            # print 'ie', ie
            
            # print shape(self.y[:,ie])
            tx=self.y[:,ie]
            tx=squeeze(tx)
            # np.save('test_z', tx)
            t1=asf(tx) # self.y[:,ie])
            tm.append(t1)
        tm=array(tm)
        
        print('shape---tm', shape(tm))
        
        tm=transpose(tm)
        sp_tm=csr_matrix(tm)
        s1=self.sp_y.T*sp_tm
        s1=0.5*(s1+s1.T)
        
        s1=s1.todense()
        s1=array(s1)
        
        s1=identity(self.ne)+s1
        svd_err=False
        
        try:
            u, w, vt=svd(s1)
        except LinAlgError:
            tm=identity(self.ne)
            svd_err=True
            print ('err_return--:', shape(tm), shape(svd))
            
            return tm, svd_err
        
        # tm=u/sqrt(1.0+w)
        tm=u/sqrt(w)
        
        tm=dot(tm, transpose(u))
        
        print ('w in tm', w[0:5])
        print ('tm in tm', tm[0:5,0])
        
        return tm,svd_err
    
    
        
    def get_inc_m(self,yobs):
        """ get the analysis increment 
                Arguments:
                    yobs    :           array of [nobs], the observations
                Returns:     
                    inc_m   :           array of [ne]
                Notes:
                    inc_m=U*W*(I+WT*W)^-1*VT*R^-0.5*(Yobs-mean_Y)
        """
    
        # v_nx,v_ny=shape(self.vy)
        dy=yobs-self.y_mean # innovation
        # solve the linear equation
        print('dy: ', shape(dy))
        ## modified by rainbow at 2020.10.09
#        inc_m=linsolve.spsolve(self.sp_s, dy)
        inc_m = spsolve(self.sp_s, dy)
        ## modification end
        
        inc_m=array(inc_m)
        
        inc_m=dot(self.y.T, inc_m)
        
        return inc_m
        # self.ana_array=transpose(ana_array)
        
    def get_x_inc(self, inc_m, x_in=None, mean_x_in=None):
        """ the analysis increment on the variables  
                Arguments:
                    yobs    :           array of [nobs], the observations
                    inc_m   :          array of [ne], the increment matrix  
                Returns:
                    xinc  :   array of [nx],  the X increments      
        """    
        if (x_in==None):
            xinc=dot(self.x, transpose(inc_m))
        else:
            if (mean_x_in==None):
                dx=array(x_in)
            else:
                dx=x_in-mean_x_in[:,newaxis]
            
            dx=dx/sqrt(self.xnorm)
            xinc=dot(dx, transpose(inc_m))
            
        return xinc
    def get_y_inc(self, y_in, inc_m, mean_y=None):
        """ the analysis increment on the variables  
                Arguments:
                    y_in    :           array of [nnew_y, ne], the 'new' observations
                    inc_m :           array of [ne],        the increment from analysis
                Returns:
                    yinc   :       array of [ny],  y increments from xinc     
        """    
        # mean_y=mean(y_in, axis=1)
        if (mean_y==None):
            dy=array(y_in[:,:])
            # -mean_y[:, newaxis]
        else:
            dy=y_in[:,:] -mean_y[:, newaxis]
        
        dy=dy/sqrt(self.xnorm)
        yinc=dot(dy, transpose(inc_m))
        return yinc
    
    def get_y_transform(self, y_in, tm, mean_y=None):
        """ the analysis increment on the variables  
        Arguments:
                        y_in    :           array of [nnew_y, ne], the 'new' observations
                        tm      :           array of [ne, ne],   the transformation matrix
                Returns:
                yt        :            array of [nnew_y, ne] by transform tm   as X'(a)   
        """    

        if (mean_y==None):
            dy=array(y_in[:,:])
            # -mean_y[:, newaxis]
        else:
            dy=y_in[:,:] -mean_y[:, newaxis]
        
        
        # mean_y=mean(y_in, axis=1)
        # dy=y_in[:,:]-mean_y[:, newaxis]
        yt=dot(dy, tm)
        return yt
        
        

 
        
if (__name__=='__main__'):
    import numpy.random as rnd
    x=array([0.3, -0.1])
    xtrue=array([0.2, 0.1])
    dx=array([1.0, 4.0])
    xen=list()
    xen.append(x)

    new_x=array(x)
    new_x[0]=x[0]+dx[0]
    xen.append(new_x)

    new_x=array(x)
    new_x[0]=x[0]-dx[0]
    xen.append(new_x)
    
    new_x=array(x)
    new_x[1]=x[1]+dx[1]
    xen.append(new_x)

    new_x=array(x)
    new_x[1]=x[1]-dx[1]
    xen.append(new_x)
    xen=array(xen)
    print(shape(xen))
    print(xen)
    y=list()
    yobs=list()
    tstep=1.0
    hx=array([1.0, 2])
    robs=list()
    R0=array([0.2, 0.4])
    nsteps=300
    ym=list()
    
    for it in range(nsteps):
        yt=xtrue*(hx*tstep)
        err=array([rnd.normal(0.0, R0[0]),rnd.normal(0.0, R0[1])])
        yt=yt + err
        yobs.append(array(yt))
        robs.append(R0*R0)
        yy=xen*(hx*tstep)
        ym.append(array(yy))
    ym=array(ym)
    #    print ym[:,1,:]
    nx, ny, nz=shape(ym)
    new_ym=zeros([nx*nz, ny], float)
    for ien in range(ny):
        yy=ym[:,ien,:]
        yy=yy.reshape(-1)
        
        new_ym[:,ien]=yy[:]

    ym=new_ym
    
    xen=transpose(xen)
    # print xen
    
    # print 'ym 0', ym[:,0]
    # print 'ym 1', ym[:,1]
    
    
    
    yobs=array(yobs)
    yobs=yobs.reshape(-1)
    
    
    robs=array(robs)
    robs=robs.reshape(-1)
    
    # print yobs
    # print robs
    
    # print robs
    
    
    # print yobs
    # print robs
    
    ett=etkf_c(xen, ym, robs)
    tm, svd_err=ett.get_tm()
    #    print array2string(tma, precision=2, suppress_small=True)
    print('='*80)
    k=ett.get_k()
    print('Kalman Gain Matrix (size)',  shape(k))
    
    # print array2string(k, precision=2, suppress_small=True)
    xinc=dot(k, transpose(yobs-ett.y_mean)) 
    xan=ett.x_mean+xinc
    inc_m=ett.get_inc_m(yobs)
    
    xinc=ett.get_x_inc(inc_m)
    xan2=ett.x_mean+xinc
    
    print('Forecast:',  x)
    print('Analysis 1:',  xan)
    print('Analysis2 :',  xan2)
    print('True Value',  xtrue)
    err_cov=diag(dx*dx)
    print('='*30+'Forcast Error Covariance'+'='*30)
    print(array2string(err_cov, precision=6, suppress_small=True))
    print(' ')
    print('-'*30+'Analysis Error Covariance'+'-'*30)
    ky=dot(k, ett.y)
    kyx=dot(ky, transpose(ett.x))
    err_cov=err_cov-kyx
    print(array2string(err_cov, precision=6, suppress_small=True))
    print(' ')
    print('-'*30+'Estimated Covariance from tma'+'-'*30)
    xen_an=dot(ett.x, tm)
    err_cov=dot(xen_an, transpose(xen_an))
    err_cov[1,1]=sqrt(err_cov[1,1])
    err_cov[0,0]=sqrt(err_cov[0,0])
    
    print(array2string(err_cov, precision=6, suppress_small=True))
    
    
    
    
        
        
        
    
        
    
        
        
