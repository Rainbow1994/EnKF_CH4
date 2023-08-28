from pylab import *
from numpy import *
import numpy.linalg as nlg
import sys
# import matrix_tool as mtool

class etkf_c:
    """ class for Ensemble Transform Kalman Filter 
    """
    def __init__(self, mean_x, mean_y, x, y, r, xref=1, xnorm=1.0, do_debug=False):
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
        self.mean_x=array(mean_x)
        print('etkf: mean_x -6, -1', self.mean_x[-6:-1])
        
        self.xnorm=xnorm
        self.x=self.x[:,:]
        self.x=self.x/sqrt(self.xnorm)
        

        dev=dot(self.x, transpose(self.x))
        
        self.y=array(y)
        self.ny=size(y[:,0])  # ---> nobs 
        self.mean_y=array(mean_y)
        self.y=self.y/sqrt(self.xnorm)
        self.yy=None
        self.r=r
        self.sqr=sqrt(r)
        # print shape(self.sqr)
        # print shape(self.y)
        # scy=self.y/self.sqr[:,newaxis]  # scale y by sqrt(r)
            
        self.scyt=transpose(self.y)/self.sqr # the scaled y=r^-0.5*y ==>scaled yt=yt*r^-0.5
        self.uy, self.wy, self.vyt=nlg.svd(self.scyt) # decomposition of tranpose(scaled_y)
        if (do_debug):
            # subplot(2,1,1)
            # plot(self.wy)
            # plot(self.scyt[0,0:50])
            # plot(self.scyt[1,0:50])
            # subplot(2,1,2)
            scx=self.x[:,:]/xref
            dev=dot(scx, transpose(scx))
            pcolor(dev)
            colorbar()
            
        
            
        
            # show()
        
       #self.wy has only one dimension with size of min[ny, ne]
        print('etkf: shape(y)',  shape(self.y))
        print('etkf: y[0:8,-6]',  self.y[0:8,-6])
        print('etkf: r',  self.r[0:8])
        print('etkf: wy', self.wy[0:8])
        
        #self.md=1.0/(1.0+self.wy*self.wy) #  1d array of min[ny, ne] 
        self.md = divide(ones_like(self.wy), ones_like(self.wy)+self.wy*self.wy)
        # in calculation of  K , we need 1.0/(1.0+w*wt)  of [nobs, nobs]
        # in calculation of transform matrix, we need 1.0/(1.0+wt*w) of [ne, ne]
        
        idx=where(self.wy>0.0)
        idx=squeeze(idx)
        self.n_wy=size(idx)
        self.uw=self.uy[:, 0:self.n_wy]*self.wy[0:self.n_wy]  # should have the shape of [ne, nobs] but    now only have [ne, self.n_wy]
    def reset(self, mean_x, mean_y, x, y, r, xref=1.0, xnorm=1.0, do_debug=False):
        """ initialize the class 
        Arguments:  
        x   :    array of dimension [nx, ne]  with nx = the number of the  variables,  ne=the number of ensemble 
                y   :   array  of dimension [nobs, ne]   ensemble of model forecasts 
                r   :   observation error covariances given in 1 dimension
            Notes:  
                    self.x=(x-mean(x))/sqrt(2) 
                   self.y=(y-mean(y))/sqrt(2) 
        """
        
        dy=array(y)
        dy=dy/sqrt(xnorm)
        sqr=sqrt(r)
        scyt=transpose(dy)/sqr # the scaled y=r^-0.5*y ==>scaled yt=yt*r^-0.5
        
        try: 
            uy, wy, vyt=nlg.svd(scyt)
        except nlg.LinAlgError:
            return -1
        
        
        self.nx, self.ne=shape(x)
        self.x=array(x)
        self.xnorm=xnorm
        # print self.nx
        self.x=self.x/sqrt(self.xnorm)
        self.y=array(y)
        self.ny=size(y[:,0])  # ---> nobs 
        self.mean_y=array(mean_y)
        self.y=self.y/sqrt(self.xnorm)
        self.yy=None
        self.r=r
        self.sqr=sqrt(r)
        # print shape(self.sqr)
        # print shape(self.y)
        # scy=self.y/self.sqr[:,newaxis]  # scale y by sqrt(r)
        
        self.scyt=transpose(self.y)/self.sqr # the scaled y=r^-0.5*y ==>scaled yt=yt*r^-0.5
        self.uy, self.wy, self.vyt=nlg.svd(self.scyt) # decomposition of tranpose(scaled_y)
        # self.wy has only one dimension with size of min[ny, ne]
        self.md=1.0/(1.0+self.wy*self.wy) #  1d array of min[ny, ne] 
        # in calculation of  K , we need 1.0/(1.0+w*wt)  of [nobs, nobs]
        # in calculation of transform matrix, we need 1.0/(1.0+wt*w) of [ne, ne]
        if (do_debug):
            # plot(self.wy)
            # plot(self.scyt[0,0:50])
            # plot(self.scyt[1,0:50])
            scx=self.x[:,:]/xref
            dev=dot(scx, transpose(scx))
            pcolor(dev)
            colorbar()
            
        idx=where(self.wy>0.0)
        idx=squeeze(idx)
        self.n_wy=size(idx)
        self.uw=self.uy[:, 0:self.n_wy]*self.wy[0:self.n_wy]  # should have the shape of [ne, nobs] but    now only have [ne, self.n_wy
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
        x=self.mean_x+dot(self.x, transpose(inc_m))
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

        inv_s=ones(self.ne, float)
        print('etkf: m_wy', self.n_wy)
        inv_s[0:self.n_wy]=sqrt(self.md[0:self.n_wy])
        print('etkf:  inv_s', inv_s[0:6])
        tm=self.uy*inv_s
        tm=dot(tm, transpose(self.uy))
        # ['tm 0:20']
        
        print('-'*8+"etkf: xtm"+'-'*8)
        print(array2string(tm[-5:-1,-5:-1],precision=5, suppress_small=True))
        svd_err = False
        return tm, svd_err
    
    def get_k(self):
        """ calculate analysis gain k matrix 
            Arguments:
                None
            Returns: 
                k matrix    :   array of [nobs, nobs]=[ny, ny]
            Notes:
                K=X'*U*W*(I+WT*W)^-1*VT*R^-0.5
        """

        inv_s=ones(self.ny, float) # filled with one 

        inv_s[0:self.n_wy]=self.md[0:self.n_wy]
        
        uw_nx, uw_ny=shape(self.uw)
        rny=self.ny
        if (uw_ny<rny):
            rny=uw_ny
        
        y_inv=self.uw[:,0:rny]*inv_s[0:rny]
        y_inv_v=dot(y_inv, self.vyt[0:rny, :]) #dimension of [ne, ny]
        k=dot(self.x, y_inv_v) # dimension [nx, ne] X [ne, ny]
        k=k/self.sqr
        return k
        
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
        print(shape(yobs), shape(self.mean_y))
        
        dy=yobs-self.mean_y # innovation
        print('etkf: dy in get_inc_m', dy[0:6], max(dy), min(dy))
        # print dy
        dy=dy/self.sqr
        inc_m=dot(self.vyt, transpose(dy))
        
        inc_m=squeeze(inc_m)
        inc_m=transpose(inc_m)
        inv_s=ones(self.ny, float) # filled with one 
        inv_s[0:self.n_wy]=self.md[0:self.n_wy]
        inc_m=inv_s*inc_m
        uw_nx, uw_ny=shape(self.uw) 
        rny=self.ny
        if (uw_ny<self.ny):
            rny=uw_ny
        
        inc_m=dot(self.uw[:,0:rny], transpose(inc_m[0:rny]))
        
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
        
        
    def get_h(self):
        """ H matrix: it could be wrong when other perturbations  is not added as required """
        h=zeros([self.ny, self.nx], float)
        for iy in range(self.ny):
            for ix in range(self.nx):
                dy=self.y[iy, :]
                dx=self.x[ix, :]
                sel_x=where(abs(dx)>100)
                if (len(sel_x)>0):
                    sel_x=squeeze(sel_x)
                    # if (ix==2):
                    #     print size(sel_x), dx[sel_x], dy[sel_x]
                    if (size(sel_x)==2):
                        dydx=dy[sel_x]/dx[sel_x]
                        h[iy, ix]=mean(dydx)
                    elif (size(sel_x)==1):
                        h[iy, ix]=dy[sel_x]/dx[sel_x]
                    else:
                        print('wrong too many', ix, sel_x)
                else:
                    print('wrong no one', ix)
                
        return h
    def get_avg_kernel(self):
        """ it could be wrong when other perturb is added """
        h=self.get_h()
        k=self.get_k()
        # print max(k[2,:])
        # print max(1.0e6*h[:, 2])
        
        avg_kernel=dot(k, h)
        print(avg_kernel[0, 0], avg_kernel[1,1], avg_kernel[2,2]) 
        return avg_kernel
    
    def get_avg_kernel_from_xa(self, xa, xf=None):
        """ it could be wrong when other perturb is added """
        import numpy.linalg as nlg
        if (xf==None):
            dxf=self.x
        else:
            dxf=array(xf)
            dxf=dxf/sqrt(self.xnorm)

        xcorf=dot(dxf, transpose(dxf))
        
        dxa=array(xa)
        dxa=dxa/sqrt(self.xnorm)
        xcora=dot(dxa, transpose(dxa))
        inv_xcorf=nlg.inv(xcorf)
        kh=dot(xcora, inv_xcorf)
        kh=identity(size(kh[0,:]))-kh
        return kh
    # print shape(xcor0)


 
        
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
    tm=ett.get_tm()
    #    print array2string(tma, precision=2, suppress_small=True)
    print('='*80)
    k=ett.get_k()
    print('Kalman Gain Matrix (size)',  shape(k))
    
    # print array2string(k, precision=2, suppress_small=True)
    xinc=dot(k, transpose(yobs-ett.mean_y)) 
    xan=ett.mean_x+xinc
    inc_m=ett.get_inc_m(yobs)
    
    xinc=ett.get_x_inc(inc_m)
    xan2=ett.mean_x+xinc
    
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
    
    
    
    
        
        
        
    
        
    
        
        
