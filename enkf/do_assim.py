from numpy import *
import time_module as tm
import etkf_cor as etc
import etkf_half as et

def do_assim(tdx, tdy, mean_x, mean_y, yobs, yerr, xref, xnorm, use_sparse=True):  
    
    ''' 
    Inversion parameters calculation
    
    Args:
        tdx: prior flux factor pertubations(delta_X)
        tdy: prior concentration pertubations(delta_Y)
        mean_x: prior flux factors(x)
        mean_y: model results with prior being input(H(x))
        yobs: observation results
        yerr: observation errors
        xref, xnorm: Calculation parameters 
        use_sparse: switch to use sparse matrix in calculation
    Returns:
        xinc: Increments to prior state
        xtm: update factors to pertubations
        inc_m: Increment parameters

    '''
    
    if (use_sparse):
        # using sparse matrix
        
        full_r = yerr.astype(float)
            
        ett = etc.etkf_c(mean_x, mean_y, tdx, tdy, full_r, \
                       xref = xref, xnorm = xnorm)
    else:
        # pass
        r = diag(yerr)
        ett=et.etkf_c(mean_x, mean_y, tdx, tdy, r, \
                          xref = xref, xnorm = xnorm)
            
    svd_err=True
    itry=0
        
    while (svd_err):
        xtm, svd_err = ett.get_tm()
        if (svd_err):
            itry=itry+1
            print('problem with svd', itry)    
            if (itry>5):
                # too many tries
                ne_try=size(xtm[0,:])
                # set inc_m as zero of (size_ne)
                inc_m=zeros(ne_try, float)
                break

                
            # reconstruct the matrix but with enlarged error covariance. 
            retry_enlarge_factor=(1.10)**itry
            
            if (use_sparse):
                
                
                full_r=retry_enlarge_factor*full_r
                
                full_r=full_r.astype(float)        
                print('shape--full-r:', shape(full_r))
                ett=etc.etkf_c(mean_x, mean_y, tdx, tdy, full_r, \
                               xref=xref, xnorm=xnorm)
            else:
                r=retry_enlarge_factor*r
                print('shape--r:', shape(full_r))       
                ett=et.etkf_c(mean_x, mean_y, tdx, tdy, r, \
                              xref=xref, xnorm=xnorm)
                    
        else:
            # if svd work
            inc_m=ett.get_inc_m(yobs)
            
            break
            
    xinc = dot(tdx, inc_m)
                
        
    return xinc, xtm, inc_m
    
    