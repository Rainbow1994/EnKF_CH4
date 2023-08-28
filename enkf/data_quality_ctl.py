import numpy as np
import xarray as xr
from scipy.stats import sem, t
import sample_ch4_gosat as scg


def data_quality_ctl(smpl_data):
    '''
    Data quality control for sampled data(data filter)
    
    '''
    
    ## step1: Acoording to obs and results

    ns = (smpl_data['obs'] > 0.0) & (smpl_data['oerr']>0.0) & (smpl_data['hm0'] > 0.0)
    dqc1 = smpl_data.sel(n=ns)
    
    nl =  np.all( dqc1['hm'] > 0, axis =1) 
    dqc2 = dqc1.sel(n=nl)
    
    ## step2: Acoording to difference of obs and model results
    
    dif = abs(dqc2['obs'] - dqc2['hm0'])
    n = len(dif.n)
    m = dif.mean()
    #std_err = sem(dif_data)
    #h = std_err * t.ppf((1 + confidence) / 2, n - 1)
    
    #nd = (dif <= m+h) & (dif>=m-h)
    nd = (dif <= 3* m )
    dqc_data = dqc2.sel(n=nd)
    
    if dqc_data['obs'].size == 0:
        return False
    else:
        dqc_data.attrs['Description'] =  'Matched sample result after data quality control'
        return dqc_data
    

def tccon_quality_ctl(smpl_data):
    '''
    Data quality control for tccon sampled data(data filter)
    
    '''
    
    ## step1: Acoording to obs and results

    ns = (smpl_data['obs'] > 0.0) & (smpl_data['oerr']>0.0) & (smpl_data['hm0'] > 0.0)
    dqc1 = smpl_data.sel(ns=ns)
    
    nl =  np.all( dqc1['hm'] > 0, axis =1) 
    dqc2 = dqc1.sel(ns=nl)
    
    ## step2: Acoording to difference of obs and model results
    
    dif = abs(dqc2['obs'] - dqc2['hm0'])
    n = len(dif.ns)
    m = dif.mean()
    #std_err = sem(dif_data)
    #h = std_err * t.ppf((1 + confidence) / 2, n - 1)
    
    #nd = (dif <= m+h) & (dif>=m-h)
    nd = (dif <= 3* m )
    dqc_data = dqc2.sel(ns=nd)
    
    dqc_data.attrs['Description'] =  'Matched sample result after data quality control'
    return dqc_data

if (__name__=="__main__"): 

    doy = 37
    yyyy = 2015
    cur_step = 0
    samp_data = scg.get_hm(cur_step, doy, yyyy) 
    result = data_quality_ctl(samp_data)
    
    ## exec time : 41 ms
    ## samping + quality control : 40.1s
    