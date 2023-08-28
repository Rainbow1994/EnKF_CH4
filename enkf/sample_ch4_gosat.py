#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
#import sys
#sys.path.append('/home/shzhu/enkf_ch4/inv_code')
import numpy as np
import numba as nb
import pandas as pd
import xarray as xr 
from xbpch import open_bpchdataset
from scipy import interpolate
from scipy.interpolate import interp1d
import geos_chem_def as gcdf
import time_module as tm
import rerun_geos as rg
import os
from global_land_mask import globe


def sat_read(filename):
    sat = xr.open_dataset(filename)
    sat = sat.assign_coords(m=(sat.m), n=(sat.n))
    ns = (sat['sensor_zenith_angle'] < 75) & (sat['xch4_quality_flag'] == 0)& (sat['xch4'] > 0)
    sat_data = sat.sel(n=ns)
    
    nl= np.all( sat_data['pressure_levels'] > 0, axis =1)
    sat_data1 = sat_data.sel(n=nl)
    
    return sat_data1

#@nb.jit(nopython=True,parallel=True) 
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def setup_daily_feild_nc(model_filename,
                      lons,\
                     lats):
    
    model_ch4 = xr.open_dataset(model_filename).load()
    model_ch4.close()

    lon =[]
    lat =[]
    for i_lon, i_lat in zip(lons, lats):
        lon.append(find_nearest(model_ch4.lon.values, i_lon))
        lat.append(find_nearest(model_ch4.lat.values, i_lat))
        
    ind_x = xr.DataArray(lon, dims='n')
    ind_y = xr.DataArray(lat, dims='n')
    model_ch4 = model_ch4.isel(lon=ind_x, lat=ind_y)
#     dr_ch4 = []
#     for var in model_ch4.data_vars:
#         #ch4 = model_ch4[var].interp(lon=lons, lat=lats)
#         ch4 = model_ch4[var].isel(lon=ind_x, lat=ind_y)
# #         ch4.values = pd.DataFrame(ch4.values).fillna(method='ffill',\
# #                                                       axis=0).fillna(method='bfill',\
# #                                                                      axis=0).fillna(method='ffill',\
# #                                                                                     axis=1).fillna(method='bfill',\
# #                                                                                                    axis=1).values
#         dr_ch4.append(ch4)
#     model_ch4 = xr.merge(dr_ch4)
    return model_ch4

def setup_daily_feild(model_filename,
                      diaginfo_filename,
                      tracerinfo_filename,
                      lons,\
                     lats):
    
    model_ch4 = open_bpchdataset(model_filename,\
                                 diaginfo_file= diaginfo_filename,\
                                 tracerinfo_file=tracerinfo_filename)


    lon =[]
    lat =[]
    for i_lon, i_lat in zip(lons, lats):
        lon.append(find_nearest(model_ch4.lon.values, i_lon))
        lat.append(find_nearest(model_ch4.lat.values, i_lat))
        
    ind_x = xr.DataArray(lon, dims='n')
    ind_y = xr.DataArray(lat, dims='n')
    model_ch4 = model_ch4.isel(lon=ind_x, lat=ind_y)
#     dr_ch4 = []
#     for var in model_ch4.data_vars:
#         #ch4 = model_ch4[var].interp(lon=lons, lat=lats)
#         ch4 = model_ch4[var].isel(lon=ind_x, lat=ind_y)
# #         ch4.values = pd.DataFrame(ch4.values).fillna(method='ffill',\
# #                                                       axis=0).fillna(method='bfill',\
# #                                                                      axis=0).fillna(method='ffill',\
# #                                                                                     axis=1).fillna(method='bfill',\
# #                                                                                                    axis=1).values
#         dr_ch4.append(ch4)
#     model_ch4 = xr.merge(dr_ch4)
    return model_ch4

def setup_nc_daily_feild(combine_en,\
                         lons,\
                         lats):
    
    

    lon =[]
    lat =[]
    for i_lon, i_lat in zip(lons, lats):
        lon.append(find_nearest(combine_en.lon.values, i_lon))
        lat.append(find_nearest(combine_en.lat.values, i_lat))
        
    ind_x = xr.DataArray(lon, dims='n')
    ind_y = xr.DataArray(lat, dims='n')
    model_ch4 = combine_en.isel(lon=ind_x, lat=ind_y)
#     dr_ch4 = []
#     for var in combine_en.data_vars:
#         ch4 = combine_en[var].isel(lon=ind_x, lat=ind_y)        
# #         ch4 = combine_en[var].interp(lon=lons, lat=lats)
# #         ch4.values = pd.DataFrame(ch4.values).fillna(method='ffill',\
# #                                                   axis=0).fillna(method='bfill',\
# #                                                                  axis=0).fillna(method='ffill',\
# #                                                                                 axis=1).fillna(method='bfill',\
# #                                                                                                axis=1).values
#         dr_ch4.append(ch4)
#     model_ch4 = xr.merge(dr_ch4)
    return model_ch4

def setup_pres_weight(pres):
    
    n = len(pres)
    p0 = np.full([n], 0.01, dtype=np.float)
    gross_p = pres[:, 0] - p0
    p_mins = np.concatenate((pres[:, 1:], p0[:, None]), axis=1)
    delta_p = np.divide((pres - p_mins) , gross_p[:, None])
    delta_p = np.array(delta_p)
    return delta_p

#@nb.jit(nopython=True,parallel=True) 
def vertical_intpl_2d(po_x, og_x, og_y):
    po_y = np.zeros_like(po_x)
    for i in range(len(og_x[:, 0])):
        #po_y[i, :] = np.interp(po_x[i, :], og_x[i, :],og_y[i, :])
        ff = interp1d(og_x[i, :], og_y[i, :], fill_value="extrapolate")
        po_y[i, :] = ff(po_x[i, :])

    return po_y


def get_xgp(obs_pres, obs_ave, model_pres, model_data, pres, delta_p, method='linear'):
    '''
    interpolate observation average kernel & model results and caculate xgp 
    
    method: 
        'linear': linear interpolation
        'log'：Logarithmic interpolation
       
    '''
    
    
    if (method == 'linear'):
        pres = pres
        obs_pres = obs_pres
        model_pres = model_pres
        
    elif (method == 'log'):
        if (np.all(pres > 0)) and (np.all(obs_pres>0)) and (np.all(model_pres>0)):
            pres = np.log10(pres)
            obs_pres = np.log10(obs_pres)
            model_pres = np.log10(model_pres)
        else:
            print(np.all(pres > 0), np.all(obs_pres>0), np.all(model_pres>0))
            raise ValueError("PRESSURE VALUES ERROR")
    else:
        raise AssertionError('ERROR WITH INTERPOLATION METHOD')
        
    
    ave_ker = vertical_intpl_2d(pres, obs_pres, obs_ave)
    md_interp = vertical_intpl_2d(pres, model_pres, model_data)
    
    xgp = (ave_ker* md_interp*delta_p).sum(axis=1)
    
    return xgp

def combine_xbpch(en1, tra_num1, en2, tra_num2):
    
    en3 = xr.Dataset()
    tracers = ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA']
    for i_tra in tracers:
        for i in range(tra_num2):
            new_name = i_tra + str(i + 1 + tra_num1)
            old_name = i_tra + str(i + 1)

            #en2 = en2.rename_vars({old_name:new_name})#[new_name].values = en1[old_name].values
            en3[new_name] = en2[old_name]
#    en1 = en1.drop(['IJ_AVG_S_', 'time_bnds'])
#    en2 = en2.drop(['IJ_AVG_S_CH4', 'PEDGE_S_PSURF'])
    en = xr.merge([en1,en3])
    en.attrs = en1.attrs
    en.astype(np.float32)
    
    return en


def get_hm(cur_step, doy, yyyy, sampl_file):
    '''
    read single and ensemble geos_chem run diagnostic ouput
    resample them to satellite observation sites
    according to layer pressure and its weight function to caculate XCH4
    
    '''
    
    ## step1: get time information
    tst=tm.doy_to_utc(doy, 0, yyyy)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '') 
    
    
    ## step2: read satellite observation data
    #sat_fn = 'ESACCI-GHG-L2-CH4-GOSAT-OCPR-'+tst[0:8]+'-fv7.0.nc'
    sat_fn = 'UoL-GHG-L2-CH4-GOSAT-OCPR-'+tst[0:8]+'-fv9.0.nc'
    sat_fn = gcdf.sat_path + str(yyyy) + '/'+ sat_fn
    print(sat_fn)
    if os.path.exists(sat_fn):
        #trying to open a file in read mode
        sat_data = sat_read(sat_fn)
    else:
        return False
    
   
    
    ### step3: get xgp0 
    sat_ave = sat_data['xch4_averaging_kernel'].values
    xgp0 = ((1 - sat_ave) * sat_data['pressure_weight'] *\
        sat_data['ch4_profile_apriori']).sum(dim='m')
    xgp0 = np.array(xgp0.values)
    
   
    ## step4: read geos-chem single tracer(rerun) output and do resample
        # sampling girds
#     lons = xr.DataArray(sat_data['longitude'].values, dims='n')
#     lats = xr.DataArray(sat_data['latitude'].values, dims='n')
    lons = sat_data['longitude'].values
    lats = sat_data['latitude'].values
    
        # output file information
    enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step, 1, 2)
    sp_fn = 'ts_satellite.'+enaf+'.'+tst[0:8]+'.bpch'
    sp_fn = gcdf.sp_diagn_path+'ND51/'+sp_fn
    sp_diag = gcdf.sp_diagn_path +'diaginfo_sp.dat'
    sp_tra = gcdf.sp_diagn_path +'tracerinfo_sp.dat'

        # resampling bpch output
    sp_mod = setup_daily_feild(sp_fn,\
                               sp_diag,\
                               sp_tra,\
                               lons,\
                              lats) 
    

    ## step5: combine observation and model pressure layers and delta pres

    model_pres = sp_mod['PEDGE_S_PSURF'].values
    sat_pres = sat_data['pressure_levels'].values
    pres = np.concatenate((model_pres, sat_pres), 1)
    pres = np.sort(pres, axis=1)[:, ::-1]
    delta_p = setup_pres_weight(pres)
    

    ## step6: caculation xgp result of single tracers[hm0]
    sgl_data = sp_mod['IJ_AVG_S_CH4'].values 
    
    sp_xgp = get_xgp(sat_pres, sat_ave, model_pres, sgl_data, pres, delta_p, method='log')
    hm0 = xgp0 + sp_xgp
    hm0 = np.array(hm0)
    
    
    ## step7: caculation xgp result of ensemble tracers[hm]
    em_step_list,em_name_list  = rg.get_ensemble_list(yyyy, doy, \
                                                    entry_table_name= gcdf.entry_table)
    tag_xgp =[]
    var_list= [''.join(['IJ_AVG_S_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    for em_step, em_name in zip(em_step_list,em_name_list):
        #em_st = int(em_name[10:12])
        #em_end = int(em_name[17:19])
        #em_enaf = r'ST%3.3d.EN%4.4d-EN%4.4d' % (em_step, em_st, em_end)
        em_fn_0 ='ts_satellite.'+ em_name +'.'+tst[0:8]+'.bpch'
        em_fn = gcdf.en_diagn_path+ '/' + str(yyyy) + '/'+em_fn_0
        if not os.path.exists(em_fn):
            em_fn = gcdf.en_diagn_path+ '/' + str(yyyy - 1) + '/'+em_fn_0
        em_diag = gcdf.en_diagn_path  + '/' 'diaginfo_em.dat'
        em_tra = gcdf.en_diagn_path  + '/''tracerinfo_em.dat'
        em_mod = setup_daily_feild(em_fn,\
                                    em_diag,\
                                    em_tra,\
                                    lons,\
                                    lats)
    
        for var in var_list:
                
            em_data = em_mod[var].values
            em_xgp = get_xgp(sat_pres, sat_ave, model_pres, em_data, pres, delta_p, method='log')
            tag_xgp.append(em_xgp + xgp0)
                
    hm = np.array(tag_xgp)
    hm= hm.transpose()
    
    
    ## step8: integration of sampling results
    ns, ntag = np.shape(hm)    
    samp_data = xr.Dataset({'hm0':(['n'], hm0),
                            'hm':(['n','t'], hm)},
                            coords= {'n': sat_data.n,
                                    't':  range(ntag),
                                    'ref_time' : tm.doy_to_utc(doy, 0, yyyy)})
                                 
        
    
    samp_data['obs']= sat_data['xch4']
    samp_data['oerr'] = sat_data['xch4_uncertainty']
    samp_data['lon'] = sat_data['longitude']
    samp_data['lat'] = sat_data['latitude']
    samp_data['time'] = sat_data['time']

    
        
    ## step9: add some description to variables    
    samp_data.attrs['Description'] = 'Matched sample result without data quality control '
        
    samp_data['hm0'].attrs['long_name'] = 'model_xch4'
    samp_data['hm0'].attrs['standard_name'] = 'sample result of GEOS-Chem single run'
    samp_data['hm0'].attrs['units'] = '1e-9'
    samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
    samp_data['hm'].attrs['standard_name'] = 'sample result of GEOS-Chem pertubation run'
    samp_data['hm'].attrs['units'] = '1e-9'
#     samp_data.to_netcdf(sampl_file)
        
    return samp_data 

def get_hm2(cur_step, doy, yyyy, sampl_file):
    '''
    read single and ensemble geos_chem run diagnostic ouput
    resample them to satellite observation sites
    according to layer pressure and its weight function to caculate XCH4
    ### note: read seperated ensemble tagrun!
    '''
    
    ## step1: get time information
    tst=tm.doy_to_utc(doy, 0, yyyy)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '') 
    
    
    ## step2: read satellite observation data
    #sat_fn = 'ESACCI-GHG-L2-CH4-GOSAT-OCPR-'+tst[0:8]+'-fv7.0.nc'
    sat_fn = 'UoL-GHG-L2-CH4-GOSAT-OCPR-'+tst[0:8]+'-fv9.0.nc'
    sat_fn = gcdf.sat_path + str(yyyy) + '/'+ sat_fn
    
    if os.path.exists(sat_fn):
        #trying to open a file in read mode
        sat_data = sat_read(sat_fn)
    else:
        return False
    
   
    
    ### step3: get xgp0 
    sat_ave = sat_data['xch4_averaging_kernel'].values
    xgp0 = ((1 - sat_ave) * sat_data['pressure_weight'] *\
        sat_data['ch4_profile_apriori']).sum(dim='m')
    xgp0 = np.array(xgp0.values)
    
   
    ## step4: read geos-chem single tracer(rerun) output and do resample
        # sampling girds
#     lons = xr.DataArray(sat_data['longitude'].values, dims='n')
#     lats = xr.DataArray(sat_data['latitude'].values, dims='n')  
    lons = sat_data['longitude'].values
    lats = sat_data['latitude'].values
        # output file information
    #enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step, 1, 2)
    #sp_fn = 'ts_satellite.'+enaf+'.'+tst[0:8]+'.bpch'
    enaf=r'EN%4.4d-EN%4.4d' % (1, 2)
    sp_fn = 'ts_satellite.4x5.'+enaf+'.'+tst[0:8]+'.bpch'
    sp_fn = gcdf.sp_diagn_path+'ND51/'+sp_fn
    sp_diag = gcdf.sp_diagn_path +'diaginfo_sp.dat'
    sp_tra = gcdf.sp_diagn_path +'tracerinfo_sp.dat'

        # resampling bpch output
    sp_mod = setup_daily_feild(sp_fn,\
                               sp_diag,\
                               sp_tra,\
                               lons,\
                              lats) 
    

    ## step5: combine observation and model pressure layers and delta pres

    model_pres = sp_mod['PEDGE_S_PSURF'].values
    sat_pres = sat_data['pressure_levels'].values
    pres = np.concatenate((model_pres, sat_pres), 1)
    pres = np.sort(pres, axis=1)[:, ::-1]
    delta_p = setup_pres_weight(pres)
    

    ## step6: caculation xgp result of single tracers[hm0]
    sgl_data = sp_mod['IJ_AVG_S_CH4'].values 
    
    sp_xgp = get_xgp(sat_pres, sat_ave, model_pres, sgl_data, pres, delta_p, method='log')
    hm0 = xgp0 + sp_xgp
    hm0 = np.array(hm0)
    
    
    ## step7: caculation xgp result of ensemble tracers[hm]
    em_step_list1,em_name_list1  = rg.get_ensemble_list(yyyy, doy, \
                                                    entry_table_name= gcdf.entry_table1)
    em_step_list2,em_name_list2  = rg.get_ensemble_list(yyyy, doy, \
                                                    entry_table_name= gcdf.entry_table2)
    tag_xgp = []
    tra_num1 = gcdf.tra_num1
    tra_num2 = gcdf.tra_num2
    var_list= [''.join(['IJ_AVG_S_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    for em_step, em_name1, em_name2 in zip(em_step_list1,em_name_list1, em_name_list2):
        
        em_fn ='ts_satellite.'+ em_name1 +'.'+tst[0:8]+'.bpch'
        em_fn1 = gcdf.en_diagn_path+ '/' + str(yyyy) + '/'+em_fn
        if not os.path.exists(em_fn1):
            em_fn1 = gcdf.en_diagn_path+ '/' + str(yyyy - 1) + '/'+em_fn
            
        em_fn ='ts_satellite.'+ em_name2 +'.'+tst[0:8]+'.bpch'
        em_fn2 = gcdf.en_diagn_path+ '/' + str(yyyy) + '/'+em_fn
        
        if not os.path.exists(em_fn2):
            em_fn2 = gcdf.en_diagn_path+ '/' + str(yyyy - 1) + '/'+em_fn
            
        en1 = open_bpchdataset(em_fn1,\
                               diaginfo_file = gcdf.em_diag,\
                               tracerinfo_file = gcdf.em_tra1)

        en2 = open_bpchdataset(em_fn2,\
                               diaginfo_file = gcdf.em_diag,\
                               tracerinfo_file = gcdf.em_tra2)
        

        combine_em = combine_xbpch(en1, tra_num1, en2, tra_num2)

        em_mod = setup_nc_daily_feild(combine_em,\
                                      lons,\
                                      lats)
    
        for var in var_list:
                
            em_data = em_mod[var].values
            em_xgp = get_xgp(sat_pres, sat_ave, model_pres, em_data, pres, delta_p, method='log')
            tag_xgp.append(em_xgp + xgp0)
                
    hm = np.array(tag_xgp)
    hm= hm.transpose()
    
    
    ## step8: integration of sampling results
    ns, ntag = np.shape(hm)    
    samp_data = xr.Dataset({'hm0':(['n'], hm0),
                            'hm':(['n','t'], hm)},
                            coords= {'n': sat_data.n,
                                    't':  range(ntag),
                                    'ref_time' : tm.doy_to_utc(doy, 0, yyyy)})
                                 
        
    
    samp_data['obs']= sat_data['xch4']
    samp_data['oerr'] = sat_data['xch4_uncertainty']
    samp_data['lon'] = sat_data['longitude']
    samp_data['lat'] = sat_data['latitude']
    samp_data['time'] = sat_data['time']

    
        
    ## step9: add some description to variables    
    samp_data.attrs['Description'] = 'Matched sample result without data quality control '
        
    samp_data['hm0'].attrs['long_name'] = 'model_xch4'
    samp_data['hm0'].attrs['standard_name'] = 'sample result of GEOS-Chem single run'
    samp_data['hm0'].attrs['units'] = '1e-9'
    samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
    samp_data['hm'].attrs['standard_name'] = 'sample result of GEOS-Chem pertubation run'
    samp_data['hm'].attrs['units'] = '1e-9'
#     samp_data.to_netcdf(sampl_file)
        
    return samp_data 

def get_hm0(cur_step, doy, yyyy, hm, sampl_file):
    '''
    read single and ensemble geos_chem run diagnostic ouput
    resample them to satellite observation sites
    according to layer pressure and its weight function to caculate XCH4
    ### note: read seperated ensemble tagrun!
    '''
    
    ## step1: get time information
    tst=tm.doy_to_utc(doy, 0, yyyy)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '') 
    
    
    ## step2: read satellite observation data
    #sat_fn = 'ESACCI-GHG-L2-CH4-GOSAT-OCPR-'+tst[0:8]+'-fv7.0.nc'
    sat_fn = 'UoL-GHG-L2-CH4-GOSAT-OCPR-'+tst[0:8]+'-fv9.0.nc'
    sat_fn = gcdf.sat_path + str(yyyy) + '/'+ sat_fn
    
    if os.path.exists(sat_fn):
        #trying to open a file in read mode
        sat_data = sat_read(sat_fn)
    else:
        return False
    
    
    
    
    ### step3: get xgp0 
    sat_ave = sat_data['xch4_averaging_kernel'].values
    xgp0 = ((1 - sat_ave) * sat_data['pressure_weight'] *\
        sat_data['ch4_profile_apriori']).sum(dim='m')
    xgp0 = np.array(xgp0.values)
    
   
    ## step4: read geos-chem single tracer(rerun) output and do resample
        # sampling girds
#     lons = xr.DataArray(sat_data['longitude'].values, dims='n')
#     lats = xr.DataArray(sat_data['latitude'].values, dims='n')  
    lons = sat_data['longitude'].values
    lats = sat_data['latitude'].values
    
    
        # output file information
    enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step, 1, 2)
    sp_fn = 'ts_satellite.'+enaf+'.'+tst[0:8]+'.bpch'
    sp_fn = gcdf.sp_diagn_path+'ND51/'+sp_fn
    sp_diag = gcdf.sp_diagn_path +'diaginfo_sp.dat'
    sp_tra = gcdf.sp_diagn_path +'tracerinfo_sp.dat'

        # resampling bpch output
    sp_mod = setup_daily_feild(sp_fn,\
                               sp_diag,\
                               sp_tra,\
                               lons,\
                              lats) 
    

    ## step5: combine observation and model pressure layers and delta pres

    model_pres = sp_mod['PEDGE_S_PSURF'].values
    sat_pres = sat_data['pressure_levels'].values
    pres = np.concatenate((model_pres, sat_pres), 1)
    pres = np.sort(pres, axis=1)[:, ::-1]
    delta_p = setup_pres_weight(pres)
    

    ## step6: caculation xgp result of single tracers[hm0]
    sgl_data = sp_mod['IJ_AVG_S_CH4'].values 
    
    sp_xgp = get_xgp(sat_pres, sat_ave, model_pres, sgl_data, pres, delta_p, method='log')
    hm0 = xgp0 + sp_xgp
    hm0 = np.array(hm0)
    
    
    
    ## step8: integration of sampling results
    ns, ntag = np.shape(hm)    
    samp_data = xr.Dataset({'hm0':(['n'], hm0),
                            'hm':(['n','t'], hm.values)},
                            coords= {'n': sat_data.n,
                                    't':  range(ntag),
                                    'ref_time' : tm.doy_to_utc(doy, 0, yyyy)})
                                 
        
    
    samp_data['obs']= sat_data['xch4']
    samp_data['oerr'] = sat_data['xch4_uncertainty']
    samp_data['lon'] = sat_data['longitude']
    samp_data['lat'] = sat_data['latitude']
    samp_data['time'] = sat_data['time']
    
    ### remove Tibet concentration 
    tibet = setup_daily_feild_nc(model_filename='/data/geos-chem/mask/Tibet.nc',
                             lons = samp_data['lon'].values,\
                             lats = samp_data['lat'].values)
    print('original_obs_number:', len(samp_data['n'].values))
    ll = tibet['MASK']==0
    samp_data1 = samp_data.sel(n = ll)
    print('obs_number_after_filter:', len(samp_data1['n'].values))
    
    ## selected land data
    lon = samp_data['lon'].values
    lat = samp_data['lat'].values
    is_on_land = globe.is_land(lat, lon)
    samp_data = samp_data.sel(n = is_on_land)
    if len(samp_data['lon'].values) == 0:
        return False   
    ## step9: add some description to variables    
    samp_data.attrs['Description'] = 'Matched sample result without data quality control '
        
    samp_data['hm0'].attrs['long_name'] = 'model_xch4'
    samp_data['hm0'].attrs['standard_name'] = 'sample result of GEOS-Chem single run'
    samp_data['hm0'].attrs['units'] = '1e-9'
    samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
    samp_data['hm'].attrs['standard_name'] = 'sample result of GEOS-Chem pertubation run'
    samp_data['hm'].attrs['units'] = '1e-9'
    #samp_data.to_netcdf(sampl_file)
        
    return samp_data1 

if (__name__=="__main__"): 

#     doy = 33
#     yyyy = 2015
    
    for yyyy in np.arange(2014, 2015):
        days = 365
        if (yyyy%4 == 0 and yyyy%400 != 0):
            days = 366
        for idoy in range(days):
            doy =  idoy + 1
            time = tm.doy_to_utc(doy, 0, yyyy)
            sampl_file = '/data/shzhu/org_ch4/global_4x5/gosat_samp/gosat_sampl' +  "." + time[:10] + ".nc"
    #samp_data = get_hm(cur_step, doy, yyyy)
            print(idoy, yyyy, time)
            samp_data = get_hm2(23, doy, yyyy,sampl_file)
    #print(samp_data)
    ## exec time : 38.7s
            