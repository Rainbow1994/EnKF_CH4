#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
import sys
import os
import numpy as np
import numba as nb
import pandas as pd
import xarray as xr 
from datetime import datetime


#@nb.jit(nopython=True,parallel=True) 
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    
def field_sampling(lons, lats, data):
    
    lon =[]
    lat =[]
    for i_lon, i_lat in zip(lons, lats):
        lon.append(find_nearest(data.lon.values, i_lon))
        lat.append(find_nearest(data.lat.values, i_lat))

    ind_x = xr.DataArray(lon, dims='n')
    ind_y = xr.DataArray(lat, dims='n')
    samp_rlt = data.isel(lon=ind_x, lat=ind_y)
    return samp_rlt

#@nb.jit(nopython=True,parallel=True) 
def setup_pres_weight(pres):

    n = len(pres)
    p0= 1
    #p0 = np.full([n], 0.01, dtype=np.float)
    gross_p = pres[0] - p0
    p_mins = np.zeros(len(pres))
    p_mins[:-1]  = pres[1:]
    p_mins[-1] = p0
    delta_p = np.divide((pres - p_mins) , gross_p)
    delta_p = np.array(delta_p)
    return delta_p
    
def restart_nc_read(mod,\
                surf_pres,\
                lons,\
                lats):


    pres_levels = (mod['hyam']+ mod['hybm']*(surf_pres/100.0))
    mod['PEDGE_S_PSURF'] = pres_levels
    lon =[]
    lat =[]
    for i_lon, i_lat in zip(lons, lats):
        lon.append(find_nearest(mod.lon.values, i_lon))
        lat.append(find_nearest(mod.lat.values, i_lat))
        
    ind_x = xr.DataArray(lon, dims='n')
    ind_y = xr.DataArray(lat, dims='n')
    mod_sel = mod.isel(time = 0, lon=ind_x, lat=ind_y)
    mod_sel = mod_sel.transpose('ilev', 'n','lev')
    
    return mod_sel

#@nb.jit(nopython=True,parallel=True) 
def get_gc_xch4(tropp, gc_pres, gc_field):
    '''
    get stratopheric xch4 using tropopause
    Args:
        tropp: 1-d tropopause[n]
        gc_pres: 2-d pressure levels of GC[lev, n]
        gc_field: 2-d gc tracer field[lev, n]
    Return:
        gc_xch4: 1-d stratospheric xch4[n]
    '''
    nobs = len(tropp)
    gc_xch4 = []

    for i in range(nobs):
        tropo = tropp[i]
        pres = gc_pres[i,:]   
        gc_mask = (pres*100<tropo)&(pres>1.)
        pres1 = pres[gc_mask]

        gc_ch4_pl =  gc_field[i,:][gc_mask]
        gc_pres_wgt = setup_pres_weight( pres1 )
        gc_xch4.append(np.nansum(gc_pres_wgt*gc_ch4_pl)*1e9)
    gc_xch4 = np.array(gc_xch4)
    return gc_xch4


#@nb.jit(nopython=True,parallel=True) 
def find_cor_fac(lat, lat_bins, cor_factors):
    '''
    Notes: lat_bins should be strictly monotonically increasing
    '''
    for lat_v, cor_f in zip(lat_bins[::-1], cor_factors[::-1]): 
        if lat > lat_v:
            return cor_f
        else:
            continue
            

#@nb.jit(nopython=True,parallel=True) 
def gc_stra_cor(res_field, tropp_layer, latitudes, cor_factors, lat_bins):
    
    nlev, nlat, nlon = np.shape(res_field)
    new_res = res_field.copy() 
    lev_st = np.min(tropp_layer)
    for ilev in np.arange(int(lev_st-1),int(nlev)):
        
        for ilat, lat in zip(np.arange(nlat), latitudes):
            pos = (tropp_layer[ilat,:] <= ilev+1)
            if np.all(~pos):
                continue
            cor_value = find_cor_fac(lat, lat_bins, cor_factors)
            
            if cor_value == None:
                continue
                
            new_res[ilev, ilat, pos] = res_field[ilev, ilat, pos] * cor_value
            
    return new_res



def get_zonal_cor(ace_xch4, gc_xch4, longitude, latitude, lat_bins):
    '''
    get gc zonal bias correction factors according latitude bins
    
    Args:
        ace_xch4: ACE-FTS measurements
        gc_xch4: GEOS-Chem sampling stratospheric xch4
        lontitude: lontitude array of measurement locations
        latitude: latitude array of measurement locations
        lat_bins: latitude bins for calculation
    Return:
        cor_factor: zontal mean bias correction factors according to latitude bins input
    '''
    
    ## step: create new dataframe
    data = pd.DataFrame(columns = ['ace_xch4', 'gc_xch4','latitude','longitude'])
    data['ace_xch4'] = np.array(ace_xch4)
    data['gc_xch4'] = np.array(gc_xch4)
    data['longitude'] = np.array(longitude)
    data['latitude'] = np.array(latitude)
    
    ## step2: filter data

    data['diff'] = data['ace_xch4']-data['gc_xch4']
    data1 = data[abs(data['diff'])< 3*data['diff'].std()].reset_index()
    data1['cor_factor'] = data1['ace_xch4']/ data1['gc_xch4']
    
    ## step3: calculate zontal averaged factor
    cor_factor = data1['cor_factor'].groupby(pd.cut(data1['latitude'], bins=lat_bins)).mean()
    
    return cor_factor



def Restart_bias_cor(org_res, \
                     ace_file, \
                     tropp_file, \
                     sf_pres_file, \
                     bins = [-90,-80, -70, -60, -50, -40, -30, -15, 0, 15, 30, 40, 50, 60, 70, 80, 90]):
    '''
        Stratopheric Bias Correction for GC restart files
    '''
    ## Step1: open orginal restart file and get timestamp
#     org_res = xr.open_dataset(restart_file).load()
    
    sea = org_res['time'].dt.season.values[0]
    year = org_res['time'].dt.year.values[0]
    
    ## Step2: get ACE-FTS climatological xch4 field & tropopause field
    
    cor_fac = xr.open_dataset(ace_file).load()

    ace_xch4 = (cor_fac['xch4'].sel(season = sea)*\
                cor_fac['xch4_factors'].sel(season=sea, year = year).values)#.isel(season =0)
    ace_xch4_stacked = ace_xch4.stack(n=['lat','lon'])
    lons = ace_xch4_stacked[ace_xch4_stacked.notnull()]['lon'].values
    lats = ace_xch4_stacked[ace_xch4_stacked.notnull()]['lat'].values
    
    nobs = len(lons)
    
    tropp = xr.open_dataset(tropp_file).sel(time = org_res.time).isel(time =0)
    
    ## Step3: Sampling for Tropopause and Restart fields
    tropp_ns =field_sampling(lons, lats, tropp)
    
    pres_surf = xr.open_dataset(sf_pres_file)['PS'].isel(time =0)
    res_ns = restart_nc_read(mod = org_res,\
                             surf_pres = pres_surf,\
                             lons = lons,\
                             lats =lats)
    
    ## Step4: get stratopheric xch4 of restart fiels
    gc_xch4 = get_gc_xch4(tropp = tropp_ns['pres'].values,\
                          gc_pres = res_ns['PEDGE_S_PSURF'].values,\
                          gc_field = res_ns['SpeciesRst_CH4'].values)
    ## Step5: get zonal bias correction factors
    
#     bins = [-90,-80, -70, -60, -50, -40, -30, -15, 0, 15, 30, 40, 50, 60, 70, 80, 90]
    cor_factors = get_zonal_cor(ace_xch4 = np.array(ace_xch4_stacked[ace_xch4_stacked.notnull()].values),\
                                gc_xch4 = np.array(gc_xch4),\
                                longitude = lons,\
                                latitude = lats,\
                                lat_bins = bins)  
    
    ## Step6: apply cor factors for restart fields high than Tropopause
    
    new_res = gc_stra_cor(res_field = org_res['SpeciesRst_CH4'].values.squeeze(), \
                      tropp_layer = org_res['Met_TropLev'].values.squeeze(), \
                      latitudes = org_res['lat'].values, \
                      cor_factors = cor_factors.values, \
                      lat_bins = bins[:-1])
    
    org_res['org_SpeciesRst_CH4'] = org_res['SpeciesRst_CH4'].copy()
    org_res['SpeciesRst_CH4'][dict(time = 0)] = new_res
    
    return org_res


if (__name__ == '__main__'):
    new_res = Restart_bias_cor(restart_file = '/data/shzhu/org_ch4/global_2x25/Restart/Restart.20100101_0000z.nc4' ,\
                               ace_file = '/data/shzhu/measurements/ACE-FTS/ACE_Xch4_Climatology.nc',\
                               tropp_file = '/data/geos-chem/tropp/Tropp.mon.2010-2020.nc',\
                               sf_pres_file='/data/geos-chem/ctm/GEOS_2x2.5/MERRA2/2010/01/MERRA2.20100101.I3.2x25.nc4')
    new_res.to_netcdf('test.nc')
