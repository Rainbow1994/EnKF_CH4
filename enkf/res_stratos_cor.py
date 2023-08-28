#! /home/shzhu/apps/miniconda3/envs/idp/bin/python


### readme####
# apply factors to ats-fts field
import sys
import os
import numpy as np
import numba as nb
import pandas as pd
import xarray as xr 
from datetime import datetime

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

def setup_pres_weight(pres):
    
    n = len(pres)
#     print(pres)
    p0= pres[-1]
    #p0 = np.full([n], 0.01, dtype=np.float)
    gross_p = pres[0] - p0
#     p_mins = np.zeros(len(pres))
#     p_mins[:-1]  = pres[1:]
#     p_mins[-1] = p0
    delta_p = np.divide((pres[:-1]- pres[1:] ) , gross_p)
    delta_p = np.array(delta_p)
    return delta_p


    
    

def cor_res_field(res_field, zonal_ch4, Troplev, lat, zonal_fac, mon_factor):
    '''
    restart field correction using ace-fts strtospheric zonal field
    '''
    nlev, nlat, nlon = np.shape(res_field)
    mod_field = res_field.copy()
    
        
    for ilat, llat in zip(range(nlat),lat):
        if np.all(np.isnan(zonal_ch4[:, ilat])):
            continue
        else:
            ace_ch4_pl = zonal_ch4[:, ilat]

        lat_bins = zonal_fac['lat_bins'].values
        fac_idx = find_closest(lat_bins, llat)
        if ~np.isnan(zonal_fac.isel(lat_bins = fac_idx).values):
            ace_ch4_pl = zonal_fac.isel(lat_bins = fac_idx).values* ace_ch4_pl
            
        fac = mon_factor.isel(lat_bins = fac_idx).values
        if np.isnan(fac):
            fac = 0.0
        
#         if llat < 45:
#             fac = np.abs(llat)/200
#         else:
#             fac = np.abs(llat)/100
        
        mean_lon = np.nanmean(res_field[33,ilat, :]-ace_ch4_pl[33])
        std_lon = np.nanstd(res_field[33,ilat, :]-ace_ch4_pl[33])
        for ilon in range(nlon):
            pres_lev = np.int(Troplev[ilat, ilon]) -1
            top = 34
            if pres_lev>=top:
                continue
            else:
                fac1 = fac 
                if abs(res_field[pres_lev,ilat, ilon]-ace_ch4_pl[pres_lev])>(mean_lon +3*std_lon):
                    fac1 = fac1/10.0
                else:
                    fac1 = fac1
                
                mod_field[pres_lev:top,ilat, ilon] = (1.0-fac1)*res_field[pres_lev:top,ilat, ilon]+ fac1*ace_ch4_pl[pres_lev:top]
   
                    


    return mod_field

def Restart_bias_cor(org_res,\
                     pres_file,\
                     ace_file,\
                     field_type):
    '''
        Stratopheric Bias Correction for GC restart files
        field_type = month or season
    '''
    ## Step1: open orginal restart file and get timestamp
#     org_res = xr.open_dataset(restart_file).load()
#     org_res.close()
    
#     year = org_res['time'].dt.year.values[0]
    
    ## Step2: get ACE-FTS climatological xch4 field 
    
    cor_fac = xr.open_dataset(ace_file).load()
    if field_type=='month':
        month = org_res['time'].dt.month.values[0]
        zonal_ch4 = cor_fac['CH4'].sel(month = month)
    elif field_type == 'season':
        sea = org_res['time'].dt.season.values[0]
        zonal_ch4 = cor_fac['CH4'].sel(season = sea)

    
    year = org_res['time'].dt.year.values[0]
#     month = np.int(org_res['time'].dt.month.values[0])
    
    zonal_fac = cor_fac['factors'].sel(season = sea, year= year)
    weight_fac =  cor_fac['weight_factor'].sel(season = sea)
#     ## Step3: calculating field factor for mass conservation
    
#     cur_zonal_fac = ace_field_factor(res_field =org_res['SpeciesRst_CH4'].values.squeeze(),\
#                                      zonal_ch4 = zonal_ch4.values,\
#                                      Troplev = org_res['Met_TropLev'].values.squeeze(),\
#                                      lat = org_res.lat.values,\
#                                      zonal_fac =zonal_fac,\
#                                      n_lev =3)
    
    
#     month_factors = np.array([0.14896956, 0.36515615, 0.71026034, 1.        , 0.99264037,
#        0.96737397, 0.81044566, 0.48053762, 0.35884427, 0.12333926,
#        0.00609595, 0.        ])
#     index = month-1
#     mon_factor = month_factors[index]
#     print(month, mon_factor)
    ## Step4: apply cor factors for restart fields high than Tropopause
    
    mod_field = cor_res_field(res_field = org_res['SpeciesRst_CH4'].values.squeeze(),\
                              zonal_ch4 = zonal_ch4.values,\
                              Troplev = org_res['Met_TropLev'].values.squeeze(),\
                              lat = org_res.lat.values,\
                              zonal_fac = zonal_fac,\
                              mon_factor =weight_fac)
    
    org_res['org_SpeciesRst_CH4'] = org_res['SpeciesRst_CH4'].copy()
    org_res['SpeciesRst_CH4'][dict(time = 0)] = mod_field
    
    return org_res

def cor_res_field1(res_field, zonal_ch4, Troplev, lat, zonal_fac, n_lev, mon_factor):
    '''
    restart field correction using ace-fts strtospheric zonal field
    '''
    nlev, nlat, nlon = np.shape(res_field)
    mod_field = res_field.copy()

    for ilat, llat in zip(range(nlat),lat):
        if np.all(np.isnan(zonal_ch4[:, ilat])):
            continue
        else:
            ace_ch4_pl = zonal_ch4[:, ilat]

        lat_bins = zonal_fac['lat_bins'].values
        fac_idx = find_closest(lat_bins, llat)
   
        if ~np.isnan(zonal_fac.isel(lat_bins = fac_idx).values):
            ace_ch4_pl = zonal_fac.isel(lat_bins = fac_idx).values* ace_ch4_pl
            
        fac = mon_factor.isel(lat_bins = fac_idx).values
        if np.isnan(fac):
            fac = 0.0
  #         if np.abs(llat) < 45:
#             fac = np.abs(llat)/200
#         else:
#             fac = np.abs(llat)/100
#         fac = np.abs(llat)/100
        mean_lev = np.int(np.nanmean(Troplev[ilat,:]))
        mean_lon = np.nanmean(res_field[mean_lev,ilat, :]-ace_ch4_pl[mean_lev])
        std_lon = np.nanstd(res_field[mean_lev,ilat, :]-ace_ch4_pl[mean_lev])
        for ilon in range(nlon):
            pres_lev = np.int(Troplev[ilat, ilon]-1) 
            fac1 = fac
            if abs(res_field[pres_lev,ilat, ilon]-ace_ch4_pl[pres_lev])>(mean_lon +3*std_lon):
                fac1 = fac1/10.0
            else:
                fac1 = fac1

            top = pres_lev+n_lev
            
            mod_field[pres_lev:top,ilat, ilon] = (1.0-fac1)*res_field[pres_lev:top,ilat, ilon]+ fac1*ace_ch4_pl[pres_lev:top]
         

    return mod_field

def ace_field_factor(res_field, zonal_ch4, Troplev, lat,  n_lev, zonal_fac, pres_levels):
    '''
    calculating field factors for mass conservation
    
    '''
    nlev, nlat, nlon = np.shape(res_field)
  
    res_stratos = np.full((n_lev, nlat, nlon), np.nan)
    ace_stratos = np.full((n_lev, nlat, nlon), np.nan)
    xch4_res = np.full(( nlat, nlon), np.nan)
    xch4_ace = np.full(( nlat, nlon), np.nan)
#     print(np.shape(zonal_ch4), np.nanmean(zonal_ch4))
    for ilat, llat in zip(range(nlat), lat):
        
        if np.all(np.isnan(zonal_ch4[:, ilat])):
            continue
        else:
            ace_ch4_pl = zonal_ch4[:, ilat]

        lat_bins = zonal_fac['lat_bins'].values
        fac_idx = find_closest(lat_bins, llat)
        if ~np.isnan(zonal_fac.isel(lat_bins = fac_idx).values):
            ace_ch4_pl = zonal_fac.isel(lat_bins = fac_idx).values* ace_ch4_pl
        
        for ilon in range(nlon):
            pres_wgt = setup_pres_weight(pres_levels[:, ilat, ilon])
            pres_lev = np.int(Troplev[ilat, ilon]) 
            
            res_stratos[:, ilat, ilon] = res_field[pres_lev:pres_lev+n_lev, ilat, ilon]
            ace_stratos[:, ilat, ilon] = ace_ch4_pl[pres_lev:pres_lev+n_lev]
            
            xch4_res[ilat, ilon] = np.nansum(pres_wgt[pres_lev:pres_lev+n_lev]*res_stratos[:, ilat, ilon])
            xch4_ace[ilat, ilon] = np.nansum(pres_wgt[pres_lev:pres_lev+n_lev]*ace_stratos[:, ilat, ilon])
            
    field_factor = np.nanmean(xch4_res)/np.nanmean(xch4_ace)
    
    zonal_fac1 = zonal_fac*field_factor 
    
    return zonal_fac1

def Restart_bias_cor1(org_res,\
                      pres_file,\
                      ace_file,\
                      field_type):
    '''
        Stratopheric Bias Correction for GC restart files
        field_type = month or season
    '''
    ## Step1: open orginal restart file and get timestamp
#     org_res = xr.open_dataset(restart_file).load()
#     org_res.close()
    
    pres = xr.open_dataset(pres_file)['Met_PEDGEDRY'].isel(time =0).load()
    pres.close()
    
#     year = org_res['time'].dt.year.values[0]
    
    ## Step2: get ACE-FTS climatological xch4 field 
    
    cor_fac = xr.open_dataset(ace_file).load()
    if field_type=='month':
        month = org_res['time'].dt.month.values[0]
        zonal_ch4 = cor_fac['CH4'].sel(month = month)
    elif field_type == 'season':
        sea = org_res['time'].dt.season.values[0]
        zonal_ch4 = cor_fac['CH4'].sel(season = sea)


    year = org_res['time'].dt.year.values[0]
#     month = np.int(org_res['time'].dt.month.values[0])
    
    zonal_fac = cor_fac['growth_factors'].sel(season = sea, year= year)
    weight_fac =  cor_fac['weight_factor'].sel(season = sea)

#     month_factors = np.array([0.14896956, 0.36515615, 0.71026034, 1.        , 0.99264037,\
#                               0.96737397, 0.81044566, 0.48053762, 0.35884427, 0.12333926,\
#                               0.00609595, 0.        ])
#     index = month-1
#     mon_factor = month_factors[index]
#     print(month, mon_factor)
    ## Step3: calculating field factor for mass conservation
    zonal_fac = ace_field_factor(res_field = org_res['SpeciesRst_CH4'].values.squeeze(),\
                                 zonal_ch4 = zonal_ch4.values,\
                                 Troplev = org_res['Met_TropLev'].values.squeeze(),\
                                 lat = org_res.lat.values,\
                                 n_lev =3,\
                                 zonal_fac = zonal_fac,\
                                 pres_levels = pres.values.squeeze())
    
    ## Step4: apply cor factors for restart fields high than Tropopause
    
    mod_field = cor_res_field1(res_field = org_res['SpeciesRst_CH4'].values.squeeze(),\
                               zonal_ch4 = zonal_ch4.values,\
                               Troplev = org_res['Met_TropLev'].values.squeeze(),\
                               lat = org_res.lat.values,\
                               zonal_fac = zonal_fac,\
                               n_lev =3,\
                               mon_factor = weight_fac)
    
     
    
    
    
    org_res['org_SpeciesRst_CH4'] = org_res['SpeciesRst_CH4'].copy()
    org_res['SpeciesRst_CH4'][dict(time = 0)] = mod_field 
    
    return org_res


if (__name__ == '__main__'):
    period = pd.date_range(start = '2010-01-01', end = '2019-12-01', freq = 'MS')
    filenm =  '/data/shzhu/org_ch4/global_2x25/Restart/Restart.YYYYMMDD_0000z.nc4'
    pres_file = '/data/shzhu/org_ch4/global_2x25/met/GEOSChem.LevelEdgeDiags.YYYYMMDD_0000z.nc4'
#     filenm =  './Restart/Restart.YYYYMMDD_0000z.nc4'
    new_file = '/home/shzhu/python/2021.07.18_new_stratos_cor_scheme/mod_res_1/Restart.YYYYMMDD_0000z.nc4'

    for date in period:
        print(date)
        filename = filenm.replace('YYYYMMDD', date.strftime('%Y%m%d'))
        presflnm = pres_file.replace('YYYYMMDD', date.strftime('%Y%m%d'))
        new_filename =new_file.replace('YYYYMMDD', date.strftime('%Y%m%d'))
        org_res = xr.open_dataset(filename).load()
        org_res.close()
        new_res = Restart_bias_cor(org_res = org_res,\
                                   pres_file = presflnm,\
                                   ace_file = '/data/shzhu/measurements/ACE-FTS/ace_ch4_sea_pres_new.nc',
                                   field_type = 'season')#zonal_daily_sea.nc
#         new_res = Restart_bias_cor1(org_res = org_res ,\
#                                     pres_file = presflnm,\
#                                     ace_file = '/data/shzhu/measurements/ACE-FTS/ace_ch4_sea_pres_new.nc',
#                                     field_type = 'season')#zonal_daily_sea.nc
    
#         new_res.to_netcdf(new_filename)
