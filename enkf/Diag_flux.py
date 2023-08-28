#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
from numpy import *
import time_module as tm
import xarray as xr
import geos_chem_def as gcdf
#import Diag_state as dgstate
from calendar import monthrange

def daily_posterior_flux(outflnm, doy, yyyy, emis_fn, mask_fn, inv_factor):
    ''' 
    Output the daily prior & posterior fluxes and the corresponding errors

    Args:
        outflnm: output file path + filename
        doy: time stamp
        eims_fn: prior flux nc file (path + filname)
        mask_fn: regional mask nc file !!! note: level dimension name(different regions) of mask must be i_reg
        inv_factor: Daily inverion state for posterior flux calculation

    Returns:
        None

    '''
    ## get variable information from config file 
    n_sous = gcdf.n_sous
    n_regs = gcdf.n_regs
    sous_name = gcdf.emis_var_name
    
    ## define the timestamp used in prior field selection
    
    time = list(date)
    time[8:10] = '01'
    timestamp =''.join(time) 
    
    ## read prior flux and mask file 
    emis_reg = xr.open_dataset(emis_fn).sel(time = timestamp)
    emis_err = emis_reg 
    mask = xr.open_dataset(mask_fn)
    ## adjustment to mask file ( delete in future use )
    
    mask = mask.rename({'lev':'i_reg'})
    mask = mask.transpose('lon', 'lat', 'i_reg')
    mask = mask.assign_coords({'lon': emis_reg['lon'].values,
                               'lat': emis_reg['lat'].values,
                               'i_reg': range(20)})
    
    ## enlarge inv_factor form 1-d[0:n_regs] to 2-d global or regional girds[nlon,nlat] by apply mask
    inv_factor = inv_factor.transpose('i_reg','i_sous')
    inv_grid =xr.Dataset()
    for var in inv_factor.data_vars:
        inv_grid[var] = (inv_factor[var]*mask['MASK']).sum(dim = 'i_reg') 
    inv_grid = inv_grid.transpose( 'i_sous','lat','lon')
    
    ## build dataset contain prior and posterior information
    post_info = xr.Dataset()
    for var in sous_name:
        flux0_nm = var + '_flux0'  # prior flux
        flux_nm = var + '_flux'    # posterior flux
        err0_nm = var + '_err0'    # prior error
        err_nm = var + '_err'      # posterior error
        
        post_info[flux0_nm] = emis_reg[var] + inv_grid['cur_mean_x0'].sel(i_sous = var) * emis_err[var]
        post_info[flux_nm] = emis_reg[var] +  inv_grid['cur_mean_x'].sel(i_sous = var) * emis_err[var]
        post_info[err0_nm] = sqrt(inv_grid['cur_err_x0'].sel(i_sous = var)) * emis_err[var]
        post_info[err_nm] = sqrt(inv_grid['cur_err_x'].sel(i_sous = var)) * emis_err[var]
    
    post_info.drop(['time','i_sous'])
    post_info.attrs['Description'] = 'Daily assimilation results'
    post_info.attrs['Vars_prefix'] = 'Emission category'
    post_info.attrs['_0'] = 'Prior information'
    post_info.attrs['_'] = 'Posterior information'
    post_info.attrs['Reference time'] = date

    post_info.to_netcdf(outflnm)
    
    return None
    
def daily_posterior_conc(outflnm,\
                         time,\
                         lons,\
                         lats,\
                         mod_conc0,\
                         mod_conc,\
                         obs_conc,\
                         inv_conc,\
                         obs_err0,\
                         obs_err,\
                         mod_err,\
                         inv_err,\
                         project
                        ):
    ''' 
    Output the daily prior & posterior fluxes and the corresponding errors

    Args:
        outflnm   : output file path + filename
        doy       : time stamp
        time      : inversion date
        lons, lats: positions of the observations
        mod_conc0 : original model result(hm0)
        mod_conc  : daily prior model result(all_mean_y)
        obs_conc  : observations(xch4)
        inv_conc  : daily inversion concentration
        obs_err0  : original observation errors
        obs_err   : adjusted observation errors
        mod_err   : model errors
        inv_err   : invesion errors
        
    Notes:
        all errors are error **2

    Returns:
        None

    '''
    
    nobs = len(lons)
    
    inv_conc = xr.Dataset({'mod_conc0': (['nobs'], mod_conc0),
                           'mod_conc': (['nobs'], mod_conc),
                           'obs_conc': (['nobs'], obs_conc),
                           'inv_conc': (['nobs'], inv_conc),
                           'obs_err': (['nobs'], sqrt(obs_err) ),
                           'obs_err0': (['nobs'], sqrt(obs_err0)),
                           'mod_err': (['nobs'], sqrt(mod_err) ),
                           'inv_err': (['nobs'], sqrt(inv_err) ),
                           'obs_lon': (['nobs'], lons),
                           'obs_lat': (['nobs'], lats),
                           'project': (['nobs'], project)
                          },
                            coords = {'nobs': arange(nobs)})
    
    inv_conc.attrs['Description'] = 'Daily inversion concentration (ppbv)'
    inv_conc.attrs['obs_lon'] = 'Longitudes of the observations'
    inv_conc.attrs['obs_lat'] = 'Latitudes of the observations'
    inv_conc.attrs['mod_conc0'] = 'Original model result(hm0)'
    inv_conc.attrs['mod_conc'] = 'Daily prior model result(all_mean_y)'
    inv_conc.attrs['obs_conc'] = 'Observations(xch4)'
    inv_conc.attrs['inv_conc'] = 'Daily inversion concentration'
    inv_conc.attrs['obs_err0'] = 'Original observation errors'
    inv_conc.attrs['obs_err'] = 'Adjusted observation errors'
    inv_conc.attrs['mod_err'] = 'Model errors(dyr)'
    inv_conc.attrs['inv_err'] = 'Daily invesion errors'
    inv_conc.attrs['Reference time'] =  time
    
    
    inv_conc.to_netcdf(outflnm)
    return None

    
def whole_posterior_flux(outflnm, periods, sst, end, nstep , emis_fn, mask_fn, inv_factors):
    '''
    Output prior & posterior fluxes and the corresponding errors during the whole period
    
    Args:
        outflnm: output file path + filename
        periods: time series of the whole periods
        
        nstep: the number the invesion steps
        emis_fn: prior emissions file(path + filname)
        mask_fn: mask file(path + filename)
        inv_factors: inverion state for posterior flux calculation (state increments & error factors) 
        
    Returns:
        None
    '''
    if len(periods) != nstep:
        print('The periods are not consistent with steps!' )
        return None
    
    ## get variable information from config file 
    n_sous = gcdf.n_sous
    n_regs = gcdf.n_regs
    sous_name = gcdf.emis_var_name
    
    ## define the time series used in prior field selection
  
    ## read prior flux and mask file 
    emis_reg = xr.open_dataset(emis_fn).sel(time = slice(sst, end))
    emis_err = emis_reg 
    mask = xr.open_dataset(mask_fn)
    ## adjustment to mask file ( delete in future use )
    mask = mask.rename({'lev':'i_reg'})
    mask = mask.transpose('lon','lat','i_reg')
    mask = mask.assign_coords({'lon': emis_reg['lon'].values,
                               'lat': emis_reg['lat'].values,
                               'i_reg': range(n_regs)})
    
    ## enlarge inv_factor form 1-d[0:n_regs] to 2-d global or regional girds[nlon,nlat] by apply mask
    inv_factors = inv_factors.transpose('time','i_reg','i_sous')
    inv_grid =xr.Dataset()
    for var in inv_factors.data_vars:
        temp = inv_factors[var].isel(i_reg = 0)*mask['MASK'].isel(i_reg = 0)
        for i_reg in range(1, n_regs):
            median= inv_factors[var].isel(i_reg = i_reg)*mask['MASK'].isel(i_reg = i_reg)
            temp = temp + median
        inv_grid[var] = temp 
        
    inv_grid = inv_grid.transpose( 'i_sous','time','lat','lon')
    ## build dataset contain prior and posterior information
    post_info = xr.Dataset()
    for var in sous_name:
        flux0_nm = var + '_flux0'  # prior flux
        flux_nm = var + '_flux'    # posterior flux
        err0_nm = var + '_err0'    # prior error
        err_nm = var + '_err'      # posterior error
        
        post_info[flux0_nm] = emis_reg[var] + inv_grid['mean_x0'].sel(i_sous = var).values * emis_err[var]
        post_info[flux_nm] = emis_reg[var] + gcdf.pri_err* inv_grid['mean_x'].sel(i_sous = var).values * emis_err[var]
        post_info[err0_nm] = sqrt(inv_grid['err_x0'].sel(i_sous = var).values) * emis_err[var]
        post_info[err_nm] = sqrt(inv_grid['err_x'].sel(i_sous = var).values) * emis_err[var]
        
        ## added 18/01/2021 by rainbow for prior flux check 
#         meanx0_nm = var + '_meanx0' # grid meanx0
#         meanx_nm = var + '_meanx'   # grid meanx
        
#         post_info[meanx0_nm] =  inv_grid['mean_x0'].sel(i_sous = var)
#         post_info[meanx_nm] =  inv_grid['mean_x'].sel(i_sous = var)
        ## end ##
    
    post_info.attrs['Description'] = 'Assimilation results during the whole inversion period'
    post_info.attrs['Vars_prefix'] = 'Emission category'
    post_info.attrs['_0'] = 'Prior information'
    post_info.attrs['_'] = 'Posterior information'
    
    post_info.to_netcdf(outflnm)
    
    return None


if (__name__=='__main__'):
    
    idoy = 60
    yyyy = 2015
    
    date = tm.doy_to_utc(idoy, 0, yyyy)
    n_sous = 2
    n_regs = 20
    sous_name =['CH4_AN' , 'CH4_NA']
    cur_mean_x0 = random.rand(n_sous,n_regs)
    cur_mean_x = random.rand(n_sous,n_regs)
    cur_err_x0 = random.rand(n_sous,n_regs)
    cur_err_x = random.rand(n_sous,n_regs)
    daily_assim_factor = xr.Dataset({'cur_mean_x0':(['i_sous','i_reg'], cur_mean_x0),
                                   'cur_mean_x':(['i_sous','i_reg'], cur_mean_x),
                                   'cur_err_x0':(['i_sous','i_reg'], cur_err_x0),
                                   'cur_err_x':(['i_sous','i_reg'], cur_err_x)},
                                  coords = {'i_sous': sous_name,
                                            'i_reg': range(n_regs)})
    daily_assim_factor.attrs['Description'] = 'Assimilated inversion state of current day'
    daily_assim_factor.attrs['cur_mean_x0'] = 'Prior increment parameters'
    daily_assim_factor.attrs['cur_mean_x'] = 'Posterior increment parameters'
    daily_assim_factor.attrs['cur_err_x0'] = 'Porior error factors'
    daily_assim_factor.attrs['cur_err_x'] = 'Posterior error factors' 
    daily_assim_factor.attrs['Reference time'] = date
    
    outflnm = gcdf.inv_daily_state_outpath + "/" + "Daily_assim_state" + "." + date[:10] + ".nc"
    daily_posterior_flux(outflnm = outflnm,\
                        date = date,\
                        emis_fn = gcdf.emis_file,\
                        mask_fn = gcdf.mask_file,\
                        inv_factor = daily_assim_factor)
    
    
    