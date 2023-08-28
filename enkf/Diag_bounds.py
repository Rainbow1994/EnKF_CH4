#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
from numpy import *
import time_module as tm
import xarray as xr
import geos_chem_def as gcdf 
from rerun_geos import get_ensemble_list

def daily_sgl_bounds_update(doy, yyyy, bounds_fn, en_bounds_fn, new_bounds_fn, cur_meanx):
    '''
    update prior single tracer boundaries for nested assimilation 
    
    Args:
        doy, yyyy: curent date
        bounds_path: original boundary files path
        en_bound_fn: ensemble boundary files
        new_bounds_path: new boundaries files path
        cur_meanx: boundary increment factors
        
    Returns:
        None
    '''
    ##  Date and filenames information ###### 
    tst=tm.doy_to_utc(doy, 0, yyyy)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '')
    
    em_step_list,em_name_list = get_ensemble_list(yyyy, doy, gcdf.entry_table)
    em_bounds_names = [''.join([i, '.BoundaryConditions.', tst[0:8] ,'_0000z.nc4']) for i in em_name_list]
    var_list= [''.join(['SpeciesBC_CH4_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    
    ## Calculation of increments ############
    sgl_file =  ''.join([bounds_fn, '/', 'ST000.EN0001-EN0002.BoundaryConditions.', tst[0:8] ,'_0000z.nc4'])
    sgl_bounds = xr.open_dataset(sgl_file)
    
    add_data = xr.zeros_like(sgl_bounds['SpeciesBC_CH4'])
    ic = 0
    for file in em_bounds_names:
        en_bounds = xr.open_dataset(en_bounds_fn + file)
        for ivar in var_list:
            add_data.values = add_data.values + en_bounds[ivar].values * cur_meanx[ic]
            ic = ic+1
            
    ## Renew Boundaries #####################
    new_bounds = xr.Dataset()
    new_bounds = sgl_bounds + add_data
    new_bounds['SpeciesBC_CH4']
    for var in new_bounds.data_vars:
        new_bounds[var].attrs = sgl_bounds[var].attrs
    new_bounds.attrs = sgl_bounds.attrs 
    
    new_file = ''.join([new_bounds_fn, '/', 'ST000.EN0001-EN0002.BoundaryConditions.', tst[0:8] ,'_0000z.nc4'])

    new_bounds.to_netcdf(new_file)
    
    return None 



def daily_en_bounds_update(doy, yyyy, en_bounds_fn, new_bounds_fn, cur_xtm):
    '''
    update prior ensemble boundaries for nested assimilation 
    
    Args:
        doy, yyyy: time stamp
        en_bounds_fn: original ensemble boundarie files path
        new_bounds_fn: new ensemble boundaries files path
        cur_xtm: ensemble boundaries factors
        
    Returns:
        None
    '''
    ##  Date and filenames information ###### 
    tst=tm.doy_to_utc(doy, 0, yyyy)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '')
    
    em_step_list,em_name_list = get_ensemble_list(yyyy, doy, gcdf.entry_table)
    em_bounds_names = [''.join([i, '.BoundaryConditions.', tst[0:8] ,'_0000z.nc4']) for i in em_name_list]
    var_list= [''.join(['SpeciesBC_CH4_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    
    ## update the ensemble boudaries
    ic = 0
    for file in em_bounds_names:
        en_bounds = xr.open_dataset(en_bounds_fn + file)
        new_bounds = xr.zeros_like(en_bounds)
        for ivar in var_list:
            new_bounds[ivar].values = en_bounds[ivar].values * cur_xtm[ic]
            ic = ic + 1
        new_bounds.to_netcdf(new_bounds_fn + file)
        
    return None
        
    
def step_sgl_bounds_update(doy_list, yyyy_list, bounds_path, en_bounds_path, new_bounds_path, step_meanx):
    '''
    update prior single tracer boundary for nested assimilation 
    
    Args:
        doy_list: [doy_sst, doy_end]
        yyyy_list: [yyyy_sst, yyyy_end]
        bounds_path: original single boundary files path
        en_bounds_path: original ensemble boundary files path
        new_bounds_path: new single boundary file path
        step_meanx: boundary increment factors of the curent period
        
    Returns:
        None
    '''
    
    ###determine the update time periods ##########
    sst = doy_list[0]
    yyyy = yyyy_list[0]
    
    if (yyyy_list[1] != yyyy_list[0]):
        days_of_year = 366 if calendar.isleap(yyyy) else 365
        end = doy_list[0] + doy_list[1]
    else:
        end = doy_list[1]
        
    for doy in range(sst, end):
        
        ###daily update loop ##########
        daily_sgl_bounds_update(doy = doy,\
                                yyyy = yyyy,\
                                bounds_fn = bounds_path,\
                                en_bounds_fn = en_bounds_path,\
                                new_bounds_fn = new_bounds_path,\
                                cur_meanx = step_meanx )
        
    return None
    
def step_en_bounds_update(doy_list, yyyy_list, en_bounds_path, new_bounds_path, step_xtm):
    
    '''
    update prior ensemble boundaries for nested assimilation 
    
    Args:
        doy_list: [doy_sst, doy_end]
        yyyy_list: [yyyy_sst, yyyy_end]
        en_bounds_path: original ensemble boundary files path  
        new_bounds_path: new ensemble boundary files path
        step_xtm: ensemble boundaries factors of the curent period
        
    Returns:
        None
    '''

    ###determine the update time periods ##########
    sst = doy_list[0]
    yyyy = yyyy_list[0]
    
    if (yyyy_list[1] != yyyy_list[0]):
        days_of_year = 366 if calendar.isleap(yyyy) else 365
        end = doy_list[0] + doy_list[1]
    else:
        end = doy_list[1]
    
    for doy in range(sst, end):
        
        ###daily update loop ##########
        daily_en_bounds_update(doy = doy,\
                               yyyy = yyyy,\
                               en_bounds_fn = en_bounds_path,\
                               new_bounds_fn = new_bounds_path,\
                               cur_xtm = step_xtm)
        
    return None

if (__name__=='__main__'):
    
    nx = 40
    cur_meanx = random.rand(nx)
    idoy = 1
    yyyy = 2015
    bounds_fn = gcdf.sgl_bounds_path
    en_bounds_fn = gcdf.en_bounds_path
    new_bounds_fn = gcdf.new_sgl_bounds_path
    daily_sgl_bounds_update(idoy, yyyy, bounds_fn, en_bounds_fn, new_bounds_fn, cur_meanx)
    
    ne = nx
    cur_xtm = random.rand(ne)
    new_bounds_fn = gcdf.new_en_bounds_path
    daily_en_bounds_update(idoy, yyyy, en_bounds_fn, new_bounds_fn, cur_xtm)
    