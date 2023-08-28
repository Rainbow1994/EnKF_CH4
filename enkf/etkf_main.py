#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
from numpy import *
import time_module as tm
import xarray as xr
import geos_chem_def as gcdf 
import numpy.linalg as alg
import rerun_geos as rg
import sample_ch4_gosat as scg
import data_quality_ctl as dqc
import construct_obs_err as ocor
import prior_err_cov as perr
import do_assim as da
import Diag_state as dgstate
import Diag_flux as dgflux
import Diag_bounds as dgbounds
from calendar import monthrange
import assim_step as assim
import os
import sys
#import assim_step_tccon as assim

if sys.platform.startswith('linux'):
    os.system('clear')
    
f = open('geos_chem_def.py', 'r')
content = f.read()
print(content)
f.close()

assim = assim.assim_step(use_time_cor = gcdf.use_time_cor,\
                         use_fixed_mod_err = gcdf.use_fixed_mod_err,\
                         use_full_r = gcdf.use_full_r,\
                         use_sparse = gcdf.use_sparse,\
                         add_rnd_err = gcdf.add_rnd_err,\
                         do_update = gcdf.do_update,\
                         full_rerun = gcdf.full_rerun,\
                         sampl_input = gcdf.sampl_input,\
                         daily_prior_output = gcdf.daily_prior_output,\
                         daily_para_output = gcdf.daily_para_output,\
                         cumulative_para_output = gcdf.cumulative_para_output,\
                         step_para_output = gcdf.step_para_output,\
                         daily_posterior_output = gcdf.daily_posterior_output,\
                         daily_updata_en_bounds = gcdf.daily_updata_en_bounds,\
                         daily_updata_sgl_bounds = gcdf.daily_updata_sgl_bounds,\
                         daily_conc_output = gcdf.daily_conc_output)


posterior_diag = gcdf.posterior_flux_diag
update_en_bounds = gcdf.update_en_bounds
update_sgl_bounds = gcdf.update_sgl_bounds

step_st = gcdf.step_start 
step_end = gcdf.step_end
nstep = step_end - step_st
print(nstep)

## define the time series used in prior field selection
time_sst = tm.doy_to_utc(gcdf.st_doy, 0, gcdf.st_yyyy)
tem_res = gcdf.temp_res
periods = [time_sst]
if tem_res == 'mon':
    mon = int(time_sst[5:7]) 
    yyyy = gcdf.st_yyyy
    doys = 0
    #mon = step_st + 1
    for i in range(nstep):
        if mon == 13:
            yyyy =yyyy +1
            mon = 1
            doys += monthrange(yyyy, mon)[1]
        else:
            doys += monthrange(yyyy, mon)[1]
        mon += 1
        periods.append(tm.doy_to_utc(gcdf.st_doy + doys, 0, gcdf.st_yyyy))
else:
    for i in range(nstep):
        periods.append(tm.doy_to_utc(gcdf.st_doy + tem_res, 0, gcdf.st_yyyy))
    

periods = periods[:-1]
print('periods:', len(periods), periods)


for i in range(step_st, step_end):
    
    print('Step:', i)
    
    ### step inversion loop
    doy_sst = assim.rerun_doy_st
    yyyy_sst = assim.rerun_yyyy_st
    
    assim.do_one_step(update_step = gcdf.update_step,\
                      rerun_step = gcdf.rerun_step,\
                      sample_step = i)
    doy_end = assim.rerun_doy_st
    yyyy_end = assim.rerun_yyyy_st

    doy_list = [doy_sst, doy_end]
    yyyy_list = [yyyy_sst, yyyy_end]
    
    ### update ensemble boundaries using step inversion parameters
    
    if (update_en_bounds):
        
        tdx = dot(assim.dx, assim.xtm)
        step_err_x = diag(dot(tdx, transpose(tdx)))
        
        dgbounds.step_en_bounds_update(doy_list = doy_list,\
                                       yyyy_list = yyyy_list,\
                                       en_bounds_path = gcdf.en_bounds_path,\
                                       new_bounds_path = gcdf.new_en_bounds_path,\
                                       step_xtm = step_err_x )
        
    ### update single tracer boundaries using step inversion parameters
    
    if (update_sgl_bounds):
       
        dgbounds.step_en_bounds_update(doy_list = doy_list,\
                                       yyyy_list = yyyy_list,\
                                       bounds_path = gcdf.sgl_bounds_path,\
                                       en_bounds_path = gcdf.en_bounds_path,\
                                       bounds_new_path = gcdf.new_sgl_bounds_path,\
                                       step_meanx = assim.mean_x)
    
if (posterior_diag):
    ### output posterior information during the whole inversion period
    
    ### the inversion parameters
    outflnm = gcdf.whole_state_path + '/' + 'whole_inv_state.' + str(step_st) + '_' + str(step_end) + '.nc' 
    
    inv_factors = dgstate.output_whole_states(outflnm = outflnm,\
                                             periods = periods,\
                                             nstep = nstep,\
                                             mean_x0 = assim.whole_mean_x0,\
                                             err_x0 = assim.whole_err_x0,\
                                             mean_x = assim.whole_mean_x,\
                                             err_x = assim.whole_err_x,\
                                             diag_whole_state = True)
    print(inv_factors)
    
    ### output prior & posterior fluxes and their error during the inverion period
    outflnm = gcdf.whole_inv_flux_path + '/' + 'post_flux.' + str(step_st) + '_' + str(step_end) + '.nc' 

    dgflux.whole_posterior_flux(outflnm = outflnm,\
                                periods = periods,\
                                sst = periods[0],\
                                end = periods[-1],\
                                nstep = nstep,\
                                emis_fn = gcdf.emis_file,\
                                mask_fn = gcdf.mask_file,\
                                inv_factors = inv_factors)
    
    