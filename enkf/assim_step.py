#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
import sys
sys.path.append('/home/shzhu/enkf_ch4/rundirs/grid_2x25_new/')
import numpy as np
from numpy import *
import time_module as tm
import xarray as xr
import geos_chem_def as gcdf 
import numpy.linalg as nlg
import rerun_geos as rg
import sample_ch4_gosat as scg
import obs.obs_operator as obs
import obs.operator_config as obs_config
import data_quality_ctl as dqc
import construct_obs_err as ocor
import prior_err_cov as perr
import do_assim as da
import Diag_state as dgstate
import Diag_flux as dgflux
import Diag_bounds as dgbounds
import calendar
import os 
from datetime import datetime
from global_land_mask import globe

class assim_step:
    def __init__(self, \
                 use_time_cor = True,\
                 use_fixed_mod_err = True,\
                 use_full_r = True,\
                 use_sparse = True,\
                 add_rnd_err = False,\
                 do_update = True,\
                 full_rerun = True,\
                 sampl_input = False,\
                 daily_prior_output = False,\
                 daily_para_output = False,\
                 cumulative_para_output = False,\
                 step_para_output = True,\
                 daily_posterior_output = False,\
                 daily_updata_en_bounds = False,\
                 daily_updata_sgl_bounds = False,\
                 daily_conc_output = False):
        
        """ self initialization  """
        
        # module switches
        self.use_time_cor = use_time_cor
        self.use_fixed_mod_err = use_fixed_mod_err
        self.use_full_r = use_full_r
        self.use_sparse = use_sparse
        self.add_rnd_err = add_rnd_err
        self.do_update = do_update
        self.full_rerun = full_rerun
        self.sampl_input = sampl_input
        # switches for paramters diagnosis
        self.daily_state_output = daily_prior_output
        self.daily_para_output = daily_para_output
        self.cumul_para_output = cumulative_para_output
        self.step_para_output = step_para_output
        # switches for daily posterior concentration & flux output
        self.daily_flux_output = daily_posterior_output
        self.daily_conc_output = daily_conc_output
        # switches for daily boundaries update
        self.daily_en_bounds = daily_updata_en_bounds
        self.daily_sgl_bounds = daily_updata_sgl_bounds
        
        
        # statistical coefficients
        self.xnorm = gcdf.xnorm
        self.xref = gcdf.xref
        self.mod_bias = array(gcdf.mod_bias)
        self.bias_err = array(gcdf.bias_err)
        
        self.xtm_enlarge_factor = gcdf.xtm_enlarge_factor
        self.tcor_factor = gcdf.tcor_factor
        self.tcor_len = gcdf.tcor_len
        
        # state matrix size
        self.nx, self.ne = 0,0
        self.ntracers = gcdf.ntracers
        self.nbias = gcdf.nbias
        
        # inversion parameters
        self.lag_window = gcdf.inv_lag_window 
        self.inv_time_step = gcdf.temp_res
        self.first_rerun = False    
        
        # initialize step information
        
        self.maxstep = gcdf.inv_lag_window + gcdf.step_start # for state update
        self.cur_step = gcdf.step_start
        self.cur_yyyy = None
        self.cur_doy = None
        
        # ETKF parameters (daily loop)
        self.mean_x = None              # increment to xf with temporal correlation 
        self.xtm = None                 # rotations current transform matrix(iterative)
        
        
        # ETKF parameters (step loop)
        self.dx = None                  # prior error covariance of each step  
        self.mean_x0 = None             # prior parameters(for future use)     
        self.cor_mean_x = None          # temporal covariance at current step 
        self.err = None                 # error matrix before current step inversion
        
        self.step_dx= list()           # store the dx of each step 
        self.x_inc = list()            # store the increment of each step
        self.mean_x_step= list()       # store the cumulative increments after current step
        
        
        # step iterative state coefficients
        self.whole_mean_x0 = None
        self.whole_mean_x = None
        self.whole_err_x0 = None
        self.whole_err_x = None

        
        # rerun settings 
        
        self.rerun_yyyy_st = gcdf.st_yyyy
        self.rerun_doy_st = gcdf.st_doy
        self.rerun_yyyy_end = gcdf.st_yyyy 
        self.rerun_doy_end = gcdf.st_doy
        
    def do_one_step(self, update_step = 0, rerun_step = 6, sample_step = 0):

        """ do one assimilation and rerun model to the next step  """

        # determine the current step time period
        self.cur_doy = self.rerun_doy_st
        self.cur_yyyy = self.rerun_yyyy_st 
        cur_st_doy  = self.rerun_doy_st
        
        yyyy, mm,dd = tm.doy_to_time_array(self.cur_doy, self.cur_yyyy)
        
        obs_operator = obs.obs_operator(obs_dict = obs_config.measurements)
        
        if (self.inv_time_step == 'mon'):
            tm_res = calendar.monthrange(yyyy,mm)[1]
        else:
            tm_res = self.inv_time_step
            
        cur_end_doy = cur_st_doy + tm_res
        self.rerun_doy_end = cur_end_doy
        
        ## -----------------------------------------------------------
        ## step1: rerun geoschem to get model result at current step
        ## -----------------------------------------------------------
        ## Args:
        ##      doy_st, doy_end: rerun period
        ##      cur_yyyy: year for rerun
        ##      cur_step: current step
        ##      ** : rerun and configuration directories are set in geos_chem_def.py
        ## Return:
        ##      None
        ## Output:
        ##      concentration fields of single tracer
        ##      coresponding boundaries
        ## -----------------------------------------------------------

        
        if (self.cur_step >= rerun_step):
            self.first_rerun = True if (self.cur_step == 0) else False
            
            rg.rerun_geos_chem(doy_st = cur_st_doy,\
                               doy_end = cur_end_doy,\
                               cur_yyyy = self.cur_yyyy,\
                               cur_step = self.cur_step,\
                               fst_rerun = self.first_rerun)
            
            
        ## -----------------------------------------------------------
        ## step2: prepare state parameters for current step 
        ## -----------------------------------------------------------
  
        all_cor_mean_x = self.prepare_step_state()
        self.prepare_whole_state()
        
        print('#########################')
        print('CURRENT STEP:', self.cur_step)
        print('BEFORE STEP INVERSTION!')
        print('=========================')
        print('mean_x', self.mean_x)
        print('=========================')
        print('whole_mean_x', self.whole_mean_x)
        print('#########################')
        ## -----------------------------------------------------------
        ## step3: initialize or update inversion parameters
        ## -----------------------------------------------------------
        # Daily inv paramters: xinc, cur_inc_m, cur_xtm
        xinc = zeros(shape(self.mean_x), float)
        cur_inc_m = xinc
        cur_xtm = zeros(shape(self.xtm), float)
        
        # Cumulative inv parameters: all_xinc, all_cur_xtm, all_cur_inc_m
        all_xinc = zeros(shape(self.mean_x), float)
        all_cur_inc_m = all_xinc
        all_cur_xtm = zeros(shape(self.xtm), float)
        
        # Step inv parameters: step_xinc, step_xtm, step_inc_m
        step_xinc = zeros(shape(self.mean_x), float)
        step_inc_m = step_xinc
        step_xtm = zeros(shape(self.xtm), float)
        
        print('Prior state preparation finished!')

        if (self.cur_step >= update_step):
            
            cur_step = self.cur_step
            usd_doy = 0
            fst_update = True

            # enter each day loop 
            for idoy in range(cur_st_doy, cur_end_doy):
                
                
                yyyy, mm,dd = tm.doy_to_time_array(idoy, self.cur_yyyy)
                print('Current date is: %4d-%02d-%02d'%(yyyy,mm,dd))
                
                ## -----------------------------------------------------------
                ## step4: resample simple and ensemble model result 
                ## -----------------------------------------------------------
                ## Args:
                ##      idoy: curent days of year
                ##      yyyy: curent year
                ##      cur_step: current step
                ##      ** : model output and observation directories are 
                ## Return:
                ##      sampling result
                ##            hm0: model result at observation sites
                ##            hm: ensemble model results
                ##            obs: observations
                ##            err: observation errors
                ## 
                ## -----------------------------------------------------------
                date0 = '%4d%02d%02d'%(yyyy,mm,dd)
                date = datetime.strptime(date0, '%Y%m%d')
                
                
                if gcdf.new_hm0:
                    sample_rlt = obs_operator.hm0_concat(date, cur_step) #### new hm0
                else:
                    sample_rlt = obs_operator.hm_concat(date) #### sampling exited hm of dif measurements
#                 print(date0)
#                 sample_rlt.to_netcdf('./test_'+date0+'.nc')   
                if not sample_rlt:
                    usd_doy = usd_doy + 1 
                    
                    continue
             
                
                print('Daily resampling finished')
                ## -----------------------------------------------------------
                ## step5: data quality control to get matched yobs, error, hm0, hm
                ## -----------------------------------------------------------


                data_ctl_rlt = dqc.data_quality_ctl(sample_rlt)
                
                if not data_ctl_rlt: ## no data
                    usd_doy = usd_doy + 1 
                    
                    continue

                all_mean_y0 = data_ctl_rlt['hm0'].values          # model results(prior information) 
                all_yobs = data_ctl_rlt['obs'].values             # observation results
                all_h = data_ctl_rlt['hm'].values                 # H: delta Yf_0/ delta Xf_0
                all_oerr = data_ctl_rlt['oerr'].values            # observation error
                all_lon = data_ctl_rlt['lon'].values              # longitudes of observation sites
                all_lat = data_ctl_rlt['lat'].values              # latitudes of observation sites
                
                print('Data filtering finished')

                ## -----------------------------------------------------------
                ## step6: adjust simulation result(hm0) by adding temporal correlation
                ## -----------------------------------------------------------

                cor_adj = dot(all_h, all_cor_mean_x)       # add temperal covariance
                nobs = size(all_lon)
                
                if (usd_doy > 0):
                    hm0_adj = dot(all_h, all_xinc)           # apply updated assim paras 
                else:
                    hm0_adj = zeros(nobs, float)


                if (self.full_rerun):
                    history_adj = hm0_adj + cor_adj

                else:
                    # the correlation adj has been the included in self.mean_x
                    history_xinc = self.mean_x
                    history_adj = dot(all_h, history_xinc)
                    
                print('maxium/min all_h:',  max(all_h.flat), min(all_h.flat))
                print('maxium/min all_cor_mean_x:',  max(all_cor_mean_x), min(all_cor_mean_x))
                print('maxium/min all_xinc:',  max(all_xinc), min(all_xinc))
                print('maxium/min history_adj:',  max(history_adj), min(history_adj))

                all_mean_y = all_mean_y0 + history_adj           # adjusted model results
                
                print('maxium/min all_mean_y:',  max(all_mean_y), min(all_mean_y))
                
                print('Model result adjustment finished!')

                ## -----------------------------------------------------------
                ## step7: determine the finial error R of observarion
                ## -----------------------------------------------------------

                dy  = dot(all_h, self.dx)                   # H(delta Xf_{k-1}) = delta Yf_(k-1)
                xdy = dot(dy, self.xtm)                     # delta Yf_k
                xdy = dot(xdy, transpose(xdy))              # delta Yf_k * (delta Yf_k)T
                dyr = diag(xdy)                             # prior concention error **2
                # ydev: deviation bewteen observations and model results
                ydev = (all_mean_y - all_yobs) **2
                
                obs_err = ocor.get_obs_err(lon = all_lon,\
                                           lat = all_lat,\
                                           oerr = all_oerr,\
                                           ydev = ydev,\
                                           model_r = dyr,\
                                           use_fixed_mod_err = self.use_fixed_mod_err,\
                                           use_full_r = self.use_full_r)


                print('Observation error construction finished!')

                ## -----------------------------------------------------------
                ## step8: do assimilations to get inv parameters at current day
                ## -----------------------------------------------------------

                # add random error to observation results (defalut is false)
                if (self.add_rnd_err):
                    sqrt_err = sqrt(all_yobs_err)
                    rnd_err = rnd.normal(sqrt_err)
                    all_yobs = all_yobs + 0.4 * rnd_err

                # prepare prior parameters befor entering ETKF 
                dx = dot(self.dx, self.xtm)                  # delta Xf_k
                dy = dot(dy, self.xtm)                       # delta Yf_k
                all_mean_x = self.mean_x
                print('Congratulations! all variables are avalible~~~')
                
                if (self.daily_state_output):
                    sstep = r'%2.2d' % (self.cur_step)
                    #time = tm.doy_to_utc(idoy, 0, self.cur_yyyy)
                    output_file = gcdf.daily_state_path + "/" + "daily_prior" + "." + time[:10] + ".nc"
                    dgstate.daily_prior_state(outflnm = output_file,\
                                              time = time,\
                                              lons = all_lon,\
                                              lats = all_lat,\
                                              delta_x = dx,\
                                              hm = all_h,\
                                              mean_x = all_mean_x,\
                                              mod_rlt = all_mean_y,\
                                              obs = all_yobs,\
                                              obs_err = diag(obs_err))
                    
                ##############################################################
                ##### Methods:
                ## xf: prior flux 
                ## delta Xf: an ensemble of pertubation states-- prior error
                ## xa: posterior flux
                ## delta Xa: posterior error
                ## --------ETKF-----------------------------------------------
                ## xa = xf + delta Xf * wa
                ## delta Xa = delta Xf * omega
                ##
                ## Purpose: caculating 'wa' and 'omega'
                ## wa = (delta Yf)T *(delta Yf * (delta Yf)T +R)^(-1)(y_obs - H(xf))
                ## omega = (I + (delta Yf)T * R^(-1) * delta Yf)^(-1/2)
                ##------------------------------------------------------------
                ##### Inputs:
                ## dx: prior flux pertubation--delta X(delta Xf)
                ## dy: prior concentration pertubation --delta Y(delta Yf)
                ## mean_x: prior estimate parameters (wa)
                ## all_mean_y: model results from CTM (H(xf))
                ## all_yobs: observation results(y_obs)
                ## obs_err: observation error(R)
                ##------------------------------------------------------------
                ##### Output:
                ## xinc: increment to prior estimates Xf (wa * delta Xf) 
                ## cur_xtm: current factors apply for delta Xf(omega)
                ## cur_inc_m: parameters for pertubations (wa)
                ##############################################################
                n_cur_fac = self.ntracers
                #print('dx:', dx[-n_cur_fac:, -n_cur_fac:])
                #print('dy:', dy[-n_cur_fac:, -n_cur_fac:])
                print('max(dy):', amax(dy[-n_cur_fac:, -n_cur_fac:]))
                print('min(dy):', amin(dy[-n_cur_fac:, -n_cur_fac:]))
                print('max(mean_x), min(mean_x):', amax(all_mean_x), amin(all_mean_x))
                print('max(mean_y), min(mean_y):', amax(all_mean_y), amin(all_mean_y))
                #print('observations:', all_yobs)
                print('max(obs):', max(all_yobs))
                print('min(obs):', min(all_yobs))
                #print('obs_errs:', obs_err)
                print('max(obs_err):', amax(obs_err))
                print('min(obs_err):', amin(obs_err))
                
                xinc, cur_xtm, cur_inc_m = da.do_assim(tdx = dx,\
                                                       tdy = dy,\
                                                       mean_x = all_mean_x,\
                                                       mean_y = all_mean_y,\
                                                       yobs = all_yobs,\
                                                       yerr = obs_err,\
                                                       xref = self.xref, \
                                                       xnorm = self.xnorm, \
                                                       use_sparse = self.use_sparse)
                
                print('Inversion calculation finished!')
                
                # update inv parameters daily after assimilation 
                print('Step prior error:',  diag(self.err))
                print('---max, ---min:',  amax(diag(self.err)), amin(diag(self.err)))
                tdx=dot(self.dx, self.xtm)
                new_err=dot(tdx, transpose(tdx))
                print('Daily prior error:', diag(new_err))
                print('---max, ---min:',  amax(diag(new_err)), amin(diag(new_err)))
                
                self.mean_x = self.mean_x + xinc
                self.xtm = dot(self.xtm, cur_xtm)
                
                tdx=dot(self.dx, self.xtm)
                new_err=dot(tdx, transpose(tdx))
                print('posterior error:', diag(new_err))
                print('---max, ---min:',  amax(diag(new_err)), amin(diag(new_err)))
                
                new_hm0_adj = dot(all_h, xinc)
                assim_conc = all_mean_y + new_hm0_adj
                inv_dev = (assim_conc - all_yobs)**2
                mod_dev = (all_mean_y0 - all_yobs)**2
                
                print('model deviation(max/min/mean):', amax(sqrt(mod_dev)),amin(sqrt(mod_dev)), mean(sqrt(mod_dev)))
                print('Daily prior deviation(max/min/mean):', amax(sqrt(ydev)), amin(sqrt(ydev)), mean(sqrt(ydev)))
                print('inv deviation(max/min/mean):', amax(sqrt(inv_dev)),amin(sqrt(inv_dev)), mean(sqrt(inv_dev)))
                
                # accumulative update on current day in current step
                if (fst_update):
                    all_xinc = array(xinc)
                    all_cur_xtm = array(cur_xtm)
                    all_cur_inc_m = array(cur_inc_m)
                    fst_update = False
                else:
                    all_xinc = all_xinc + array(xinc)
                    all_cur_inc_m = all_cur_inc_m + dot(all_cur_xtm, cur_inc_m)
                    all_cur_xtm = dot(all_cur_xtm, cur_xtm)
                    
                print('Step cum xinc max/min:', amax(all_xinc), amin(all_xinc))
                print('Step cum xtm max/min:', amax(diag(all_cur_xtm)), amin(diag(all_cur_xtm)))   
                
                usd_doy = usd_doy + 1 
                
                ## -----------------------------------------------------------
                ## step9: Daily Diagnosis (inv_parameters, inv_state)
                ## -----------------------------------------------------------
                ############ Daily posterior concentration output############################
                if (self.daily_conc_output):
                    
                    tdy = dot(all_h, tdx)
                    tdy = dot(tdy, transpose(tdy))
                    time = tm.doy_to_utc(idoy, 0, self.cur_yyyy)
                    outflnm = gcdf.inv_daily_conc_outpath + "/" + "Daily_conc" + "." + time[:10] + ".nc"
                    
                   
                    dgflux.daily_posterior_conc(outflnm = outflnm,\
                                                time = time,\
                                                lons = all_lon,\
                                                lats = all_lat,\
                                                mod_conc0 = all_mean_y0,\
                                                mod_conc = all_mean_y,\
                                                obs_conc = all_yobs,\
                                                inv_conc = assim_conc,\
                                                obs_err0 = all_oerr,\
                                                obs_err = diag(obs_err),\
                                                mod_err = dyr,\
                                                inv_err = diag(tdy),\
                                                project = data_ctl_rlt['project'].values
                                               )                                              
                    
                    
                
                ############ Daily inv_paras output############################
                if (self.daily_para_output):
                    sstep = r'%2.2d' % (self.cur_step)
                    #time = tm.doy_to_utc(idoy, 0, self.cur_yyyy)
                    outflnm = gcdf.inv_daily_outpath + "/" + "daily_assim" + "." + time[:10] + ".nc"
                    dgstate.etkf_paras(outflnm = outflnm,\
                                       istep = sstep,\
                                       period = time,\
                                       diag_mode = 'Daily',\
                                       xinc = xinc,\
                                       xtm = cur_xtm,\
                                       inc_m = cur_inc_m)
                    
                    print('Daily inversion parameter output finished!')
                    
                
                ############ Cumulative inv_paras output########################
                if (self.cumul_para_output):
                    sstep = r'%2.2d' % (self.cur_step)
                    sst = tm.doy_to_utc(cur_st_doy, 0, self.cur_yyyy)
                    send = tm.doy_to_utc(idoy, 0, self.cur_yyyy)
                    period = list([sst, send])
                    outflnm = gcdf.inv_cum_outpath + "/" + "cumul_assim" + "." + time[:10] + ".nc"
                    dgstate.etkf_paras(outflnm = outflnm,\
                                       istep = sstep,\
                                       period = period,\
                                       diag_mode = 'Cumulative',\
                                       xinc = all_xinc,\
                                       xtm = all_cur_xtm,\
                                       inc_m = all_cur_inc_m)
                 
                    print('Step cumulative parameter output finished!')
                    
                ############ Output Daily inv_state output #####################
                ######## for Daily posterior flux and its error calculatin######
                if ((self.cur_step > self.lag_window) & self.daily_flux_output):
                    print('current step:', self.cur_step)
                    print('lag_window:', self.lag_window)
                    print('daily_flux_output:', daily_flux_output)
                    ## calculating the state of the current time
                    cur_mean_x0 = self.mean_x0
                    cur_mean_x = self.mean_x
                    cur_err_x0 = diag(dot(self.dx, transpose(self.dx)))
                    tdx = self.dx * self.xtm
                    cur_err_x = diag(dot(tdx, transpose(tdx)))
                    
                    ## inv_state used to calculate daily posterior flux
                    time = tm.doy_to_utc(idoy, 0, self.cur_yyyy)
                    outflnm = gcdf.inv_daily_state_outpath + "/" + "Daily_assim_state" + "." + time[:10] + ".nc"
                    inv_paras = dgstate.output_daily_state(outflnm = outflnm,\
                                                           date = time,\
                                                           mean_x0 = cur_mean_x0,\
                                                           err_x0 = cur_err_x0,\
                                                           mean_x = cur_mean_x,\
                                                           err_x = cur_err_x,\
                                                           nc_output = False)
                    
                    ## Daily posterior flux and error calculation and output
                    outflnm = gcdf.inv_daily_flux_outpath + "/" + "Daily_flux" + "." + time[:10] + ".nc"
                    dgflux.daily_posterior_flux(outflnm = outflnm,\
                                                doy = idoy,\
                                                yyyy = self.cur_yyyy,\
                                                emis_fn = gcdf.emis_file,\
                                                mask_fn = gcdf.mask_file,\
                                                inv_factor = inv_paras)
                    
                    print('Daily posterior flux output finished!')
                    
                ######## update daily boundaries with single tracer ######    
                if (self.daily_sgl_bounds):
                    
                    dgbounds.daily_sgl_bounds_update(doy = idoy,\
                                                     yyyy = self.cur_yyyy,\
                                                     bounds_fn = gcdf.sgl_bounds_path,\
                                                     en_bounds_fn = gcdf.en_bounds_path,\
                                                     new_bounds_fn = gcdf.new_sgl_bounds_path,\
                                                     cur_meanx = self.mean_x)
                    
                    print('Daily single boundary update finished!')
                    
                ######## update daily boundaries with ensemble tracers ######   
                if (self.daily_en_bounds):    
                    tdx = self.dx * self.xtm  ## tdy = dy * self.xtm
                    cur_err_x = diag(dot(tdx, transpose(tdx)))
                    cur_err_x = sqrt(cur_err_x)
                    
                    dgbounds.daily_en_bounds_update(doy = idoy,\
                                                    yyyy = self.cur_yyyy,\
                                                    en_bounds_fn = gcdf.en_bounds_path,\
                                                    new_bounds_fn = gcdf.new_en_bounds_path,\
                                                    cur_xtm = cur_err_x )
                    
                    print('Daily ensemble boundaries update finished!')
                    
            ## -----------------------------------------------------------
            ## step10: output inversion parameters of current step 
            ## -----------------------------------------------------------

            step_xinc = all_xinc            
            step_xtm = all_cur_xtm
            step_inc_m = all_cur_inc_m
            ##cur_inc_m: accumulative increment of current step
           
        
            ############ Step inv_paras output########################
            if (self.step_para_output):
                sstep = r'%2.2d' % (self.cur_step)
                outflnm = gcdf.inv_step_outpath + "/" + "step_assim" + "." + sstep + ".nc"
            
                sst = tm.doy_to_utc(cur_st_doy, 0, self.cur_yyyy)
                send = tm.doy_to_utc(cur_end_doy, 0, self.cur_yyyy)
                period = list([sst, send])
                print('Step xinc max/min:', amax(step_xinc), amin(step_xinc))
                print('Step xtm max/min:', amax(diag(step_xtm)), amin(diag(step_xtm)))
                dgstate.etkf_paras(outflnm = outflnm,\
                                   istep = sstep,\
                                   period = period,\
                                   diag_mode = 'Step',\
                                   xinc = step_xinc,\
                                   xtm = step_xtm,\
                                   inc_m = step_inc_m)
                
            print('Step inversion parameters output finished!')
      
            
        else:

            sstep = r'%2.2d' % (self.cur_step)

            incflnm = gcdf.inv_step_outpath + "/" + "step_assim" + "." + sstep + ".nc"

            etkf_pars = xr.open_dataset(incflnm)

            step_xinc = etkf_pars['xinc'].values
            step_xtm = etkf_pars['xtm'].values
            step_inc_m = etkf_pars['inc_m'].values

            # update inv parameters of current step
            self.mean_x = self.mean_x + step_xinc
            self.xtm = dot(self.xtm, step_xtm)
            
            print('Step inversion parameters input finished!')
            
        # update inv state during whole process
        self.update_whole_state()
        print('#########################')
        print('AFTER STEP INVERSTION!')
        print('=========================')
        print('mean_x', self.mean_x)
        print('=========================')
        print('whole_mean_x', self.whole_mean_x)
        print('#########################')

        print('Whole state update finished!')
        
        # update inv_paramters generated in current step

 #       adx = dot(self.dx, self.xtm)
 #       self.step_dx.append(adx)
 #       self.x_inc.append(step_xinc)
 #       self.mean_x_step.append(self.mean_x)
        
        # update time paramters for the next step
        
        days_of_year = 366 if calendar.isleap(self.cur_yyyy) else 365
            
        if (cur_end_doy > days_of_year):
            self.rerun_doy_st = cur_end_doy - days_of_year
            self.rerun_yyyy_st += 1
        else:
            self.rerun_doy_st = self.rerun_doy_end
        
        self.cur_step = self.cur_step + 1

        ## -----------------------------------------------------------
        ## step11: prepare restart files for next step
        ## -----------------------------------------------------------
        
        if (self.cur_step >= rerun_step and self.cur_step != gcdf.step_end):
            if (self.full_rerun):
                # the correlation xinc will be added to the system 
                restart_xinc_m = step_xinc + all_cor_mean_x
            else:
                # self.mean_x include the correlation inc
                restart_xinc_m = self.mean_x

            print('restart_xinc_m: ', restart_xinc_m)
            
            if (gcdf.sep_en):
                rg.prepare_restart3(cur_step = self.cur_step,\
                                   doy_st = cur_st_doy,\
                                   doy_end = cur_end_doy,\
                                   yyyy = self.cur_yyyy,\
                                   inc_m = restart_xinc_m,\
                                   fst_rst = self.first_rerun)
            else:
                rg.prepare_restart(cur_step = self.cur_step,\
                                   doy_st = cur_st_doy,\
                                   doy_end = cur_end_doy,\
                                   yyyy = self.cur_yyyy,\
                                   inc_m = restart_xinc_m,\
                                   fst_rst = self.first_rerun)
            
            print('Restart files preparation finished!')
            

    
    def prepare_step_state(self):
        '''
        initialize or update transform parameters of current step and calculate 
        temporal and spatial error covariance 
        
        Args:
            self.dx          
            self.xtm
            self.mean_x0
            self.mean_x
            self.err
            self.cor_mean_x
            
        Returns:
            full_cor_mean_x: temporal correlations of the passed periods 
            
        '''
        if gcdf.prior_err_input:
            file  =gcdf.err_file.replace('YYYY', str(self.cur_yyyy)).replace('ddd', str(self.cur_doy).zfill(3))
            err_cor = xr.open_dataset(file)['err_cov'].values
        else:
            err_cor = perr.spatial_err_cov(emis_fn = gcdf.emis_file,\
                                           mask_fn = gcdf.mask_file,\
                                           yyyy = self.cur_yyyy,\
                                           doy = self.cur_doy,\
                                           cor_length = gcdf.cor_length)
        pos = self.nbias # default = 0
        nx  = self.ntracers
        ne  = nx
        
        # initialize the parameters of curent step
        dx = ones(nx, float)
        dx = dx/sqrt(self.xnorm)
        dx = diag(dx)
        
        dx = dot(dx, err_cor)            # pertubation parameters (omega)
        dx = gcdf.pri_err * dx
        mean_x = zeros(nx, float)        # increment parameters (wa)
        #print('err_cor', shape(err_cor))
        
        
        if (self.cur_step == gcdf.step_start):
            # --------------------------------------------  
            # structure of x / whole_x
            #     0--nbias-1:  bias  in observations
            #     nbias--nbias+nregion: the coefficients for regional perturbations in period  1
            #     nbias+nregion---nbias+nregion+nregion: coefficients  for regional perturbations 2
            #
            # ------------------------------------------
            # error covariances
            #    prior x 
            #        P=dot(dx, transpose(dx))
            #           
            #    posteriori x
            #        pdx=dot(dx, xtm)
            #        P=dot(pdx, transpose(pdx))
            #            
            #----------------------------------------------
            
            # determine the size of parameter matrix
            nx = pos + nx
            ne = nx
            self.nx, self.ne = nx, ne
            
            # initialize increment paras for frist assim
            self.mean_x0 = zeros(nx,float)
            self.mean_x = zeros(nx, float)
            self.mean_x[0:pos] = self.mod_bias
            
            # the error state paras at current step
            all_dx = identity(ne, float)
            dx_bias = array([self.bias_err]*pos)
            dx_bias = diag(dx_bias)
            all_dx[0:pos, 0:pos] = dx_bias[:,:]
            all_dx[pos:nx, pos:ne] = dx[:,:]
            
            # initialize pertucation paras for frist assim
            self.dx = array(all_dx)
            self.xtm = identity(ne,float)
            
            # initialize temp-spatial covariance 
            cor_mean_x = zeros(nx, float)
            self.cor_mean_x = cor_mean_x
            
            # initialize error matrix
            new_err=dot(self.dx, transpose(self.dx))
            self.err=array(new_err)
               
        else:
            
            if(self.cur_step < self.maxstep):
                
                # determine the size of paras matrix 
                old_nx, old_ne = shape(self.dx)
                new_nx = old_nx + nx
                new_ne = old_ne + ne
                nxs = old_nx    # tmp_old_nx
                nxe = old_ne    # tmp_old_ne
                
                # determine retained index of paras matrix
                retain_x_idx = list(range(0, old_nx))
                retain_x_idx = array(retain_x_idx)
                
                
            else:
                # determine the size of paras matrix 
                old_nx, old_ne = shape(self.dx)
                new_nx = old_nx 
                new_ne = old_ne 
                nxs = old_nx - nx
                nxe = old_ne - ne
                
                # determine retained index of paras matrix
                retain_x_idx = list(range(0, pos)) + list(range(pos + nx, old_nx))
                retain_x_idx = array(retain_x_idx)
        
        
            # update self.dx for current step 
            # (add spatial correlation)
            new_dx = zeros([new_nx, new_ne],float)
            bdx = self.dx[retain_x_idx, :]
            bdx = bdx[:, retain_x_idx]
            new_dx[0:nxs, 0:nxe] = bdx[0:nxs, 0:nxe]
            new_dx[nxs:new_nx, nxe:new_ne] = dx[0:nx, 0:ne]
            self.dx = array(new_dx)    
            self.nx, self.ne = shape(self.dx)
            
            # prepare self.xtm for current step 
            #     |  xtm_old     0  |
            # xtm=|                 |* enlarge_factor
            #     |    0         I  |
            
            self.xtm = self.xtm_enlarge_factor * self.xtm
            bxtm = self.xtm[retain_x_idx, :]
            bxtm = bxtm[:, retain_x_idx]
            xtm = identity(new_ne, float)
            xtm[0:nxs, 0:nxe] = bxtm[0:nxs,  0:nxe]
            self.xtm = array(xtm)
            
            # add temporal correlation  into prior error matrix
            #
            # temporal correlated DX is generated as  
            # 
            # DX=
            #    |  DX0     0   |       | DX0  0   |  | I       0  |
            #    |  cor_dx  DX1 |   ==  |  0   DX1 |  | tc_xtm, I  |
            #
            #
            # as a results, 
            # the error covariance for newly included x1 will be
            #   B=DX1*(I+tc_xtm*tc_xtm^T)*DX1^T, so a rescaling is needed 
            
            tc_xtm = identity(new_ne, float)
            
            nperiod = int((nxs - pos)/nx)
            sum_row_factor = zeros(nx, float)
            sum_row_factor[:] = 1.0
            ipos = pos
            
            # tc_xtm is inform of 
            # |a, 0, 0, ..., 0_nx,a2, 0, 0, 0, ...,   0_nx, ..., ..., |  
            # |0, a, 0, ..., 0_nx, 0, a2, 0, 0, ...,  0_nx, ... , ... |
            # |0, 0, a, ..., 0_nx, 0. 0,  a2, 0, ..., 0_nx, ... , ... |
            # |--- period -2 ---| -------period -1  -------|---period 0---|
            for iperiod in range(nperiod):
                
                cor_factor = self.tcor_factor * exp(-(nperiod-iperiod)/self.tcor_len)
                
                for irow in range(nx):
                    tc_xtm[nxs + irow, ipos + irow] = cor_factor
                    sum_row_factor[irow] = sum_row_factor[irow] + cor_factor * cor_factor
                ipos = ipos + nx
                
            # weight function to aviod enlarging self.dx
            
            sum_row_factor = divide(ones_like(sum_row_factor), sqrt(sum_row_factor))

            
            # update self.dx for current step once more 
            # (add temporal correlation)
            for irow in range(nx):
                self.dx[nxe + irow, :] = sum_row_factor[irow] * self.dx[nxe + irow, :]
                
            self.dx = dot(self.dx, tc_xtm)
            
            
            # update error matrix
            tdx=dot(self.dx, self.xtm)
            new_err=dot(tdx, transpose(tdx))
            
            self.err=array(new_err)
            
            # caculate temp-spatial covariance for self.mean_x
            base_mean_x = self.mean_x[retain_x_idx]
            #base_mean_x = base_mean_x[0:nxs]
            
            base_dx = self.dx[0:nxs, 0:nxs]
            cor_dx = self.dx[nxs:, 0:nxs]
            
            cor_mean_x = nlg.solve(base_dx, base_mean_x)
            cor_mean_x = dot(cor_dx, cor_mean_x)
            self.cor_mean_x = cor_mean_x
            
            # prepare self.mean_x for current step (before assim)
            new_meanx = zeros(new_nx, float)
            new_meanx[0:nxs] = self.mean_x[retain_x_idx]
            new_meanx[nxs:new_nx] = mean_x[0:nx] + cor_mean_x[0:nx]
            self.mean_x = array(new_meanx)
            print('############################################')
            print('cor_mean_x: ', amax(cor_mean_x), amin(cor_mean_x), mean(cor_mean_x))
            print('############################################')
            # self.mean_x0: zero (for future use)
            new_meanx0 = zeros(new_nx, float)
            new_meanx0[0:nxs] = self.mean_x0[retain_x_idx]
            new_meanx0[nxs:new_nx] = mean_x[0:nx]
            self.mean_x0=array(new_meanx0)

            
        # return value: the temp-spatial correlation to mean_x in currentstep
    
        full_cor_mean_x = zeros(shape(self.mean_x), float)
        full_cor_mean_x[-nx:] = cor_mean_x[:]
        #print('full_cor_mean_x', shape(full_cor_mean_x))
        #print('cor_mean_x', shape(cor_mean_x))
        return full_cor_mean_x
                
      
    def prepare_whole_state(self):
        '''
        prepare accumulated mean_x and err for diagnosis
              
        '''
        n_tra = self.ntracers
        
        if (self.cur_step < self.maxstep):
            
            self.whole_mean_x0 = array(self.mean_x0)
            self.whole_mean_x = array(self.mean_x)
                
            new_err = dot(self.dx, transpose(self.dx))
            self.whole_err_x0 = diag(new_err)
            self.whole_err_x = diag(new_err)
            
        else:
            ### xxxx[retain_nwhole:]  
            
            old_nwhole = size(self.whole_mean_x)
            new_whole_mean_x = zeros(old_nwhole + n_tra, float)
            new_whole_mean_x[0:old_nwhole] = self.whole_mean_x0[:]
            new_whole_mean_x[old_nwhole:] = self.mean_x0[-n_tra:]
            
            self.whole_mean_x0 = array(new_whole_mean_x)
            
            new_whole_mean_x[0:old_nwhole] = self.whole_mean_x[:]
            new_whole_mean_x[old_nwhole:] = self.mean_x[-n_tra:]
            
            self.whole_mean_x = array(new_whole_mean_x)

            # err
            new_err = dot(self.dx, transpose(self.dx))
            new_err = diag(new_err)
            new_whole_err_x = zeros(old_nwhole +  n_tra, float)
            new_whole_err_x[0:old_nwhole] = self.whole_err_x0[:]
            new_whole_err_x[old_nwhole:] = new_err[-n_tra:]
            
            self.whole_err_x0 = array(new_whole_err_x)
         
            new_err = diag(self.err)
            new_whole_err_x[0:old_nwhole] = self.whole_err_x[:]
            new_whole_err_x[old_nwhole:] = new_err[-n_tra:]
            
            self.whole_err_x = array(new_whole_err_x)
            
#        del new_whole_mean_flux
#        del new_whole_mean_x
        
        return None
            
            
    def update_whole_state(self):
        '''
        update accumulated mean_x and err for diagnosis
        
        notes: _x0: prior parameters
               _x : posterior parameters
        
        '''
        ncur_x = size(self.mean_x)
        nwhole_x = size(self.whole_mean_x)
        nold_x = nwhole_x-ncur_x
        
        post_err = dot(self.dx, self.xtm)
        post_err = dot(post_err, transpose(post_err))
        post_err = diag(post_err)

        
        if (nold_x == 0): # self.cur_step< self.maxstep
            
            self.whole_mean_x = array(self.mean_x)
            self.whole_err_x = array(post_err)
            
        else:  # current step larger than lagwinow length            

        
            self.whole_mean_x[nold_x:nwhole_x] = self.mean_x[:]
            self.whole_err_x[nold_x:nwhole_x] = post_err[:]

if (__name__=='__main__'):
    
    main = assim_step(use_time_cor = True,\
                      use_fixed_mod_err = True,\
                      use_full_r = True,\
                      use_sparse = True,\
                      add_rnd_err = False,\
                      do_update = True,\
                      full_rerun = True,\
                      daily_prior_output = True,\
                      daily_para_output = False,\
                      cumulative_para_output = False,\
                      step_para_output = True,\
                      daily_posterior_output = False,\
                      daily_updata_en_bounds = False,\
                      daily_updata_sgl_bounds = False)
    for i in range(4):
        main.do_one_step(update_step = 3, rerun_step = 3)