#!/exports/csce/datastore/geos/users/s2008420/miniconda/base/envs/geo/bin/python
from numpy import *
import time_module as tm
import xarray as xr
import geos_chem_def as gcdf 


def daily_prior_state(outflnm, time, lons, lats, delta_x, hm, mean_x, mod_rlt, obs, obs_err):
    ''' 
    Output the daily etkf states before inversion 

    Args:
        outflnm: output file path + filname
        time: inversion date
        lons, lats: positions of the states
        delta_x: daily ensemble pertubations of x(delta_X)
        hm: linear relationship of pertubations H: delta Yf_0/ delta Xf_0
        mean_x: parameters to prior states
        mod_rlt: model results after adjustment
        obs: observations use in inversion
        obs_err: observation errors

    Returns:
        None

    '''
    ## nobs: the number of observations
    ## ne: the number of ensemble pertubations, generally equals to lagwindows * ntracers
    ## nx: the length of increments matrix, equals to ne
    
    nobs = len(lons)
    nx = len(mean_x)
    ne = nx


    prior_state = xr.Dataset({'delta_x': (['nx', 'ne'], delta_x),
                             'hm': (['nobs', 'ne'], hm),
                             'mean_x': (['nx'], mean_x),
                             'mod_rlt': (['nobs'], mod_rlt),
                             'obs': (['nobs'], obs),
                             'obs_err': (['nobs'], obs_err),
                             'obs_lon': (['nobs'], lons),
                             'obs_lat': (['nobs'], lats)},
                            coords = {'nobs': arange(nobs),
                                     'nx': arange(nx),
                                     'ne': arange(ne),
                                     'reference_time': time})
    
    prior_state.attrs['Description'] = 'Daily prior states before inversion'
    prior_state.attrs['obs_lon'] = 'Longitudes of the observations'
    prior_state.attrs['obs_lat'] = 'Latitudes of the observations'
    prior_state.attrs['delta_x'] = 'Daily ensemble pertubations of x(delta_X)'
    prior_state.attrs['hm'] = 'H: delta Yf_0/ delta Xf_0'
    prior_state.attrs['mean_x'] = 'Increment parameters to prior states'
    prior_state.attrs['mod_rlt'] = 'Model results on observation sites'
    prior_state.attrs['obs'] = 'Daily observations use in inversion'
    prior_state.attrs['obs_err'] = 'Observation errors'
    
    prior_state.to_netcdf(outflnm)
    
    
def etkf_paras(outflnm, istep, period, diag_mode, xinc, xtm, inc_m):
    ''' 
    Output the daily, cumulative or step inversion parameters after inversion

    Args:
        outflnm: output file path + filename
        istep: current step
        period: inversion period
        diag_mode: diagnostic mode: 'Daily', 'Cumulative', 'Step'
        ------------------------
        if output daily results:
            xinc: current daily increment to prior estimates Xf (wa * delta Xf_0)  
            xtm: current factors apply for delta Xf(omega)
            inc_m: parameters for pertubations (wa) 
        ------------------------    
        if output daily cumulative results:
            all_xinc: cumulative increment to prior estimates Xf_0 (wa * delta Xf_0)  
            all_cur_xtm: cumulative factors apply for delta Xf_0(omega)
            all_cur_inc_m: accumulative increment in current step
        ------------------------    
        if output step results:
            step_xinc: cumulative increment of current step  
            step_xtm: cumulative factors of current step  
            step_inc_m: accumulative increment of current step (wa * delta Xf) 

    Returns:
        None

    '''
    ## ne: the length of pertubations([ne,ne])
    ## nx: the size the states, generally equals to lagwindows * ntracers
    
    ne, ne = shape(xtm)
    nx = size(inc_m)
    
    xnx, xne = arange(nx), arange(ne)
    
    etkf_pars = xr.Dataset({'xinc':(['nx'], xinc),
                            'xtm':(['nx','ne'], xtm),
                            'inc_m':(['nx'], inc_m)},
                           coords = {'nx': xnx,
                                     'ne': xne,
                                     'Step': istep,
                                     'Diag_mod': diag_mode,
                                     'Reference_time': period})
    
    
    etkf_pars.attrs['Description'] = 'Assimilation inversion parameters'
    etkf_pars.attrs['Reference_time'] = 'Time stamp'
    etkf_pars.attrs['Step'] = 'Step stamp'
    etkf_pars.attrs['Diag_mod'] = 'Diagnostic mode'
    etkf_pars.attrs['xinc'] = 'Increment to prior estimates Xf_0(wa * delta Xf)'
    etkf_pars.attrs['xtm'] = 'factors apply for delta Xf(omega)'
    if diag_mode is 'Daily':
        etkf_pars.attrs['inc_m'] = 'Parameters for pertubations (wa)'
    elif  diag_mode is 'Cumulative': 
        etkf_pars.attrs['inc_m'] = 'accumulative increment on current day form the begining '
    else:
        etkf_pars.attrs['inc_m'] = 'accumulative increment of current step (wa * delta Xf)'
    
    etkf_pars.to_netcdf(outflnm)
    
    

def output_daily_state(outflnm, date, mean_x0, err_x0, mean_x, err_x, nc_output):
    '''
    output cumulative daily assimilation parameters 
    
    Args:
        outflnm: output file path + filename
        date: Time stamp
        mean_x0: prior increment parameters 
        mean_x: the finial increment to prior state at current time
        err_x0: prior error factors
        err_x: posterior error factors after iteration at current time
        nc_output: True or False switch for output 
    Returns:
        daily_assim_factor
        
    '''
    n_step = gcdf.inv_lag_window 
    n_sous = gcdf.n_sous
    n_regs = gcdf.n_regs
    sous_name = gcdf.emis_var_name
    
    ## reshape 1-D parameters to 3-D (nstep, n_sous, n_regs)
    mean_x0 = reshape(mean_x0, (n_step, n_sous, n_regs)) 
    mean_x = reshape(mean_x, (n_step, n_sous, n_regs)) 
    err_x0 = reshape(err_x0, (n_step, n_sous, n_regs))
    err_x = reshape(err_x, (n_step, n_sous, n_regs))
    
    ## collect the completed inversion state at current time
    cur_mean_x0 = squeeze(mean_x0[0, :, :]) 
    cur_mean_x = squeeze(mean_x[0, :, :])
    cur_err_x0 = squeeze(err_x0[0, :, :])
    cur_err_x = squeeze(err_x[0, :, :])
    
    
    ## write into a dataset
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
    
    if (nc_output):
        daily_assim_factor.to_netcdf(outflnm)
    
    return daily_assim_factor
    

def output_whole_states(outflnm, periods, nstep, mean_x0, err_x0, mean_x, err_x, diag_whole_state):
    '''
    output inversion state of the whole periods
    
    Args:
        outflnm: output file path + filename
        time: inversion period
        nstep: inversion steps
        mean_x0: prior increment parameters during the whole periods 
        mean_x: the finial increment to prior state after iterative assimilation
        err_x0: prior error factors
        err_x: posterior error factors after iteration
        
    Returns:
        assim_factor
    '''
    n_step = len(periods)
    n_sous = gcdf.n_sous
    n_regs = gcdf.n_regs
    sous_name = gcdf.emis_var_name
    ## reshape 1-D parameters to 3-D (nstep, n_sous, n_regs)
    mean_x0 = reshape(mean_x0, (n_step, n_sous, n_regs)) 
    mean_x = reshape(mean_x, (n_step, n_sous, n_regs)) 
    err_x0 = reshape(err_x0, (n_step, n_sous, n_regs))
    err_x = reshape(err_x, (n_step, n_sous, n_regs))
    
    ## write into a dataset
    assim_factor = xr.Dataset({'mean_x0':(['time','i_sous','i_reg'], mean_x0),
                               'mean_x':(['time','i_sous','i_reg'], mean_x),
                               'err_x0':(['time','i_sous','i_reg'], err_x0),
                               'err_x':(['time','i_sous','i_reg'], err_x)},
                              coords = {'time': periods,
                                        'i_sous': sous_name,
                                        'i_reg': range(n_regs)})
    
    assim_factor.attrs['Description'] = 'Assimilated inversion state of the whole periods'
    assim_factor.attrs['mean_x0'] = 'Prior increment parameters'
    assim_factor.attrs['mean_x'] = 'Posterior increment parameters'
    assim_factor.attrs['err_x0'] = 'Porior error factors'
    assim_factor.attrs['err_x'] = 'Posterior error factors'    
    
    
    if (diag_whole_state):
        assim_factor.to_netcdf(outflnm)
        
    return assim_factor
    


if (__name__=='__main__'):
    
    # daily_prior_state
    nobs = 1060
    nx = 40
    ne = nx
    delta_x = random.rand(nx, ne)
    hm = random.rand(nobs, ne)
    mean_x = random.rand(nx)
    mod_rlt = random.rand(nobs)
    obs = random.rand(nobs)
    obs_err = random.rand(nobs)
    lons = random.rand(nobs)
    lats = random.rand(nobs)
    time = tm.doy_to_utc(28, 0, 2015)
    outflnm = gcdf.daily_state_path + "/" + "daily_prior" + "." + time[:10] + ".nc"
    daily_prior_state(outflnm, time, lons, lats, delta_x, hm, mean_x, mod_rlt, obs, obs_err)
    
    # etkf_paras
    xinc = random.rand(nx)
    xtm = random.rand(ne, ne)
    inc_m = random.rand(nx)
    diag_mode = 'Daily'
    istep = 1
    
    outflnm = gcdf.inv_daily_outpath + "/" + "daily_assim" + "." + time[:10] + ".nc"
    etkf_paras(outflnm, istep, time, diag_mode, xinc, xtm, inc_m)
    
    
    

