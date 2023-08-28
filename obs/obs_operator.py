#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
import sys
sys.path.append('/home/shzhu/enkf_ch4/rundirs/grid_2x25_new/obs/')
import numpy as np
import numba as nb
import pandas as pd
import xarray as xr 
from xbpch import open_bpchdataset
from scipy import interpolate
from scipy.interpolate import interp1d
import operator_config as config
import tropomi
import gosat 
import noaa 
import os

class obs_operator:
    '''
    integration of measurements sampling results 
    '''
    def __init__(self,\
                 obs_dict):
        '''
        Args:
            obs_dict: measurements used in assimilation
            
        '''
        
        self.obs_dict = obs_dict
  
        
    def hm_concat(self, date):
        '''
        Concatenate sampling result from different measurements
        note: just concatenate existed sampling results
        
        Args:
            date: current date in datetime format
        Return:
            
        '''
 
        meas = []
        for mea in self.obs_dict:

            mea_dic =eval('config.'+mea+'_info')
            file =  mea_dic['sampl_path'].replace('YYYY', date.strftime('%Y'))\
            +mea_dic['sampl_file'].replace('YYYYMMDD', date.strftime('%Y%m%d'))
            if not os.path.exists(file):
                continue
            data = xr.open_dataset(file)
            if mea == 'gosat':
                data = data.drop_vars('time')
            if mea == 'tropomi':
                data = data.drop_vars('obs_raw')  
            if mea == 'noaa':
                projects = config.noaa_info['project']
                noaa_list = list()
                for item in projects:
                    data0 = data.where(data['project'] == item, drop =True)
                    noaa_list.append(data0)
                del data
                data = xr.concat(noaa_list, 'n')
                data = data.drop_vars(['altitude', 'long_name'])
            meas.append(data)
        if len(meas) == 0:
            print('NO MEASUREMENTS ON:', date.strftime('%Y-%m-%d'))
            return False
        daily_mea = xr.concat(meas, dim = 'n')
        daily_mea['site'] = daily_mea['n'].values.astype(np.str)
        daily_mea = daily_mea.assign_coords({'n': range(len(daily_mea['n']))})
        daily_mea.attrs['Description'] = 'Concatenation result from'+ '-'.join(self.obs_dict)
        
        return daily_mea
    
    
    def hm0_concat(self, date, cur_step):
        '''
        Concatenate sampling result from different measurements
        note: Concatenate updated rerun output
        
        '''
        meas = []
       
        for mea in config.measurements:
            
            mea_obj = eval(mea + '.'+ mea+'_obs')(date = date,\
                                                  info= eval('config.'+mea+'_info'),\
                                                  sp_diagn_info = config.sp_diagn_info,\
                                                  en_info_list = config.en_info_list,\
                                                  en_item_table = config.entry_table)
            
            
            mea_dic = eval('config.'+mea+'_info')
            
            hm_file = mea_dic['sampl_path'].replace('YYYY', date.strftime('%Y'))\
            +mea_dic['sampl_file'].replace('YYYYMMDD', date.strftime('%Y%m%d'))
            if not os.path.exists(hm_file):
                continue
            data = mea_obj.gen_hm0(hm_file, False, cur_step)

            if mea == 'gosat':
                data = data.drop_vars('time')
            if mea == 'tropomi':
                data = data.drop_vars('obs_raw')  
            if mea == 'noaa':
                projects = config.noaa_info['project']
                noaa_list = list()
#                 print(data)
                for item in projects:
                    data0 = data.where(data['project'] == item, drop =True)
                    noaa_list.append(data0)
                del data
                data = xr.concat(noaa_list, 'n')
                data = data.drop_vars(['altitude', 'long_name'])
            meas.append(data)
        if len(meas) == 0:
            print('NO MEASUREMENTS ON:', date.strftime('%Y-%m-%d'))
            return False  
        daily_mea = xr.concat(meas, dim = 'n')
        daily_mea['site'] = daily_mea['n'].values.astype(np.str)
        daily_mea = daily_mea.assign_coords({'n': range(len(daily_mea['n']))})
        daily_mea.attrs['Description'] = 'Concatenation result from '+ '-'.join(self.obs_dict)

        return daily_mea


if __name__ == '__main__':
    
    date_st = config.date_st
    date_end = config.date_end

    period = pd.date_range(start = date_st, end = date_end, freq= 'D')
    ch4_hm = obs_operator(obs_dict = config.measurements)
    for date in period:
        print(date)
        out_file = 'test.'+ date.strftime('%Y%m%d')+'.nc'
        conc_data = ch4_hm.hm0_concat(date)
        conc_data.to_netcdf(out_file)