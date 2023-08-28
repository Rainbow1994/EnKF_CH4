#! /home/shzhu/apps/miniconda3/envs/idp/bin/python

import numpy as np
import numba as nb
import pandas as pd
import xarray as xr 
from xbpch import open_bpchdataset
from scipy import interpolate
from scipy.interpolate import interp1d
import operator_config as config
from datetime import datetime
import os

class noaa_obs:
    '''
    assimilation sampling for noaa observations
    '''
    def __init__(self,\
                 date,\
                 info,\
                 sp_diagn_info,\
                 en_info_list,\
                 en_item_table):
        '''
        Args:
            date: current date
            noaa_info: noaa observations information
            sp_diagn_info: single diagnosis information
            en_info_list: ensemble tagrun file list(contains ensemble outputs information)
            en_item_table: ensemble step information at current day
        '''
        self.date = date
        self.noaa_info = info
        self.sp_info = sp_diagn_info
        self.en_info = en_info_list
        self.en_item = en_item_table

    #@nb.jit(nopython=True,parallel=True) 
    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    
    def setup_nc_daily_feild(self,
                             combine_en,\
                             lons,\
                             lats):
    
    

        lon =[]
        lat =[]
        for i_lon, i_lat in zip(lons, lats):
            lon.append(self.find_nearest(combine_en.lon.values, i_lon))
            lat.append(self.find_nearest(combine_en.lat.values, i_lat))

        ind_x = xr.DataArray(lon, dims='ns')
        ind_y = xr.DataArray(lat, dims='ns')
        model_ch4 = combine_en.isel(lon=ind_x, lat=ind_y)

        return model_ch4
    
    def get_ensemble_list(self):


        """ read into entry table """

        fl_entry = open(self.en_item, "r")
        lines=fl_entry.readlines()
        fl_entry.close()
        yyyy = int(self.date.strftime('%Y'))
        doy = int(self.date.strftime('%j'))
        ist=0
        table_entry_list=list()
        
        for line in lines:
            if ('mem_st' in line):
                break
            ist=ist+1

        for line in lines[ist+1:]:
            line=line.strip()
            if (len(line)>0):
                terms=line.split(',')
                yyyyddd_st=terms[3].strip()+terms[5].strip()
                yyyyddd_st=int(yyyyddd_st)
                yyyyddd_end=terms[4].strip()+terms[6].strip()
                yyyyddd_end=int(yyyyddd_end)
                vals=[yyyyddd_st, yyyyddd_end, int(terms[1]), int(terms[2]), int(terms[0]), terms[8]]
                table_entry_list.append(vals)

        cur_time=yyyy*1000+doy
#         em_st_list=list()
#         em_end_list=list()
#         em_yst_list=list()
        em_step_list=list()
        em_name_list=list()

        """ check the ensemble member to available """
        for vals in table_entry_list:
            if (cur_time>=vals[0] and cur_time <vals[1]):
#                 em_st_list.append(vals[2])
#                 em_end_list.append(vals[3])
#                 em_yyyy=int(vals[0]/1000.0)
#                 em_yst_list.append(em_yyyy)
                em_step=vals[4]
                em_step_list.append(em_step)
                em_name=vals[5]
                em_name_list.append(em_name)
#         return em_st_list, em_end_list, em_yst_list, em_step_list
        return em_name_list
    
    def combine_xbpch(self, step):
        '''
        combine ensemble run result with seperated emission sources& sources
        Args:
            file_info_list: tag run output files information
               eg, file1_info = {
               'filepath': '/data/shzhu/tagch4_2x25/ND51/2010/',
               'filename': 'ts_satellite.ST000.EN0001-EN0037.20100101.bpch',
               'diaginfo_file' : '/data/shzhu/tagch4_2x25/ND51/diaginfo_em.dat',
               'tracerinfo_file' : '/data/shzhu/tagch4_2x25/ND51/tracerinfo_EN0001-EN0037.dat',
               'var_prefix': ['IJ_AVG_S_TAG_AN','IJ_AVG_S_TAG_FF'],
               'trac_num':37,
               'trac_st':1}
        Return:
            en: combined result
        '''
        
        date = self.date.strftime('%Y%m%d')
        yyyy = self.date.strftime('%Y')
        
        var_list = []
        print(step)
        for file in self.en_info:
            
            filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
            filename = file['filepath']+ str(yyyy) +'/'+filename
            if not os.path.exists(filename):
                cur_yyyy = str(int(yyyy) -1)
                filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
                filename = file['filepath']+ str(cur_yyyy) +'/'+filename
            
            tag_file = open_bpchdataset(filename,
                                   diaginfo_file=file['diaginfo_file'],
                                   tracerinfo_file=file['tracerinfo_file'])
            file_new = xr.Dataset()
            trac_st = file['trac_st']
            trac_num = file['trac_num']

            for i_tra in file['var_prefix']:
                for i in range(trac_num):
                    new_name = i_tra + str(i + trac_st)
                    old_name = i_tra + str(i + 1)

                    file_new[new_name] = tag_file[old_name]
            var_list.append(file_new)
        en = xr.merge(var_list)
        en.attrs = tag_file.attrs
        en.astype(np.float32)

        return en
    
    def gen_hm(self, sampl_file):
        '''
        read single and ensemble geos_chem run diagnostic ouput
        resample them to satellite observation sites
        according to layer pressure and its weight function to caculate XCH4
        ### note: read seperated ensemble tagrun!
        '''

        ## step1: get time information
        date = self.date.strftime('%Y%m%d')
        yyyy = self.date.strftime('%Y')

        ## step2: read satellite observation data
        
        obs_fn = self.noaa_info['filepath'] +self.noaa_info['filename']

        if os.path.exists(obs_fn):
            #trying to open a file in read mode
            obs = xr.open_dataset(obs_fn).load()
            obs = obs.sel(time =self.date.strftime('%Y-%m-%d'))
            if len(obs.where(obs['ch4'] >0, drop = True).ns) == 0:
                return False
            else:
                obs_sel = obs.where(obs['ch4'] >0, drop = True)
   
        else:
            return False


#         print( obs_sel.ns)
        
        ## step3: read geos-chem single tracer(rerun) output and do resample
            # sampling girds 
#         lons = np.array([site['lon'] for site in self.noaa_info['sites_list']])
#         lats = np.array([site['lat'] for site in self.noaa_info['sites_list']])
        lons = obs_sel['longitude'].values
        lats = obs_sel['latitude'].values

            # resampling bpch output
        filename = self.sp_info['filepath'] + self.sp_info['filename'].replace('YYYYMMDD', date)

        sp_ch4 = open_bpchdataset(filename,\
                                  diaginfo_file=self.sp_info['diaginfo_file'],\
                                  tracerinfo_file=self.sp_info['tracerinfo_file'])
        
        sp_mod = self.setup_nc_daily_feild(combine_en = sp_ch4,\
                                           lons = lons,\
                                           lats = lats) 

        

        sgl_data = sp_mod[self.sp_info['var_prefix']].values 

        hm0 = sp_mod[self.sp_info['var_prefix']].isel(lev = slice(None, 6)).mean(dim = 'lev').values


        ## step7: caculation xgp result of ensemble tracers[hm]
        em_step_list = self.get_ensemble_list()
        
        tag_xgp = []
        
        var_list= [''.join([config.emis_prefix,sous[4:],str(i_reg+1)]) for sous in config.emis_vars for i_reg in range(config.n_regs) ]
        for em_step in em_step_list:
            
            combine_em = self.combine_xbpch(em_step[:5])
            em_mod = self.setup_nc_daily_feild(combine_em,\
                                               lons,\
                                               lats)

            for var in var_list:

                em_data = em_mod[var].isel(lev = slice(None, 6)).mean(dim = 'lev').values
               
                tag_xgp.append(em_data)

        if np.array(tag_xgp).ndim == 1:
            tag_xgp = tag_xgp[None,:]
        hm = np.array(tag_xgp)
        hm = hm.transpose()


        ## step8: integration of sampling results
        ns, ntag = np.shape(hm) 
#         print(ns, ntag)
        samp_data = xr.Dataset(
            {
                'hm0':(['n'], hm0),
                'hm':(['n','t'], hm),
                'obs':(['n'], obs_sel['ch4'].values*1e9),
                'oerr':(['n'], obs_sel['ch4_unc'].values*1e9),
                'lon':(['n'],  obs_sel['longitude'].values),
                'lat':(['n'], obs_sel['latitude'].values),
                'project':(['n'], obs_sel['project'].values),
                'long_name':(['n'], obs_sel['long_name'].values),
                'altitude' :(['n'], obs_sel['altitude'].values)
            },
            coords= {
                'n': obs_sel['ns'].values,
                't':  np.arange(ntag),
#                                         'ref_time' : self.date.strftime('%Y-%m-%d')
            }
        )



#         samp_data['obs']= obs_sel['ch4']
#         samp_data['oerr'] = obs_sel['ch4_unc']
#         samp_data['lon'] = obs_sel['longitude']
#         samp_data['lat'] = obs_sel['latitude']
#         samp_data['time'] = self.date.strftime('%Y-%m-%d')



        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Matched noaa sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'noaa sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'noaa sample result of GEOS-Chem pertubation run'
        samp_data['hm'].attrs['units'] = '1e-9'
        
        samp_data.to_netcdf(sampl_file)

        return samp_data
    
    def gen_hm0(self, hm_file, sampl_file, cur_step):
        '''
        update geos_chem run diagnostic ouput
        resample them to satellite observation sites
        according to layer pressure and its weight function to caculate XCH4
        ### note: read existed ensemble tagrun sampling result!
        '''

         ## step1: get time information
        date = self.date.strftime('%Y%m%d')
        yyyy = self.date.strftime('%Y')

        ## step2: read satellite observation data
        
        obs_fn = self.noaa_info['filepath'] +self.noaa_info['filename']

        if os.path.exists(obs_fn):
            #trying to open a file in read mode
            obs = xr.open_dataset(obs_fn).load()
            obs = obs.sel(time =self.date.strftime('%Y-%m-%d'))
            if len(obs.where(obs['ch4'] >0, drop = True).ns) == 0:
                return False
            else:
                obs_sel = obs.where(obs['ch4'] >0, drop = True)
   
        else:
            return False



        
        ## step3: read geos-chem single tracer(rerun) output and do resample
            # sampling girds 
#         lons = np.array([site['lon'] for site in self.noaa_info['sites_list']])
#         lats = np.array([site['lat'] for site in self.noaa_info['sites_list']])
        lons = obs_sel['longitude'].values
        lats = obs_sel['latitude'].values

            # resampling bpch output
        filename = self.sp_info['filepath'] + self.sp_info['filename'].replace('YYYYMMDD', date).replace('ST000', 'ST'+ str(cur_step).zfill(3))

        sp_ch4 = open_bpchdataset(filename,\
                                  diaginfo_file=self.sp_info['diaginfo_file'],\
                                  tracerinfo_file=self.sp_info['tracerinfo_file'])
        
        sp_mod = self.setup_nc_daily_feild(combine_en = sp_ch4,\
                                           lons = lons,\
                                           lats = lats) 

        

        sgl_data = sp_mod[self.sp_info['var_prefix']].values 

        hm0 = sp_mod[self.sp_info['var_prefix']].isel(lev = slice(None, 6)).mean(dim = 'lev').values


        ## step7: integration of sampling results
        hm = xr.open_dataset(hm_file)['hm']
        ns, ntag = np.shape(hm)
       

        samp_data = xr.Dataset(
            {
                'hm0':(['n'], hm0),
                'hm':(['n','t'], hm),
                'obs':(['n'], obs_sel['ch4'].values*1e9),
                'oerr':(['n'], obs_sel['ch4_unc'].values*1e9),
                'lon':(['n'],  obs_sel['longitude'].values),
                'lat':(['n'], obs_sel['latitude'].values),
                'project':(['n'], obs_sel['project'].values),
                'long_name':(['n'], obs_sel['long_name'].values),
                'altitude' :(['n'], obs_sel['altitude'].values)
            },
            coords= {
                'n': obs_sel['ns'].values,
                't':  np.arange(ntag),
#                                         'ref_time' : self.date.strftime('%Y-%m-%d')
            }
        )

#         samp_data['obs']= obs_sel['ch4']
#         samp_data['oerr'] = obs_sel['ch4_unc']
#         samp_data['lon'] = obs_sel['longitude']
#         samp_data['lat'] = obs_sel['latitude']
#         samp_data['time'] = self.date.strftime('%Y-%m-%d')



        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Matched noaa sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'noaa sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'noaa sample result of GEOS-Chem pertubation run'
        samp_data['hm'].attrs['units'] = '1e-9'
        
        if sampl_file:  
            samp_data.to_netcdf(sampl_file)

        return samp_data
    
    
    
    def add_project_tag(self, hm_file, sampl_file):
        
        '''
        add project information for existed sample results
        '''
        ## step1: get time information
        date = self.date.strftime('%Y%m%d')
        yyyy = self.date.strftime('%Y')

        ## step2: read satellite observation data
        
        obs_fn = self.noaa_info['filepath'] +self.noaa_info['filename']

        if os.path.exists(obs_fn):
            #trying to open a file in read mode
            obs = xr.open_dataset(obs_fn).load()
            obs.close()
            obs = obs.sel(time =self.date.strftime('%Y-%m-%d'))
            if len(obs.where(obs['ch4'] >0, drop = True).ns) == 0:
                return False
            else:
                obs_sel = obs.where(obs['ch4'] >0, drop = True)
   
        else:
            return False




        ## step3: integration of sampling results
        samp = xr.open_dataset(hm_file).load()
        samp.close()
        hm = samp['hm'].values
        hm0 = samp['hm0'].values
        ns, ntag = np.shape(hm)
       

        samp_data = xr.Dataset(
            {
                'hm0':(['n'], hm0),
                'hm':(['n','t'], hm),
                'obs':(['n'], obs_sel['ch4'].values*1e9),
                'oerr':(['n'], obs_sel['ch4_unc'].values*1e9),
                'lon':(['n'],  obs_sel['longitude'].values),
                'lat':(['n'], obs_sel['latitude'].values),
                'project':(['n'], obs_sel['project'].values),
                'long_name':(['n'], obs_sel['long_name'].values),
                'altitude' :(['n'], obs_sel['altitude'].values)
            },
            coords= {
                'n': obs_sel['ns'].values,
                't':  np.arange(ntag),
#                                         'ref_time' : self.date.strftime('%Y-%m-%d')
            }
        )

#         samp_data['obs']= obs_sel['ch4']
#         samp_data['oerr'] = obs_sel['ch4_unc']
#         samp_data['lon'] = obs_sel['longitude']
#         samp_data['lat'] = obs_sel['latitude']
#         samp_data['time'] = self.date.strftime('%Y-%m-%d')



        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Matched noaa sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'noaa sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'noaa sample result of GEOS-Chem pertubation run'
        samp_data['hm'].attrs['units'] = '1e-9'
        
        if sampl_file:  
            samp_data.to_netcdf(sampl_file)

        return samp_data 
        

if(__name__=='__main__'):
    
    date_st = config.date_st
    date_end = config.date_end
    period = pd.date_range(start = date_st, end = date_end, freq= 'D')
    for day in period:
        print(day.strftime('%Y-%m-%d'))
        noaa_sampl = noaa_obs(date = day,\
                                info = config.noaa_info,\
                                sp_diagn_info = config.sp_diagn_info ,\
                                en_info_list = config.en_info_list,\
                                en_item_table = config.entry_table)
        
        sampl_path = config.noaa_info['sampl_path'].replace('YYYY', day.strftime('%Y'))
#         if not os.path.exists(sampl_path):
#             mkdir_cmd = 'mkdir ' + sampl_path
#             os.system(mkdir_cmd)
        
        sampl_file = sampl_path + config.noaa_info['sampl_file'].replace('YYYYMMDD', day.strftime('%Y%m%d'))
#         noaa_sampl.gen_hm(sampl_file)
        hm_file = sampl_file
        noaa_sampl.add_project_tag(hm_file, sampl_file)
   