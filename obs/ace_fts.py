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

class ace_obs:
    '''
    assimilation sampling for GOSAT observations
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
            gosat_info: gosat observations information
            sp_diagn_info: single diagnosis information
            en_info_list: ensemble tagrun file list(contains ensemble outputs information)
            en_item_table: ensemble step information at current day
        '''
        self.date = date
        self.gosat_info = info
        self.sp_info = sp_diagn_info
        self.en_info = en_info_list
        self.en_item = en_item_table
        
#     def sat_read(self, filename):
        
#         sat = xr.open_dataset(filename)
#         sat = sat.assign_coords(m=(sat.m), n=(sat.n))
#         ns = (sat['sensor_zenith_angle'] < 75) & (sat['xch4_quality_flag'] == 0)& (sat['xch4'] > 0)
#         sat_data = sat.sel(n=ns)

#         nl= np.all( sat_data['pressure_levels'] > 0, axis =1)
#         sat_data1 = sat_data.sel(n=nl)
        
#         return sat_data1
        
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

        ind_x = xr.DataArray(lon, dims='n')
        ind_y = xr.DataArray(lat, dims='n')
        model_ch4 = combine_en.isel(lon=ind_x, lat=ind_y)

        return model_ch4
    
#     def setup_pres_weight(self, pres):
    
#         n = len(pres)
#         p0 = np.full([n], 0.01, dtype=np.float)
#         gross_p = pres[:, 0] - p0
#         p_mins = np.concatenate((pres[:, 1:], p0[:, None]), axis=1)
#         delta_p = np.divide((pres - p_mins) , gross_p[:, None])
#         delta_p = np.array(delta_p)
#         return delta_p
    def setup_pres_weight(self, pres):

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
    
    #@nb.jit(nopython=True,parallel=True) 
    def vertical_intpl_2d(self, po_x, og_x, og_y):
        po_y = np.zeros_like(po_x)
        for i in range(len(og_x[:, 0])):
            #po_y[i, :] = np.interp(po_x[i, :], og_x[i, :],og_y[i, :])
            ff = interp1d(og_x[i, :], og_y[i, :], fill_value="extrapolate")
            po_y[i, :] = ff(po_x[i, :])

        return po_y
    
    def get_xgp(self, obs_pres, obs_ave, model_pres, model_data, pres, delta_p, method='linear'):
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


        ave_ker = self.vertical_intpl_2d(pres, obs_pres, obs_ave)
        md_interp = self.vertical_intpl_2d(pres, model_pres, model_data)

        xgp = (ave_ker* md_interp*delta_p).sum(axis=1)

        return xgp
    
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
#         cur_step = 'ST'+ str(step).zfill(3)
        print(step)
        for file in self.en_info:
            
            filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
            filename = file['filepath']+ str(yyyy) +'/'+filename
            if not os.path.exists(filename):
                cur_yyyy = str(int(yyyy) -1)
                filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
                filename = file['filepath']+ str(cur_yyyy) +'/'+filename
            print(filename)
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

        ## step2: read satellite observation data --- ACE-FTS
        
        sat_fn = xr.open_dataset(self.ace_info['filename'])
        #gosat_info['filepath'].replace('YYYY', str(yyyy)) +self.gosat_info['filename'].replace('YYYYMMDD', date)
        ll = (sat_fn['time'].dt.year == self.date.year) & \
        (sat_fn['time'].dt.month == self.date.month) & \
        (sat_fn['time'].dt.day == self.date.day)
        ace_daily = sat_fn.sel(index = ll)

        if ace_daily.index.size == 0:  ## no observations
            return False, False




        ## step3: read geos-chem single tracer(rerun) output and do resample
            # sampling girds 
        lons = ace_daily['lon'].values
        lats = ace_daily['lat'].values
        tropp_Lev = ace_daily['tropp_lev'].values

            # resampling bpch output
        filename = self.sp_info['filepath'] + self.sp_info['filename'].replace('YYYYMMDD', date)

        sp_ch4 = open_bpchdataset(filename,\
                                  diaginfo_file=self.sp_info['diaginfo_file'],\
                                  tracerinfo_file=self.sp_info['tracerinfo_file'])
        
        sp_mod = self.setup_nc_daily_feild(combine_en = sp_ch4,\
                                           lons = lons,\
                                           lats = lats) 
        
        ## step4: calculate pressure weight function on measurements
        gc_pres_dry = ace_daily['gc_pres_dry'].values

        ## step5: combine observation and model pressure layers and delta pres

        model_pres = sp_mod['PEDGE_S_PSURF'].values
        
        ####  get tropopause levels 
    tropp_lev = field_sampling_3d(lons = lons,\
                                lats = lats,\
                                times = times,\
                                data = gc_field['Met_TropLev'])
        sat_pres = sat_data['pressure_levels'].values
        pres = np.concatenate((model_pres, sat_pres), 1)
        pres = np.sort(pres, axis=1)[:, ::-1]
        delta_p = self.setup_pres_weight(pres)


        ## step6: caculation xgp result of single tracers[hm0]
        sgl_data = sp_mod[self.sp_info['var_prefix']].values 

        sp_xgp = self.get_xgp(sat_pres, sat_ave, model_pres, sgl_data, pres, delta_p, method='log')
        hm0 = xgp0 + sp_xgp
        hm0 = np.array(hm0)


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

                em_data = em_mod[var].values
                em_xgp = self.get_xgp(sat_pres, sat_ave, model_pres, em_data, pres, delta_p, method='log')
                tag_xgp.append(em_xgp + xgp0)

        hm = np.array(tag_xgp)
        hm= hm.transpose()


        ## step8: integration of sampling results
        ns, ntag = np.shape(hm)    
        samp_data = xr.Dataset({'hm0':(['n'], hm0),
                                'hm':(['n','t'], hm)},
                                coords= {'n': sat_data.n,
                                        't':  range(ntag),
                                        'ref_time' : self.date.strftime('%Y-%m-%d')})



        samp_data['obs']= sat_data['xch4']
        samp_data['oerr'] = sat_data['xch4_uncertainty']
        samp_data['lon'] = sat_data['longitude']
        samp_data['lat'] = sat_data['latitude']
        samp_data['time'] = sat_data['time']



        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Matched ACE-FTS sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'ACE-FTS sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'ACE-FTS sample result of GEOS-Chem pertubation run'
        samp_data['hm'].attrs['units'] = '1e-9'
        samp_data.to_netcdf(sampl_file)

        return samp_data 
    
    def gen_hm0(self, hm_file, sampl_file):
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
        
        sat_fn = self.gosat_info['filepath'].replace('YYYY', str(yyyy)) +self.gosat_info['filename'].replace('YYYYMMDD', date)

        if os.path.exists(sat_fn):
            #trying to open a file in read mode
            sat_data = self.sat_read(sat_fn)
        else:
            return False



        ### step3: get xgp0 
        sat_ave = sat_data['xch4_averaging_kernel'].values
        xgp0 = ((1 - sat_ave) * sat_data['pressure_weight'] *\
            sat_data['ch4_profile_apriori']).sum(dim='m')
        xgp0 = np.array(xgp0.values)


        ## step4: read geos-chem single tracer(rerun) output and do resample
            # sampling girds 
        lons = sat_data['longitude'].values
        lats = sat_data['latitude'].values

            # resampling bpch output
        filename = self.sp_info['filepath'] + self.sp_info['filename'].replace('YYYYMMDD', date)

        sp_ch4 = open_bpchdataset(filename,\
                                  diaginfo_file=self.sp_info['diaginfo_file'],\
                                  tracerinfo_file=self.sp_info['tracerinfo_file'])
        
        sp_mod = self.setup_nc_daily_feild(combine_en = sp_ch4,\
                                           lons = lons,\
                                           lats = lats) 


        ## step5: combine observation and model pressure layers and delta pres

        model_pres = sp_mod['PEDGE_S_PSURF'].values
        sat_pres = sat_data['pressure_levels'].values
        pres = np.concatenate((model_pres, sat_pres), 1)
        pres = np.sort(pres, axis=1)[:, ::-1]
        delta_p = self.setup_pres_weight(pres)


        ## step6: caculation xgp result of single tracers[hm0]
        sgl_data = sp_mod[self.sp_info['var_prefix']].values 

        sp_xgp = self.get_xgp(sat_pres, sat_ave, model_pres, sgl_data, pres, delta_p, method='log')
        hm0 = xgp0 + sp_xgp
        hm0 = np.array(hm0)


        ## step7: integration of sampling results
        hm = xr.open_dataset(hm_file)['hm']
        ns, ntag = np.shape(hm)
   
        samp_data = xr.Dataset(
            {
                'hm0':(['n'], hm0),
                'hm':(['n','t'], hm)},
            coords= {'n': sat_data.n,
                     't':  range(ntag),
                     'ref_time' : self.date.strftime('%Y-%m-%d')})



        samp_data['obs']= sat_data['xch4']
        samp_data['oerr'] = sat_data['xch4_uncertainty']
        samp_data['lon'] = sat_data['longitude']
        samp_data['lat'] = sat_data['latitude']
        samp_data['time'] = sat_data['time']



        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Matched gosat sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'gosat sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'gosat sample result of GEOS-Chem pertubation run'
        samp_data['hm'].attrs['units'] = '1e-9'
        if sampl_file:  
            samp_data.to_netcdf(sampl_file)

        return samp_data 



    
if (__name__=="__main__"):
    
    date_st = config.date_st
    date_end = config.date_end
    period = pd.date_range(start = date_st, end = date_end, freq= 'D')
    for day in period:
        print(day.strftime('%Y-%m-%d'))
        gosat_sampl = gosat_obs(date = day,\
                                info = config.gosat_info,\
                                sp_diagn_info = config.sp_diagn_info ,\
                                en_info_list = config.en_info_list,\
                                en_item_table = config.entry_table)
        
        sampl_path = config.gosat_info['sampl_path'].replace('YYYY', day.strftime('%Y'))
        if not os.path.exists(sampl_path):
            mkdir_cmd = 'mkdir ' + sampl_path
            os.system(mkdir_cmd)
        
        sampl_file = sampl_path + config.gosat_info['sampl_file'].replace('YYYYMMDD', day.strftime('%Y%m%d'))
        gosat_sampl.gen_hm(sampl_file)