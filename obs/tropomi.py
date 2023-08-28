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

class tropomi_obs:
    '''
    assimilation sampling for TROPOMI observations
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
            tropomi_info: tropomi observations information
            sp_diagn_info: single diagnosis information
            en_info_list: ensemble tagrun file list(contains ensemble outputs information)
            en_item_table: ensemble step information at current day
        '''
        self.date = date
        self.tropomi_info = info
        self.sp_info = sp_diagn_info
        self.en_info = en_info_list
        self.en_item = en_item_table
        
    def sat_read( self, filename):
        with xr.open_dataset(filename) as sat:
            sat.load()
        #sat = xr.open_dataset(filename)
        #sat = sat.assign_coords(lon=sat.lon, lat= sat.lat)
        sel_obs = (sat['solar_zenith_angle'] < 75) & (sat['qa_value'] > 0.005) & (sat['xch4_corrected'] > 0)\
        &(sat['xch4_precision'] > 0) #& (sat['surface_pressure'] > 0)
        sat_data = sat.sel(nobs=sel_obs)
#         sat_pres = np.linspace(sat_data['surface_pressure'], sat_data['surface_pressure']-(sat_data['dp']*12.0),13,\
#                                endpoint=True, dtype=np.float32, axis=1)
#         nl=np.all( sat_pres > 0, axis =1)
# #         print(sat_data, np.shape(nl))
#         sat_data1 = sat_data.sel(nobs = nl)
        return sat_data
        
        
    #@nb.jit(nopython=True,parallel=True) 
    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    
 
    
    #@nb.jit(nopython=True,parallel=True) 
    def vertical_intpl_2d(self, po_x, og_x, og_y):
        po_y = np.zeros_like(po_x)
        for i in range(len(og_x[:, 0])):
            #po_y[i, :] = np.interp(po_x[i, :], og_x[i, :],og_y[i, :])
            ff = interp1d(og_x[i, :], og_y[i, :], fill_value="extrapolate")
            po_y[i, :] = ff(po_x[i, :])

        return po_y
    
    
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
    
    def setup_pres_weight(self, pres):
    
        n = len(pres)
        p0 = np.full([n], 0.01, dtype=np.float)
        gross_p = pres[:, 0] - p0
        p_mins = np.concatenate((pres[:, 1:], p0[:, None]), axis=1)
        delta_p = np.divide((pres - p_mins) , gross_p[:, None])
        delta_p = np.array(delta_p)
        return delta_p
    
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
        em_step_list=list()
        em_name_list=list()

        """ check the ensemble member to available """
        for vals in table_entry_list:
            if (cur_time>=vals[0] and cur_time <vals[1]):

                em_step=vals[4]
                em_step_list.append(em_step)
                em_name=vals[5]
                em_name_list.append(em_name)
        return em_name_list
    
    def grid_flatten(self, model_lon, model_lat, dlat, dlon,  data):
        '''
        interpolate tropomi measurements to model grid mean

        Args:
            model_lon: longitude array in model grid
            model_lat: latitude array in model grid
            dlat     : delta latitude
            dlon     : delta longitude
            data     : tropomi L2 point retrievals
        Returns:
            sat_grid_data: tropomi grid flatten dataset 
        '''

        ### variables initialization
        n_lon = len(model_lon)
        n_lat = len(model_lat)
        n_lev = len(data['nlayer'])
        sat_lat = data['lat'].values
        sat_lon = data['lon'].values
        xch4_grid = np.zeros((n_lat,n_lon))+np.nan
        xch4_corrected_grid = np.zeros((n_lat,n_lon))+np.nan
        xch4_precision_grid = np.zeros((n_lat,n_lon))+np.nan
        ak_grid = np.zeros((n_lat,n_lon,n_lev))+np.nan
        qa_value = np.zeros((n_lat,n_lon))+np.nan
        ch4_ap_grid = np.zeros((n_lat,n_lon,n_lev))+np.nan

        psurf_grid = np.zeros((n_lat,n_lon))+np.nan
        dp_grid = np.zeros((n_lat,n_lon))+np.nan

        grid_lon = np.full((n_lat, n_lon), np.nan)
        grid_lat = np.full((n_lat, n_lon), np.nan)

        ### assign values to grid variables
        for ilat, lat in enumerate(model_lat):
#             if np.abs(lat) < 68.:
            for ilon, lon in enumerate(model_lon):

                pos = np.where(np.logical_and(np.abs(sat_lat - lat) <= dlat/2, np.abs(sat_lon - lon) <= dlon/2))[0]
                
                if len(pos)>0:
                    oerr = data['xch4_precision'].values[pos]
                    xch4_grid[ilat,ilon] = np.mean(data['xch4'].values[pos])
                    xch4_corrected_grid[ilat,ilon] = np.mean(data['xch4_corrected'].values[pos])
                    xch4_precision_grid[ilat,ilon] = np.mean(data['xch4_precision'].values[pos])
                    ch4_ap_grid[ilat,ilon,:] = np.mean(data['ch4_profile_apriori'].values[pos,:]/data['dry_air_subcolumns'].values[pos,:] *1e9,axis=0)
                    ak_grid[ilat,ilon,:] = np.mean(data['xch4_column_averaging_kernel'].values[pos,:],axis=0)
                    qa_value[ilat,ilon] = np.mean(data['qa_value'].values[pos])
                    dp_grid[ilat,ilon] = np.mean(data['dp'].values[pos])
                    psurf_grid[ilat,ilon] = np.mean(data['surface_pressure'].values[pos])
                    grid_lon[ilat,ilon] = np.mean(data['lon'].values[pos])
                    grid_lat[ilat,ilon] = np.mean(data['lat'].values[pos])
#                     print(np.max(oerr), np.min(oerr))
       
#                     xch4_grid[ilat,ilon] = np.sum(np.divide(data['xch4'].values[pos],oerr[:]))/np.sum(np.divide(1.0,oerr[:]))
#                     xch4_corrected_grid[ilat,ilon] = np.sum(data['xch4_corrected'].values[pos]/oerr[:])/np.sum(1.0/oerr[:])

#                     xch4_precision_grid[ilat,ilon] = np.sum(data['xch4_precision'].values[pos]/oerr[:])/np.sum(1.0/oerr[:])
#                     ch4_ap_grid[ilat,ilon,:] = np.sum((data['ch4_profile_apriori'].values[pos,:]/data['dry_air_subcolumns'].values[pos,:]*1e9)/oerr[:, None],axis=0)/np.sum(1.0/oerr[: ,None],axis=0)
#                     ak_grid[ilat,ilon,:] = np.sum(data['xch4_column_averaging_kernel'].values[pos,:]/oerr[:, None],axis=0)/np.sum(1.0/oerr[: ,None],axis=0)
#                     qa_value[ilat,ilon] = np.sum(data['qa_value'].values[pos]/oerr[:])/np.sum(1.0/oerr[:])

#                     psurf_grid[ilat,ilon] = np.sum(data['surface_pressure'].values[pos]/oerr[:])/np.sum(1.0/oerr[:])
#                     dp_grid[ilat,ilon] = np.sum(data['dp'].values[pos]/oerr[:])/np.sum(1.0/oerr[:])

#                     grid_lon[ilat,ilon] = np.sum(data['lon'].values[pos]/oerr[:])/np.sum(1.0/oerr[:])
#                     grid_lat[ilat,ilon] = np.sum(data['lat'].values[pos]/oerr[:])/np.sum(1.0/oerr[:])

        ### create a new dataset
        field_nc = xr.Dataset(
            {'xch4': (['lat', 'lon'], xch4_grid),
             'xch4_corrected': (['lat', 'lon'], xch4_corrected_grid),
             'xch4_precision':(['lat', 'lon'], xch4_precision_grid),
             'ch4_profile_apriori':([ 'lat', 'lon', 'lev'], ch4_ap_grid),
             'xch4_column_averaging_kernel': ([ 'lat', 'lon', 'lev'], ak_grid),
             'qa_value':(['lat', 'lon'], qa_value),

             'surface_pressure':(['lat', 'lon'],psurf_grid),
             'dp':(['lat', 'lon'],dp_grid),
             'grid_lon':(['lat', 'lon'],grid_lon),
             'grid_lat':(['lat', 'lon'],grid_lat)

            },
            coords={
                'lon': model_lon,
                'lat': model_lat,
                'lev': data['nlayer'].values
            }  
        )

        #### flatten 2-d field
        sat_flatten = field_nc.stack(n=['lat', 'lon']).transpose('n', 'lev')

        sat_flatten = sat_flatten.assign_coords({'n': range(len(sat_flatten['n']))})

        sat_flatten = sat_flatten.where(~np.isnan(sat_flatten['xch4']), drop= True)

        sat_flatten = sat_flatten.assign_coords({'n': range(len(sat_flatten['n']))})


        return sat_flatten
    
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
        mm = self.date.strftime('%m')
        ## step2: read satellite observation data
        
        sat_fn = self.tropomi_info['filepath'].replace('YYYY', str(yyyy)).replace('MM',str(mm))\
        +self.tropomi_info['filename'].replace('YYYYMMDD', date)

        if os.path.exists(sat_fn):
            #trying to open a file in read mode
            sat = self.sat_read(sat_fn)
#             sat = xr.open_dataset(sat_fn).load()
#             sat.close()
        else:
            return False
        

        grid = xr.open_dataset(self.tropomi_info['grid_file']).load()
        grid.close()
        
        sat_data = self.grid_flatten(model_lon = grid['lon'].values,\
                                model_lat = grid['lat'].values,\
                                dlat = self.tropomi_info['dlat'],\
                                dlon = self.tropomi_info['dlon'],\
                                data = sat)
 
        if sat_data.n.size ==0:
            return False

        ### step3: get xgp0  #
        
        
        sat_ave = sat_data['xch4_column_averaging_kernel'].values
        xgp0 = ((1 - sat_data['xch4_column_averaging_kernel']) * sat_data['ch4_profile_apriori']*(1./12.)).sum(dim ='lev')
#         xgp0 = ((1 -  sat_ave) *sat_data['ch4_profile_apriori']*(1./12.)).sum(dim ='lev')
        xgp0 = np.array(xgp0.values)
        print(xgp0)


        ## step4: read geos-chem single tracer(rerun) output and do resample
            # sampling girds 
        lons = sat_data['grid_lon'].values
        lats = sat_data['grid_lat'].values

            # resampling bpch output

        filename = self.sp_info['filepath'] + self.sp_info['filename'].replace('YYYYMMDD', date)

           
        sp_ch4 = open_bpchdataset(filename,\
                                  diaginfo_file= self.sp_info['diaginfo_file'],\
                                  tracerinfo_file=self.sp_info['tracerinfo_file'])
        
        sp_mod = self.setup_nc_daily_feild(combine_en = sp_ch4,\
                                           lons = lons,\
                                           lats = lats) 


        ## step5: combine observation and model pressure layers and delta pres

        model_pres = sp_mod['PEDGE_S_PSURF'].values # hPa
        surf_pres = sat_data['surface_pressure'].values/100.0 # Pa-->hPa 
        dp = sat_data['dp'].values/100.0 # Pa-->hPa 
        level_pres = np.linspace(surf_pres, surf_pres-(dp*12.0),13,\
                               endpoint=True, dtype=np.float32, axis=1)
        sat_pres = (level_pres[:, :-1]+ level_pres[:,1:])/2.0 ## used for interpolation of average kernel and aprior profile
        pres = np.concatenate((model_pres, level_pres), 1)
        pres = np.sort(pres, axis=1)[:, ::-1]
        mid_pres = (pres[:,:-1]+ pres[:,1:])/2.0
        delta_p = (pres[:,:-1]-pres[:,1:])/(pres[:,0] - pres[:,-1] )[:, None]
#         delta_p = self.setup_pres_weight_1(pres)


        ## step6: caculation xgp result of single tracers[hm0]
        sgl_data = sp_mod[self.sp_info['var_prefix']].values 
        sp_xgp = self.get_xgp(sat_pres, sat_ave, model_pres, sgl_data, mid_pres, delta_p, method='linear')
        hm0 = xgp0 + sp_xgp
        hm0 = np.array(hm0)
        
#         albedo = sat_data['surface_albedo'].values
#         c1 = 1.0173
#         c2 = -0.1538
#         c3 = 0.2036
#         hm0_cor = hm0*(c1 +c2*albedo + c3*albedo**2)
#         print(hm0)
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
                em_xgp = self.get_xgp(sat_pres, sat_ave, model_pres, em_data, mid_pres, delta_p, method='linear')
                tag_xgp.append(em_xgp + xgp0)

        hm = np.array(tag_xgp)
        hm= hm.transpose()


        ## step8: integration of sampling results
        ns, ntag = np.shape(hm)    
        samp_data = xr.Dataset({'hm0':(['n'], hm0),
                                'hm':(['n','t'], hm),
                                'obs':(['n'],sat_data['xch4_corrected'].values),
                                'obs_raw':(['n'],sat_data['xch4'].values),
                                'oerr':(['n'],sat_data['xch4_precision'].values),
                                'lon':(['n'],lons),
                                'lat':(['n'],lats),
                                'project':(['n'], list(['tropomi']* ns)),
                               },
                               
                                coords= {'n': sat_data.n,
                                        't':  range(ntag),
                                        'ref_time' : self.date.strftime('%Y-%m-%d')})






        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Matched tropomi sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'Tropomi sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'Tropomi sample result of GEOS-Chem pertubation run'
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
        mm = self.date.strftime('%m')
        ## step2: read satellite observation data
        
        sat_fn = self.tropomi_info['filepath'].replace('YYYY', str(yyyy)).replace('MM',str(mm))\
        +self.tropomi_info['filename'].replace('YYYYMMDD', date)

        if os.path.exists(sat_fn):
            #trying to open a file in read mode
            sat = self.sat_read(sat_fn)
#             sat = xr.open_dataset(sat_fn).load()
#             sat.close()
        else:
            return False
        

        grid = xr.open_dataset(self.tropomi_info['grid_file']).load()
        grid.close()
        
        sat_data = self.grid_flatten(model_lon = grid['lon'].values,\
                                model_lat = grid['lat'].values,\
                                dlat = self.tropomi_info['dlat'],\
                                dlon = self.tropomi_info['dlon'],\
                                data = sat)
        
        if sat_data.n.size ==0:
            return False

        ### step3: get xgp0  #
        
        
        sat_ave = sat_data['xch4_column_averaging_kernel'].values
        xgp0 = ((1 - sat_data['xch4_column_averaging_kernel']) * sat_data['ch4_profile_apriori']*(1./12.)).sum(dim ='lev')
#         xgp0 = ((1 -  sat_ave) *sat_data['ch4_profile_apriori']*(1./12.)).sum(dim ='lev')
        xgp0 = np.array(xgp0.values)
        print(xgp0)


        ## step4: read geos-chem single tracer(rerun) output and do resample
            # sampling girds 
        lons = sat_data['grid_lon'].values
        lats = sat_data['grid_lat'].values

            # resampling bpch output

#         filename = self.sp_info['filepath'] + self.sp_info['filename'].replace('YYYYMMDD', date)
        filename = self.sp_info['filepath'] + self.sp_info['filename'].replace('YYYYMMDD', date).replace('ST000', 'ST'+ str(cur_step).zfill(3))
           
        sp_ch4 = open_bpchdataset(filename,\
                                  diaginfo_file= self.sp_info['diaginfo_file'],\
                                  tracerinfo_file=self.sp_info['tracerinfo_file'])
        
        sp_mod = self.setup_nc_daily_feild(combine_en = sp_ch4,\
                                           lons = lons,\
                                           lats = lats) 


        ## step5: combine observation and model pressure layers and delta pres

        model_pres = sp_mod['PEDGE_S_PSURF'].values # hPa
        surf_pres = sat_data['surface_pressure'].values/100.0 # Pa-->hPa 
        dp = sat_data['dp'].values/100.0 # Pa-->hPa 
        level_pres = np.linspace(surf_pres, surf_pres-(dp*12.0),13,\
                               endpoint=True, dtype=np.float32, axis=1)
        sat_pres = (level_pres[:, :-1]+ level_pres[:,1:])/2.0 ## used for interpolation of average kernel and aprior profile
        pres = np.concatenate((model_pres, level_pres), 1)
        pres = np.sort(pres, axis=1)[:, ::-1]
        mid_pres = (pres[:,:-1]+ pres[:,1:])/2.0
        delta_p = (pres[:,:-1]-pres[:,1:])/(pres[:,0] - pres[:,-1] )[:, None]
#         delta_p = self.setup_pres_weight_1(pres)


        ## step6: caculation xgp result of single tracers[hm0]
        sgl_data = sp_mod[self.sp_info['var_prefix']].values 
        sp_xgp = self.get_xgp(sat_pres, sat_ave, model_pres, sgl_data, mid_pres, delta_p, method='linear')
        hm0 = xgp0 + sp_xgp
        hm0 = np.array(hm0)
        
#         albedo = sat_data['surface_albedo'].values
#         c1 = 1.0173
#         c2 = -0.1538
#         c3 = 0.2036
#         hm0_cor = hm0*(c1 +c2*albedo + c3*albedo**2)
#         print(hm0)
        ## step7: caculation xgp result of ensemble tracers[hm]
        hm = xr.open_dataset(hm_file)['hm']
        ns, ntag = np.shape(hm)
        
        samp_data = xr.Dataset({'hm0':(['n'], hm0),
                                'hm':(['n','t'], hm),
                                'obs':(['n'],sat_data['xch4_corrected'].values),
                                'obs_raw':(['n'],sat_data['xch4'].values),
                                'oerr':(['n'],sat_data['xch4_precision'].values),
                                'lon':(['n'],lons),
                                'lat':(['n'],lats),
                                'project':(['n'], list(['tropomi']* ns)),
                               },
                               
                                coords= {'n': sat_data.n,
                                        't':  range(ntag),
                                        'ref_time' : self.date.strftime('%Y-%m-%d')})






        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Matched tropomi sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'Tropomi sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'Tropomi sample result of GEOS-Chem pertubation run'
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
        tropomi_sampl = tropomi_obs(date = day,\
                                tropomi_info = config.tropomi_info,\
                                sp_diagn_info = config.sp_diagn_info ,\
                                en_info_list = config.en_info_list,\
                                en_item_table = config.entry_table)
        
        sampl_path = config.tropomi_info['sampl_path'].replace('YYYY', day.strftime('%Y'))
        if not os.path.exists(sampl_path):
            mkdir_cmd = 'mkdir ' + sampl_path
            os.system(mkdir_cmd)
        
        sampl_file = sampl_path + config.tropomi_info['sampl_file'].replace('YYYYMMDD', day.strftime('%Y%m%d'))
        
        tropomi_sampl.gen_hm(sampl_file)