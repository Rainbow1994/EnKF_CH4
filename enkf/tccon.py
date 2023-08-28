#! /home/shzhu/apps/miniconda3/envs/idp/bin/python

import numpy as np
import numba as nb
import pandas as pd
import xarray as xr 
from xbpch import open_bpchdataset
from scipy import interpolate
from scipy.interpolate import interp1d
import time_module as tm
import rerun_geos as rg
import geos_chem_def as gcdf 
import os

class tccon_obs:
    '''
    assimilation & validation procedure using TCCON dataset
    '''
    def __init__(self, sel_id):
        '''
        sel_id: sites id arrat selected for assimilations or validation
        '''
        self.sel_id = sel_id

    
    
    def obs_read(self,\
                 tccon_obs,\
                 tccon_ids
                ):
        '''
        read tccon daily observations for assimilation

        Args:
            tccon_obs: daily mean tccon observation results

        Return:
            obs_sel: select valid site column concentration
            site_id: a list contains ids of the valid site name 
        '''
        
        obs = xr.open_dataset(tccon_obs).sel(ns = tccon_ids)
        #obs_sel = obs.where(obs['xch4_ppm']!= -999.99, drop=True)
        sel_id =obs['ns'].where(obs['xch4_ppm']!= -999.99).dropna(dim = 'ns')
        obs_sel = obs.sel(ns = sel_id)
        #print(obs_sel)
        site_id = obs_sel['ns'].values

        return obs_sel, site_id

    #@nb.jit(nopython=True,parallel=True) 
    def find_nearest(self,\
                     array,\
                     value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def mod_xbpch_read(self,\
                       mod_file,\
                       diaginfo_file,\
                       tracerinfo_file,\
                       lons,\
                       lats):
        '''
        read model daily xbpch output
        '''

        mod = open_bpchdataset(mod_file,\
                               diaginfo_file= diaginfo_file,\
                               tracerinfo_file=tracerinfo_file)

        lon =[]
        lat =[]
        for i_lon, i_lat in zip(lons, lats):
            lon.append(self.find_nearest(mod.lon.values, i_lon))
            lat.append(self.find_nearest(mod.lat.values, i_lat))

        ind_x = xr.DataArray(lon, dims='n')
        ind_y = xr.DataArray(lat, dims='n')
        mod_sel = mod.isel(lon=ind_x, lat=ind_y)



        return mod_sel




    def mod_nc_read(self,\
                    mod_file,\
                    time,\
                    lons,\
                    lats):
        '''
        read model netcdf output and interpolate to observation sites
        '''

        mod = xr.open_dataset(mod_file)

        mod_sel = mod.sel(time = time).interp(lon=lons, lat=lats)

        for var in mod_sel.data_vars:
            mod_sel[var].values = pd.DataFrame(mod_sel[var].values).fillna(method='ffill',\
                                                          axis=0).fillna(method='bfill',\
                                                                         axis=0).fillna(method='ffill',\
                                                                                        axis=1).fillna(method='bfill',\
                                                                                                       axis=1).values

        return mod_sel


    #@nb.jit(nopython=True,parallel=True) 
    def setup_comb_pres(self,\
                        obs_pres,\
                        mod_pres):

        '''
        combine model and obs profile pressure levels and 
        Args:
            obs_pres: observation pressure levels
            mod_pres: model pressure levels

        return:
            comb_pres: combined result

        '''

        pres = np.concatenate((obs_pres, mod_pres), 1)
        comb_pres = np.sort(pres, axis=1)[:, ::-1]

        return comb_pres

    #@nb.jit(nopython=True,parallel=True) 
    def setup_pres_weight(self,\
                          pres,\
                          method):

        '''
        calculate pressure weight of each level
        Args:
            pres: pressure levels
            method: 'linear' or 'log'

        return:

            pres_wgt: the pressure weight of the level(delta(P_k)/ P_surf)
        '''
        if (method == 'linear'):
            pres = pres

        elif (method == 'log'):
            if np.all(pres > 0):
                pres = np.log10(pres)
            else:
                raise ValueError("PRESSURE VALUES ERROR")
        else:
            raise AssertionError('ERROR WITH INTERPOLATION METHOD')

        n = len(pres)
        p0 = np.full([n], 0.01, dtype=np.float)
        gross_p = pres[:, 0] - p0
        p_mins = np.concatenate((pres[:, 1:], p0[:, None]), axis=1)
        delta_p = (pres - p_mins) / gross_p[:, None]
        pres_wgt = np.array(delta_p)

        return pres_wgt


    #@nb.jit(nopython=True,parallel=True) 
    def vertical_intpl_2d(self, po_x, og_x, og_y):
        '''
        2-d array interpolation

        '''
        po_y = np.zeros_like(po_x)
        for i in range(len(og_x[:, 0])):
            ff = interp1d(og_x[i, :], og_y[i, :], fill_value="extrapolate")
            po_y[i, :] = ff(po_x[i, :])
        return po_y


    #@nb.jit(nopython=True,parallel=True) 
    def prof_intpl_2d(self, in_pres, out_pres, var, method):
        '''
        2-d interpolation using pressure using different method
        Args:
            in_pres: orginal pressure levels
            out_pres: target pressure levels
            var: interpolation variables
            method: 'linear' or 'log'
        Return:
            intpl_var
        '''
        if (method == 'linear'):
            pres = in_pres

        elif (method == 'log'):
            if np.all(in_pres > 0):
                pres = np.log(in_pres)
                out_pres = np.log(out_pres)
            else:
                raise ValueError("PRESSURE VALUES ERROR")
        else:
            raise AssertionError('ERROR WITH INTERPOLATION METHOD')

        intpl_var = self.vertical_intpl_2d(out_pres, pres, var)

        return intpl_var

    #@nb.jit(nopython=True,parallel=True) 
    def prof_intpl_1d(self, in_pres, out_pres, var, method):
        '''
        1-d interpolation using pressure using different method
        Args:
            in_pres: orginal pressure levels
            out_pres: target pressure levels
            var: interpolation variables
            method: 'linear' or 'log'
        Return:
            intpl_var
        '''
        if (method == 'linear'):
            pres = in_pres

        elif (method == 'log'):
            if np.all(in_pres > 0):
                pres = np.log(in_pres)
                out_pres = np.log(out_pres)
            else:
                raise ValueError("PRESSURE VALUES ERROR")
        else:
            raise AssertionError('ERROR WITH INTERPOLATION METHOD')

        ff = interp1d(pres, var, fill_value="extrapolate")
        intpl_var = ff(out_pres)

        return intpl_var  



    #@nb.jit(nopython=True,parallel=True) 
    def get_xgp0(self,\
                 ave_ker,\
                 pres_wgt,\
                 apri_prof):
        '''
        get daily xgp0

        method:
            Xgas = sum(apri_prof * pres_wgt) + (mod_prof - apri_prof)* ave_ker * pres_wgt)
                 = sum((1 - ave_ker) * apri_prof * pres_wgt) + sum(mod_prof * ave_ker * pres_wgt)
                   -----------------------------------------   ----------------------------------
                 =                    xgp0                                   xgp

        Notes: 
            No need interpolation, just use observation information
        Args:
            ave_ker: average kernel
            pres_wgt: pressure weight function
            apri_prof: a prior file for Xgas calculation

        Return:
            xgp0 : seen in above method

        '''

        xgp0 = ((1 - ave_ker) * apri_prof * pres_wgt).sum(axis = 1)

        return xgp0

    #@nb.jit(nopython=True,parallel=True) 
    def get_xgp(self,\
                mod_prof,\
                ave_ker,\
                pres_wgt):
        '''
        get daily xgp

        method:
            Xgas = sum(apri_prof * pres_wgt) + (mod_prof - apri_prof)* ave_ker * pres_wgt)
                 = sum((1 - ave_ker) * apri_prof * pres_wgt) + sum(mod_prof * ave_ker * pres_wgt)
                   -----------------------------------------   ----------------------------------
                 =                    xgp0                                   xgp

        Notes: 
            all input variable should after interpolation according new combined pressure levels
        Args:
            ave_ker: average kernel
            pres_wgt: pressure weight function
            apri_prof: a prior file for Xgas calculation

        Return:
            xgp : seen in above method

        '''

        xgp = (mod_prof * ave_ker * pres_wgt).sum(axis=1)

        return xgp






    def tccon_hm(self,\
                 cur_step,\
                 doy,\
                 yyyy,\
                 sampl_file):
        '''
        read single and ensemble geos_chem run diagnostic ouput
        resample them to TCCON observation sites
        according to layer pressure and its weight function to caculate XCH4

        notes:
            read no seperated ensemble tagrun

        Args:
            cur_step: current step
            doy: days of current year
            yyyy: current year
            sampl_file: output filename

        Return:
            samp_data: sampling result
        '''
        method = 'log'

        ## step1: get time information
        tst= tm.doy_to_utc(doy, 0, yyyy)
        time = tst[:10]
        tst=tst.replace('-', '')
        tst=tst.replace(':', '') 


        ## step2: read tccon observation data
        tccon_fln = 'tccon.' + time +'.nc'
        tccon_fln = gcdf.tccon_path + '/' + str(yyyy) + '/' + tccon_fln
        if (os.path.exists(tccon_fln) and self.sel_id.size != 0):
            #trying to open a file in read mode
            tccon_obs, site_id = self.obs_read(tccon_obs = tccon_fln,\
                                               tccon_ids = self.sel_id)
        else:
            return False
        lons = tccon_obs['longitude'].values
        lats = tccon_obs['latitude'].values
        
        ### step3: setup single diagnosis output filename
        enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step, 1, 2)
        sp_fn = 'ts_satellite.'+enaf+'.'+tst[0:8]+'.bpch'
        sp_fn = gcdf.sp_diagn_path+'ND51/'+sp_fn
        sp_diag = gcdf.sp_diagn_path +'diaginfo_sp.dat'
        sp_tra = gcdf.sp_diagn_path +'tracerinfo_sp.dat'

        sp_mod = self.mod_xbpch_read(mod_file = sp_fn,\
                                diaginfo_file = sp_diag ,\
                                tracerinfo_file = sp_tra,\
                                lons = lons,\
                                lats = lats
                               )

        ### step4: get xgp0 
        mod_pres = sp_mod['PEDGE_S_PSURF'].values
        n_obs = len(mod_pres)
        obs_nlev = len(tccon_obs['ak_P_hPa'])
        obs_pres = np.zeros((n_obs, obs_nlev))
        for i in range(n_obs):
            obs_pres[i,:] = tccon_obs['ak_P_hPa'].values
        pres_wgt = self.setup_pres_weight(obs_pres, method = method)

        xgp0 = self.get_xgp0(ave_ker = tccon_obs['ak_ch4'].values,\
                        pres_wgt = pres_wgt[0,:],\
                        apri_prof = tccon_obs['prior_ch4'].values)

        
        
        ## step5: calculate intergration parameters(press_weight, average_kernel, aprior_profiles)
        ## combine observation and model pressure levels
        

        comb_pres = self.setup_comb_pres(obs_pres = obs_pres,\
                                    mod_pres = mod_pres)

        ## calculate pressure weight factor& column dry air factor
        pres_wgt = self.setup_pres_weight(comb_pres, method = method)

        ## interpolation apply for average kernel, a-prior profile and model profiles 
        obs_aprior = tccon_obs['prior_ch4'].values
        obs_ave_ker = tccon_obs['ak_ch4'].values
        mod_prof = sp_mod['IJ_AVG_S_CH4'].values
        mod_prof = self.prof_intpl_2d(mod_pres, comb_pres, mod_prof, method = method)
        ave_ker = self.prof_intpl_2d(obs_pres, comb_pres, obs_ave_ker, method = method)

        ## step6: caculation xgp result of single tracers[hm0]
        sp_xgp = self.get_xgp(mod_prof = mod_prof,\
                         ave_ker = ave_ker,\
                         pres_wgt = pres_wgt)

        hm0 = np.array(xgp0 + sp_xgp)


        ## step7: caculation xgp result of ensemble tracers[hm]
        em_step_list,em_name_list  = rg.get_ensemble_list(yyyy, doy, \
                                                        entry_table_name= gcdf.entry_table)
        tag_xgp =[]
        var_list= [''.join(['IJ_AVG_S_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
        for em_step, em_name in zip(em_step_list,em_name_list):

            em_fn_0 ='ts_satellite.'+ em_name +'.'+tst[0:8]+'.bpch'
            em_fn = gcdf.en_diagn_path+ '/' + str(yyyy) + '/'+em_fn_0
            if not os.path.exists(em_fn):
                em_fn = gcdf.en_diagn_path+ '/' + str(yyyy - 1) + '/'+em_fn_0
            em_diag = gcdf.em_diag
            em_tra = gcdf.em_tra
            em_mod = self.mod_xbpch_read(mod_file = em_fn,\
                                    diaginfo_file = em_diag ,\
                                    tracerinfo_file = em_tra,\
                                    lons = lons,\
                                    lats = lats
                                   )

            for var in var_list:

                mod_prof = em_mod[var].values
                mod_prof = self.prof_intpl_2d(mod_pres, comb_pres, mod_prof, method = method)
                em_xgp =  self.get_xgp(mod_prof = mod_prof,\
                                  ave_ker = ave_ker,\
                                  pres_wgt = pres_wgt)
                tag_xgp.append(em_xgp + xgp0)

        hm = np.array(tag_xgp)
        hm= hm.transpose()

        ## step6: integration of sampling results
        ns, ntag = np.shape(hm)    
        samp_data = xr.Dataset({'hm0':(['ns'], hm0),
                                'hm':(['ns','t'], hm)},
                                coords= {'ns': tccon_obs['ns'].values,
                                        't':  range(ntag)})
                                        #'ref_time' : tm.doy_to_utc(doy, 0, yyyy)



        samp_data['obs']= tccon_obs['xch4_ppm'] * 1000
        samp_data['oerr'] = tccon_obs['xch4_ppm_error'] * 1000
        samp_data['lon'] = tccon_obs['longitude']
        samp_data['lat'] = tccon_obs['latitude']
        samp_data['time'] = tccon_obs['date']



        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Tccon sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'sample result of GEOS-Chem pertubation run'
        samp_data['hm'].attrs['units'] = '1e-9'
        #samp_data.to_netcdf(sampl_file)

        return samp_data 
    
    def combine_em(self, en1, tra_num1, en2, tra_num2):
    
        en3 = xr.Dataset()
        tracers = ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA']
        for i_tra in tracers:
            for i in range(tra_num2):
                new_name = i_tra + str(i + 1 + tra_num1)
                old_name = i_tra + str(i + 1)

                #en2 = en2.rename_vars({old_name:new_name})#[new_name].values = en1[old_name].values
                en3[new_name] = en2[old_name]
    #    en1 = en1.drop(['IJ_AVG_S_', 'time_bnds'])
    #    en2 = en2.drop(['IJ_AVG_S_CH4', 'PEDGE_S_PSURF'])
        en = xr.merge([en1,en3])
        en.attrs = en1.attrs
        en.astype(np.float32)

        return en
    
    def tccon_hm2(self,\
                 cur_step,\
                 doy,\
                 yyyy,\
                 sampl_file):
        '''
        read single and ensemble geos_chem run diagnostic ouput
        resample them to TCCON observation sites
        according to layer pressure and its weight function to caculate XCH4

        notes:
             read seperated ensemble tagrun

        Args:
            cur_step: current step
            doy: days of current year
            yyyy: current year
            sampl_file: output filename

        Return:
            samp_data: sampling result
        '''
        method = 'linear'

        ## step1: get time information
        tst= tm.doy_to_utc(doy, 0, yyyy)
        time = tst[:10]
        tst=tst.replace('-', '')
        tst=tst.replace(':', '') 


        ## step2: read tccon observation data
        tccon_fln = 'tccon.' + time +'.nc'
        tccon_fln = gcdf.tccon_path + '/' + str(yyyy) + '/' + tccon_fln
        if (os.path.exists(tccon_fln) and self.sel_id.size != 0):
            #trying to open a file in read mode
            tccon_obs, site_id = self.obs_read(tccon_obs = tccon_fln,\
                                               tccon_ids = self.sel_id)
        else:
            return False
        lons = tccon_obs['longitude'].values
        lats = tccon_obs['latitude'].values
        
        ### step3: setup single diagnosis output filename
        enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step, 1, 2)
        sp_fn = 'ts_satellite.'+enaf+'.'+tst[0:8]+'.bpch'
        sp_fn = gcdf.sp_diagn_path+'ND51/'+sp_fn
        sp_diag = gcdf.sp_diagn_path +'diaginfo_sp.dat'
        sp_tra = gcdf.sp_diagn_path +'tracerinfo_sp.dat'

        sp_mod = self.mod_xbpch_read(mod_file = sp_fn,\
                                diaginfo_file = sp_diag ,\
                                tracerinfo_file = sp_tra,\
                                lons = lons,\
                                lats = lats
                               )

        ### step4: get xgp0 
        mod_pres = sp_mod['PEDGE_S_PSURF'].values
        n_obs = len(mod_pres)
        obs_nlev = len(tccon_obs['ak_P_hPa'])
        obs_pres = np.zeros((n_obs, obs_nlev))
        for i in range(n_obs):
            obs_pres[i,:] = tccon_obs['ak_P_hPa'].values
        pres_wgt = self.setup_pres_weight(obs_pres, method = method)
        #print(pres_wgt[0,:])
        xgp0 = self.get_xgp0(ave_ker = tccon_obs['ak_ch4'].values,\
                        pres_wgt = pres_wgt[0,:],\
                        apri_prof = tccon_obs['prior_ch4'].values)

        

        ## step5: calculate intergration parameters(press_weight, average_kernel, aprior_profiles)
        ## combine observation and model pressure levels
        

        comb_pres = self.setup_comb_pres(obs_pres = obs_pres,\
                                    mod_pres = mod_pres)

        ## calculate pressure weight factor& column dry air factor
        pres_wgt = self.setup_pres_weight(comb_pres, method = method)

        ## interpolation apply for average kernel, a-prior profile and model profiles 
        obs_aprior = tccon_obs['prior_ch4'].values
        obs_ave_ker = tccon_obs['ak_ch4'].values
        mod_prof = sp_mod['IJ_AVG_S_CH4'].values
        mod_prof = self.prof_intpl_2d(mod_pres, comb_pres, mod_prof, method = method)
        ave_ker = self.prof_intpl_2d(obs_pres, comb_pres, obs_ave_ker, method = method)
       

        ## step6: caculation xgp result of single tracers[hm0]
        sp_xgp = self.get_xgp(mod_prof = mod_prof,\
                         ave_ker = ave_ker,\
                         pres_wgt = pres_wgt)
        #print(xgp0)
        hm0 = np.array(xgp0 + sp_xgp)


        ## step7: caculation xgp result of ensemble tracers[hm]
        em_step_list1,em_name_list1  = rg.get_ensemble_list(yyyy, doy, \
                                                    entry_table_name= gcdf.entry_table1)
        em_step_list2,em_name_list2  = rg.get_ensemble_list(yyyy, doy, \
                                                    entry_table_name= gcdf.entry_table2)
        tra_num1 = gcdf.tra_num1
        tra_num2 = gcdf.tra_num2
        
        tag_xgp =[]
        var_list= [''.join(['IJ_AVG_S_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
        for em_step, em_name1, em_name2 in zip(em_step_list1,em_name_list1, em_name_list2):
        
            em_fn ='ts_satellite.'+ em_name1 +'.'+tst[0:8]+'.bpch'
            em_fn1 = gcdf.en_diagn_path+ '/' + str(yyyy) + '/'+em_fn
            if not os.path.exists(em_fn1):
                em_fn1 = gcdf.en_diagn_path+ '/' + str(yyyy - 1) + '/'+em_fn

            em_fn ='ts_satellite.'+ em_name2 +'.'+tst[0:8]+'.bpch'
            em_fn2 = gcdf.en_diagn_path+ '/' + str(yyyy) + '/'+em_fn

            if not os.path.exists(em_fn2):
                em_fn2 = gcdf.en_diagn_path+ '/' + str(yyyy - 1) + '/'+em_fn
            
            en1 = self.mod_xbpch_read(mod_file = em_fn1,\
                                      diaginfo_file = gcdf.em_diag,\
                                      tracerinfo_file = gcdf.em_tra1,\
                                      lons = lons,\
                                      lats = lats
                                     )
            
            en2 = self.mod_xbpch_read(mod_file = em_fn2,\
                                    diaginfo_file = gcdf.em_diag,\
                                    tracerinfo_file = gcdf.em_tra2,\
                                    lons = lons,\
                                    lats = lats
                                   )
            
            em_mod = self.combine_em(en1, tra_num1, en2, tra_num2)

            for var in var_list:

                mod_prof = em_mod[var].values
                mod_prof = self.prof_intpl_2d(mod_pres, comb_pres, mod_prof, method = method)
                em_xgp =  self.get_xgp(mod_prof = mod_prof,\
                                  ave_ker = ave_ker,\
                                  pres_wgt = pres_wgt)
                tag_xgp.append(em_xgp + xgp0)

        hm = np.array(tag_xgp)
        hm= hm.transpose()

        ## step6: integration of sampling results
        ns, ntag = np.shape(hm)    
        samp_data = xr.Dataset({'hm0':(['ns'], hm0),
                                'hm':(['ns','t'], hm)},
                                coords= {'ns': tccon_obs['ns'].values,
                                        't':  range(ntag)})
                                        #'ref_time' : tm.doy_to_utc(doy, 0, yyyy)



        samp_data['obs']= tccon_obs['xch4_ppm'] * 1000
        samp_data['oerr'] = tccon_obs['xch4_ppm_error'] * 1000
        samp_data['lon'] = tccon_obs['longitude']
        samp_data['lat'] = tccon_obs['latitude']
        samp_data['time'] = tm.doy_to_utc(doy, 0, yyyy)#tccon_obs['date']



        ## step9: add some description to variables    
        samp_data.attrs['Description'] = 'Tccon sample result without data quality control '

        samp_data['hm0'].attrs['long_name'] = 'model_xch4'
        samp_data['hm0'].attrs['standard_name'] = 'sample result of GEOS-Chem single run'
        samp_data['hm0'].attrs['units'] = '1e-9'
        samp_data['hm'].attrs['long_name'] = 'pertubation_xch4'
        samp_data['hm'].attrs['standard_name'] = 'sample result of GEOS-Chem pertubation run'
        samp_data['hm'].attrs['units'] = '1e-9'
        samp_data.to_netcdf(sampl_file)

        return samp_data

    # def tccon_diag(self,\
    #                doy,\
    #                yyyy,\
    #                diag_output):
    
    
if (__name__=="__main__"):
    
    import tccon_sites as tc
    import geos_chem_def as gcdf
    import rerun_geos as rg
    import calendar

    assim_tc = tccon_obs(sel_id = tc.assim_id)
  #  valid_tc = tccon_obs(sel_id = tc.valid_id)
    doy = 0
    yyyy = 2016
    #cur_step = 0
    for cur_step in range(12, 24):
        mon_len = calendar.monthrange(yyyy,cur_step-11)[1]
        for idoy in range(mon_len):
            doy +=1
            time = tm.doy_to_utc(doy, 0, yyyy)
            print(time)
            sampl_file = gcdf.sampl_file +  "." + time[:10] + ".nc"
            assim_data = assim_tc.tccon_hm2(cur_step = cur_step,\
                                        doy = doy,\
                                        yyyy = yyyy,\
                                        sampl_file = sampl_file)
    
#     valid_data = valid_tc.tccon_hm2(cur_step = cur_step,\
#                                 doy = doy,\
#                                 yyyy = yyyy,\
#                                 sampl_file = None)
   # print(assim_data)
  #  print(valid_data)
    