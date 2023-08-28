#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
# %load rerun_geos.py
import sys
#sys.path.append('/home/shzhu/enkf_ch4/inv_code')
import os
import numpy as np
import numba as nb
import xarray as xr
import time_module as tm
import geos_chem_def as gcdf
from collections import defaultdict
import restart_cor as res_cor
# import res_stratos_cor_height as res_stra_cor
import res_stratos_cor as res_stra_cor
import calendar
from datetime import datetime, timedelta
import restart_config as res_config

def get_ensemble_list(yyyy, doy, entry_table_name):
    

    """ read into entry table """
    
    fl_entry = open(entry_table_name, "r")
    lines=fl_entry.readlines()
    fl_entry.close()
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
    em_st_list=list()
    em_end_list=list()
    em_yst_list=list()
    em_step_list=list()
    em_name_list=list()
    
    """ check the ensemble member to available """
    for vals in table_entry_list:
        if (cur_time>=vals[0] and cur_time <vals[1]):
            em_st_list.append(vals[2])
            em_end_list.append(vals[3])
            em_yyyy=int(vals[0]/1000.0)
            em_yst_list.append(em_yyyy)
            em_step=vals[4]
            em_step_list.append(em_step)
            em_name=vals[5]
            em_name_list.append(em_name)
        
    
#    return em_st_list, em_end_list, em_yst_list, em_step_list
    return em_step_list,em_name_list

def combine_res_nc(res1, tra_num1, res2, tra_num2):
    '''
    combine seperate ensemble tagrun restart files
    
    '''
    tracers = ['SpeciesRst_CH4_TAG_AN', 'SpeciesRst_CH4_TAG_NA']
    res3 = xr.Dataset()
#     drop_vars= []
#     for var in res2.data_vars:
#         if 'TAG' not in var:
#             drop_vars.append(var)
    for i_tra in tracers:
        for i in range(tra_num2):
            new_name = i_tra + str(i + 1 + tra_num1)
            old_name = i_tra + str(i + 1)

            #res2 = res2.rename_vars({old_name:new_name})#[new_name].values = en1[old_name].values
            res3[new_name] = res2[old_name] 
    #res2 = res2.drop(drop_vars)
    res = xr.merge([res1,res3])
    res.attrs = res1.attrs
    
    return res

def combine_res_nc2(step, year, doy):
    '''
    combine seperate ensemble tagrun restart files
    
    '''
    date = datetime(year, 1, 1) + timedelta(doy - 1)
    yyyy = date.strftime('%Y')
    date = date.strftime('%Y%m%d')
    var_list = []
    for file in gcdf.res_info_list:
        filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
        filename = file['filepath']+ str(yyyy) +'/'+filename
        if not os.path.exists(filename):
            cur_yyyy = str(int(yyyy) -1)
            filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
            filename = file['filepath']+ str(cur_yyyy) +'/'+filename
        print(filename)    
        tag_file = xr.open_dataset(filename)
        file_new = xr.Dataset()
        trac_st = file['trac_st']
        trac_num = file['trac_num']
        for i_tra in file['var_prefix']:
            for i in range(trac_num):
                new_name = i_tra + str(i + trac_st)
                old_name = i_tra + str(i + 1)

                file_new[new_name] = tag_file[old_name]
        var_list.append(file_new)
    en_res = xr.merge(var_list)
    en_res.attrs = tag_file.attrs
        
    
    return en_res

def combine_res_nc3(step, year, doy):
    '''
    combine seperate ensemble tagrun restart files
    
    '''
    date = datetime(year, 1, 1) + timedelta(doy - 1)
    yyyy = date.strftime('%Y')
    date = date.strftime('%Y%m%d')
    var_list = []
    for file in res_config.res_info_list:
        filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
        filename = file['filepath']+ str(yyyy) +'/'+filename
        if not os.path.exists(filename):
            cur_yyyy = str(int(yyyy) -1)
            filename = file['filename'].replace('YYYYMMDD', date).replace('ST000', step)
            filename = file['filepath']+ str(cur_yyyy) +'/'+filename
        print(filename)    
        tag_file = xr.open_dataset(filename)
        file_new = xr.Dataset()
        trac_old = file['old_item']
        trac_new = file['new_item']
        for i_tra in file['var_prefix']:
            for old, new in zip(trac_old, trac_new):
                new_name = i_tra + str(new)
                old_name = i_tra + str(old)

                file_new[new_name] = tag_file[old_name]
        var_list.append(file_new)
    en_res = xr.merge(var_list)
    en_res.attrs = tag_file.attrs
    return en_res

def prepare_restart( cur_step, doy_st, doy_end, yyyy, inc_m, fst_rst):
    ''' 
    perepare restart files for geos-chem rerun in the next step
    
    '''
   
    # read base restart file generated by rerunning  geos in current step 
    og_enaf = r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step -1, 1, 2)
    enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step , 1, 2)
    
    tst = tm.doy_to_utc(doy_end, 0, yyyy)
    tst = tst.replace('-', '')
    tst = tst.replace(':', '')
    days_of_year = 366 if calendar.isleap(yyyy) else 365
            
    if (doy_end > days_of_year):
 
        cur_yyyy = yyyy + 1
    else:
        cur_yyyy = yyyy
    #en_tst = tm.doy_to_utc(doy_st, 0, yyyy)
    #en_tst = en_tst.replace('-', '')
    #en_tst = en_tst.replace(':', '')
    
    full_restart_name='Restart.'+enaf+'.'+tst[0:8]+'_0000z.nc4'
    newfile=gcdf.new_res_path+"/"+full_restart_name
    og_restart_name = 'Restart.'+og_enaf+'.'+tst[0:8]+'_0000z.nc4'
    og_restart_name = gcdf.sgl_res_path+"/"+ og_restart_name
    print('og_restart_name: ', og_restart_name)
    print('newfile  : ',  newfile)
    
    # if first step: there is no adjustment
#     if fst_rst:
        
#         os_cp_cmd = 'cp '+ og_restart_name +' '+ newfile
#         os.system(os_cp_cmd)
        
#         return None
    
    new_data0 = xr.open_dataset(og_restart_name)
    
    
    # read base restart file generated by rerunning  geos in current step 
    ic = 0
    tracers_list= [sous +str(i_reg+1) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    ens_pos = gcdf.entry_table
    em_step_list,em_name_list = get_ensemble_list(yyyy, doy_st, entry_table_name = ens_pos)
    var_list= [''.join(['SpeciesRst_CH4_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    add_data = xr.zeros_like(new_data0['SpeciesRst_CH4'])
    print(em_step_list)
    if len(em_name_list) ==0:
        print('There is something wrong in ensemble restart files!')
        print('Orginal restart files', og_restart_name)
        print('Modified restart file:', newfile)
        return None
    for em_step, em_name in zip(em_step_list,em_name_list):

        em_restart_name = 'Restart.'+em_name+'.'+ tst[0:8] +'_0000z.nc4'
        
        em_restart_name = gcdf.en_res_path+"/"+ str(cur_yyyy)+ '/'+ em_restart_name
        print('em_restart:' ,em_restart_name)
        em_restart = xr.open_dataset(em_restart_name)
        
        for ivar in var_list:
            
            add_data.values = add_data.values + gcdf.pri_err*em_restart[ivar].values * inc_m[ic]
            ic = ic+1
            
    # combine base res and adjustment together 
    new_data = new_data0
    new_data['SpeciesRst_CH4'].values = new_data0['SpeciesRst_CH4'].values + add_data.values
    new_data.to_netcdf(newfile)
    
    
def prepare_restart2( cur_step, doy_st, doy_end, yyyy, inc_m, fst_rst):
    ''' 
    perepare restart files for geos-chem rerun in the next step
    
    '''
   
    # read base restart file generated by rerunning  geos in current step 
    og_enaf = r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step -1, 1, 2)
    enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step , 1, 2)
    
    tst = tm.doy_to_utc(doy_end, 0, yyyy)
    tst = tst.replace('-', '')
    tst = tst.replace(':', '')
    days_of_year = 366 if calendar.isleap(yyyy) else 365
            
    if (doy_end > days_of_year):
 
        cur_yyyy = yyyy + 1
    else:
        cur_yyyy = yyyy
    #en_tst = tm.doy_to_utc(doy_st, 0, yyyy)
    #en_tst = en_tst.replace('-', '')
    #en_tst = en_tst.replace(':', '')
    
    full_restart_name='Restart.'+enaf+'.'+tst[0:8]+'_0000z.nc4'
    newfile=gcdf.new_res_path+"/"+full_restart_name
    og_restart_name = 'Restart.'+og_enaf+'.'+tst[0:8]+'_0000z.nc4'
    og_restart_name = gcdf.sgl_res_path+"/"+ og_restart_name
    print('og_restart_name: ', og_restart_name)
    print('newfile  : ',  newfile)
    
    # if first step: there is no adjustment
#     if fst_rst:
        
#         os_cp_cmd = 'cp '+ og_restart_name +' '+ newfile
#         os.system(os_cp_cmd)
        
#         return None
    
    new_data0 = xr.open_dataset(og_restart_name)
    
    
    # read base restart file generated by rerunning  geos in current step 
    ic = 0
    #tracers_list= [sous +str(i_reg+1) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    tra_num1 = gcdf.tra_num1
    tra_num2 = gcdf.tra_num2
    em_step_list1,em_name_list1  = get_ensemble_list(yyyy, doy_st, \
                                                    entry_table_name= gcdf.entry_table1)
    em_step_list2,em_name_list2  = get_ensemble_list(yyyy, doy_st, \
                                                    entry_table_name= gcdf.entry_table2)
    
    var_list= [''.join(['SpeciesRst_CH4_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    add_data = xr.zeros_like(new_data0['SpeciesRst_CH4'])
    #print(em_step_list1)
    if len(em_name_list) ==0:
        print('There is something wrong in ensemble restart files!')
        print('Orginal restart files', og_restart_name)
        print('Modified restart file:', newfile)
        return None
    for em_step, em_name1, em_name2 in zip(em_step_list1,em_name_list1, em_name_list2):

        em_restart_name1 = 'Restart.'+em_name1+'.'+ tst[0:8] +'_0000z.nc4'
        em_restart_name2 = 'Restart.'+em_name2+'.'+ tst[0:8] +'_0000z.nc4'
        
        em_restart1 = gcdf.en_res_path+"/"+ str(cur_yyyy)+ '/'+ em_restart_name1
       
        em_restart2 = gcdf.en_res_path+"/"+ str(cur_yyyy)+ '/'+ em_restart_name2
        if not os.path.exists(em_restart1):
            em_restart1 =  gcdf.en_res_path+ '/' + str(cur_yyyy - 1) + '/'+em_restart_name1
            
        if not os.path.exists(em_restart2):    
            em_restart2 =  gcdf.en_res_path+ '/' + str(cur_yyyy - 1) + '/'+em_restart_name2
            
        print('em_restart1:' ,em_restart1)
        print('em_restart2:' ,em_restart2)
        res1 = xr.open_dataset(em_restart1)
        res2 = xr.open_dataset(em_restart2)
        em_restart = combine_res_nc(res1, tra_num1, res2, tra_num2)
        for ivar in var_list:

            add_data.values = add_data.values + gcdf.pri_err*em_restart[ivar].values * inc_m[ic]
            ic = ic+1
            
    # combine base res and adjustment together 
    new_data = new_data0
    new_data['SpeciesRst_CH4'].values = new_data0['SpeciesRst_CH4'].values + add_data.values
    new_data.to_netcdf(newfile)

def prepare_restart3( cur_step, doy_st, doy_end, yyyy, inc_m, fst_rst):
    ''' 
    perepare restart files for geos-chem rerun in the next step
    
    '''
   
    # read base restart file generated by rerunning  geos in current step 
    og_enaf = r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step -1, 1, 2)
    enaf=r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step , 1, 2)
    
    tst = tm.doy_to_utc(doy_end, 0, yyyy)
    tst = tst.replace('-', '')
    tst = tst.replace(':', '')

    
    full_restart_name='Restart.'+enaf+'.'+tst[0:8]+'_0000z.nc4'
    newfile=gcdf.new_res_path+"/"+full_restart_name
    og_restart_name = 'Restart.'+og_enaf+'.'+tst[0:8]+'_0000z.nc4'
    og_restart_name = gcdf.sgl_res_path+"/"+ og_restart_name
    print('og_restart_name: ', og_restart_name)
    print('newfile  : ',  newfile)

    
    new_data0 = xr.open_dataset(og_restart_name)
    
    
    # read base restart file generated by rerunning  geos in current step 
    ic = 0
    #tracers_list= [sous +str(i_reg+1) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
#     tra_num1 = gcdf.tra_num1
#     tra_num2 = gcdf.tra_num2
    em_step_list,em_name_list  = get_ensemble_list(yyyy, doy_st, \
                                                    entry_table_name= gcdf.entry_table)
    
    var_list= [''.join(['SpeciesRst_CH4_TAG_',sous[4:],str(i_reg+1)]) for sous in gcdf.emis_var_name for i_reg in range(gcdf.n_regs) ]
    add_data = xr.zeros_like(new_data0['SpeciesRst_CH4'])
    print(em_name_list)
    if len(em_name_list) ==0:
        print('There is something wrong in ensemble restart files!')
        print('Orginal restart files', og_restart_name)
        print('Modified restart file:', newfile)
        return None
    for em_step in em_name_list:
#         print(em_step[:5])
        em_restart = combine_res_nc3(em_step[:5], yyyy, doy_end)
        en_var_list = list(em_restart.keys())
#         if not var_list.sort() == en_var_list():
#             print('Something wrong with ensemble restart files combination!!!')
#             print( tst)
        for ivar in var_list:

            add_data.values = add_data.values + gcdf.pri_err*em_restart[ivar].values * inc_m[ic]
            ic = ic+1
            
    # combine base res and adjustment together 
    new_data = new_data0
    new_data['SpeciesRst_CH4_rerun'] = new_data0['SpeciesRst_CH4'].copy()
    new_data['SpeciesRst_CH4'].values = new_data0['SpeciesRst_CH4'].values + add_data.values
    
    # bias correction for new restart file
    if (gcdf.lat_bias_cor):
        psl_file = gcdf.pres_file.replace('YYYYMMDD', tst[0:8])
#         new_data = res_cor.Restart_bias_cor(new_data, \
#                                    ace_file = gcdf.ace_file, \
#                                    tropp_file = gcdf.tropp_file, \
#                                    sf_pres_file = ps_file, \
#                                    bins = gcdf.lat_bins)
#         new_data = res_stra_cor.Restart_bias_cor(org_res = new_data,\
#                                                  met_file =ps_file,\
#                                                  ace_file = gcdf.ace_file,\
#                                                  field_type = 'season')
        new_data = res_stra_cor.Restart_bias_cor1(org_res = new_data,\
                                                 pres_file = psl_file,\
                                                 ace_file = gcdf.ace_file,\
                                                 field_type = 'season')
    
    new_data.to_netcdf(newfile)

    
def rerun_geos_chem(doy_st, doy_end, cur_yyyy, cur_step, fst_rerun):
    '''
     rerun geos-chem with single tracer to get model simulation result in current step
     
    '''
    
    ## step1: determine geos-chem simulation period
        
    tst=tm.doy_to_utc(doy_st, 0, cur_yyyy)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '')
    
    tend=tm.doy_to_utc(doy_end, 0, cur_yyyy)
    tend=tend.replace('-', '')
    tend=tend.replace(':', '')
    
    rerun_path = gcdf.rerun_path
    ## step2: prepare restart file and HEMCO_Diagn.rc for current step
    enaf = r'ST%3.3d.EN%4.4d-EN%4.4d' % (cur_step, 1, 2)
    full_restart_name = 'Restart.'+enaf+'.'+tst[0:8]+'_0000z'+'.nc4'
    
    if (fst_rerun):
        
        ## copy defalut restart file at the first step  
        
        default_restart_file = gcdf.default_res_file
        os_cp_cmd = 'cp '+ default_restart_file+' '+ gcdf.new_res_path  + full_restart_name
        print('COPYING RESTART FILE',os_cp_cmd)
        os.system(os_cp_cmd)
        print('Finish copying restart file to res_file')
        
        ##  copy HEMCO_Diagn.rc for whole rerun periods (doesn't change)
        temp_Diagn = gcdf.rerun_temp_path + 'HEMCO_Diagn.temp'
        new_Diagn = rerun_path + 'HEMCO_Diagn.rc'
        os_cp_cmd = 'cp '+temp_Diagn+' '+ new_Diagn
        os.system(os_cp_cmd)
        
    
    ## step3: copy input.geos, HEMCO_Config.rc, HISTORY.rc from /Config_file 
    gen_config_file(tst, tend, enaf)
    
    input_new_file = 'input_'+enaf+'.geos'
    Config_new_file = 'HEMCO_Config_'+enaf +'.rc'
    His_new_file = 'HISTORY_' + enaf + '.rc'
    
   
    os.system('cp '+ gcdf.rerun_config_path+input_new_file+' '+rerun_path+'input.geos' )
    os.system('cp '+ gcdf.rerun_config_path+Config_new_file+' '+rerun_path+'HEMCO_Config.rc' )
    os.system('cp '+ gcdf.rerun_config_path+His_new_file+' '+rerun_path+'HISTORY.rc' )
    
        
    print('Finish copying configure files')
     
    ## step4: rerun geos.chem
    os.chdir( rerun_path )
    os.system('pwd')
#   os_rerun_cmd =  rerun_path+'geos.mp'  
#    os.system(os_rerun_cmd)
#    if not fst_rerun:
        
    os.system('./geos.mp' )
    
    print('Finish geos_chem rerun')
    run_dir = gcdf.root 
    os.chdir(run_dir)
    os.system('pwd')
    
    
def gen_config_file(tst, tend, enaf):    
    '''
    generate configure files
    
    '''
   # enafix = r'ST%3.3d.EN%4.4d-EN%4.4d' % (1, 1, 2)
    
    ### input.geos
    
    diag_species = '1 512'
    input_new_file = 'input_'+enaf+'.geos'
    temp_input = gcdf.rerun_temp_path + 'input.geos.'+ gcdf.grid+'.temp'
    new_input = gcdf.rerun_config_path + input_new_file
    
    tm_res = gcdf.temp_res
    if (tm_res =='mon'):
        tm_res_str = '00000100 000000'
    else:
        tm_res_str ='%08d'%(tm_res) + ' ' + '0'*6  
    
    # assign emis parameters
    emis_path = gcdf.emis_path
    sous_name = gcdf.sous_name
    emis_file_name = gcdf.emis_file_name
    emis_var_name = gcdf.emis_var_name
    emis_time = gcdf.emis_time
    emis_unit = gcdf.emis_unit
    emis_rank = gcdf.emis_rank
    n_sous=gcdf.n_sous
    
    input_file = open(temp_input, 'r')
    output_file = open(new_input, 'w')
    ts = input_file.read().split('------------------------+------------------------------------------------------')
    for t in ts:
        lines = t.strip('\n').split('\n')
        if lines[0] == '%%% SIMULATION MENU %%% :':
            for line in lines: 
                line=line.replace('$RUNPATH', gcdf.rerun_path)
                line=line.replace('$DATAPATH', gcdf.data_path)
                if (line.split(':')[0].strip(' ') == 'Start YYYYMMDD, hhmmss'):
                    output_file.write(line.split(':')[0]+': '+ tst +'\n')
                elif (line.split(':')[0].strip(' ') == 'End   YYYYMMDD, hhmmss'):
                    output_file.write(line.split(':')[0]+': '+ tend +'\n')
                else:
                    output_file.write(line+'\n')
        elif 'ND51' in lines[0]:
            for line in lines:
                if (line.split(':')[0].strip(' ') == 'Tracers to include'):
                    output_file.write(line.split(':')[0]+': '+ diag_species + '\n')
                else:
                    line=line.replace('$DATAPATH', gcdf.rerun_diagn_path+'/ND51')
                    line=line.replace('STYYY.ENXXXX-ENXXXX',  enaf)
#                     line=line.replace('STYYY.ENXXXX-ENXXXX',  r'ST%3.3d.EN%4.4d-EN%4.4d' % (0, 1, 2))
                    output_file.write(line+'\n')
                    
        elif lines[0] == '%%% CH4 MENU %%%        :':
            output_file.write(lines[0]+'\n') 
            for line in lines: 
                if (line.split(':')[0].strip(' ') == 'Use TCCON obs operator?'):
                    output_file.write(line+'\n')
                if (line.split(':')[0].strip(' ') == 'Use GOSAT obs operator?'):
                    output_file.write(line+'\n')
                if (line.split(':')[0].strip(' ') == 'Ensemble Temporal Res'):
                    output_file.write(line.split(':')[0]+': '+ tm_res_str+'\n')
                if (line.split(':')[0].strip(' ') == 'The Number of Regions'):
                    output_file.write(line.split(':')[0]+': '+ str(1)+'\n')    
                if (line.split(':')[0].strip(' ') == 'The Number of Sources'):
                    output_file.write(line.split(':')[0]+': '+ str(n_sous)+'\n')
                if (line.split(':')[0].strip(' ') == 'CH4 Sources Entries --->'):
                    output_file.write(line+'\n')
                    if (len(sous_name) == n_sous):
                        for i_sous in range(n_sous):
                            output_file.write('Sources name            : '+ sous_name[i_sous]+'\n')        
                #output_file.write(line+'\n')
        else:
            for line in lines:
                line=line.replace('$DATAPATH', gcdf.rerun_diagn_path)
                line=line.replace('STYYY.ENXXXX-ENXXXX',  enaf)
                
                output_file.write(line+'\n')
        output_file.write('------------------------+------------------------------------------------------'+'\n')  

    output_file.close()
    input_file.close()      
    print(enaf + ' input.geos done')
    
    
    ### HEMCO_Config.rc
    full_restart_name = 'Restart.'+enaf+'.'+tst[0:8]+'_0000z'+'.nc4'
    input_resflnm =  gcdf.new_res_path + full_restart_name
    Config_new_file = 'HEMCO_Config_'+enaf +'.rc'
    temp_config = gcdf.rerun_temp_path+ 'HEMCO_Config.temp'
    new_config = gcdf.rerun_config_path+ Config_new_file
    config_input = open(temp_config, "r")  #temp_path
    config_output = open(new_config, "w")
    
    
    
    # build a dict contains emission information
    info = defaultdict(dict)
    for iemis in range (n_sous):
        info['ch4emis{0}'.format(iemis)]['var_name']  = emis_var_name[iemis]
        info['ch4emis{0}'.format(iemis)]['sous_name'] = sous_name[iemis]
        info['ch4emis{0}'.format(iemis)]['file_path'] = emis_path + emis_file_name[iemis]
        info['ch4emis{0}'.format(iemis)]['emis_time'] = emis_time[iemis]
        info['ch4emis{0}'.format(iemis)]['emis_unit'] = emis_unit
        info['ch4emis{0}'.format(iemis)]['emis_rank'] = emis_rank
        
    emi_sec =[]
    cat = 1
    for key, value in info.items():
        line ='### CH4 EMISSION by ' + info[key]['var_name']+'\n' 
        emi_sec.append(line)
        emi_sec.append(' '+'\n')
        emi_sec.append('0  ' + info[key]['sous_name'] + ' '*6 + info[key]['file_path'] + ' '* 3 + info[key]['var_name']\
                            + ' '*3  + info[key]['emis_time'] +' C '+ info[key]['emis_rank'] + ' '\
                            +  info[key]['emis_unit']+ ' '*2 + 'CH4' + ' - '\
                            +  str(cat) + '  1' +'\n')
        cat = cat+1
    lines = config_input.readlines() 
    for index, line in enumerate(lines):
        if '(((SP_CH4 ' in line:
            for iemi in range(len(emi_sec)):
                lines.insert(index + iemi + 1, emi_sec[iemi])
                
    for line in lines:
        
        line = line.replace('$HEM_PATH', gcdf.hem_path) 
        line = line.replace('$MET_PATH', gcdf.met_path) 
        line = line.replace('$DiagnPrefix', gcdf.rerun_diagn_path + '/HEMCO_diagnostics' ) 
        
        line = line.replace('$RESPATH', input_resflnm)
        line = line.replace('$MASKFILE', gcdf.sgl_mask_file)
        config_output.writelines(line)
        
    config_input.close()    
    config_output.close() 
    print(enaf + ' HEMCO_Config.rc done')
    
    
    ### HISTORY.rc
    temp_His = gcdf.rerun_temp_path+ 'HISTORY.temp'
    His_new_file = 'HISTORY_' + enaf + '.rc'
    new_His  = gcdf.rerun_config_path+ His_new_file
     
    
    hist_input = open(temp_His, "r")
    hist_output = open(new_His, "w")
    lines = hist_input.readlines()
    
    resflnm = gcdf.sgl_res_path +'Restart.'+ enaf+'.%y4%m2%d2_%h2%n2z.nc4'
    for line in lines:
        line = line.replace('$Restart_Path', resflnm )
        line = line.replace('$Diagn_Path', gcdf.rerun_diagn_path+'/Boundary/'+ enaf)
        line = line.replace('$Restart_Resolution', tm_res_str) 
        line = line.replace('$Output_resolution', tm_res_str)
        hist_output.writelines(line)
    hist_input.close()
    hist_output.close()
    
    print(enaf + ' HISTORY.rc done')
    print('Config files done!')
    
    
    
if (__name__=="__main__"): 
#    tst= '20150601 000000'
#    tend = '20150701 000000'
#    tm_res= gcdf.temp_res 
#    enaf = 'ST005.EN0001-EN0002'
#    gen_config_file(tst, tend, tm_res, enaf) 

#     doy_st =1
#     doy_end = 32
#     cur_yyyy = 2015
#     cur_step = 0
#     fst_rst = False
#     rerun_geos_chem(doy_st, doy_end, cur_yyyy, cur_step, fst_rst)
    cur_step = 1
    doy_st = 1 
    doy_end = 32 
    yyyy = 2015
    inc_m = np.random.rand(104)
    fst_rst = True
    prepare_restart2( cur_step, doy_st, doy_end, yyyy, inc_m, fst_rst)
    