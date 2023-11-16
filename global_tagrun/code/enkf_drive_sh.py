#!/home/shzhu/apps/miniconda3/envs/geo/bin/python
# -*- coding: utf-8 -*-
"""
this program is used to pre-calculate the response functions for the pre-defined perturbations (basis functions)
Created on Thu Dec 12 11:13:17 2019

@author: rainbow s2008420
"""

# includes modules   
import sys                        # system commands 
import os                         #  operational system commands 
from calendar import monthrange
import time_module as tm          #  time conversion
import gen_GC_Config as gcc       #  generate config files for GEOS-Chem
import time as systime 
import geos_chem_def as gcdf      # deflaut settings
import read_em_cfg as emcfg       # configurations for pre-define perturbation functions
from numpy import *
import signal                     # use for kill execution

print('*'*80)
print(''*30+'CH4 ENSEMBLE RUN DRIVER'+'*'*30)
print('*'*80)
print(' '*30)    

#>/ps/read in configuration for ensemble run

#re_run_now = gcdf.rerun_now
##/c/restart date
#re_run_date = gcdf.rerun_date

##/c/starting year/mm/dd

#yyyy = gcdf.st_yyyy  # the year  
#mm = gcdf.st_mm     # the month 
#dd = gcdf.st_dd     #  the day

temp_res = gcdf.temp_res
nstep = gcdf.nstep #  the last day of geos-chem simulation is given by ntime * temp_res 


##/c/lag window

lag_window = gcdf.inv_lag_window


##/c/on to reset the initial distribitribution to fixed values

new_restart = gcdf.new_restart

##/c/perturbation configuration files (BF)

em_cfg = gcdf.cfg_pb_file
##/c/read in starting doy/yyyy/number of pb functions/stored files
em_doy_st, em_yyyy_st, em_npb, em_ch4flnm =emcfg.read_em_conf(em_cfg)

print(em_doy_st)
print(em_yyyy_st)
print(em_npb)


##/c/desk clock for jobs starting

gmt = systime.gmtime()

##/c/output file for configurations to be used in inversion
ftt = open('./ens_pos.dat', "w")

##/c/the head lines  

line =\
   'geos_chem run at %4.4d%2.2d%2.2d, %2.2d:%2.2d:%2.2d' %\
    (gmt[0], gmt[1], gmt[2], gmt[3], gmt[4], gmt[5])

ftt.write(line+'\n')
line = r'temp_res: %4.4d  nstep: %4.4d' % (30, nstep)
ftt.write(line+'\n')
##/c/time step/member start/member end/year start/year end/doy start/doy end/BF file/extension of model output/method_for_select_member_for_inversion/
line =\
   'step,  mem_st,  mem_end,  year_st,  year_end,\
    doy_st,  doy_end,  ch4flnm, output_name, index_method'
ftt.write(line+'\n')

#>/ps/do the ensemble simulations


select_runs = gcdf.select_runs
default_restart_file = gcdf.restart_file
sub_sous = gcdf.sub_sous
sub_tracers = gcdf.sub_tracers

print('nstep',nstep)
print('select-runs',select_runs)

print('step range:',list(range(gcdf.nstep_start, nstep)))

#gcdf.restart_file

#>>/ps/loop0/over each step

for istep in range(gcdf.nstep_start, nstep): # run through the emission time period 
    npbf = em_npb[istep]
    pbflnm = em_ch4flnm[istep]
    yyyy_st = em_yyyy_st[istep]
    doy_st = em_doy_st[istep]
    yyyy, mm, dd = tm.doy_to_time_array(doy_st, yyyy_st)
    ##/c/time resolution  (the exact temporal resolutions is controled by the emission
    if (gcdf.temp_res == 'mon'):    
        temp_res = monthrange(yyyy, mm)[1]
        rst_res = 'mon'
    else:
        temp_res = gcdf.temp_res
        rst_res = temp_res
    ##c/tau_st in seconds since 1985.1.1.0.0
    tau_st = tm.get_tau(yyyy,mm, dd)

    print('-'*10+'starting year, month, day:', yyyy, mm, dd)
    
    #istep_end = min(istep+lag_window, len(em_doy_st)-1)
    istep_end = istep+lag_window
    print('-'*10+'istep_start, istep_end:', istep, istep_end)
    
    if(istep == istep_end):
        gmt = systime.gmtime()
        line =\
        'geos_chem end at %4.4d%2.2d%2.2d, %2.2d:%2.2d:%2.2d' %\
        (gmt[0], gmt[1], gmt[2], gmt[3], gmt[4], gmt[5])
        ftt.write(line+'\n')
        ##c/ close ensemble run file
        #ftt.close()
        os.kill(os.getpid(), signal.SIGKILL)
        
    
    if (sub_sous):
        nsub_sous = gcdf.n_sous_loop 
    # sub emission sources loop
    else:
        nsub_sous = 1
        
    if (sub_tracers):
        nsub_run = gcdf.n_tracer_loop
        a_mst = gcdf.tracer_start
        a_mend = gcdf.tracer_end 
        input_n_regs = gcdf.tracer_num
    else:
        nsub_run = 1
        a_mst = [1]
        a_mend = [gcdf.n_regs]
        input_n_regs = [gcdf.n_regs]
        
    

    print('range(n_subsous),range(nsub_run):',nsub_sous, nsub_run )

    

    #>>>/ps/loop0/loop1/sub_run to cover each emission sources

    ##/c/1/relative positition of the first used perturbation ensemble
    ##/c/2/(i.e, the second tracer) in the pb functions
    pbst = 0
    #exit()
    
    for isub_sous in range(nsub_sous):
        sous_name =[]
        emis_file_name = []
        emis_var_name=[]
        emis_time = []
        ## define sous_name used in input.geos
        if (sub_sous):
            sous_name.append(gcdf.sous_name[isub_sous])
            emis_file_name.append(gcdf.emis_file_name[isub_sous])
            emis_var_name.append(gcdf.emis_var_name[isub_sous])
            emis_time.append(gcdf.emis_time[isub_sous])
            input_n_sous=1
            
        else:
            sous_name = gcdf.sous_name
            emis_file_name = gcdf.emis_file_name
            emis_var_name = gcdf.emis_var_name
            emis_time = gcdf.emis_time
            input_n_sous=gcdf.n_sous
        #>>>/ps/loop0/loop1/sub_run to cover sub_region species	
        for isub_run  in range(nsub_run): 
            mst = a_mst[isub_run]
            mend = a_mend[isub_run]
            line =\
				r'%3.3d, %4.4d, %4.4d, %4.4d, %3.3d, %3.3d, %3.3d' %\
				(istep, \
				mst, mend, \
				em_yyyy_st[istep], em_yyyy_st[istep_end],
				em_doy_st[istep], \
				em_doy_st[istep_end])
       
    
            ##c/print time_start
            ##c/generate file extension for restart files
            if (sub_sous):
                enaf =r'ST%3.3d.SR%2.2d-%2.2d.EN%4.4d-EN%4.4d' % (istep,input_n_sous,isub_sous+1, mst, mend)
              
            else:
                enaf =r'ST%3.3d.EN%4.4d-EN%4.4d' % (istep, mst, mend)
            input_new_file = 'input_'+enaf+'.geos'
            Config_new_file = 'HEMCO_Config_'+enaf +'.rc'
            Diag_new_file = 'HEMCO_Diagn_' + enaf + '.rc'
            His_new_file = 'HISTORY_' + enaf + '.rc'
            ##c finalize the information line for outputs/0 means the 
            ##default indexing method

            line = line+ ','+pbflnm+','+enaf+',0\n'
            #print('line',line)
            ftt.write(line+'\n')
        
            
            #>>>>/ps/loop0/loop1/loop2/run each step
            
            print('em_doy_st[istep:istep_end]',istep,istep_end,em_doy_st[istep:istep_end])
        
            #>>>>/ps/loop0/loop1/loop2/if select run start 
            if (istep in gcdf.select_runs):
                yst = em_yyyy_st[istep]
                dst = em_doy_st[istep]
                print('yst:', yst)
                    
                tst = tm.doy_to_utc(dst, sec=0, yyyy=yst)
                print('tst:', tst)
                    
                tst = tst.replace('-', '')
                tst = tst.replace(':', '')
                time_start = r'%4.4d%3.3d' % (yst, dst)
                                
                #>>>>/ps/generate restart file
               
                
                output_resflnm = gcdf.res_path+'Restart.'+ enaf+'.%y4%m2%d2_%h2%n2z.nc4'
                
                ##c/generate configuration file for ensemble run
                print('======>Step 1: Generate input file <======')
                gcc.gen_input_geos(run_step = istep,\
                                   YYYY =[em_yyyy_st[istep], em_yyyy_st[istep_end]],\
                                   DOY  =[em_doy_st[istep], em_doy_st[istep_end]],\
                                   member_start=mst,\
                                   member_end=mend,\
                                   tmpfile = 'input.geos.' + gcdf.grid + '.temp', \
                                   newfile = input_new_file, \
                                   temp_path = gcdf.temp_path,\
                                   run_path= './',\
                                   data_path= gcdf.data_path,\
                                   out_path = gcdf.config_path,\
                                   diagn_path= gcdf.diagn_path,\
                                   n_regs = input_n_regs[isub_run],\
                                   n_sous = input_n_sous,\
                                   sous_name = sous_name,\
                                   temporal_res = temp_res)
                    
        #            os.system("mv input.geos.new input.geos")
                print('======>Step 2: Generate HEMCO_Config file <======')
                gcc.gen_hemco_config(run_step = istep,\
                                    temp_path= gcdf.temp_path,\
                                    data_path = gcdf.data_path,\
                                    met_path =gcdf.met_path,\
                                    emis_path =gcdf.emis_path,\
                                    out_path = gcdf.config_path,\
                                    mask_file = gcdf.mask_path+gcdf.mask_file[isub_run],\
                                    diagn_path= gcdf.diagn_path,\
                                    emis_file_name = emis_file_name,\
                                    emis_var_name = emis_var_name,\
                                    emis_time = emis_time,\
                                    emis_unit = gcdf.emis_unit,\
                                    emis_rank = gcdf.emis_rank,\
                                    member_start=mst, \
                                    member_end=mend,\
                                    tmpfile = 'HEMCO_Config.temp',\
                                    newfile = Config_new_file,\
                                    resfile = default_restart_file,\
                                    n_regs = input_n_regs[isub_run],\
                                    n_sous = input_n_sous,\
                                    sous_name = sous_name)
                    
                print('======>Step 3: Generate HEMCO_Diagn file <======')
                gcc.gen_hemco_diagn(out_path = gcdf.config_path,\
                                    diagn_name = Diag_new_file,\
                                    member_start = mst,\
                                    member_end = mend,\
                                    n_regs = input_n_regs[isub_run],\
                                    n_sous = input_n_sous,\
                                    sous_name = sous_name)
                    
                print('======>Step 4: History file<=====')
                gcc.gen_history(temp_path= gcdf.temp_path,\
                                out_path = gcdf.config_path,\
                                res_path = gcdf.res_path,\
                                tmpfile  = 'HISTORY.temp',\
                                newfile  = His_new_file,\
                                tm_res   = rst_res,\
                                resfile = output_resflnm,\
                                diagfile = gcdf.bound_path+'/'+ enaf)    
                                
                #>>>>>/ps/loop0/loop1/loop2/print information
                print('ISTEP, GCDF.SELECT_RUNS',istep, gcdf.select_runs)

                print('original input file, new input file:', gcdf.config_path+input_new_file, './'+'input.geos' )
                print('original HEMCO_Config file, new HEMCO_Config file:',gcdf.config_path+Config_new_file, './'+'HEMCO_Config.rc' )
                print('original HEMCO_Diagn file, new HEMCO_Diagn file:',gcdf.config_path+Diag_new_file, './'+'HEMCO_Diagn.rc' )
                print('original HISTORT file, new HISTORY file:',gcdf.config_path+His_new_file, './'+'HISTORY.rc' )
                #>>>>> copy config files to run path
                print('start to copy files!')
                os.system('cp '+ gcdf.config_path+input_new_file+' '+'./input.geos' )
                os.system('cp '+ gcdf.config_path+Config_new_file+' '+'./HEMCO_Config.rc' )
                os.system('cp '+ gcdf.config_path+Diag_new_file+' '+'./HEMCO_Diagn.rc' )
                os.system('cp '+ gcdf.config_path+His_new_file+' '+'./HISTORY.rc' )
                
                print('finish copying files!')
                #>>>>> lauch the CTM 
                #os.system('nohup ./geos.mp > '+ enaf+'.file 2>&1 &')
                os.system('./geos.mp > '+ enaf+'.file')
#                os.system('./geos.mp')
            
            #>>>>/ps/loop0/loop1/loop2/if select run end        

        
        #>>>>/ps/loop0/loop1/loop2/ sub species(regions) end
        

    #>>>/ps/loop0/loop1/sub sources end

    
#>>/ps/loop0/each step end

gmt = systime.gmtime()
line =\
   'geos_chem end at %4.4d%2.2d%2.2d, %2.2d:%2.2d:%2.2d' %\
    (gmt[0], gmt[1], gmt[2], gmt[3], gmt[4], gmt[5])
ftt.write(line+'\n')
##c/ close ensemble run file
ftt.close()
