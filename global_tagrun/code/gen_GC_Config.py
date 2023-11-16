# -*- coding: utf-8 -*-
"""
 generate GEOS-Chem outer files(input.geos; HEMCO_Config.rc; HEMCO_Diagn.rc) 
 
 by Rainbow

 09/12/2019
"""
from numpy import *
import time_module as tm
import geos_chem_def as gcdf
from collections import defaultdict

def gen_input_geos(run_step,\
                   YYYY,\
                   DOY,\
                   member_start =1,\
                   member_end = 40,\
                   tmpfile="input.geos.temp",\
                   newfile="input.geos" ,\
                   temp_path = './temp/',\
                   run_path='./',\
                   data_path='./',\
                   out_path='./Config_files',\
                   diagn_path='./Diagn_files',\
                   n_regs = 20,\
                   n_sous = 2,\
                   sous_name = ['CH4_TAG_AN','CH4_TAG_NA'],\
                   temporal_res = 30,\
                   **keywords \
                          ):
    """ create the new input to drive the ensemble run
    YYYY: years of the time serier
    DOY:  doys of the time serier
    em_doy: the date need emission
    temp_path: the path of input.geos.temp
    n_regs: the number of regions
    n_sous: the number of emission sources(defualt incloude anthr and nature)
    sous_name: the name of emission sources     
    temporal_res: ensemble run temperal resolution (unit: days)
    
    """
    
    yst=YYYY[0]
    dst=DOY[0]
    yed=YYYY[1] 
    ded=DOY[1]
    # list for start and ending simulation time 
    
    tst=tm.doy_to_utc(dst, sec=0, yyyy=yst)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '')
    
    tend=tm.doy_to_utc(ded, sec=0, yyyy=yed)
    tend=tend.replace('-', '')
    tend=tend.replace(':', '')
    diag_species = '  '
    for i_spc in range(n_sous * n_regs + 1):
        diag_species = diag_species + str(i_spc + 1) + ' '*2
    diag_species = diag_species + '512'
    enafix=r'ST%3.3d.EN%4.4d-EN%4.4d' % (run_step, member_start, member_end)
    enfix =r'EN%4.4d-EN%4.4d' % (member_start, member_end)
    temp_res_d = r'%08d'%(temporal_res)
    
    input_file = open(temp_path+tmpfile, 'r')
    output_file = open(out_path+newfile, 'w')
    ts = input_file.read().split('------------------------+------------------------------------------------------')

    for t in ts:
        lines = t.strip('\n').split('\n')
        if lines[0] == '%%% SIMULATION MENU %%% :':
            for line in lines: 
                line=line.replace('$RUNPATH', run_path)
                line=line.replace('$DATAPATH', data_path)
                if (line.split(':')[0].strip(' ') == 'Start YYYYMMDD, hhmmss'):
                    output_file.write(line.split(':')[0]+': '+ tst +'\n')
                elif (line.split(':')[0].strip(' ') == 'End   YYYYMMDD, hhmmss'):
                    output_file.write(line.split(':')[0]+': '+ tend +'\n')
                else:
                    output_file.write(line+'\n') 
                
        elif lines[0] == '%%% CH4 MENU %%%        :':
            output_file.write(lines[0]+'\n') 
            for line in lines: 
                if (line.split(':')[0].strip(' ') == 'Use TCCON obs operator?'):
                    output_file.write(line+'\n')
                if (line.split(':')[0].strip(' ') == 'Use GOSAT obs operator?'):
                    output_file.write(line+'\n')
                if (line.split(':')[0].strip(' ') == 'Ensemble Temporal Res'):
                    output_file.write(line.split(':')[0]+': '+ str(temp_res_d)+' ' + '0'*6+'\n')
                if (line.split(':')[0].strip(' ') == 'The Number of Regions'):
                    output_file.write(line.split(':')[0]+': '+ str(n_regs)+'\n')    
                if (line.split(':')[0].strip(' ') == 'The Number of Sources'):
                    output_file.write(line.split(':')[0]+': '+ str(n_sous)+'\n')
                if (line.split(':')[0].strip(' ') == 'CH4 Sources Entries --->'):
                    output_file.write(line+'\n')
                    if (len(sous_name) == n_sous):
                        for i_sous in range(n_sous):
                            output_file.write('Sources name            : '+ sous_name[i_sous]+'\n')
                        
        elif lines[0] == '%%% ADVECTED SPECIES MENU %%%:':  
            output_file.write(lines[0]+'\n')
            output_file.write('Type of simulation      : 9' + '\n')
            output_file.write('Species Entries ------->: Name' + '\n')
            output_file.write('Species name            : CH4' + '\n')
            for i_sous in range(n_sous):
                for i_regs in range(n_regs):
                    output_file.write('Species name            : '+ sous_name[i_sous]+str(i_regs+1)+'\n')
        

        elif 'ND51' in lines[0]:
            for line in lines:
                if (line.split(':')[0].strip(' ') == 'Tracers to include'):
                    output_file.write(line.split(':')[0]+': '+ diag_species + '\n')
                else:
                    line=line.replace('$DATAPATH', diagn_path)
                    line=line.replace('STYYY.ENXXXX-ENXXXX',  enafix)
                    output_file.write(line+'\n')	
                
                #output_file.write(line+'\n')
        else:
            for line in lines:
                line=line.replace('$DATAPATH', diagn_path)
                line=line.replace('STYYY.ENXXXX-ENXXXX',  enafix)
                line=line.replace('ENXXXX-ENXXXX',  enfix)
                output_file.write(line+'\n')
        output_file.write('------------------------+------------------------------------------------------'+'\n')  
    print('reach the end')
    output_file.close()
    input_file.close()




def gen_hemco_diagn(out_path = './Config_files',\
                    diagn_name = 'HEMCO_Diagn.rc',\
                    member_start = 1,\
                    member_end = 40,\
                    n_regs = 20,\
                    n_sous = 2,\
                    sous_name = ['CH4_TAG_AN','CH4_TAG_NA'],\
                    **keywords\
                    ):
    """ create the new HEMCO_Diagn.rc to drive the ensemble run
   
    n_regs: the number of regions
    n_sous: the number of emission sources(defualt incloude anthr and nature)
    sous_name: the name of emission sources     
    
    """

    filename = out_path + diagn_name
    diag_file = open(filename,'w')
    diag_file.write('# Name                    Spec   ExtNr  Cat Hier Dim OutUnit\n')
    diag_file.write('EMIS_CH4_TOTAL            CH4   -1     -1   -1   2   molec/cm2/s\n')
    for i_sous in range(n_sous):
        for i_regs in range(n_regs):
            name = 'EMIS_'+sous_name[i_sous]+str(i_regs+1)
            name_size = 17
            name = name.ljust(name_size)
            cat  = str(i_sous*n_regs+i_regs+1)
            cat  = cat.rjust(2)
     #   print(cat)
            diag_file.write(name + '         CH4    0     '+ cat +'   -1   2   molec/cm2/s'+'\n')
    diag_file.close()   
 


    

def gen_hemco_config(run_step,\
                    member_start =1,\
                    member_end = 40,\
                    temp_path ='./temp/',\
                    data_path = '/data/geos-chem/ctm',\
                    met_path ='/scratch/local/shzhu/GC_DATA/ctm/GEOS_4x5/GEOS_FP',\
                    emis_path ='/exports/csce/datastore/geos/users/s2008420/',\
                    out_path = './Config_files',\
                    mask_file = '/exports/csce/datastore/geos/users/s2008420/global_mask_20.nc',\
                    diagn_path='./Diagn_files',\
                    emis_file_name = ['emi_anthro_s.nc','emi_nature_s.nc'],\
                    emis_var_name = ['Anthr' , 'Nature'],\
                    emis_time = ['2012/1/1/0','2015/1-12/1/0'],\
                    emis_unit = 'molec/cm2/s',\
                    emis_rank = 'xy',\
                    tmpfile = 'HEMCO_Config.temp',\
                    newfile = 'HEMCO_Config.rc',\
                    resfile = './OutputDir/GEOSChem.Restart.%y4%m2%d2_%h2%n2z.nc4',\
                    n_regs = 20,\
                    n_sous = 2,\
                    sous_name = ['CH4_TAG_AN','CH4_TAG_NA'],\
                    **keywords):
    
    """ create the new HEMCO_Config.rc to drive the ensemble run
   
    temp_path: the path of input.geos.temp
    met_path: the path of meterological fields
    emis_path: the path of emission files
    emis_file_name: the filename of emission nc files
    emis_var_name: the vaiable names in emission.nc 
    emis_time: the timeseries of data available ()
    emis_unit: the unit of emissions data
    emis_rank: the rank of emissions data
    ########## number of parameters in above lists should be same as n_sous##################
    tmpfile: temperal history file name
    n_regs: the number of regions
    n_sous: the number of emission sources(defualt incloude anthr and nature)
    sous_name: the name of emission sources
    
    """
    print(emis_var_name,n_sous,n_regs)
    info = defaultdict(dict)
    for iemis in range (n_sous):
        info['ch4emis{0}'.format(iemis)]['var_name'] = emis_var_name[iemis]
        info['ch4emis{0}'.format(iemis)]['sous_name'] = sous_name[iemis]
        info['ch4emis{0}'.format(iemis)]['file_path'] = emis_path + emis_file_name[iemis]
        info['ch4emis{0}'.format(iemis)]['emis_time'] = emis_time[iemis]
        info['ch4emis{0}'.format(iemis)]['emis_unit'] = emis_unit
        info['ch4emis{0}'.format(iemis)]['emis_rank'] = emis_rank
    tagemi =[]
    i_sous = 0
    for key, value in info.items():
        line ='### CH4 EMISSION by ' + info[key]['var_name']+'\n' 
        tagemi.append(line)
        tagemi.append(' '+'\n')  
        for i_regs in range(n_regs):
            tag_name =  info[key]['sous_name']+str(i_regs+1)
            tag_name_size = 12
            tag_name = tag_name.ljust(tag_name_size)
            cat  = str(i_sous * n_regs + i_regs + 1)
            cat  = cat.rjust(2)
            tagemi.append('0  ' + tag_name + ' '*6 + info[key]['file_path'] + ' '* 3 + info[key]['var_name']\
                            + ' '*3  + info[key]['emis_time'] +' C '+ info[key]['emis_rank'] + ' '\
                            +  info[key]['emis_unit']+ ' '*2 + tag_name + ' - '\
                            +  cat + '  1' +'\n')
        i_sous += 1
        tagemi.append(' '+'\n') 
        
    config_input = open(temp_path + tmpfile, "r")  #temp_path
    config_output = open(out_path + newfile, "w")
    lines = config_input.readlines()
    config_input.close()
    
    for index, line in enumerate(lines):
        if '(((TAGCH4' in line:
            for iemi in range(len(tagemi)):
                lines.insert(index + iemi + 1, tagemi[iemi])
    enafix=r'ST%3.3d.EN%4.4d-EN%4.4d' % (run_step, member_start, member_end)   
    hemco_path = data_path+'/HEMCO' 
    for line in lines:
        line=line.replace('$HEMCOPATH', hemco_path) 
        line=line.replace('$METPATH', met_path) 
        line=line.replace('$Diagn_Path', diagn_path)
        line=line.replace('$MASK_FILE',  mask_file)
        line=line.replace('$RESPATH', resfile) 
    
    
        config_output.writelines(line)
        
    config_output.close() 
    
    
def gen_history(temp_path='./temp/',\
                out_path = './Config_files',\
                tmpfile  = 'HISTORY.temp',\
                newfile  = 'HISTORY.rc',\
                tm_res= 30,\
                resfile = './OutputDir/GEOSChem.Restart.%y4%m2%d2_%h2%n2z.nc4',\
                diagfile = './Res',\
                **keywords
                ):
    """ create the new HEMCO_Diagn.rc to drive the ensemble run
   
    temp_path: path of history.temp
    run_path: path of GEOS-Chem rundirs
    tmpfile: name of history.temp
    newfile: name of history.rc
    tm_res= 30: termporal resolution of ensemble run (unit: day)
    res_name: name of restart files
    
    """
    if (tm_res =='mon'):
        tm_res_str = '00000100 000000'
    else:
        tm_res_str ='%08d'%(tm_res) + ' ' + '0'*6    
    
    print(tm_res_str)
    hist_input = open(temp_path + tmpfile, "r")
    hist_output = open(out_path + newfile, "w")
    lines = hist_input.readlines()
    hist_input.close()

    for line in lines:
        line = line.replace('$Restart_Path', resfile )
        line = line.replace('$Diagn_Path', diagfile)
        line = line.replace('$Restart_Resolution', tm_res_str) 
        line = line.replace('$Output_resolution', tm_res_str)
        hist_output.writelines(line)
    hist_output.close()        

if (__name__=="__main__"):
    
    
    data_path='/scratch/local/shzhu/GC_DATA/ctm'
    YYYY=[2016, 2016]
    DOY=[32, 32+28]
    
    
    gen_input_geos(YYYY,\
                    DOY,\
                    member_start =1,\
                    member_end = 40,\
                    tmpfile="input.geos.temp", \
                    newfile="input.geos" , \
                    temp_path = './temp/',\
                    run_path= './',\
                    data_path=data_path,\
                    n_regs = 20,\
                    n_sous = 2,\
                    sous_name = ['CH4_TAG_AN','CH4_TAG_NA'],\
                    temporal_res = 30)
    gen_hemco_diagn(run_path = './',\
                    diagn_name = 'HEMCO_Diagn.rc',\
                    member_start =1,\
                    member_end = 40,\
                    n_regs = 20,\
                    n_sous = 2,\
                    sous_name = ['CH4_TAG_AN','CH4_TAG_NA'])  
                    
    gen_hemco_config(temp_path='./temp/',\
                     run_path ='./',\
                     met_path ='/scratch/local/shzhu/GC_DATA/ctm/GEOS_4x5/GEOS_FP',\
                     emis_path ='/exports/csce/datastore/geos/users/s2008420/',\
                     emis_file_name = ['emi_anthro_s.nc','emi_nature_s.nc'],\
                     emis_var_name = ['Anthr' , 'Nature'],\
                     emis_time = ['2012/1/1/0','2015/1-12/1/0'],\
                     emis_unit = 'molec/cm2/s',\
                     emis_rank = 'xy',\
                     member_start =1,\
                     member_end = 40,\
                     tmpfile = 'HEMCO_Config.temp',\
                     newfile = 'HEMCO_Config.rc',\
                     n_regs = 20,\
                     n_sous = 2,\
                     sous_name = ['CH4_TAG_AN','CH4_TAG_NA'])
    gen_history(temp_path='./temp/',\
                run_path= './',\
                tmpfile  = 'HISTORY.temp',\
                newfile  = 'HISTORY.rc',\
                tm_res= 30,\
                res_name = 'GEOSChem.Restart.%y4%m2%d2_%h2%n2z.nc4')