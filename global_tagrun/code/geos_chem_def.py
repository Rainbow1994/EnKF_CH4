# default directories
import time_module as tm

en_path ='./'
#root = '/home/shzhu/global_tagrun/'
root = '_root_'
temp_path = root + 'temp/'
config_path = root+ 'gc_config/'
#met_path = '/data/geos-chem/ctm/GEOS_2x2.5/MERRA2'  #HEMCO_Config.rc
grid = '_grid_'
met_path = '/data/GC_data/ExtData/GEOS_' + grid + '/MERRA2'  #HEMCO_Config.rc
data_path = '/data/GC_data/ExtData' # input.geos
emis_path ='/data/geos-chem/emis_ch4/' #path of emission datas 
mask_path = '/data/geos-chem/mask/'
#en_diag_path = '/data/shzhu/tagch4/'
en_path = '_diag_'
year = '_year_'
res_path = en_path+ 'Restart/'+ year +'/'
diagn_path = en_path + 'ND51/'+ year +'/'
bound_path =  en_path + '/Boundary/'+ year +'/'
#out_file ='/home/shzhu/global_tagrun/2015/enkf_output/ens_pos.dat' # the model output
#log_file = './enkf.log'
##c/time resultion in day

temp_res='mon' # it can be fixed or flexible(should modified main code)
inv_lag_window=3
nstep=12
nstep_start=0

cfg_pb_file= root + '/ens_cfg/ens_ch4_cfg.'+ year

##c/force the system to re-initialize to be a fixed value if True
new_restart= True
fixed_init_val=1.0e-8
restart_file = temp_path + 'GEOSChem.Restart.20160701_0000z.nc4'
#rerun_now=True
#rerun_date='20130101'
ntracers = 186  # not in use
sub_sous = False
n_sous_loop  = 1


sub_tracers = True

# n_tracer_loop = 1
# mask_file = ['asian_mask_37.01x01.nc']
# tracer_start = [1] #[1, maxtracer+1]
# tracer_end =[37] #[1, maxtracer+1]
# tracer_num = [37]

# n_tracer_loop = 1
# mask_file =['global_res_mask_46.2x25.nc']#['asian_mask_47.2x25.nc'] #['global_res_mask_46.2x25.nc']
# tracer_start =[48]#[1] #[48] #[1, maxtracer+1]
# tracer_end = [93]#[47]#[93] #[1, maxtracer+1]
# tracer_num = [46]#[47]#[46]

# n_tracer_loop = 2
# mask_file = ['asian_mask_37.01x01.nc','global_res_mask_45.01x01.nc']
# tracer_start = [1,38] #[1, maxtracer]
# tracer_end =[37,82] #[1, maxtracer]
# tracer_num = [37,45]

n_tracer_loop = 2
mask_file = ['asian_mask_47.2x25.nc','global_res_mask_46.2x25.nc']
tracer_start = [1,48] #[1, maxtracer]
tracer_end =[47,93] #[1, maxtracer]
tracer_num = [47,46]


##c selected run
select_runs=list(range(_st,end_))
#select_sous=

#def_tracerinfo_file=run_path+"/tracerinfo.dat"
#def_diaginfo_file=run_path+"/diaginfo.dat"


# added by rainbow
n_regs = 22 # used in single run ---- no sub_tracers
n_sous = 1  # used in no sub_sous run, if sub_sous == True n_sous =1
# sous_name = ['CH4_TAG_AN','CH4_TAG_NA'] # should be line with the n_sous
# emis_var_name = ['CH4_AN' , 'CH4_NA']  # should be line with the n_sous
# emis_file_name = ['CH4_EMIS_2x25_v2021.nc','CH4_EMIS_2x25_v2021.nc']  # should be line with the n_sous
# emis_time = ['2009-2025/1-12/1/0','2009-2025/1-12/1/0'] 

# sous_name = ['CH4_TAG_AN','CH4_TAG_NA'] # should be line with the n_sous
# emis_var_name = ['CH4_AN' , 'CH4_NA']  # should be line with the n_sous
# emis_file_name = ['CH4_EMIS_2x25_Liang_2010_2020.nc','CH4_EMIS_2x25_Liang_2010_2020.nc']
# emis_time = ['2009-2025/1-12/1/0','2009-2025/1-12/1/0'] 
# sous_name = ['CH4_TAG_AN','CH4_TAG_FF']
# ## for HEMCO_Config.RC
# emis_file_name = ['CH4_EMIS_AN_NA_FF_2X25_2010s.nc','CH4_EMIS_AN_NA_FF_2X25_2010s.nc'] 
# emis_var_name = ['CH4_AN' , 'CH4_FF']

sous_name = ['CH4_TAG_NA'] # should be line with the n_sous
emis_var_name = ['CH4_NA']  # should be line with the n_sous
emis_file_name = ['CH4_EMIS_2x25_v2021.nc']  # should be line with the n_sous

emis_time = ['2009-2025/1-12/1/0']#,'2009-2025/1-12/1/0']  # should be line with the n_sous
# #emis_time = ['2019/1-12/1/0','2019/1-12/1/0']
emis_unit = 'kg/m2/s'
emis_rank = 'xy'
