##############################
## description: gosat & siberia in 2020 using climatological OH field

#----------------------------
## inversion switches 
#----------------------------
use_time_cor = False              # temporal correlation using in prior pertubations(assim_step.py)
use_fixed_mod_err = True         # fixed model error using in calculating observation errors(construct_obs_err.py)
use_full_r = False               # spatial correlation using in observation errors(construct_obs_err.py)
use_sparse = False                # sparse matrix using in calculating inversion parameters(do_assim.py)
add_rnd_err = False              # add random error to observation results(assim_step.py)
do_update = True                 # inversion calculation 
full_rerun = True                # rerun ctm at every step
sampl_input = True               # exited sampling result
#tccon_input = True               # exited tccon sampling result
#----------------------------
## diagnosis switches
#----------------------------
daily_prior_output = False       
daily_para_output = False
cumulative_para_output = False
step_para_output = True
daily_posterior_output = False
daily_conc_output = True
daily_updata_en_bounds = False
daily_updata_sgl_bounds = False
posterior_flux_diag = True
update_en_bounds = False
update_sgl_bounds = False



#----------------------------
## observation and model output path
#----------------------------
root = '/home/shzhu/enkf_ch4/'
sat_path = '/data/GOSAT/GOSAT-XCH4/v9.0/gosat_v9/'
tccon_path = '/data/shzhu/tccon_daily/'
#sp_diagn_path = '/data/shzhu/rerun_ch4/grid_4x5/'
#sp_diagn_path = '/data/shzhu/rerun_ch4/grid_2x25/AN_NA/' # NA_new_land, NA_new, NA_new_siberia, NA_land_siberia
sp_diagn_path = '/data/shzhu/rerun_ch4/grid_2x25/OH_test/two_step_OH/' # clim_OH #full_OH # two_step_OH
#sp_diagn_path = '/data/shzhu/org_ch4/global/'

en_diag = '/data/shzhu/tagch4_lf/'
en_diagn_path = en_diag + 'ND51/'
#sampl_file = '/data/shzhu/org_ch4/global/gosat_samp/gosat_sampl'
#sampl_filepath = '/data/shzhu/org_ch4/global_2x25/sampl_rlt/GOSAT_AN_NA/YYYY/'
sampl_filepath = '/data/shzhu/org_ch4/global_total_lf/sampl_rlt/GOSAT/YYYY/'
#tccon_file = '/data/shzhu/global_inv/sampl_rlt/TCCON/tccon_sampl'
new_hm0 = True
new_sampl_filepath = '/data/shzhu/global_inv/grid_2x25/flask_insitu/gosat_samp_rlt/YYYY/'
sampl_filename = 'gosat_sampl.YYYYMMDD.nc'

sep_en = True
############################ important ########################################!!!!!!!!!!!
entry_table = '/home/shzhu/enkf_ch4/ens_pos/ens_pos_2020.dat'#'ens_pos_2020.dat
# entry_table = '/home/shzhu/enkf_ch4/ens_pos/ens_pos_2010_2020.EN001-EN037.dat'
# entry_table = '/home/shzhu/enkf_ch4/ens_pos/ens_pos_2018.5-2020.dat'
############################ important ########################################!!!!!!!

# em_diag = en_diagn_path + 'diaginfo_em.dat'
# em_tra1 = en_diagn_path + 'tracerinfo_EN0001-EN0037.dat'
# em_tra2 = en_diagn_path + 'tracerinfo_EN0038-EN0082.dat'
# tra_num1 = 37
# tra_num2 = 45


prior_err_input = True
err_file = '/data/shzhu/org_ch4/global_2x25/prior_cov_liang_new/prior_err_cor.YYYYDddd.nc'
#----------------------------
## Boundaries path
#----------------------------
sgl_bounds_path = sp_diagn_path + 'Boundary/'
new_sgl_bounds_path = sp_diagn_path + 'Mod_Boundary/'
en_bounds_path = en_diagn_path +'Boundary/'
new_en_bounds_path = en_diagn_path + '/Mod_Boundary/'

#----------------------------
## diagnosis output path
#----------------------------

inv_diag_path = '/data/shzhu/global_inv/grid_2x25/OH_test/two_step_OH/'
# Daily prior state before inversion
daily_state_path = inv_diag_path + 'Daily_prior_state/'
# Daily inversion parameters
inv_daily_outpath = inv_diag_path + 'Daily_inv_pars/'
# cumulative inversion parameters 
inv_cum_outpath = inv_diag_path + 'Step_cum_inv_pars/'
# Daily posterior flux and its error
inv_daily_state_outpath = inv_diag_path + 'Daily_state/'
inv_daily_flux_outpath = inv_diag_path + 'Daily_flux/'
# Daily posterior concentration and its error
inv_daily_conc_outpath = inv_diag_path + 'Daily_conc/'
# Step parameters
inv_step_outpath = inv_diag_path + 'Step_state/'
# inversion parameters during the whole period
whole_state_path = inv_diag_path + 'Whole_state/'
# posterior flux during the whole period
whole_inv_flux_path = inv_diag_path + 'Whole_flux/'
# Daily tccon validation result 
tccon_valid_outpath = inv_diag_path + 'tc_valid'

#----------------------------
## assimilation settings
#----------------------------

# initial information
st_yyyy = 2020                      # start year
st_doy =  1                        # start day for the year
step_start = 0                      # start step
step_end = 12                       # end step

# update and rerun information
rerun_step = 1                       # rerun ctm step
update_step = 0                      # inversion calculation step 

# inversion parameters
temp_res = 'mon'  # it can be fixed or flexible(with modify main code)
inv_lag_window = 3
ntracers = 186
nbias = 0

#----------------------------
## prior pertubations settings
#----------------------------
xtm_enlarge_factor = 1.0
tcor_factor = 1.0                   # temporal correlation parameters
tcor_len = 1.0                      # temporal correlation parameters
pri_err = 0.7
#----------------------------
## Obs error settings
#----------------------------
#use_fixed_mod_err = True
#use_full_r = True
model_err = 0.0
cor_length = 600
max_cor = 10
obs_err_scale = 0.5

#----------------------------
## Model error settings(CTM)
#----------------------------
# inversion error statistics
mod_bias = 0
bias_err = 0
xnorm = 1.0
xref = 1.0


#----------------------------
## restart stratospheric zonal bias correction
#----------------------------
lat_bias_cor = False
# ace_file =  '/data/shzhu/measurements/ACE-FTS/ace_ch4_sea_height.nc'
# ace_file =  '/data/shzhu/measurements/ACE-FTS/ace_ch4_sea_fac_height.nc'
# ace_file =  '/data/shzhu/measurements/ACE-FTS/ace_ch4_sea_pres_new.nc'
# ace_file = '/data/shzhu/measurements/ACE-FTS/v4.1/gc_cor/ace_ch4_clim_test.nc' ## old fields with new growth factors
ace_file = '/data/shzhu/measurements/ACE-FTS/v4.1/gc_cor/ace_ch4_clim.nc'## new fields
pres_file = '/data/shzhu/org_ch4/global_2x25/met/GEOSChem.LevelEdgeDiags.YYYYMMDD_0000z.nc4'

#----------------------------
## rerun geos-chem settings
#----------------------------
grid = '2x25' #'4x5'
# Config settings for geos-chem
data_path = '/data/geos-chem/ctm'  # input.geos
met_path =  data_path  + '/GEOS_2x2.5/MERRA2'  #HEMCO_Config.rc
#met_path =  data_path  + '/GEOS_4x5/MERRA2'
hem_path =  data_path + '/HEMCO' #path for HEMCO FILES 
emis_path = '/data/geos-chem/emis_ch4/'
emis_file = '/data/geos-chem/emis_ch4/CH4_EMIS_2x25_Liang_2010_2020.nc'
mask_file = '/data/geos-chem/mask/global_mask_93.2x25.nc'

# input and output file path for rerun 
rerun_path = root + 'rerun_geos_2/'
rerun_temp_path = rerun_path +'temp/' 
rerun_config_path = rerun_path +'Config_files/'
rerun_diagn_path = sp_diagn_path
sgl_mask_file = rerun_temp_path + 'sgl_mask.nc'

# default_res_file = rerun_temp_path+'GEOSChem.Restart.20100101_0000z.nc4'
default_res_file = '/data/fuyu/testdata/ForZhu/GCout.restart.CH4_OH.20200101_0000z.nc'
# default_res_file = '/data/shzhu/rerun_ch4/grid_2x25/Prior_lf_2/Mod_Restart/Restart.ST100.EN0001-EN0002.20180501_0000z.nc4'
# default_res_file = '/data/shzhu/org_ch4/global_2x25/Restart/Restart.20100101_0000z.nc4'
sgl_res_path = rerun_diagn_path +'Restart/'
#en_res_path = '/data/shzhu/tagch4/Restart/'
new_res_path = rerun_diagn_path +'Mod_Restart/'


    

# for HEMCO_Config.RC
n_sous = 2
n_regs = 93
sous_name = ['CH4_TAG_AN','CH4_TAG_NA']
#emis_file_name = ['CH4_EMIS_new_version_2010s.nc','CH4_EMIS_new_version_2010s.nc']
emis_file_name = ['CH4_EMIS_2x25_Liang_2010_2020.nc', 'CH4_EMIS_2x25_Liang_2010_2020.nc']
emis_var_name = ['CH4_AN' , 'CH4_NA']
emis_time = ['2009-2025/1-12/1/0','2009-2025/1-12/1/0']
emis_unit = 'kg/m2/s'
emis_rank = 'xy'

# update restart files settings
#gout_res_path = rerun_path + 'Geos_Res_files/'
#mdf_res_path = rerun_path + 'Modif_Res_files/'