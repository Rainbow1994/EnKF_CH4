#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
'''
Configuration for gosat observations sampling!!
'''
from datetime import datetime
import sys
sys.path.append('/home/shzhu/enkf_ch4/rundirs/grid_2x25_new/')
import geos_chem_def as gcdf
# import siberia as si
### 
year = 2010
# mon = 2
date_st = datetime(2010,1,1) 
date_end = datetime(2010,1,5) 
measurements = ['gosat', 'noaa']
# measurements = ['tropomi']
n_regs = 93
# emis_vars = ['CH4_AN' , 'CH4_FF', 'CH4_NA']
emis_vars = ['CH4_AN' , 'CH4_NA']
emis_prefix = 'IJ_AVG_S_TAG_'
########################################################################
### gosat observation information output sampling output information
gosat_info = {
    'filepath': '/data01/observation/GOSAT/GOSAT-XCH4/v9.0/gosat_v9/YYYY/',
    'filename': 'UoL-GHG-L2-CH4-GOSAT-OCPR-YYYYMMDD-fv9.0.nc',
    'sampl_path': '/data/shzhu/org_ch4/global_total_lf/sampl_rlt/GOSAT/YYYY/', 
#    'sampl_path': '/data/shzhu/global_inv/sampl_rlt/GOSAT/',
#     'sampl_path': '/data/shzhu/org_ch4/global_4x5/sampl_rlt/GOSAT/YYYY/',

    'sampl_file': 'gosat_sampl.YYYYMMDD.nc'
}
#     'xch4': 'xch4',
#     'sza': 'sensor_zenith_angle',
#     'q_value': 'xch4_quality_flag',
#     'pres_lev': 'pressure_levels',
#     'ave_kernel': 'xch4_averaging_kernel',
#     'pres_weight': 'pressure_weight',
#     'apri_profile': 'ch4_profile_apriori',
#     'dims': ['n', 'm'],
#     'lon': 'longitude',
#     'lat': 'latitude'
# }
###########################################################################
### Tropomi observation information output sampling output information
tropomi_info = {
    'filepath': '/data/TROPOMI/CH4/L2_Daily/YYYY/MM/',
    'filename': 'TROPOMI_XCH4_Daily_YYYYMMDD.nc',
    'sampl_path': '/data/shzhu/org_ch4/global_total_lf/sampl_rlt/Tropomi/YYYY/',
    'sampl_file': 'tropomi_sampl.YYYYMMDD.nc',
    'grid_file': '/data/geos-chem/ctm/GEOS_2x2.5/MERRA2/2016/01/MERRA2.20160101.CN.2x25.nc4',
    'dlon' : 2.0,
    'dlat' : 2.5
}

###########################################################################
### NOAA observation information output sampling output information
noaa_info = {
    'filepath': '/data/shzhu/NOAA_CH4/202110/obspack_ch4_1_GLOBALVIEWplus_v4.0_2021-10-14/field/',
#     'filename': 'siberia_ch4_obs.nc',
    'filename': 'noaa_siberia_surface_inversion.nc',#'noaa_surface_inversion.nc',
#     'sites_list': si.sites_list,
    'sampl_path': '/data/shzhu/org_ch4/global_2x25/sampl_rlt/NOAA/YYYY/',
    'sampl_file': 'noaa_surface_sampl.YYYYMMDD.nc',
    'project': ['flask', 'insitu','pfp','tower','siberia']#['insitu']#['flask','pfp','tower']#['flask', 'insitu','pfp','tower']
}

###########################################################################
### TCCON observation information output sampling output information
tccon_info = {
    'filepath': '/data/shzhu/NOAA_CH4/siberia_ch4/field/',
    'filename': 'siberia_ch4_obs.nc',
    'sampl_path': '/data/shzhu/org_ch4/global_2x25/sampl_rlt/TCCON_AN_NA/YYYY/',
    'sampl_file': 'noaa_sampl.YYYYMMDD.nc'
}



##########################################################################
### single run output path
# sp_diagn_info = {
#     'filepath': '/data/shzhu/org_ch4/global_2x25/ND51/',
#     'filename': 'ts_satellite.YYYYMMDD.bpch',
#     'diaginfo_file' : '/data/shzhu/org_ch4/global_2x25/diaginfo_sp.dat',
#     'tracerinfo_file' : '/data/shzhu/org_ch4/global_2x25/tracerinfo_sp.dat',
#     'var_prefix': 'IJ_AVG_S_CH4'

# }
sp_diagn_info = {
    'filepath': gcdf.sp_diagn_path + '/ND51/',
#     'filename': 'ts_satellite.YYYYMMDD.bpch',
    'filename': 'ts_satellite.ST000.EN0001-EN0002.YYYYMMDD.bpch',
    'diaginfo_file' : '/data/shzhu/rerun_ch4/grid_2x25/Prior_lf/diaginfo_sp.dat',
    'tracerinfo_file' : '/data/shzhu/rerun_ch4/grid_2x25/Prior_lf/tracerinfo_sp.dat',
    'var_prefix': 'IJ_AVG_S_CH4'

}
### ensemble items table
# entry_table = '/home/shzhu/enkf_ch4/ens_pos/ens_pos_2018.5-2020.dat'
entry_table = '/home/shzhu/enkf_ch4/ens_pos/ens_pos_2010_2020.EN001-EN037.dat'#ens_pos_2014_2015.EN001-EN037.dat'

## ensemble run information
en_file1_info = {
    'filepath': '/data/shzhu/tagch4/ND51/',
    'filename': 'ts_satellite.ST000.EN0001-EN0047.YYYYMMDD.bpch',
    'diaginfo_file' : '/data/shzhu/tagch4/ND51/diaginfo_em.dat',
    'tracerinfo_file' : '/data/shzhu/tagch4/ND51/tracerinfo_EN0001-EN0047.dat',
    'var_prefix': ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA'],
    'trac_num': 47,
    'trac_st':1
}


en_file2_info = {
    'filepath': '/data/shzhu/tagch4/ND51/',
    'filename': 'ts_satellite.ST000.EN0048-EN0093.YYYYMMDD.bpch',
    'diaginfo_file' : '/data/shzhu/tagch4/ND51/diaginfo_em.dat',
    'tracerinfo_file' : '/data/shzhu/tagch4/ND51/tracerinfo_EN0048-EN0093.dat',
    'var_prefix': ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA'],
    'trac_num': 46,
    'trac_st':48
}

en_info_list = [en_file1_info, en_file2_info]
# en_file1_info = {
#     'filepath': '/data/shzhu/tagch4/ND51/',
#     'filename': 'ts_satellite.ST000.EN0001-EN0037.YYYYMMDD.bpch',
#     'diaginfo_file' : '/data/shzhu/tagch4/ND51/diaginfo_em.dat',
#     'tracerinfo_file' : '/data/shzhu/tagch4/ND51/AN/tracerinfo_EN0001-EN0037.dat',
# #    'var_prefix': ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA'],
#     'var_prefix': ['IJ_AVG_S_TAG_AN'],
#     'trac_num': 37,
#     'trac_st':1
# }

# en_file3_info = {
#     'filepath': '/data/shzhu/tagch4_2x25_AN/ND51/',
#     'filename': 'ts_satellite.ST000.EN0001-EN0037.YYYYMMDD.bpch',
#     'diaginfo_file' : '/data/shzhu/tagch4/ND51/diaginfo_em.dat',
# #    'tracerinfo_file' : '/data/shzhu/tagch4/ND51/tracerinfo_EN0001-EN0037.dat',
#     'tracerinfo_file' : '/data/shzhu/tagch4_2x25_AN/ND51/tracerinfo_EN0001-EN0037.dat',
# #    'var_prefix': ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA'],
#      'var_prefix': ['IJ_AVG_S_TAG_NA'],
#     'trac_num': 37,
#     'trac_st':1
# }

# en_file2_info = {
#     'filepath': '/data/shzhu/tagch4/ND51/',
#     'filename': 'ts_satellite.ST000.EN0038-EN0082.YYYYMMDD.bpch',
#     'diaginfo_file' : '/data/shzhu/tagch4/ND51/diaginfo_em.dat',
#     'tracerinfo_file' : '/data/shzhu/tagch4/ND51/AN/tracerinfo_EN0038-EN0082.dat',
# #     'var_prefix': ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA'],
#     'var_prefix': ['IJ_AVG_S_TAG_AN'],
#     'trac_num': 45,
#     'trac_st':38
# }

# en_file4_info = {
#     'filepath': '/data/shzhu/tagch4_2x25_AN/ND51/',
#     'filename': 'ts_satellite.ST000.EN0038-EN0082.YYYYMMDD.bpch',
#     'diaginfo_file' : '/data/shzhu/tagch4/ND51/diaginfo_em.dat',
# #     'tracerinfo_file' : '/data/shzhu/tagch4/ND51/tracerinfo_EN0038-EN0082.dat',
#     'tracerinfo_file' : '/data/shzhu/tagch4_2x25_AN/ND51/tracerinfo_EN0038-EN0082.dat',
# #     'var_prefix': ['IJ_AVG_S_TAG_AN', 'IJ_AVG_S_TAG_NA'],
#     'var_prefix': ['IJ_AVG_S_TAG_NA'],
#     'trac_num': 45,
#     'trac_st':38
# }

# en_info_list = [en_file1_info, en_file2_info, en_file3_info, en_file4_info]
#en_info_list = [en_file3_info, en_file4_info]