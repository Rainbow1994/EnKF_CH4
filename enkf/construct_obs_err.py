# %load construct_obs_err.py
from math import radians, cos, sin, asin, sqrt
from numpy import *
import numba as nb
import geos_chem_def as gcdf
import time_module as tm
import sample_ch4_gosat as scg
import data_quality_ctl as dqc

@nb.jit(nopython=True) 
def haversine(lon1, lat1, lon2, lat2): # (Decimal Degrees)

    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    
    # Caculate the mean radius of earth between two grids
    lat0 = [radians(lat1), radians(lat1)]
    r_e = 6378.137   # Radius at sea level at equator
    r_p = 6356.752   # Radius at poles
    r_lat = []
    for lat in lat0:
        
        r_lat.append(sqrt(((r_e**2.0*cos(lat))**2.0+(r_p**2.0*sin(lat))**2.0)/((r_e*cos(lat))**2.0+(r_p*sin(lat))**2.0)) )
    r = (r_lat[0] +r_lat[1])*0.5  # mean radius
    
    #r = float32(6371.0088)
 
    # Convert decimal degrees to radians
 #   lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    
    # haversine formula
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat*0.5)**2 + cos(lat1) * cos(lat2) * sin(dlon*0.5)**2
    b = sqrt(a)
    c = 2.0 * asin(b) 
    
    # return radius (unit:km)
    radius = c * r 
    
    return radius


@nb.njit
def get_circle_distance_xy(lon, lat):
    
    """
    Calculate the circle distance of two points in lon_lat grids 
    matrics on the earth (specified in decimal degrees)
    """
    
    row = len(lon)
    col = len(lat)
    dist=zeros((row, col), float32)
    for irow, (lon1, lat1) in enumerate(zip(lon, lat)):
        for icol, (lon2 ,lat2) in enumerate(zip(lon, lat)):
            dist[irow,icol] = haversine(lon1, lat1, lon2, lat2)
    return dist


def get_mod_err(olon, olat, devflnm='obs_dev.nc', def_err=1.0):
    '''
     observation deviation(inversion model error)
     
    '''
    nobs=size(olon)
    
    varnames=['longitude', 'latitude', 'dev', 'std', 'count']
    grdlon, grdlat, obs_dev, obs_std, obs_count=ofb.ncf_read(devflnm, varnames)
    lon_idx=searchsorted(grdlon, olon)
    ngrdlon=size(grdlon)
    lon_idx=where(lon_idx>=ngrdlon, ngrdlon-1, lon_idx)
    
    lat_idx=searchsorted(grdlat, olat)
    ngrdlat=size(grdlat)
    lat_idx=where(lat_idx>=ngrdlat, ngrdlat-1, lat_idx)
    mod_obs=zeros(nobs, float)

    for iobs in range(nobs):
        std_dev=obs_std[lon_idx[iobs],lat_idx[iobs]]
        mean_dev=obs_dev[lon_idx[iobs],lat_idx[iobs]]
        if (std_dev<0.0):
            mod_obs[iobs]=def_err
        else:
            mod_obs[iobs]=2.5*sqrt(std_dev**2+mean_dev**2)

    return mod_obs

def construct_obs_err_cov(lon, lat, r, cor_length, max_cor, do_debug=False):
    """ 
    generat the error correlation between observations 
        
    """
    
    nobs = len(r)
    err_cov = zeros([nobs, nobs], float)
    dist = get_circle_distance_xy(lon, lat)
    factor = dist/cor_length
    usd_idx = where(factor < max_cor)
    err_cov[usd_idx] = exp(-factor[usd_idx])
    err = diag(r)
#    err = squeeze(err)
#    err = sqrt(err)
    
    obs_r = dot(err_cov, err)
    obs_r = dot(obs_r, err_cov)
#    obs_r = dot(obs_r, transpose(obs_r))
    #obs_r = 0.5 *(transpose(obs_r) + obs_r)
    
    return obs_r

def get_obs_err(lon,\
                lat,\
                oerr,\
                ydev,\
                model_r,\
                use_fixed_mod_err = True,\
                use_full_r = True):
    '''
    adjust original observation error
     
    '''
    obs_num = len(oerr)
    print('obs_nums', obs_num)
    if (use_fixed_mod_err):
        if (obs_num>1):
            all_mod_err=array(oerr)
            all_mod_err[:]=gcdf.model_err*gcdf.model_err
        else:
            all_mod_err=gcdf.model_err*gcdf.model_err
                    
    else:
        all_mod_err = get_mod_err(olon, olat, devflnm, def_err=1.0)
        # print 'max min model_err', max(model_err), min(model_err)
        all_mod_err=all_mod_err*all_mod_err
            
    obs_err_scale =gcdf.obs_err_scale                
    all_yobs_err=(obs_err_scale* oerr)**2 + all_mod_err
    
#    ydev = where(ydev > 25.0 * (all_yobs_err + model_r), 10 * ydev, ydev)
    
    all_yobs_err = all_yobs_err + 0.001*ydev
    
    
    if (use_full_r):
#        err = sqrt(all_yobs_err)
        obs_r = construct_obs_err_cov(lon = lon,\
                                      lat = lon,\
                                      r = all_yobs_err,\
                                      cor_length = gcdf.cor_length,\
                                      max_cor = gcdf.max_cor,\
                                      do_debug = False)
    else:
        obs_r = diag(all_yobs_err)
        
 #   print('obs_r', obs_r)
    
    return obs_r

    
if (__name__=="__main__"):
    
    
    doy = 28
    yyyy = 2015
    cur_step = 0
    samp_data = scg.get_hm(cur_step, doy, yyyy) 
    result = dqc.data_quality_ctl(samp_data)
    
    ydev = 0.0
    oerr = result['oerr'].values
    lons = result['lon'].values
    lats = result['lat'].values
    r = get_obs_err(lons,lats,oerr, ydev,\
                model_r = 1.2,\
                use_fixed_mod_err = True,\
                use_full_r = True)
    
    
    ## exec time : 42.4s