
from math import radians, cos, sin, asin, sqrt
from numpy import *
import pandas as pd
from datetime import datetime
import xarray as xr
import time_module as tm
import numba as nb
import geos_chem_def as gcdf

@nb.jit(nopython=True) 
def a_haversine(lon1, lat1, lon2, lat2): # (Decimal Degrees)

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
    if b >1:
        b = 1
    c = 2.0 * asin(b) 
    
    # return radius (unit:km)
    radius = c * r 
    
    return radius

@nb.njit#(nopython=True) 
def get_circle_distance_xy(lonx, latx, lony, laty):
    
    """
    Calculate the circle distance of two points in lon_lat grids 
    matrics on the earth (specified in decimal degrees)
    """
    
    row = len(lonx)
    col = len(lony)
    dist=zeros((row, col), float32)
    for irow, (lon1, lat1) in enumerate(zip(lonx, latx)):
        for icol, (lon2 ,lat2) in enumerate(zip(lony, laty)):
            dist[irow,icol] = a_haversine(lon1, lat1, lon2, lat2)
    return dist

def spatial_err_cov(emis_fn, mask_fn, yyyy, doy,\
                    cor_length = 500.0):
    """ 
    generate the error correlation between regions
        lon, lat: longitude and latitude list of regions 
        reg_idx_lst: list of cells in each regions
        reg_err_lst: list of flux error in each regions
        cor_length: the parameter to construct the regional
                    error correlation 
    """
    emis_reg =  xr.open_dataset(emis_fn)
    mask = xr.open_dataset(mask_fn)
    
    time = tm.doy_to_utc(doy, 0, yyyy)
    emis_reg = emis_reg.transpose('lon', 'lat', 'time')
    emis_reg = emis_reg.sel(time = time)
    
    mask = mask.rename({'lev':'reg'})
    mask = mask.transpose('lon','lat','reg')
    emis_reg = emis_reg.assign_coords({'lon': mask['lon'].values,
                                       'lat': mask['lat'].values})
    
    emis_reg = xr.where(emis_reg == -inf, 0,emis_reg )
    emis_err = emis_reg* mask
    emis_err = emis_err.fillna(0)
    n_sous = gcdf.n_sous 
    n_regs = gcdf.n_regs 
    sous_name =gcdf.emis_var_name
    reg_idx_lst = list()
    reg_err_lst = list()

    for i_sous in sous_name:
        emis_err = emis_reg[i_sous] * mask['MASK']
        for i_reg in range(n_regs):
            used_idx = where( mask['MASK'].isel(reg = i_reg) > 0.0)
            reg_idx_lst.append(used_idx)
            data_sel = emis_err.isel(reg = i_reg).values
            reg_err_lst.append(data_sel[used_idx])

    
    lon_m, lat_m=meshgrid(emis_reg['lon'].values, emis_reg['lat'].values)
    lon_m = transpose(lon_m)
    lat_m=transpose(lat_m)
    nreg=len(reg_idx_lst)
    err_cov=zeros((nreg, nreg), dtype = float32)
    
    for i in range(nreg):
         # index for each cell in region i
        px = reg_idx_lst[i]
        wx = reg_err_lst[i]
        
        # flux error at each cell in region i
        lonx=lon_m[px]
        latx=lat_m[px]
        
        err_cov[i, i]=1.0
        for j in range(i+1, nreg):
            # index for each cell in region j
            py=reg_idx_lst[j]
            wy=reg_err_lst[j]
            # flux error at each cell in region j
            lony=lon_m[py]
            laty=lat_m[py]

            # distance (km) between cells in range i and j
            dist = get_circle_distance_xy(lonx, latx, lony, laty)



            schur_factor=zeros(shape(dist), float)
            dist = divide(dist , cor_length)
     
            usd_idx=where(dist<4.0) # only close points are used  
            #print('hhhh', usd_idx)
            schur_factor[usd_idx]=exp(-dist[usd_idx])
            schur_factor=dot(schur_factor, wy)
            schur_factor=dot(wx, schur_factor)

            errx=sum(wx)
            erry=sum(wy)
            if abs(errx *erry) > 0:
                err_cov[i, j]=0.6*sqrt(abs(schur_factor/(errx*erry)))
                err_cov[j, i]=err_cov[i, j]
                
        
    return err_cov


@nb.jit
def temporal_err_cov(n_reg = 20,\
                     n_period = 3,\
                     tcor_factor = 1.0,\
                     tcor_length = 1.0):
    
    """ 
    generate temporal correlation  into prior error matrix
    
    temporal correlated DX is generated as  
     
         DX=
            |  DX0     0   |       | DX0  0   |  | I       0  |
            |  cor_dx  DX1 |   ==  |  0   DX1 |  | tc_xtm, I  |
        
    
    as a results, the error covariance for newly included x1 will be
    B=DX1*(I+tc_xtm*tc_xtm^T)*DX1^T, so a rescaling is needed 
    """ 
    pos = 0
    sum_row_factor=zeros(n_reg, float)
    sum_row_factor[:]=1.0
    nx = n_period * n_reg
    tc_xtm =zeros(nx , float)
    
    # tc_xtm is inform of 
    # |a, 0, 0, ..., 0_nx,a2, 0, 0, 0, ...,   0_nx, ..., ..., |  
    # |0, a, 0, ..., 0_nx, 0, a2, 0, 0, ...,  0_nx, ... , ... |
    # |0, 0, a, ..., 0_nx, 0. 0,  a2, 0, ..., 0_nx, ... , ... |
    # |--- period 1 ---| -------period 2  -------|---period---|
            
    for iperiod in range(n_period):
        cor_factor=tcor_factor*exp(-(nperiod-iperiod)/tcor_len)
                
        for irow in range(n_reg):
                tc_xtm[irow,pos+irow]=cor_factor
                sum_row_factor[irow]=sum_row_factor[irow]+cor_factor*cor_factor
            
                
        pos = pos + n_reg
        
    sum_row_factor=1.0/sqrt(sum_row_factor)
    
    return sum_row_factor, tc_xtm
    

    
if (__name__=="__main__"):
    
    doy = 1
    yyyy = 2015
    emis_fn = '/exports/csce/datastore/geos/users/s2008420/enkf_ch4_sh/surface_flux/CH4_EMIS_2015.nc'
    mask_fn = gcdf.mask_file
    err_cov = spatial_err_cov(emis_fn, mask_fn, yyyy, doy,cor_length = 500.0)
    print(err_cov)
    
    ## exec time : 2 min 50 s