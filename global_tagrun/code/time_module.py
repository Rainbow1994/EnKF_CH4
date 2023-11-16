""" routine for time convertion
"""
import numpy
from pylab import *

import datetime
# second_leap=33 # there is a 32 seconds in 2000
second_leap=0
TAI0=datetime.datetime(1985, 1,1,0,0,0)
def_yyyy=2004
def day_of_year(yyyy=2005, mm=8, dd=15):
    """ get the day in a given year 
        Arguments: 
            yyyy, mm, dd: year, month and day 
        Returns:
            days in year
    """
        
    date=datetime.date(yyyy, mm, dd)
    date_ref=datetime.date(yyyy, 1, 1)
    delta=date-date_ref
    ndays=delta.days+1
    del date, date_ref, delta
    return ndays

def tai85_to_utc(tai85):
    """ convert a time in tai85 format to UTC string 
        Arguments: 
            tai85: the seconds since 1985-01-01 00:00:00
        Notes: 
            the time_leap is done in a hard (hand) way 
    """
    
    date_new=TAI0+datetime.timedelta(seconds=(tai85-second_leap))
    utc=str(date_new)
    del date_new
    return utc
def utc_to_tai85(utc):
    """ convert utc string to tai85 format (seconds since 1993-01-01 00:00:00)
        Arguments:
            utc: utc time in yyyy-mm-dd hh:mi:ss
        returns:
            tai85 in seconds
    """

    t1,t2=utc.split()
    syyyy, smm, sdd=t1.split('-')
    shh, smi, ssec=t2.split(':')
    yyyy=int(syyyy)
    mm=int(smm)
    dd=int(sdd)
    hh=int(shh)
    mi=int(smi)
    fsec=float(ssec)
    sec=int(fsec)
    iv_time=datetime.datetime(yyyy, mm, dd, hh, mi, sec)-TAI0
    tai85=3600.0*24.0*iv_time.days+iv_time.seconds+second_leap
    return tai85
def doy_to_tai85(doy, sec=0, yyyy=def_yyyy):
    """ convert day of year to tai85 format (seconds since 1993-01-01 00:00:00)
        Arguments:
            doy, sec, yyyy: day of year, seconds, and year
        returns:
            tai85 in seconds
    """

    date0=datetime.datetime(yyyy, 1, 1, 0,0,0)
    date0=date0+datetime.timedelta(days=doy-1, seconds=sec)
    iv_time=date0-TAI0
    tai85=3600.0*24.0*iv_time.days+iv_time.seconds+second_leap
    del date0
    return tai85
def doy_to_utc(doy, sec=0, yyyy=def_yyyy):
    """ convert day of year to utc string
        Arguments:
            doy, sec, yyyy: day of year, seconds, and year
        returns:
            utc: the time in utc format yyyy-mm-dd hh:mm:ss
    """
    
    date0=datetime.datetime(yyyy, 1, 1, 0,0,0)
    date0=date0+datetime.timedelta(days=doy-1, seconds=sec)
    utc=str(date0)
    del date0
    return utc
def doy_to_time_array(doy, yyyy=2005):
    utc=doy_to_utc(doy, 0, yyyy)
    yyyy, mm,dd, hh, mi, sec=utc_to_time_array(utc)
    return yyyy, mm, dd

def utc_to_time_array(utc):
    """ convert the utc string to yyyy, mm, dd, hh, mi, sec
        Arguments:
            utc: the time in utc format yyyy-mm-dd hh:mm:ss
        returns:
            yyyy, mm, dd, hh, mi, sec
    """
    sd, sh=utc.split(' ')
    syyyy, smm, sdd=sd.split('-')
    shh, smi, ssec=sh.split(':')
    yyyy=int(syyyy)
    mm=int(smm)
    dd=int(sdd)
    hh=int(shh)
    mi=int(smi)
    sec=float(ssec)
    return yyyy, mm, dd, hh, mi, sec
def time_array_to_utc(yyyy,mm, dd, hh=0, mi=0, sec=0):
    """ convert yyyy, mm, dd, hh, mi, sec to utc string
        Arguments:
            yyyy, mm, dd, hh, mi, sec: year, month, day, hour, minute and second
        returns:
            utc: the time in utc format yyyy-mm-dd hh:mm:ss
    """
    sec=int(sec)
    
    syyyy_mm_dd=r'%4.4d-%2.2d-%2.2d' %(yyyy, mm, dd)
    if (sec>=10):
        shh_mi_sec=r'%2.2d:%2.2d:%5.2f' %(hh, mi, sec)
    else:
        shh_mi_sec=r'%2.2d:%2.2d:%4.2f' %(hh, mi, sec)
    return syyyy_mm_dd+' '+shh_mi_sec

    
def get_tau(yyyy, mm, dd, hh=0, mi=0, sec=0):
    iv_time=datetime.datetime(yyyy, mm, dd, hh, mi, sec)-TAI0
    tai85=3600.0*24.0*iv_time.days+iv_time.seconds+second_leap
    return tai85


def next_doy(yyyy_in, doy_in, days=1, return_ymd=False):
    utc=doy_to_utc(doy_in, 0, yyyy_in)
    yyyy, mm, dd, hh, mi, sec=utc_to_time_array(utc)
    date0=datetime.datetime(yyyy, mm, dd, hh,mi,sec)
    date0=date0+datetime.timedelta(days=days, seconds=sec)
    utc=str(date0)
    yyyy, mm, dd, hh, mi, sec=utc_to_time_array(utc)
    doy=day_of_year(yyyy, mm, dd)
    if (return_ymd):
        return yyyy, mm, dd
    else:
        return yyyy, doy
    
def get_ut_time_slot(lt_st, lt_end, day_time_grid, lon, day_length=24.0*3600):

    tshift=day_length*lon/360.
    ut_st=lt_st-tshift
    ut_end=lt_end-tshift
    
    if (ut_end<0.0): # whole in earlier day 
        ut_st=ut_st+day_length
        ut_end=ut_end+day_length

        usd_idx=where((day_time_grid>=ut_st) & (day_time_grid<ut_end))
        usd_idx=squeeze(usd_idx)
        
    elif (ut_st>day_length): # whole in later day 
        ut_st=ut_st-day_length
        ut_end=ut_end+day_length
    
        usd_idx=where((day_time_grid>=ut_st) & (day_time_grid<ut_end))
        usd_idx=squeeze(usd_idx)
    elif (ut_st<0.0): # cross earlier day

        ut_st=ut_st+day_length
        chose1=(day_time_grid>=ut_st) & (day_time_grid<=day_length)
        chose2=(day_time_grid>=0.0) & (day_time_grid<ut_end)
        usd_idx=where(chose1 | chose2)
        usd_idx=squeeze(usd_idx)
    elif (ut_end>day_length): # cross later day
        
        ut_end=ut_end-day_length
        chose1=(day_time_grid>=ut_st) & (day_time_grid<=day_length)
        chose2=(day_time_grid>=0.0) & (day_time_grid<ut_end)
        usd_idx=where(chose1 | chose2)
        usd_idx=squeeze(usd_idx)
    else:
        usd_idx=where((day_time_grid>=ut_st) & (day_time_grid<ut_end))
        usd_idx=squeeze(usd_idx)

        
    return usd_idx
