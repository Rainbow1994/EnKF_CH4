#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
'''
NOAA observation sites information in Siberia 
'''
import numpy as np
import pandas as pd
import os 
from datetime import datetime


AZV = dict({
    'SITE_NAME': 'Azovo',
    'name': 'AZV',
    'lat': 54.7050,
    'lon': 73.0292,
    'height': 110
})

BRZ = dict({
    'SITE_NAME': 'Berezorechka',
    'name': 'BRZ',
    'lat': 56.1458,
    'lon': 84.3319,
    'height': 168
})

DEM = dict({
    'SITE_NAME': 'Demyanskoe',
    'name': 'DEM',
    'lat': 59.7914,
    'lon': 70.8711,
    'height': 75
})

IGR = dict({
    'SITE_NAME': 'Igrim',
    'name': 'IGR',
    'lat': 63.1917,
    'lon': 64.4139,
    'height': 9
})

KRS = dict({
    'SITE_NAME': 'Karasevoe',
    'name': 'KRS',
    'lat': 58.2456,
    'lon': 82.4244,
    'height': 76
})

NOY = dict({
    'SITE_NAME': 'Noyabrsk',
    'name': 'NOY',
    'lat': 63.4292,
    'lon': 75.7800,
    'height': 108
})

SVV = dict({
    'SITE_NAME': 'Savvushka',
    'name': 'SVV',
    'lat': 51.3253,
    'lon': 82.1283,
    'height': 495
})

VGN = dict({
    'SITE_NAME': 'Vaganovo',
    'name': 'VGN',
    'lat': 54.4972,
    'lon': 62.3247,
    'height': 192
})

YAK = dict({
    'SITE_NAME': 'Yakutsk',
    'name': 'YAK',
    'lat': 62.0886,
    'lon': 129.3558,
    'height': 264
})

sites_list = [AZV, BRZ, DEM, IGR, KRS, NOY, SVV, VGN, YAK ]