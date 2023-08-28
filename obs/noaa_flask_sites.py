#! /home/shzhu/apps/miniconda3/envs/idp/bin/python
'''
NOAA observation sites information in Siberia 
'''
import numpy as np
import pandas as pd
import os 
from datetime import datetime

### surface sites location ####

####################################################
###################### flask #######################

ABP = dict({
	'SITE_NAME': 'Arembepe, Bahia',
	'name': 'ABP',
	'lat': -12.77,
	'lon': -38.17,
	'height': 1.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ALT = dict({
	'SITE_NAME': 'Alert, Nunavut',
	'name': 'ALT',
	'lat': 82.4508,
	'lon': -62.5072,
	'height': 185.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ALT = dict({
	'SITE_NAME': 'Alert, Nunavut',
	'name': 'ALT',
	'lat': 82.4508,
	'lon': -62.5072,
	'height': 185.0,
	'intake_hgt': 10.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


AMS = dict({
	'SITE_NAME': 'Amsterdam Island',
	'name': 'AMS',
	'lat': -37.7983,
	'lon': 77.5378,
	'height': 55.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


AMT = dict({
	'SITE_NAME': 'Argyle, Maine',
	'name': 'AMT',
	'lat': 45.0345,
	'lon': -68.6821,
	'height': 53.0,
	'intake_hgt': 107.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


AMY = dict({
	'SITE_NAME': 'Anmyeon-do',
	'name': 'AMY',
	'lat': 36.5389,
	'lon': 126.3295,
	'height': 47.0,
	'intake_hgt': 40.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ARA = dict({
	'SITE_NAME': 'Arcturus, Queensland',
	'name': 'ARA',
	'lat': -23.8587,
	'lon': 148.4746,
	'height': 175.0,
	'intake_hgt': 10.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ASC = dict({
	'SITE_NAME': 'Ascension Island',
	'name': 'ASC',
	'lat': -7.9667,
	'lon': -14.4,
	'height': 85.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ASK = dict({
	'SITE_NAME': 'Assekrem',
	'name': 'ASK',
	'lat': 23.2625,
	'lon': 5.6322,
	'height': 2710.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


AVI = dict({
	'SITE_NAME': 'St. Croix, Virgin Islands',
	'name': 'AVI',
	'lat': 17.75,
	'lon': -64.75,
	'height': 3.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


AZR = dict({
	'SITE_NAME': 'Terceira Island, Azores',
	'name': 'AZR',
	'lat': 38.766,
	'lon': -27.375,
	'height': 19.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


BAL = dict({
	'SITE_NAME': 'Baltic Sea',
	'name': 'BAL',
	'lat': 55.35,
	'lon': 17.22,
	'height': 3.0,
	'intake_hgt': 25.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


BHD = dict({
	'SITE_NAME': 'Baring Head Station',
	'name': 'BHD',
	'lat': -41.4083,
	'lon': 174.871,
	'height': 85.0,
	'intake_hgt': 10.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


BKT = dict({
	'SITE_NAME': 'Bukit Kototabang',
	'name': 'BKT',
	'lat': -0.202,
	'lon': 100.318,
	'height': 845.0,
	'intake_hgt': 30.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


BME = dict({
	'SITE_NAME': 'St. Davids Head, Bermuda',
	'name': 'BME',
	'lat': 32.368,
	'lon': -64.6476,
	'height': 12.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


BMW = dict({
	'SITE_NAME': 'Tudor Hill, Bermuda',
	'name': 'BMW',
	'lat': 32.2647,
	'lon': -64.8788,
	'height': 30.0,
	'intake_hgt': 30.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


BRW = dict({
	'SITE_NAME': 'Barrow Atmospheric Baseline Observatory',
	'name': 'BRW',
	'lat': 71.323,
	'lon': -156.6114,
	'height': 11.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


BSC = dict({
	'SITE_NAME': 'Black Sea, Constanta',
	'name': 'BSC',
	'lat': 44.1776,
	'lon': 28.6647,
	'height': 0.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CBA = dict({
	'SITE_NAME': 'Cold Bay, Alaska',
	'name': 'CBA',
	'lat': 55.21,
	'lon': -162.72,
	'height': 21.34,
	'intake_hgt': 3.66,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CFA = dict({
	'SITE_NAME': 'Cape Ferguson, Queensland',
	'name': 'CFA',
	'lat': -19.28,
	'lon': 147.057,
	'height': 2.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CGO = dict({
	'SITE_NAME': 'Cape Grim, Tasmania',
	'name': 'CGO',
	'lat': -40.683,
	'lon': 144.69,
	'height': 94.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CGO = dict({
	'SITE_NAME': 'Cape Grim, Tasmania',
	'name': 'CGO',
	'lat': -40.683,
	'lon': 144.69,
	'height': 94.0,
	'intake_hgt': 70.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CHR = dict({
	'SITE_NAME': 'Christmas Island',
	'name': 'CHR',
	'lat': 1.7,
	'lon': -157.1518,
	'height': 0.0,
	'intake_hgt': 2.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CIB = dict({
	'SITE_NAME': 'Centro de Investigacion de la Baja Atmosfera (CIBA)',
	'name': 'CIB',
	'lat': 41.81,
	'lon': -4.93,
	'height': 845.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CMO = dict({
	'SITE_NAME': 'Cape Meares, Oregon',
	'name': 'CMO',
	'lat': 45.478,
	'lon': -123.9689,
	'height': 30.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CPA = dict({
	'SITE_NAME': 'Charles Point, Darwin',
	'name': 'CPA',
	'lat': -12.417,
	'lon': 130.567,
	'height': 3.0,
	'intake_hgt': 6.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CPT = dict({
	'SITE_NAME': 'Cape Point',
	'name': 'CPT',
	'lat': -34.3523,
	'lon': 18.4891,
	'height': 230.0,
	'intake_hgt': 30.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CRI = dict({
	'SITE_NAME': 'Cape Rama',
	'name': 'CRI',
	'lat': 15.08,
	'lon': 73.83,
	'height': 60.0,
	'intake_hgt': 6.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CRZ = dict({
	'SITE_NAME': 'Crozet Island',
	'name': 'CRZ',
	'lat': -46.4337,
	'lon': 51.8478,
	'height': 197.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


CYA = dict({
	'SITE_NAME': 'Casey, Antarctica',
	'name': 'CYA',
	'lat': -66.283,
	'lon': 110.517,
	'height': 47.0,
	'intake_hgt': 8.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


DSI = dict({
	'SITE_NAME': 'Dongsha Island',
	'name': 'DSI',
	'lat': 20.6992,
	'lon': 116.7297,
	'height': 3.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


EIC = dict({
	'SITE_NAME': 'Easter Island',
	'name': 'EIC',
	'lat': -27.1597,
	'lon': -109.4284,
	'height': 47.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ESP = dict({
	'SITE_NAME': 'Estevan Point,  British Columbia',
	'name': 'ESP',
	'lat': 49.383,
	'lon': -126.5441,
	'height': 7.0,
	'intake_hgt': 40.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


GMI = dict({
	'SITE_NAME': 'Mariana Islands',
	'name': 'GMI',
	'lat': 13.386,
	'lon': 144.656,
	'height': 0.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


GOZ = dict({
	'SITE_NAME': 'Dwejra Point, Gozo',
	'name': 'GOZ',
	'lat': 36.048,
	'lon': 14.889,
	'height': 1.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


GPA = dict({
	'SITE_NAME': 'Gunn Point',
	'name': 'GPA',
	'lat': -12.2488,
	'lon': 131.0453,
	'height': 25.0,
	'intake_hgt': 12.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


HBA = dict({
	'SITE_NAME': 'Halley Station, Antarctica',
	'name': 'HBA',
	'lat': -75.605,
	'lon': -26.21,
	'height': 30.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


HPB = dict({
	'SITE_NAME': 'Hohenpeissenberg',
	'name': 'HPB',
	'lat': 47.8011,
	'lon': 11.0245,
	'height': 936.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


HSU = dict({
	'SITE_NAME': 'Humboldt State University',
	'name': 'HSU',
	'lat': 41.0588,
	'lon': -124.75,
	'height': 0.0,
	'intake_hgt': 7.6,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


HUN = dict({
	'SITE_NAME': 'Hegyhatsal',
	'name': 'HUN',
	'lat': 46.95,
	'lon': 16.65,
	'height': 248.0,
	'intake_hgt': 96.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ICE = dict({
	'SITE_NAME': 'Storhofdi, Vestmannaeyjar',
	'name': 'ICE',
	'lat': 63.3998,
	'lon': -20.2884,
	'height': 118.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


IZO = dict({
	'SITE_NAME': 'Izana, Tenerife, Canary Islands',
	'name': 'IZO',
	'lat': 28.309,
	'lon': -16.499,
	'height': 2372.9,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


KEY = dict({
	'SITE_NAME': 'Key Biscayne, Florida',
	'name': 'KEY',
	'lat': 25.6654,
	'lon': -80.158,
	'height': 1.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


KUM = dict({
	'SITE_NAME': 'Cape Kumukahi, Hawaii',
	'name': 'KUM',
	'lat': 19.5608,
	'lon': -154.8883,
	'height': 8.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


KZD = dict({
	'SITE_NAME': 'Sary Taukum',
	'name': 'KZD',
	'lat': 44.0839,
	'lon': 76.8712,
	'height': 595.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


KZM = dict({
	'SITE_NAME': 'Plateau Assy',
	'name': 'KZM',
	'lat': 43.25,
	'lon': 77.88,
	'height': 2519.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


LEF = dict({
	'SITE_NAME': 'Park Falls, Wisconsin',
	'name': 'LEF',
	'lat': 45.9451,
	'lon': -90.2732,
	'height': 472.0,
	'intake_hgt': 396.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


LLB = dict({
	'SITE_NAME': 'Lac La Biche, Alberta',
	'name': 'LLB',
	'lat': 54.9538,
	'lon': -112.4666,
	'height': 540.0,
	'intake_hgt': 6.1,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


LLN = dict({
	'SITE_NAME': 'Lulin',
	'name': 'LLN',
	'lat': 23.47,
	'lon': 120.87,
	'height': 2862.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


LMP = dict({
	'SITE_NAME': 'Lampedusa',
	'name': 'LMP',
	'lat': 35.5181,
	'lon': 12.6322,
	'height': 45.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MAA = dict({
	'SITE_NAME': 'Mawson Station, Antarctica',
	'name': 'MAA',
	'lat': -67.617,
	'lon': 62.867,
	'height': 32.0,
	'intake_hgt': 10.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MBC = dict({
	'SITE_NAME': 'Mould Bay, Northwest Territories',
	'name': 'MBC',
	'lat': 76.247,
	'lon': -119.3533,
	'height': 30.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MEX = dict({
	'SITE_NAME': 'High Altitude Global Climate Observation Center',
	'name': 'MEX',
	'lat': 18.9841,
	'lon': -97.311,
	'height': 4464.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MHD = dict({
	'SITE_NAME': 'Mace Head, County Galway',
	'name': 'MHD',
	'lat': 53.326,
	'lon': -9.899,
	'height': 5.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MID = dict({
	'SITE_NAME': 'Sand Island, Midway',
	'name': 'MID',
	'lat': 28.2186,
	'lon': -177.3678,
	'height': 4.6,
	'intake_hgt': 0.3,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MKN = dict({
	'SITE_NAME': 'Mt. Kenya',
	'name': 'MKN',
	'lat': -0.0622,
	'lon': 37.2972,
	'height': 3644.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MLO = dict({
	'SITE_NAME': 'Mauna Loa, Hawaii',
	'name': 'MLO',
	'lat': 19.5362,
	'lon': -155.5763,
	'height': 3397.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MLO = dict({
	'SITE_NAME': 'Mauna Loa, Hawaii',
	'name': 'MLO',
	'lat': 19.5362,
	'lon': -155.5763,
	'height': 3397.0,
	'intake_hgt': 38.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


MQA = dict({
	'SITE_NAME': 'Macquarie Island',
	'name': 'MQA',
	'lat': -54.483,
	'lon': 158.967,
	'height': 6.0,
	'intake_hgt': 7.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


NAT = dict({
	'SITE_NAME': 'Farol De Mae Luiza Lighthouse',
	'name': 'NAT',
	'lat': -5.7952,
	'lon': -35.1853,
	'height': 50.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


NMB = dict({
	'SITE_NAME': 'Gobabeb',
	'name': 'NMB',
	'lat': -23.58,
	'lon': 15.03,
	'height': 456.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


NWR = dict({
	'SITE_NAME': 'Niwot Ridge, Colorado',
	'name': 'NWR',
	'lat': 40.0531,
	'lon': -105.5864,
	'height': 3523.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


OBN = dict({
	'SITE_NAME': 'Obninsk',
	'name': 'OBN',
	'lat': 55.11,
	'lon': 36.6,
	'height': 183.0,
	'intake_hgt': 301.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


OPW = dict({
	'SITE_NAME': 'Olympic Peninsula, Washington',
	'name': 'OPW',
	'lat': 48.2997,
	'lon': -124.6276,
	'height': 486.0,
	'intake_hgt': 2.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


OTA = dict({
	'SITE_NAME': 'Otway, Victoria',
	'name': 'OTA',
	'lat': -38.522,
	'lon': 142.817,
	'height': 40.0,
	'intake_hgt': 10.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


OXK = dict({
	'SITE_NAME': 'Ochsenkopf',
	'name': 'OXK',
	'lat': 50.0301,
	'lon': 11.8084,
	'height': 1022.0,
	'intake_hgt': 163.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


PAL = dict({
	'SITE_NAME': 'Pallas-Sammaltunturi, GAW Station',
	'name': 'PAL',
	'lat': 67.9733,
	'lon': 24.1157,
	'height': 565.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


PSA = dict({
	'SITE_NAME': 'Palmer Station, Antarctica',
	'name': 'PSA',
	'lat': -64.7742,
	'lon': -64.0527,
	'height': 10.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


PTA = dict({
	'SITE_NAME': 'Point Arena, California',
	'name': 'PTA',
	'lat': 38.9546,
	'lon': -123.7408,
	'height': 17.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


RPB = dict({
	'SITE_NAME': 'Ragged Point',
	'name': 'RPB',
	'lat': 13.165,
	'lon': -59.432,
	'height': 15.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SDZ = dict({
	'SITE_NAME': 'Shangdianzi',
	'name': 'SDZ',
	'lat': 40.65,
	'lon': 117.117,
	'height': 293.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SEY = dict({
	'SITE_NAME': 'Mahe Island',
	'name': 'SEY',
	'lat': -4.6824,
	'lon': 55.5325,
	'height': 2.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SGI = dict({
	'SITE_NAME': 'Bird Island, South Georgia',
	'name': 'SGI',
	'lat': -54.0,
	'lon': -38.05,
	'height': 30.0,
	'intake_hgt': 0.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SGP = dict({
	'SITE_NAME': 'Southern Great Plains, Oklahoma',
	'name': 'SGP',
	'lat': 36.607,
	'lon': -97.489,
	'height': 314.0,
	'intake_hgt': 60.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SHM = dict({
	'SITE_NAME': 'Shemya Island, Alaska',
	'name': 'SHM',
	'lat': 52.7112,
	'lon': 174.126,
	'height': 23.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SIS = dict({
	'SITE_NAME': 'Shetland Islands',
	'name': 'SIS',
	'lat': 60.089,
	'lon': -1.255,
	'height': 30.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SMO = dict({
	'SITE_NAME': 'Tutuila',
	'name': 'SMO',
	'lat': -14.2474,
	'lon': -170.5644,
	'height': 42.0,
	'intake_hgt': 18.3,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SPO = dict({
	'SITE_NAME': 'South Pole, Antarctica',
	'name': 'SPO',
	'lat': -89.98,
	'lon': -24.8,
	'height': 2810.0,
	'intake_hgt': 11.3,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SPO = dict({
	'SITE_NAME': 'South Pole, Antarctica',
	'name': 'SPO',
	'lat': -89.98,
	'lon': -24.8,
	'height': 2810.0,
	'intake_hgt': 10.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


STM = dict({
	'SITE_NAME': 'Ocean Station M',
	'name': 'STM',
	'lat': 66.0,
	'lon': 2.0,
	'height': 0.0,
	'intake_hgt': 7.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SUM = dict({
	'SITE_NAME': 'Summit',
	'name': 'SUM',
	'lat': 72.5962,
	'lon': -38.422,
	'height': 3209.54,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SYO = dict({
	'SITE_NAME': 'Syowa Station, Antarctica',
	'name': 'SYO',
	'lat': -69.0125,
	'lon': 39.59,
	'height': 14.0,
	'intake_hgt': 3.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


TAC = dict({
	'SITE_NAME': 'Tacolneston',
	'name': 'TAC',
	'lat': 52.5177,
	'lon': 1.1386,
	'height': 56.0,
	'intake_hgt': 180.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


TAP = dict({
	'SITE_NAME': 'Tae-ahn Peninsula',
	'name': 'TAP',
	'lat': 36.7376,
	'lon': 126.1328,
	'height': 16.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


THD = dict({
	'SITE_NAME': 'Trinidad Head, California',
	'name': 'THD',
	'lat': 41.0541,
	'lon': -124.151,
	'height': 107.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


TIK = dict({
	'SITE_NAME': 'Hydrometeorological Observatory of Tiksi',
	'name': 'TIK',
	'lat': 71.5965,
	'lon': 128.8887,
	'height': 19.0,
	'intake_hgt': 10.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


USH = dict({
	'SITE_NAME': 'Ushuaia',
	'name': 'USH',
	'lat': -54.8484,
	'lon': -68.3106,
	'height': 12.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


UTA = dict({
	'SITE_NAME': 'Wendover, Utah',
	'name': 'UTA',
	'lat': 39.9018,
	'lon': -113.7181,
	'height': 1327.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


UUM = dict({
	'SITE_NAME': 'Ulaan Uul',
	'name': 'UUM',
	'lat': 44.4516,
	'lon': 111.0956,
	'height': 1007.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


WIS = dict({
	'SITE_NAME': 'Weizmann Institute of Science at the Arava Institute, Ketura',
	'name': 'WIS',
	'lat': 29.9646,
	'lon': 35.0605,
	'height': 151.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


WKT = dict({
	'SITE_NAME': 'Moody, Texas',
	'name': 'WKT',
	'lat': 31.3149,
	'lon': -97.3269,
	'height': 251.0,
	'intake_hgt': 457.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


WLG = dict({
	'SITE_NAME': 'Mt. Waliguan',
	'name': 'WLG',
	'lat': 36.2879,
	'lon': 100.8964,
	'height': 3810.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ZEP = dict({
	'SITE_NAME': 'Ny-Alesund, Svalbard',
	'name': 'ZEP',
	'lat': 78.9067,
	'lon': 11.8883,
	'height': 474.0,
	'intake_hgt': 5.0,
	'project': 'surface-flask',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})

flask_list = [ABP, ALT, ALT, AMS, AMT, AMY, ARA, ASC, ASK, AVI, AZR, BAL, BHD, BKT, BME, BMW, BRW, BSC, CBA, CFA, CGO, CGO, CHR, CIB, CMO, CPA, CPT, CRI, CRZ, CYA, DSI, EIC, ESP, GMI, GOZ, GPA, HBA, HPB, HSU, HUN, ICE, IZO, KEY, KUM, KZD, KZM, LEF, LLB, LLN, LMP, MAA, MBC, MEX, MHD, MID, MKN, MLO, MLO, MQA, NAT, NMB, NWR, OBN, OPW, OTA, OXK, PAL, PSA, PTA, RPB, SDZ, SEY, SGI, SGP, SHM, SIS, SMO, SPO, SPO, STM, SUM, SYO, TAC, TAP, THD, TIK, USH, UTA, UUM, WIS, WKT, WLG, ZEP]

####################################################



####################################################
###################### insitu #######################

ABT = dict({
	'SITE_NAME': 'Abbotsford, British Columbia',
	'name': 'ABT',
	'lat': 49.0114,
	'lon': -122.3353,
	'height': 60.0,
	'intake_hgt': 33.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


ALT = dict({
	'SITE_NAME': 'Alert, Nunavut',
	'name': 'ALT',
	'lat': 82.4508,
	'lon': -62.5072,
	'height': 185.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


AMY = dict({
	'SITE_NAME': 'Anmyeon-do',
	'name': 'AMY',
	'lat': 36.5389,
	'lon': 126.3295,
	'height': 47.0,
	'intake_hgt': 20.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


BCK = dict({
	'SITE_NAME': 'Behchoko, Northwest Territories',
	'name': 'BCK',
	'lat': 62.7981,
	'lon': -115.9194,
	'height': 160.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


BHD = dict({
	'SITE_NAME': 'Baring Head Station',
	'name': 'BHD',
	'lat': -41.4083,
	'lon': 174.871,
	'height': 85.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'baseline',
	'calibration': 'WMO CH4 X2004A'
})


BLK = dict({
	'SITE_NAME': 'Baker Lake, Nunavut',
	'name': 'BLK',
	'lat': 64.3317,
	'lon': -96.0104,
	'height': 51.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


BRA = dict({
	'SITE_NAME': 'Bratt's Lake Saskatchewan',
	'name': 'BRA',
	'lat': 50.2017,
	'lon': -104.7113,
	'height': 595.0,
	'intake_hgt': 35.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


BRW = dict({
	'SITE_NAME': 'Barrow Atmospheric Baseline Observatory',
	'name': 'BRW',
	'lat': 71.323,
	'lon': -156.6114,
	'height': 11.0,
	'intake_hgt': 16.46,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


CBY = dict({
	'SITE_NAME': 'Cambridge Bay, Nunavut Territory',
	'name': 'CBY',
	'lat': 69.1284,
	'lon': -105.0577,
	'height': 35.0,
	'intake_hgt': 12.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


CDL = dict({
	'SITE_NAME': 'Candle Lake, Saskatchewan',
	'name': 'CDL',
	'lat': 53.9871,
	'lon': -105.1179,
	'height': 600.0,
	'intake_hgt': 30.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


CHL = dict({
	'SITE_NAME': 'Churchill, Manitoba',
	'name': 'CHL',
	'lat': 58.7379,
	'lon': -93.8194,
	'height': 29.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


CHM = dict({
	'SITE_NAME': 'Chibougamau, Quebec',
	'name': 'CHM',
	'lat': 49.6925,
	'lon': -74.3432,
	'height': 393.0,
	'intake_hgt': 30.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


CMN = dict({
	'SITE_NAME': 'Mt. Cimone Station',
	'name': 'CMN',
	'lat': 44.1936,
	'lon': 10.6999,
	'height': 2165.0,
	'intake_hgt': 7.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


CMN = dict({
	'SITE_NAME': 'Mt. Cimone Station',
	'name': 'CMN',
	'lat': 44.1936,
	'lon': 10.6999,
	'height': 2165.0,
	'intake_hgt': 8.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-8magl',
	'calibration': 'WMO CH4 X2004A'
})


CPS = dict({
	'SITE_NAME': 'Chapais,Quebec',
	'name': 'CPS',
	'lat': 49.8223,
	'lon': -74.9753,
	'height': 391.0,
	'intake_hgt': 8.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


CPT = dict({
	'SITE_NAME': 'Cape Point',
	'name': 'CPT',
	'lat': -34.3523,
	'lon': 18.4891,
	'height': 230.0,
	'intake_hgt': 30.0,
	'project': 'surface-insitu',
	'sel_tag': 'marine',
	'calibration': 'WMO CH4 X2004 (1983-2014), WMO CH4 X2004A (2015-present)'
})


CVO = dict({
	'SITE_NAME': 'Cape Verde Observatory',
	'name': 'CVO',
	'lat': 16.864,
	'lon': -24.8675,
	'height': 10.0,
	'intake_hgt': 30.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


EGB = dict({
	'SITE_NAME': 'Egbert, Ontario',
	'name': 'EGB',
	'lat': 44.231,
	'lon': -79.7838,
	'height': 251.0,
	'intake_hgt': 3.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


ENA = dict({
	'SITE_NAME': 'Eastern North Atlantic, Graciosa, Azores',
	'name': 'ENA',
	'lat': 39.0916,
	'lon': -28.0257,
	'height': 30.48,
	'intake_hgt': 0.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-10magl',
	'calibration': 'WMO CH4 X2004'
})


ESP = dict({
	'SITE_NAME': 'Estevan Point,  British Columbia',
	'name': 'ESP',
	'lat': 49.383,
	'lon': -126.5441,
	'height': 7.0,
	'intake_hgt': 40.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


EST = dict({
	'SITE_NAME': 'Esther, Alberta',
	'name': 'EST',
	'lat': 51.6707,
	'lon': -110.206,
	'height': 707.0,
	'intake_hgt': 3.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


ETL = dict({
	'SITE_NAME': 'East Trout Lake, Saskatchewan',
	'name': 'ETL',
	'lat': 54.3541,
	'lon': -104.9868,
	'height': 493.0,
	'intake_hgt': 105.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


FNE = dict({
	'SITE_NAME': 'Fort Nelson, British Columbia',
	'name': 'FNE',
	'lat': 58.8412,
	'lon': -122.5737,
	'height': 361.0,
	'intake_hgt': 15.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


FSD = dict({
	'SITE_NAME': 'Fraserdale',
	'name': 'FSD',
	'lat': 49.8752,
	'lon': -81.57,
	'height': 210.0,
	'intake_hgt': 40.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


GAT = dict({
	'SITE_NAME': 'Gartow',
	'name': 'GAT',
	'lat': 53.0657,
	'lon': 11.4429,
	'height': 70.0,
	'intake_hgt': 132.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-132magl',
	'calibration': 'WMO CH4 X2004A'
})


GAT = dict({
	'SITE_NAME': 'Gartow',
	'name': 'GAT',
	'lat': 53.0657,
	'lon': 11.4429,
	'height': 70.0,
	'intake_hgt': 216.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-216magl',
	'calibration': 'WMO CH4 X2004A'
})


GAT = dict({
	'SITE_NAME': 'Gartow',
	'name': 'GAT',
	'lat': 53.0657,
	'lon': 11.4429,
	'height': 70.0,
	'intake_hgt': 30.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-30magl',
	'calibration': 'WMO CH4 X2004A'
})


GAT = dict({
	'SITE_NAME': 'Gartow',
	'name': 'GAT',
	'lat': 53.0657,
	'lon': 11.4429,
	'height': 70.0,
	'intake_hgt': 341.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-341magl',
	'calibration': 'WMO CH4 X2004A'
})


GAT = dict({
	'SITE_NAME': 'Gartow',
	'name': 'GAT',
	'lat': 53.0657,
	'lon': 11.4429,
	'height': 70.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-60magl',
	'calibration': 'WMO CH4 X2004A'
})


HEL = dict({
	'SITE_NAME': 'Helgoland',
	'name': 'HEL',
	'lat': 54.1804,
	'lon': 7.8832,
	'height': 43.0,
	'intake_hgt': 110.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-110magl',
	'calibration': 'WMO CH4 X2004A'
})


HNP = dict({
	'SITE_NAME': 'Hanlan's Point, Ontario',
	'name': 'HNP',
	'lat': 43.6122,
	'lon': -79.3887,
	'height': 87.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


HPB = dict({
	'SITE_NAME': 'Hohenpeissenberg',
	'name': 'HPB',
	'lat': 47.8011,
	'lon': 11.0245,
	'height': 936.0,
	'intake_hgt': 131.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-131magl',
	'calibration': 'WMO CH4 X2004A'
})


HPB = dict({
	'SITE_NAME': 'Hohenpeissenberg',
	'name': 'HPB',
	'lat': 47.8011,
	'lon': 11.0245,
	'height': 936.0,
	'intake_hgt': 50.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-50magl',
	'calibration': 'WMO CH4 X2004A'
})


HPB = dict({
	'SITE_NAME': 'Hohenpeissenberg',
	'name': 'HPB',
	'lat': 47.8011,
	'lon': 11.0245,
	'height': 936.0,
	'intake_hgt': 93.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-93magl',
	'calibration': 'WMO CH4 X2004A'
})


HTM = dict({
	'SITE_NAME': 'Hyltemossa',
	'name': 'HTM',
	'lat': 56.0976,
	'lon': 13.4189,
	'height': 115.0,
	'intake_hgt': 150.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-150magl',
	'calibration': 'WMO CH4 X2004A'
})


HTM = dict({
	'SITE_NAME': 'Hyltemossa',
	'name': 'HTM',
	'lat': 56.0976,
	'lon': 13.4189,
	'height': 115.0,
	'intake_hgt': 30.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-30magl',
	'calibration': 'WMO CH4 X2004A'
})


HTM = dict({
	'SITE_NAME': 'Hyltemossa',
	'name': 'HTM',
	'lat': 56.0976,
	'lon': 13.4189,
	'height': 115.0,
	'intake_hgt': 70.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-70magl',
	'calibration': 'WMO CH4 X2004A'
})


INU = dict({
	'SITE_NAME': 'Inuvik,Northwest Territories',
	'name': 'INU',
	'lat': 68.3178,
	'lon': -133.5342,
	'height': 113.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


INX01 = dict({
	'SITE_NAME': 'INFLUX Site - tower 1',
	'name': 'INX01',
	'lat': 39.5805,
	'lon': -86.4207,
	'height': 256.3,
	'intake_hgt': 121.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX02 = dict({
	'SITE_NAME': 'INFLUX Site - tower 2',
	'name': 'INX02',
	'lat': 39.7978,
	'lon': -86.0183,
	'height': 266.9,
	'intake_hgt': 136.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX03 = dict({
	'SITE_NAME': 'INFLUX Site - tower 3',
	'name': 'INX03',
	'lat': 39.7833,
	'lon': -86.1651,
	'height': 226.0,
	'intake_hgt': 54.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX04 = dict({
	'SITE_NAME': 'INFLUX Site - tower 4',
	'name': 'INX04',
	'lat': 39.5926,
	'lon': -86.1099,
	'height': 249.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX05 = dict({
	'SITE_NAME': 'INFLUX Site - tower 5',
	'name': 'INX05',
	'lat': 39.8949,
	'lon': -86.2028,
	'height': 251.0,
	'intake_hgt': 125.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX07 = dict({
	'SITE_NAME': 'INFLUX Site - tower 7',
	'name': 'INX07',
	'lat': 39.7739,
	'lon': -86.2724,
	'height': 241.0,
	'intake_hgt': 58.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX08 = dict({
	'SITE_NAME': 'INFLUX Site - tower 8',
	'name': 'INX08',
	'lat': 40.0411,
	'lon': -85.9734,
	'height': 245.0,
	'intake_hgt': 41.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX09 = dict({
	'SITE_NAME': 'INFLUX Site - tower 9',
	'name': 'INX09',
	'lat': 39.8627,
	'lon': -85.7448,
	'height': 277.0,
	'intake_hgt': 130.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX10 = dict({
	'SITE_NAME': 'INFLUX Site - tower 10',
	'name': 'INX10',
	'lat': 39.7181,
	'lon': -86.1436,
	'height': 223.0,
	'intake_hgt': 40.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX11 = dict({
	'SITE_NAME': 'INFLUX Site - tower 11',
	'name': 'INX11',
	'lat': 39.8403,
	'lon': -86.1763,
	'height': 214.0,
	'intake_hgt': 130.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX13 = dict({
	'SITE_NAME': 'INFLUX Site - tower 13',
	'name': 'INX13',
	'lat': 39.7173,
	'lon': -85.9417,
	'height': 248.0,
	'intake_hgt': 87.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX14 = dict({
	'SITE_NAME': 'INFLUX Site - tower 14',
	'name': 'INX14',
	'lat': 39.9971,
	'lon': -86.7396,
	'height': 264.0,
	'intake_hgt': 76.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


INX15 = dict({
	'SITE_NAME': 'INFLUX Site - tower 15',
	'name': 'INX15',
	'lat': 39.6376,
	'lon': 86.4779,
	'height': 258.0,
	'intake_hgt': 75.0,
	'project': 'surface-insitu',
	'sel_tag': 'allhours',
	'calibration': 'WMO CH4 X2004A'
})


IPR = dict({
	'SITE_NAME': 'Ispra',
	'name': 'IPR',
	'lat': 45.807,
	'lon': 8.63,
	'height': 223.0,
	'intake_hgt': 100.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-100magl',
	'calibration': 'WMO CH4 X2004A'
})


IPR = dict({
	'SITE_NAME': 'Ispra',
	'name': 'IPR',
	'lat': 45.807,
	'lon': 8.63,
	'height': 223.0,
	'intake_hgt': 40.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-40magl',
	'calibration': 'WMO CH4 X2004A'
})


IPR = dict({
	'SITE_NAME': 'Ispra',
	'name': 'IPR',
	'lat': 45.807,
	'lon': 8.63,
	'height': 223.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-60magl',
	'calibration': 'WMO CH4 X2004A'
})


IZO = dict({
	'SITE_NAME': 'Izana, Tenerife, Canary Islands',
	'name': 'IZO',
	'lat': 28.309,
	'lon': -16.499,
	'height': 2372.9,
	'intake_hgt': 8.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


JFJ = dict({
	'SITE_NAME': 'Jungfraujoch',
	'name': 'JFJ',
	'lat': 46.55,
	'lon': 7.987,
	'height': 3570.0,
	'intake_hgt': 5.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-5magl',
	'calibration': 'WMO CH4 X2004A'
})


JFJ = dict({
	'SITE_NAME': 'Jungfraujoch',
	'name': 'JFJ',
	'lat': 46.55,
	'lon': 7.987,
	'height': 3570.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


JGS = dict({
	'SITE_NAME': 'Jejudo Gosan Suwolbong',
	'name': 'JGS',
	'lat': 33.3,
	'lon': 126.16,
	'height': 71.47,
	'intake_hgt': 6.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


KIT = dict({
	'SITE_NAME': 'Karlsruhe',
	'name': 'KIT',
	'lat': 49.0915,
	'lon': 8.4249,
	'height': 110.0,
	'intake_hgt': 100.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-100magl',
	'calibration': 'WMO CH4 X2004A'
})


KIT = dict({
	'SITE_NAME': 'Karlsruhe',
	'name': 'KIT',
	'lat': 49.0915,
	'lon': 8.4249,
	'height': 110.0,
	'intake_hgt': 200.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-200magl',
	'calibration': 'WMO CH4 X2004A'
})


KIT = dict({
	'SITE_NAME': 'Karlsruhe',
	'name': 'KIT',
	'lat': 49.0915,
	'lon': 8.4249,
	'height': 110.0,
	'intake_hgt': 30.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-30magl',
	'calibration': 'WMO CH4 X2004A'
})


KIT = dict({
	'SITE_NAME': 'Karlsruhe',
	'name': 'KIT',
	'lat': 49.0915,
	'lon': 8.4249,
	'height': 110.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-60magl',
	'calibration': 'WMO CH4 X2004A'
})


KRE = dict({
	'SITE_NAME': 'Kresin u Pacova',
	'name': 'KRE',
	'lat': 49.583,
	'lon': 15.083,
	'height': 534.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-10magl',
	'calibration': 'WMO CH4 X2004A'
})


KRE = dict({
	'SITE_NAME': 'Kresin u Pacova',
	'name': 'KRE',
	'lat': 49.583,
	'lon': 15.083,
	'height': 534.0,
	'intake_hgt': 125.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-125magl',
	'calibration': 'WMO CH4 X2004A'
})


KRE = dict({
	'SITE_NAME': 'Kresin u Pacova',
	'name': 'KRE',
	'lat': 49.583,
	'lon': 15.083,
	'height': 534.0,
	'intake_hgt': 250.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-250magl',
	'calibration': 'WMO CH4 X2004A'
})


KRE = dict({
	'SITE_NAME': 'Kresin u Pacova',
	'name': 'KRE',
	'lat': 49.583,
	'lon': 15.083,
	'height': 534.0,
	'intake_hgt': 50.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-50magl',
	'calibration': 'WMO CH4 X2004A'
})


LIN = dict({
	'SITE_NAME': 'Lindenberg',
	'name': 'LIN',
	'lat': 52.1663,
	'lon': 14.1226,
	'height': 73.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-10magl',
	'calibration': 'WMO CH4 X2004A'
})


LIN = dict({
	'SITE_NAME': 'Lindenberg',
	'name': 'LIN',
	'lat': 52.1663,
	'lon': 14.1226,
	'height': 73.0,
	'intake_hgt': 2.5,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-2.5magl',
	'calibration': 'WMO CH4 X2004A'
})


LIN = dict({
	'SITE_NAME': 'Lindenberg',
	'name': 'LIN',
	'lat': 52.1663,
	'lon': 14.1226,
	'height': 73.0,
	'intake_hgt': 40.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-40magl',
	'calibration': 'WMO CH4 X2004A'
})


LIN = dict({
	'SITE_NAME': 'Lindenberg',
	'name': 'LIN',
	'lat': 52.1663,
	'lon': 14.1226,
	'height': 73.0,
	'intake_hgt': 98.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-98magl',
	'calibration': 'WMO CH4 X2004A'
})


LLB = dict({
	'SITE_NAME': 'Lac La Biche, Alberta',
	'name': 'LLB',
	'lat': 54.9538,
	'lon': -112.4666,
	'height': 540.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


LMP = dict({
	'SITE_NAME': 'Lampedusa',
	'name': 'LMP',
	'lat': 35.5181,
	'lon': 12.6322,
	'height': 45.0,
	'intake_hgt': 8.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-8magl',
	'calibration': 'WMO CH4 X2004A'
})


LUT = dict({
	'SITE_NAME': 'Lutjewad',
	'name': 'LUT',
	'lat': 53.4037,
	'lon': 6.3529,
	'height': 1.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-60magl',
	'calibration': 'WMO CH4 X2004A'
})


MLO = dict({
	'SITE_NAME': 'Mauna Loa, Hawaii',
	'name': 'MLO',
	'lat': 19.5362,
	'lon': -155.5763,
	'height': 3397.0,
	'intake_hgt': 40.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


MNM = dict({
	'SITE_NAME': 'Minamitorishima',
	'name': 'MNM',
	'lat': 24.28,
	'lon': 153.98,
	'height': 8.0,
	'intake_hgt': 20.0,
	'project': 'surface-insitu',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


NOR = dict({
	'SITE_NAME': 'Norunda',
	'name': 'NOR',
	'lat': 60.0864,
	'lon': 17.4794,
	'height': 46.0,
	'intake_hgt': 100.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-100magl',
	'calibration': 'WMO CH4 X2004A'
})


NOR = dict({
	'SITE_NAME': 'Norunda',
	'name': 'NOR',
	'lat': 60.0864,
	'lon': 17.4794,
	'height': 46.0,
	'intake_hgt': 32.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-32magl',
	'calibration': 'WMO CH4 X2004A'
})


NOR = dict({
	'SITE_NAME': 'Norunda',
	'name': 'NOR',
	'lat': 60.0864,
	'lon': 17.4794,
	'height': 46.0,
	'intake_hgt': 58.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-58magl',
	'calibration': 'WMO CH4 X2004A'
})


OLI = dict({
	'SITE_NAME': 'Oliktok Point, Alaska',
	'name': 'OLI',
	'lat': 70.4953,
	'lon': -149.8869,
	'height': 2.0,
	'intake_hgt': 0.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-10magl',
	'calibration': 'WMO CH4 X2004'
})


OPE = dict({
	'SITE_NAME': 'Observatoire perenne de l'environnement',
	'name': 'OPE',
	'lat': 48.5619,
	'lon': 5.5036,
	'height': 390.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-10magl',
	'calibration': 'WMO CH4 X2004A'
})


OPE = dict({
	'SITE_NAME': 'Observatoire perenne de l'environnement',
	'name': 'OPE',
	'lat': 48.5619,
	'lon': 5.5036,
	'height': 390.0,
	'intake_hgt': 120.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-120magl',
	'calibration': 'WMO CH4 X2004A'
})


OPE = dict({
	'SITE_NAME': 'Observatoire perenne de l'environnement',
	'name': 'OPE',
	'lat': 48.5619,
	'lon': 5.5036,
	'height': 390.0,
	'intake_hgt': 50.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-50magl',
	'calibration': 'WMO CH4 X2004A'
})


OXK = dict({
	'SITE_NAME': 'Ochsenkopf',
	'name': 'OXK',
	'lat': 50.0301,
	'lon': 11.8084,
	'height': 1022.0,
	'intake_hgt': 163.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-163magl',
	'calibration': 'WMO CH4 X2004A'
})


OXK = dict({
	'SITE_NAME': 'Ochsenkopf',
	'name': 'OXK',
	'lat': 50.0301,
	'lon': 11.8084,
	'height': 1022.0,
	'intake_hgt': 23.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-23magl',
	'calibration': 'WMO CH4 X2004A'
})


OXK = dict({
	'SITE_NAME': 'Ochsenkopf',
	'name': 'OXK',
	'lat': 50.0301,
	'lon': 11.8084,
	'height': 1022.0,
	'intake_hgt': 90.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-90magl',
	'calibration': 'WMO CH4 X2004A'
})


PAL = dict({
	'SITE_NAME': 'Pallas-Sammaltunturi, GAW Station',
	'name': 'PAL',
	'lat': 67.9733,
	'lon': 24.1157,
	'height': 565.0,
	'intake_hgt': 5.0,
	'project': 'surface-insitu',
	'sel_tag': 'continental',
	'calibration': 'WMO CH4 X2004A'
})


PAL = dict({
	'SITE_NAME': 'Pallas-Sammaltunturi, GAW Station',
	'name': 'PAL',
	'lat': 67.9733,
	'lon': 24.1157,
	'height': 565.0,
	'intake_hgt': 5.0,
	'project': 'surface-insitu',
	'sel_tag': 'marine',
	'calibration': 'WMO CH4 X2004A'
})


PAL = dict({
	'SITE_NAME': 'Pallas-Sammaltunturi, GAW Station',
	'name': 'PAL',
	'lat': 67.9733,
	'lon': 24.1157,
	'height': 565.0,
	'intake_hgt': 5.0,
	'project': 'surface-insitu',
	'sel_tag': 'nonlocal',
	'calibration': 'WMO CH4 X2004A'
})


PAL = dict({
	'SITE_NAME': 'Pallas-Sammaltunturi, GAW Station',
	'name': 'PAL',
	'lat': 67.9733,
	'lon': 24.1157,
	'height': 565.0,
	'intake_hgt': 12.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-12magl',
	'calibration': 'WMO CH4 X2004A'
})


PUY = dict({
	'SITE_NAME': 'Puy de Dome',
	'name': 'PUY',
	'lat': 45.7719,
	'lon': 2.9658,
	'height': 1465.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-10magl',
	'calibration': 'WMO CH4 X2004A'
})


RUN = dict({
	'SITE_NAME': 'La Réunion',
	'name': 'RUN',
	'lat': -21.0796,
	'lon': 55.3841,
	'height': 2154.0,
	'intake_hgt': 6.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-6magl',
	'calibration': 'WMO CH4 X2004A'
})


RYO = dict({
	'SITE_NAME': 'Ryori',
	'name': 'RYO',
	'lat': 39.03,
	'lon': 141.82,
	'height': 260.0,
	'intake_hgt': 20.0,
	'project': 'surface-insitu',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


SAC = dict({
	'SITE_NAME': 'Saclay',
	'name': 'SAC',
	'lat': 48.7227,
	'lon': 2.142,
	'height': 160.0,
	'intake_hgt': 100.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-100magl',
	'calibration': 'WMO CH4 X2004A'
})


SAC = dict({
	'SITE_NAME': 'Saclay',
	'name': 'SAC',
	'lat': 48.7227,
	'lon': 2.142,
	'height': 160.0,
	'intake_hgt': 15.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-15magl',
	'calibration': 'WMO CH4 X2004A'
})


SAC = dict({
	'SITE_NAME': 'Saclay',
	'name': 'SAC',
	'lat': 48.7227,
	'lon': 2.142,
	'height': 160.0,
	'intake_hgt': 60.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-50magl',
	'calibration': 'WMO CH4 X2004A'
})


SGP = dict({
	'SITE_NAME': 'Southern Great Plains, Oklahoma',
	'name': 'SGP',
	'lat': 36.607,
	'lon': -97.489,
	'height': 314.0,
	'intake_hgt': 4.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-60magl',
	'calibration': 'WMO CH4 X2004'
})


SMR = dict({
	'SITE_NAME': 'Hyytiala',
	'name': 'SMR',
	'lat': 61.8474,
	'lon': 24.2948,
	'height': 181.0,
	'intake_hgt': 125.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-125magl',
	'calibration': 'WMO CH4 X2004A'
})


SMR = dict({
	'SITE_NAME': 'Hyytiala',
	'name': 'SMR',
	'lat': 61.8474,
	'lon': 24.2948,
	'height': 181.0,
	'intake_hgt': 16.8,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-17magl',
	'calibration': 'WMO CH4 X2004A'
})


SMR = dict({
	'SITE_NAME': 'Hyytiala',
	'name': 'SMR',
	'lat': 61.8474,
	'lon': 24.2948,
	'height': 181.0,
	'intake_hgt': 67.2,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-67magl',
	'calibration': 'WMO CH4 X2004A'
})


STE = dict({
	'SITE_NAME': 'Steinkimmen',
	'name': 'STE',
	'lat': 53.0431,
	'lon': 8.4588,
	'height': 29.0,
	'intake_hgt': 127.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-127magl',
	'calibration': 'WMO CH4 X2004A'
})


STE = dict({
	'SITE_NAME': 'Steinkimmen',
	'name': 'STE',
	'lat': 53.0431,
	'lon': 8.4588,
	'height': 29.0,
	'intake_hgt': 187.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-187magl',
	'calibration': 'WMO CH4 X2004A'
})


STE = dict({
	'SITE_NAME': 'Steinkimmen',
	'name': 'STE',
	'lat': 53.0431,
	'lon': 8.4588,
	'height': 29.0,
	'intake_hgt': 252.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-252magl',
	'calibration': 'WMO CH4 X2004A'
})


STE = dict({
	'SITE_NAME': 'Steinkimmen',
	'name': 'STE',
	'lat': 53.0431,
	'lon': 8.4588,
	'height': 29.0,
	'intake_hgt': 32.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-32magl',
	'calibration': 'WMO CH4 X2004A'
})


STE = dict({
	'SITE_NAME': 'Steinkimmen',
	'name': 'STE',
	'lat': 53.0431,
	'lon': 8.4588,
	'height': 29.0,
	'intake_hgt': 82.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-82magl',
	'calibration': 'WMO CH4 X2004A'
})


SVB = dict({
	'SITE_NAME': 'Svartberget',
	'name': 'SVB',
	'lat': 64.256,
	'lon': 19.775,
	'height': 235.0,
	'intake_hgt': 150.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-150magl',
	'calibration': 'WMO CH4 X2004A'
})


SVB = dict({
	'SITE_NAME': 'Svartberget',
	'name': 'SVB',
	'lat': 64.256,
	'lon': 19.775,
	'height': 235.0,
	'intake_hgt': 35.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-35magl',
	'calibration': 'WMO CH4 X2004A'
})


SVB = dict({
	'SITE_NAME': 'Svartberget',
	'name': 'SVB',
	'lat': 64.256,
	'lon': 19.775,
	'height': 235.0,
	'intake_hgt': 85.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-85magl',
	'calibration': 'WMO CH4 X2004A'
})


SYO = dict({
	'SITE_NAME': 'Syowa Station, Antarctica',
	'name': 'SYO',
	'lat': -69.0125,
	'lon': 39.59,
	'height': 14.0,
	'intake_hgt': 8.0,
	'project': 'surface-insitu',
	'sel_tag': 'representative',
	'calibration': 'Tohoku University/NIPR mole fraction scale 1987'
})


TAO = dict({
	'SITE_NAME': 'Toronto Atmospheric Observatory, University of Toronto, Ontario',
	'name': 'TAO',
	'lat': 43.66,
	'lon': -79.4,
	'height': 100.0,
	'intake_hgt': 174.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


TIK = dict({
	'SITE_NAME': 'Hydrometeorological Observatory of Tiksi',
	'name': 'TIK',
	'lat': 71.5965,
	'lon': 128.8887,
	'height': 19.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


TOH = dict({
	'SITE_NAME': 'Torfhaus',
	'name': 'TOH',
	'lat': 51.8088,
	'lon': 10.535,
	'height': 801.0,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-10magl',
	'calibration': 'WMO CH4 X2004A'
})


TOH = dict({
	'SITE_NAME': 'Torfhaus',
	'name': 'TOH',
	'lat': 51.8088,
	'lon': 10.535,
	'height': 801.0,
	'intake_hgt': 110.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-110magl',
	'calibration': 'WMO CH4 X2004A'
})


TOH = dict({
	'SITE_NAME': 'Torfhaus',
	'name': 'TOH',
	'lat': 51.8088,
	'lon': 10.535,
	'height': 801.0,
	'intake_hgt': 147.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-147magl',
	'calibration': 'WMO CH4 X2004A'
})


TOH = dict({
	'SITE_NAME': 'Torfhaus',
	'name': 'TOH',
	'lat': 51.8088,
	'lon': 10.535,
	'height': 801.0,
	'intake_hgt': 76.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-76magl',
	'calibration': 'WMO CH4 X2004A'
})


TPD = dict({
	'SITE_NAME': 'Turkey Point, Ontario',
	'name': 'TPD',
	'lat': 42.6354,
	'lon': -80.5577,
	'height': 231.0,
	'intake_hgt': 35.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


TRN = dict({
	'SITE_NAME': 'Trainou',
	'name': 'TRN',
	'lat': 47.9647,
	'lon': 2.1125,
	'height': 131.0,
	'intake_hgt': 100.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-100magl',
	'calibration': 'WMO CH4 X2004A'
})


TRN = dict({
	'SITE_NAME': 'Trainou',
	'name': 'TRN',
	'lat': 47.9647,
	'lon': 2.1125,
	'height': 131.0,
	'intake_hgt': 180.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-180magl',
	'calibration': 'WMO CH4 X2004A'
})


TRN = dict({
	'SITE_NAME': 'Trainou',
	'name': 'TRN',
	'lat': 47.9647,
	'lon': 2.1125,
	'height': 131.0,
	'intake_hgt': 50.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-50magl',
	'calibration': 'WMO CH4 X2004A'
})


TRN = dict({
	'SITE_NAME': 'Trainou',
	'name': 'TRN',
	'lat': 47.9647,
	'lon': 2.1125,
	'height': 131.0,
	'intake_hgt': 5.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-5magl',
	'calibration': 'WMO CH4 X2004A'
})


ULD = dict({
	'SITE_NAME': 'Ulleungdo',
	'name': 'ULD',
	'lat': 37.48,
	'lon': 130.9,
	'height': 220.9,
	'intake_hgt': 10.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


UTO = dict({
	'SITE_NAME': 'Uto',
	'name': 'UTO',
	'lat': 59.7839,
	'lon': 21.3672,
	'height': 8.0,
	'intake_hgt': 57.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-57magl',
	'calibration': 'WMO CH4 X2004A'
})


WSA = dict({
	'SITE_NAME': 'Sable Island, Nova Scotia',
	'name': 'WSA',
	'lat': 43.9322,
	'lon': -60.0093,
	'height': 5.0,
	'intake_hgt': 25.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid',
	'calibration': 'WMO CH4 X2004A'
})


YON = dict({
	'SITE_NAME': 'Yonagunijima',
	'name': 'YON',
	'lat': 24.47,
	'lon': 123.02,
	'height': 30.0,
	'intake_hgt': 20.0,
	'project': 'surface-insitu',
	'sel_tag': 'representative',
	'calibration': 'WMO CH4 X2004A'
})


ZEP = dict({
	'SITE_NAME': 'Ny-Alesund, Svalbard',
	'name': 'ZEP',
	'lat': 78.9067,
	'lon': 11.8883,
	'height': 474.0,
	'intake_hgt': 15.0,
	'project': 'surface-insitu',
	'sel_tag': 'allvalid-15magl',
	'calibration': 'WMO CH4 X2004A'
})

insitu_list = [ABT, ALT, AMY, BCK, BHD, BLK, BRA, BRW, CBY, CDL, CHL, CHM, CMN, CMN, CPS, CPT, CVO, EGB, ENA, ESP, EST, ETL, FNE, FSD, GAT, GAT, GAT, GAT, GAT, HEL, HNP, HPB, HPB, HPB, HTM, HTM, HTM, INU, INX01, INX02, INX03, INX04, INX05, INX07, INX08, INX09, INX10, INX11, INX13, INX14, INX15, IPR, IPR, IPR, IZO, JFJ, JFJ, JGS, KIT, KIT, KIT, KIT, KRE, KRE, KRE, KRE, LIN, LIN, LIN, LIN, LLB, LMP, LUT, MLO, MNM, NOR, NOR, NOR, OLI, OPE, OPE, OPE, OXK, OXK, OXK, PAL, PAL, PAL, PAL, PUY, RUN, RYO, SAC, SAC, SAC, SGP, SMR, SMR, SMR, STE, STE, STE, STE, STE, SVB, SVB, SVB, SYO, TAO, TIK, TOH, TOH, TOH, TOH, TPD, TRN, TRN, TRN, TRN, ULD, UTO, WSA, YON, ZEP, ]
####################################################





