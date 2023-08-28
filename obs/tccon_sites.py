import numpy as np
import pandas as pd

## tccon sites information#######

## all tccon site ########
tccon_site = {'site_id': np.array(['ae', 'an', 'bi', 'br', 'bu', 'ci', 'db', 'df', 'et', 'eu', 'fc',\
                                   'gm', 'hf', 'if', 'iz', 'jc', 'jf', 'js', 'ka', 'lh', 'll', 'lr',\
                                   'ma', 'ni', 'oc', 'or', 'pa', 'pr', 'ra', 'rj', 'so', 'sp', 'tk',\
                                   'wg', 'zs'], dtype=object),
              'long_name': np.array(['ascension01', 'anmeyondo01', 'bialystok01', 'bremen01',\
                                     'burgos01', 'pasadena01', 'darwin01', 'edwards01',\
                                     'easttroutlake01', 'eureka01', 'fourcorners01', 'garmisch01',\
                                     'hefei01', 'indianapolis01', 'izana01', 'jpl01', 'jpl02', 'saga01',\
                                     'karlsruhe01', 'lauder01', 'lauder02', 'lauder03', 'manaus01',\
                                     'nicosia01', 'lamont01', 'orleans01', 'parkfalls01', 'paris01',\
                                     'reunion01', 'rikubetsu01', 'sodankyla01', 'nyalesund01',\
                                     'tsukuba02', 'wollongong01', 'zugspitze01'], dtype=object),
              'longitude': np.array([ -14.33,  126.33,   23.02,    8.85,  120.65, -118.13,  130.89,\
                                     -117.88, -104.99,  -86.42, -108.48,   11.06,  117.17,  -86.  ,\
                                     -16.48, -118.18, -118.18,  130.29,    8.44,  169.68,  169.68,\
                                     169.68,  -60.6 ,   33.38,  -97.49,    2.11,  -90.27,    2.36,\
                                     55.49,  143.77,   26.63,   11.92,  140.12,  150.88, \
                                     10.98],dtype=np.float32),
              'latitude': np.array([ -7.92,  36.54,  53.23,  53.1 ,  18.53,  34.14, -12.43,  34.96,\
                                    54.36,  80.05,  36.8 ,  47.48,  31.9 ,  39.86,  28.3 ,  34.2 ,\
                                    34.2 ,  33.24,  49.1 , -45.05, -45.04, -45.04,  -3.21,  35.14,\
                                    36.6 ,  47.97,  45.94,  48.85, -20.9 ,  43.46,  67.37,  78.92,\
                                    36.05, -34.41,  47.42], dtype=np.float32)}

tccon_info = pd.DataFrame(tccon_site, columns = ['site_id','long_name','longitude','latitude'])
tccon_info = tccon_info.set_index('site_id')

## use for assimilation
assim_id = tccon_info.index.values
assim_lon = tccon_info.loc[assim_id, 'longitude'].values
assim_lat = tccon_info.loc[assim_id, 'latitude'].values


## use for validation
#valid_id = np.array([], dtype=object)
valid_id = tccon_info.index.values
valid_lon = tccon_info.loc[valid_id, 'longitude'].values
valid_lat = tccon_info.loc[valid_id, 'latitude'].values