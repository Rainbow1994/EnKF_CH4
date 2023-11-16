"""
read configuration file for pre-define flux perturbations 

"""

from numpy import *
#import decode_line as decl

def read_em_conf(flnm):
    # flnm:/str/in/name of perturbation flux file
    # doy: /list/integer/out/starting doy of the pulse-like perturbations
    # year:/list/integer/out/starting year of the pulse-like perturbations
    # npb: /list/integer/out/number of perturbation function in each file
    # ch4flnm:/list/str/out/binary files store the flux
    
    #ps0 local variables:
    # lines/str; fl_cfg/file handle; 

    
    #ps1/open config file/
    
    try:
        
        fl_cfg = open(flnm, 'r')
    except IOError:
        #br2 if file not exists, exit  
        print('Can not open '+flnm.strip())
        return None, None
    
    #ps2/read in lines
    
    lines = fl_cfg.readlines()
    fl_cfg.close()

    #ps3/fill outputs
    

    doy, year, npb, ch4flnm =\
       list(), list(), list(), list()
    
    for line in lines[2:]:
        line = line.strip()
        line.replace('\n', '')
        if (len(line)>0):
            terms = line.split(',')
            doy.append(int(terms[0]))
            year.append(int(terms[1]))
            npb.append(int(terms[2]))
            ch4flnm.append(terms[3].strip())

    #ps4 return the values
 #   print(doy)
    return doy, year, npb, ch4flnm

            
                        
    
    
    
    
