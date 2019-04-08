#!/usr/bin/env python
"""

"""

import os
import sys
import glob

sys.path.append(os.path.abspath('../'))
import ResolveVectorsLatitude as rvl

##############################

def usage():
    print("usage: ", sys.argv[0] + "YYYYMMDD.XXX seconds")
    print("\t YYYYMMDD.XXX: experiment directory to process")
    print("\t seconds: number of seconds to integrate data for")

    sys.exit(2)


def main():
    
    try:
        # expName = glob.glob(os.path.abspath(sys.argv[1]))[0]
        # pd = sys.argv[2]
        expName = glob.glob(os.path.abspath('20190328.006'))[0]
        pd = 'lp1-fitcal-60'
    except:
        usage()
        
    dn = os.path.dirname(expName)
    iniF = os.path.join(dn,'config_vvelsLat.ini')        
    
    vf = rvl.vvelsLat(iniF + ','+expName+'/config_exp.ini',pd)
    vf.dovels()
    vf.replotFromOutput()


if __name__ == '__main__':

    main()	
