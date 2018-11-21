#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 10:31:32 2018

@author: u300740

Script to run batch jobs.
One job is created for every year of CDR.
Important: specify version_comment, instruments and satellites!
"""

import numpy as np
import os
#from create_CDRs_v2_0 import create_CDRs

# Comment on this CDR version
version_comment = 'Second test version for own use.'
# True, if existing CDR files with the same name should be overwritten
overwrite = True
# specifiy instruments, satellites, start and end year of the mission
instruments = ['AMSUB']
satellites = dict.fromkeys(instruments)
#satellites['MHS'] = ['Noaa18']
satellites['AMSUB'] = ['Noaa15', 'Noaa16']
start_year = {i: {} for i in instruments}
end_year = {i: {} for i in instruments}
#start_year['MHS']['Metopb'] = 2013
#end_year['MHS']['Metopb'] = 2017
#start_year['MHS']['Noaa18'] = 2007
#end_year['MHS']['Noaa18'] = 2017
#start_year['AMSUB']['Noaa17'] = 2002
#end_year['AMSUB']['Noaa17'] = 2013
start_year['AMSUB']['Noaa15'] = 2003
end_year['AMSUB']['Noaa15'] = 2003
start_year['AMSUB']['Noaa16'] = 2003
end_year['AMSUB']['Noaa16'] = 2003

for i in instruments:
    for s in satellites[i]:
        for year in np.arange(start_year[i][s], end_year[i][s]+1):
            jobname = i+s+str(year)
            cmd = 'sbatch_simple -x "ctc[132-134]" {n} 16 python create_CDR_year.py {i} {s} {y} {c} {o}'.format(n=jobname, i=i, s=s, y=year, c=version_comment, o=int(overwrite))
            os.system(cmd)

