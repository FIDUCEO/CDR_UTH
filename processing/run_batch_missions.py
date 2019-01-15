#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 10:31:32 2018

@author: u300740

Script to run batch jobs.
One job is created for every year of CDR. The script create_CDR_year.py is called.
Important: 
    specify version, FCDR and CDR paths, instruments, satellites as well as
    start and end year for each satellite!
    For further CDR preferences, like grid broundaries and resolution see 
    create_CDR_year.py.
"""

import numpy as np
import os

# ----------------------------PLEASE SPECIFY-------------------------- #

# CDR version (for CDR filename)
version = '1_1'
# comment on this CDR version (will be included in each CDR file)
version_comment = 'Version_1.1_of_UTH_CDR_based_on_Microwave_FCDR_Version_4.1.'
# path to FCDR input
fcdr_path = '/scratch/uni/u237/user_data/ihans/FCDR/easy/v4_1fv2_0_1/'
# path for CDR output
cdr_path = '/scratch/uni/u237/users/tlang/CDR/CDR_UTH/CDR_files/fromFCDRv4_1'
# True, if existing CDR files with the same name should be overwritten
overwrite = True
# specifiy instruments, satellites, start and end year of the mission
instruments = ['SSMT2']
satellites = dict.fromkeys(instruments)
#satellites['MHS'] = ['Noaa18', 'Noaa19', 'Metopa', 'Metopb']
#satellites['AMSUB'] = #['Noaa15', 'Noaa16', 'Noaa17']
satellites['SSMT2'] = ['f11', 'f12', 'f14', 'f15']
#satellites['MHS'] = ['Noaa19']
start_year = {i: {} for i in instruments}
end_year = {i: {} for i in instruments}
#start_year['MHS']['Metopa'] = 2007
#end_year['MHS']['Metopa'] = 2018
#start_year['MHS']['Metopb'] = 2013
#end_year['MHS']['Metopb'] = 2018
#start_year['MHS']['Noaa19'] = 2009
#end_year['MHS']['Noaa19'] = 2018
#start_year['MHS']['Noaa18'] = 2005
#end_year['MHS']['Noaa18'] = 2018
#start_year['AMSUB']['Noaa15'] = 1999
#end_year['AMSUB']['Noaa15'] = 2011
#start_year['AMSUB']['Noaa16'] = 2001
#end_year['AMSUB']['Noaa16'] = 2014
#start_year['AMSUB']['Noaa17'] = 2002
#end_year['AMSUB']['Noaa17'] = 2014
start_year['SSMT2']['f11'] = 1994
end_year['SSMT2']['f11'] = 1995
start_year['SSMT2']['f12'] = 1994
end_year['SSMT2']['f12'] = 2001
start_year['SSMT2']['f14'] = 1997
end_year['SSMT2']['f14'] = 2005
start_year['SSMT2']['f15'] = 2000
end_year['SSMT2']['f15'] = 2005

# -----------------------------BATCH JOBS----------------------------- #
# create one batch job for each year of each satellite mission
for i in instruments:
    for s in satellites[i]:
        for year in np.arange(start_year[i][s], end_year[i][s]+1):
            jobname = i+s+str(year)
            cmd = 'sbatch_simple -x "ctc[132-134]" {n} 16 python create_CDR_year.py {i} {s} {y} {v} {vc} {fp} {cp} {o}'.format(
                    n=jobname, 
                    i=i, 
                    s=s, 
                    y=year, 
                    v=version,
                    vc=version_comment,
                    fp=fcdr_path,
                    cp=cdr_path,
                    o=int(overwrite))
            os.system(cmd)

