#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 13:21:09 2018

@author:Theresa Lang

Creates one year of CDRs.
This script is meant to be run from the command line. The following parameters
have to be specified in the call:
    1. instrument: instrument
    2. satellite: satellite/platform carriing the instrument
    3. year: year of data to be processed
    4. version comment: a comment on the version of the CDR
    5. overwrite: 1 if existing CDR files should be overwritten. 
"""
import sys
import numpy as np
import processing.utils as utils
from os.path import join, exists
from os import makedirs
import datetime
from calendar import monthrange
import matplotlib.pyplot as plt
import glob
import processing.CDRs as CDRs
import processing.FCDRs as FCDRs

instrument = sys.argv[1]
satellite = sys.argv[2]
year = int(sys.argv[3])
version_comment = sys.argv[4]
overwrite = bool(sys.argv[5])

print(overwrite)

fcdr_path='/scratch/uni/u237/user_data/ihans/FCDR/easy/v2_0fv1_1_4/'
regr_params_path='/scratch/uni/u237/users/tlang/UTH/radiative_transfer/uth_definition'
cdr_path='/scratch/uni/u237/users/tlang/CDR/CDR_UTH/processing/CDR_files'

start_date = datetime.datetime(year=year, month=1, day=1)
end_date = datetime.datetime(year=year, month=12, day=31)
cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), 
                str(start_date.year))
#cdr_path = '.'
# create directory for CDRs if it does not exist
if not exists(cdr_fullpath):
        makedirs(cdr_fullpath)

if instrument == 'HIRS':
    viewing_angles = [i for i in range(17, 41)] #+/- 12 around Nadir #center: 28,29
else:
    viewing_angles= [i for i in range(33, 58)] #+/- 13 around Nadir #center: 45,46
        
lat_boundaries = [-30, 30]
lon_boundaries = [-179, 180]

# resolution of lat/lon grid
resolution = 1.0

regr_params_file = 'regr_params_from_fitted_thres_{}.xml'.format(instrument.lower())
regr_slopes, regr_intercepts = utils.regression_params_from_xml(join(regr_params_path, regr_params_file))

date = start_date
day_delta = datetime.timedelta(days=1)
last_day_of_month = datetime.datetime(start_date.year, start_date.month, monthrange(start_date.year, start_date.month)[1])

while date <= end_date:
    print(date.month)
    gridded_days = []
    while date <= last_day_of_month and date <= end_date:
        print(date.strftime("%Y-%m-%d"))
        path = join(fcdr_path, satellite.lower(), str(date.year), '{m:02}'.format(m=date.month), '{d:02}'.format(d=date.day))
        filenames_new = glob.glob(join(path, '*.nc'))
        
        fcdrs = []
        
        if filenames_new:
            for file in filenames_new:
                f = FCDRs.FCDR.fromNetcdf(
                         fcdr_path, file, viewing_angles)
                f.calcUTH(regr_slopes, regr_intercepts)
                fcdrs.append(f)
        
            # gridded daily CDRs
            gridded_days.append(CDRs.CDR.GriddedCDRFromFCDRs(fcdrs, lat_boundaries=lat_boundaries, 
                                                        lon_boundaries=lon_boundaries, resolution=resolution))
        
        date += day_delta

    # monthly_average
    if gridded_days:
        monthly_mean = CDRs.CDR.AveragedCDRFromCDRs(gridded_days)
    # write to NetCDF
        if not monthly_mean.is_empty:
            ds = monthly_mean.toNetCDF(cdr_fullpath, comment_on_version=version_comment, overwrite=overwrite)
    
    last_day_of_month = datetime.date(date.year, date.month, monthrange(date.year, date.month)[1])
    cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), str(date.year))

    if not exists(cdr_fullpath):
        makedirs(cdr_fullpath)