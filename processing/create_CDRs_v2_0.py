#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:16:09 2018

@author: Theresa Lang

This script reads in FCDR data and uses it to create monthly averaged UTH CDRs.
"""
import numpy as np
import processing.utils as utils
import pickle
from os.path import join
import datetime
from calendar import monthrange
import matplotlib.pyplot as plt
import glob
import processing.CDRs as CDRs
import processing.FCDRs as FCDRs
import processing.plots as plots

version_comment = 'Test run'

#fcdr_path = '/scratch/uni/u237/users/tlang/CDR/FCDRs'
fcdr_path = '/scratch/uni/u237/user_data/ihans/FCDR/easy/v2_0fv1_1_4/'
#fcdr_path = '/scratch/uni/u237/user_data/ihans/FCDR/HIRS'
regr_params_path = '/scratch/uni/u237/users/tlang/UTH/radiative_transfer/uth_definition'

start_date = datetime.date(2007, 1, 1)
end_date = datetime.date(2007, 1, 31)
instrument = 'MHS'
satellite = 'Noaa18'

#cdr_path = join('CDR_files', str(instrument.upper()), str(satellite.upper()), 
#                str(start_date.year))
cdr_path = '.'

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

filenames = []
dates = []
date = start_date
day_delta = datetime.timedelta(days=1)
last_day_of_month = datetime.date(start_date.year, start_date.month, monthrange(start_date.year, start_date.month)[1])

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
    monthly_mean = CDRs.CDR.AveragedCDRFromCDRs(gridded_days)
    # write to NetCDF
    if not monthly_mean.is_empty:
        ds = monthly_mean.toNetCDF(cdr_path, comment_on_version=version_comment)
#    # save with pickle
    outfile = 'CDR_{}_{}_{}.pkl'.format(monthly_mean.time_coverage_start.strftime('%Y-%m'), monthly_mean.instrument, monthly_mean.satellite)
    with open(outfile, 'wb') as output:
        pickle.dump(monthly_mean, output, pickle.HIGHEST_PROTOCOL)
    
    last_day_of_month = datetime.date(date.year, date.month, monthrange(date.year, date.month)[1])
#    cdr_path = join('CDR_files', str(instrument.upper()), str(satellite.upper()), str(date.year))
    cdr_path = '.'