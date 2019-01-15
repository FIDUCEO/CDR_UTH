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
version = sys.argv[4]
version_comment = sys.argv[5]
version_comment = version_comment.replace('_', ' ')
fcdr_path = sys.argv[6]
cdr_path = sys.argv[7]
overwrite = bool(sys.argv[8])

# ---------------------------SCALING COEFFICIENTS ------------------------- #
# location and names of files containing the UTH scaling coefficients:
regr_params_path = '/scratch/uni/u237/users/tlang/UTH/radiative_transfer/uth_definition'
regr_params_file = 'regr_params_from_fitted_thres_{}.xml'.format(instrument.lower()) #'regr_params_old_definition_amsu.xml'
# read scaling coefficients
regr_slopes, regr_intercepts = utils.regression_params_from_xml(
        join(regr_params_path, regr_params_file)
        )
# -----------------------------START/ END DATES ----------------------------#
start_date = datetime.date(year=year, month=1, day=1)
end_date = datetime.date(year=year, month=12, day=31)

# ----------------------------PIXEL SELECTION ------------------------------#
# Satellite viewing angles used for the CDR
# only pixels close to nadir are chosen (see CDR product user guide for more information)
if instrument == 'HIRS':
    # HIRS: +/- 12 pixels around Nadir are used, center pixels: 28,29
    viewing_angles = [i for i in range(17, 41)] 
elif instrument == 'SSMT2':
    # SSMT2: +/- 5 pixels around Nadir are used, center pixels: 13,14
    viewing_angles = [i for i in range(9, 19)] 
elif instrument == 'AMSUB' or instrument == 'MHS':
    # AMSU-B: +/- 13 pixels around Nadir are used, center pixels: 45,46
    viewing_angles = [i for i in range(33, 58)] 
        
# -----------------------------LAT/LON GRID --------------------------------#
# grid boundaries: use only latitudes from 30°S to 30°N (tropical region)
lat_boundaries = [-30, 30]
lon_boundaries = [-179, 180]
# resolution of lat/lon grid
resolution = 1.0

# ----------------------------CDR PROCESSING-------------------------------#
date = start_date
day_delta = datetime.timedelta(days=1)
last_day_of_month = datetime.date(
        start_date.year, 
        start_date.month, 
        monthrange(start_date.year, start_date.month)[1]
        )
cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), 
                str(start_date.year))

# go through all days of the year
print('Year: {}'.format(year))
while date <= end_date:
    print('Month: {}'.format(date.month))
    gridded_days = []
    while date <= last_day_of_month and date <= end_date:
        print('Date: {}'.format(date.strftime("%Y-%m-%d")))
        # Full FCDR path:
        path = join(fcdr_path, 
                    satellite.lower(), 
                    str(date.year), 
                    '{m:02}'.format(m=date.month), 
                    '{d:02}'.format(d=date.day)
                    )
        # Names of all FCDR files of this day
        filenames_new = glob.glob(join(path, '*.nc')) 
        # go through all files of this day (if exisiting) and read in FCDRs
        fcdrs = []  
        print('FCDR files:')
        if filenames_new:
            for file in filenames_new:
                print(file)
                f = FCDRs.FCDR.fromNetcdf(
                         fcdr_path, 
                         file, 
                         viewing_angles
                         )
                # calculate UTH from brightness temperature using regression coefficients
                f.calcUTH(regr_slopes, regr_intercepts)
                fcdrs.append(f)
        
            # Perform pixel aggregation and averaging to get gridded daily CDRs
            gridded_days.append(CDRs.CDR.GriddedCDRFromFCDRs(
                    fcdrs, 
                    lat_boundaries=lat_boundaries, 
                    lon_boundaries=lon_boundaries, 
                    resolution=resolution))
        
        date += day_delta

    # Average all days of one month to get gridded monthly CDRs
    if gridded_days:
        monthly_mean = CDRs.CDR.AveragedCDRFromCDRs(gridded_days)
    # Write gridded monthly CDR fields to NetCDF file
        if not monthly_mean.is_empty:
            if not exists(cdr_fullpath):
                makedirs(cdr_fullpath)
            monthly_mean.toNetCDF(
                    cdr_fullpath, 
                    version=version, 
                    comment_on_version=version_comment, 
                    overwrite=overwrite)
    
    last_day_of_month = datetime.date(date.year, date.month, monthrange(date.year, date.month)[1])

