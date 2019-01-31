#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:16:09 2018

@author: Theresa Lang

This script reads in FCDR data and uses it to create gridded and averaged UTH 
CDRs. It is intended for testing and was not used in the final CDR processing.
"""
import numpy as np
import processing.utils as utils
from os.path import join, exists
from os import makedirs
import datetime
from calendar import monthrange
import glob
import processing.CDRs as CDRs
import processing.FCDRs as FCDRs

def create_CDRs(instrument, satellite, start_date, end_date, version, version_comment,
                fcdr_path, cdr_path, regr_params_path, overwrite=False):
    """ Processes all monthly CDRs for the specified instrument on the specified
    satellite in the time period from start_date to end_date.
    
    Parameters:
        instrument (str): instrument name (e.g. 'MHS', 'AMSUB')
        satellite (str): satellite name (e.g. 'MetopA', 'Noaa18')
        start_date (datetime.datetime object): start of processing time period
        end_date (datetime.datetime object): end of processing time period
        version (str): CDR version (for filename, e.g. 2_0)
        version comment (str): comment on this version of CDRs
        fcdr_path (str): full path to FCDRs
        cdr_path (str): full path to CDR output directory 
        regr_params_path (str): full path to UTH regression parameters 
        overwrite (boolean): if True, existing CDR files are overwritten, if False
            the processing ends with an error message, if a CDR file already exists.
        
    """
    
    cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), 
                    str(start_date.year))
    #cdr_path = '.'u300740
    # create directory for CDRs if it does not exist
    if not exists(cdr_fullpath):
            makedirs(cdr_fullpath)
    
    if instrument == 'HIRS':
        viewing_angles = [i for i in range(19, 39)] #+/- 10 around Nadir #center: 28,29
    elif instrument == 'SSMT2':
        viewing_angles = [i for i in range(10, 20)] # +/- 5 around Nadir #center: 13,14
    else:
        viewing_angles = [i for i in range(32, 60)] #+/- 14 around Nadir #center: 45,46
            
    lat_boundaries = [-30, 30]
    lon_boundaries = [-179, 180]
    
    # resolution of lat/lon grid
    resolution = 1.0
    
    #regr_params_file = 'regr_params_from_fitted_thres_{}.xml'.format(instrument.lower())
    regr_slopes, regr_intercepts = utils.regression_params_from_xml(regr_params_path)
    
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
            print(path)
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
                ds = monthly_mean.toNetCDF(cdr_fullpath, version, comment_on_version=version_comment, overwrite=overwrite)
        
        last_day_of_month = datetime.date(date.year, date.month, monthrange(date.year, date.month)[1])
        cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), str(date.year))

        if not exists(cdr_fullpath):
            makedirs(cdr_fullpath)


if __name__ == '__main__':
    
    cdr_path = '../CDR_files/Test/'
    fcdr_path = '/scratch/uni/u237/data/fiduceo-fcdr/easy/v4_1fv2_0_1'
    version='0_0'
    version_comment = 'Test run for own use'
    
    instruments = ['SSMT2']
    satellites = dict.fromkeys(instruments)
#    satellites['MHS'] = ['Noaa18']
#    satellites['AMSUB'] = ['Noaa19']
#    satellites['HIRS'] = ['Noaa11']
    satellites['SSMT2'] = ['f15']
    start_date = {i: {} for i in instruments}
    end_date = {i: {} for i in instruments}
#    start_date['MHS']['Metopb'] = datetime.date(2013, 2, 1)
#    end_date['MHS']['Metopb'] = datetime.date(2017, 10, 31)
#    start_date['MHS']['Noaa18'] = datetime.date(2012, 7, 1)
#    end_date['MHS']['Noaa18'] = datetime.date(2012, 7, 31)    
#    start_date['AMSUB']['Noaa19'] = datetime.date(2001, 1, 1)
#    end_date['AMSUB']['Noaa19'] = datetime.date(2001, 1, 2)
    start_date['SSMT2']['f15'] = datetime.date(2003, 1, 1)
    end_date['SSMT2']['f15'] = datetime.date(2003, 1, 30)
    
    for i in instruments:
        regr_params_file = 'regression_parameters/regr_params_{}.xml'.format(i.lower())
        for s in satellites[i]:
            create_CDRs(
                    instrument=i, 
                    satellite=s, 
                    start_date=start_date[i][s],
                    end_date=end_date[i][s], 
                    version=version, 
                    version_comment=version_comment,
                    fcdr_path=fcdr_path,
                    cdr_path=cdr_path,
                    regr_params_path=regr_params_file,
                    overwrite=True, 
                    )
