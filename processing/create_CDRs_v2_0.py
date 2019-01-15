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
from os.path import join, exists
from os import makedirs
import datetime
from calendar import monthrange
import matplotlib.pyplot as plt
import glob
import processing.CDRs as CDRs
import processing.FCDRs as FCDRs

def create_CDRs(instrument, satellite, start_date, end_date, version, version_comment,
                fcdr_path='/scratch/uni/u237/user_data/ihans/FCDR/easy/v2_0fv1_1_4/',
                regr_params_path='/scratch/uni/u237/users/tlang/UTH/radiative_transfer/uth_definition',
                regr_params_file='regr_params_from_fitted_thres_amsu.xml',
                cdr_path='/scratch/uni/u237/users/tlang/CDR/CDR_UTH/CDR_files',
                overwrite=False):

    """ Processes all monthly CDRs for the specified instrument on the specified
    satellite in the time period from start_date to end_date.
    
    Parameters:
        instrument (str): instrument name (e.g. 'MHS', 'AMSUB')
        satellite (str): satellite name (e.g. 'MetopA', 'Noaa18')
        start_date (datetime.datetime object): start of processing time period
        end_date (datetime.datetime object): end of processing time period
        version (str): CDR version (for filename, e.g. 2_0)
        version comment (str): comment on this version of CDRs
        fcdr_path (str): full path to FCDRs (optional)
        regr_params_path (str): full path to UTH regression parameters (optional)
        cdr_path (str): full path to CDR output directory (optional)
        overwrite (boolean): if True, existing CDR files are overwritten, if False
            the processing ends with an error message, if a CDR file already exists.
        
    """
    
    cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), 
                    str(start_date.year))
    #cdr_path = '.'
    # create directory for CDRs if it does not exist
    if not exists(cdr_fullpath):
            makedirs(cdr_fullpath)
    
    if instrument == 'HIRS':
        viewing_angles = [i for i in range(17, 41)] #+/- 12 around Nadir #center: 28,29
    elif instrument == 'SSMT2':
        viewing_angles = [i for i in range(9, 19)] # +/- 5 around Nadir #center: 13,14
    else:
        viewing_angles = [i for i in range(33, 58)] #+/- 13 around Nadir #center: 45,46
            
    lat_boundaries = [-30, 30]
    lon_boundaries = [-179, 180]
    
    # resolution of lat/lon grid
    resolution = 1.0
    
    #regr_params_file = 'regr_params_from_fitted_thres_{}.xml'.format(instrument.lower())
    regr_slopes, regr_intercepts = utils.regression_params_from_xml(join(regr_params_path, regr_params_file))
    
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
    #    # save with pickle
    #    outfile = 'CDR_{}_{}_{}.pkl'.format(monthly_mean.time_coverage_start.strftime('%Y-%m'), monthly_mean.instrument, monthly_mean.satellite)
    #    with open(outfile, 'wb') as output:
    #        pickle.dump(monthly_mean, output, pickle.HIGHEST_PROTOCOL)
        
        last_day_of_month = datetime.date(date.year, date.month, monthrange(date.year, date.month)[1])
        cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), str(date.year))
    #    cdr_path = '.'
        if not exists(cdr_fullpath):
            makedirs(cdr_fullpath)


if __name__ == '__main__':
    instruments = ['SSMT2']
    satellites = dict.fromkeys(instruments)
#    satellites['MHS'] = ['Noaa18']
#    satellites['AMSUB'] = ['Noaa19']
#    satellites['HIRS'] = ['Noaa11']
    satellites['SSMT2'] = ['f12']
    start_date = {i: {} for i in instruments}
    end_date = {i: {} for i in instruments}
#    start_date['MHS']['Metopb'] = datetime.date(2013, 2, 1)
#    end_date['MHS']['Metopb'] = datetime.date(2017, 10, 31)
#    start_date['MHS']['Noaa18'] = datetime.date(2012, 7, 1)
#    end_date['MHS']['Noaa18'] = datetime.date(2012, 7, 31)    
#    start_date['AMSUB']['Noaa19'] = datetime.date(2001, 1, 1)
#    end_date['AMSUB']['Noaa19'] = datetime.date(2001, 1, 2)
    start_date['SSMT2']['f12'] = datetime.date(1997, 7, 1)
    end_date['SSMT2']['f12'] = datetime.date(1997, 8, 1)
    
    for i in instruments:
        regr_params_file = 'regr_params_from_fitted_thres_{}.xml'.format(i.lower())
        for s in satellites[i]:
            create_CDRs(instrument=i, satellite=s, start_date=start_date[i][s],
                        end_date=end_date[i][s], version='0_0', version_comment='Test run for own use',
                        cdr_path='../CDR_files/Test/',
                        overwrite=True, regr_params_file=regr_params_file,
                        fcdr_path='/scratch/uni/u237/user_data/ihans/FCDR/easy/v4_1fv2_0_1/')
#            create_CDRs_year(instrument=i, satellite=s, year=2016, version_comment='Test', 
#                             cdr_path='.', overwrite=False)
        

#def create_CDRs_year(instrument, satellite, year, version_comment,
#                     fcdr_path='/scratch/uni/u237/user_data/ihans/FCDR/easy/v2_0fv1_1_4/',
#                     regr_params_path='/scratch/uni/u237/users/tlang/UTH/radiative_transfer/uth_definition',
#                     cdr_path='/scratch/uni/u237/users/tlang/CDR/CDR_UTH/processing/CDR_files',
#                     overwrite=False):
#    """ Processes all monthly CDRs for the specified instrument on the specified
#    satellite in one year.
#    
#    Parameters:
#        instrument (str): instrument name (e.g. 'MHS', 'AMSUB')
#        satellite (str): satellite name (e.g. 'MetopA', 'Noaa18')
#        year (int): year to be processed
#        version comment (str): comment on this version of CDRs
#        fcdr_path (str): full path to FCDRs (optional)
#        regr_params_path (str): full path to UTH regression parameters (optional)
#        cdr_path (str): full path to CDR output directory (optional)
#        overwrite (boolean): if True, existing CDR files are overwritten, if False
#            the processing ends with an error message, if a CDR file already exists.        
#    """
#        
#    start_date = datetime.datetime(year=year, month=1, day=1)
#    end_date = datetime.datetime(year=year, month=12, day=31)
#    cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), 
#                    str(start_date.year))
#    #cdr_path = '.'
#    # create directory for CDRs if it does not exist
#    if not exists(cdr_fullpath):
#            makedirs(cdr_fullpath)
#    
#    if instrument == 'HIRS':
#        viewing_angles = [i for i in range(17, 41)] #+/- 12 around Nadir #center: 28,29
#    else:
#        viewing_angles= [i for i in range(33, 58)] #+/- 13 around Nadir #center: 45,46
#            
#    lat_boundaries = [-30, 30]
#    lon_boundaries = [-179, 180]
#    
#    # resolution of lat/lon grid
#    resolution = 1.0
#    
#    regr_params_file = 'regr_params_from_fitted_thres_{}.xml'.format(instrument.lower())
#    regr_slopes, regr_intercepts = utils.regression_params_from_xml(join(regr_params_path, regr_params_file))
#    
#    date = start_date
#    day_delta = datetime.timedelta(days=1)
#    last_day_of_month = datetime.datetime(start_date.year, start_date.month, monthrange(start_date.year, start_date.month)[1])
#    
#    while date <= end_date:
#        print(date.month)
#        gridded_days = []
#        while date <= last_day_of_month and date <= end_date:
#            print(date.strftime("%Y-%m-%d"))
#            path = join(fcdr_path, satellite.lower(), str(date.year), '{m:02}'.format(m=date.month), '{d:02}'.format(d=date.day))
#            filenames_new = glob.glob(join(path, '*.nc'))
#            
#            fcdrs = []
#            
#            if filenames_new:
#                for file in filenames_new:
#                    f = FCDRs.FCDR.fromNetcdf(
#                             fcdr_path, file, viewing_angles)
#                    f.calcUTH(regr_slopes, regr_intercepts)
#                    fcdrs.append(f)
#            
#                # gridded daily CDRs
#                gridded_days.append(CDRs.CDR.GriddedCDRFromFCDRs(fcdrs, lat_boundaries=lat_boundaries, 
#                                                            lon_boundaries=lon_boundaries, resolution=resolution))
#            
#            date += day_delta
#    
#        # monthly_average
#        if gridded_days:
#            monthly_mean = CDRs.CDR.AveragedCDRFromCDRs(gridded_days)
#        # write to NetCDF
#            if not monthly_mean.is_empty:
#                ds = monthly_mean.toNetCDF(cdr_fullpath, comment_on_version=version_comment, overwrite=overwrite)
#        
#        last_day_of_month = datetime.date(date.year, date.month, monthrange(date.year, date.month)[1])
#        cdr_fullpath = join(cdr_path, str(instrument.upper()), str(satellite.upper()), str(date.year))
#
#        if not exists(cdr_fullpath):
#            makedirs(cdr_fullpath)