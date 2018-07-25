#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 09:20:34 2018

@author: Theresa Lang
"""
import numpy as np
from os.path import join
from typhon.arts import xml
import glob

def get_filenames_daily(path, instrument, satellite, year, month, day):
    """ Returns all filenames of files containing FCDRs for a certain instrument,
    satellite and date.
    
    Parameters:
        path (str): path to directory with FCDRs 
        instrument (str): name of instrument, e.g. MHS or Amsub 
        satellite (str): name of satellite, e.g. Metopb
        year (int): year 
        month (int): month, e.g. 1 for january  
        day (int): day 
    """
    
    files = 'FIDUCEO_FCDR_L1C_{i}_{s}_{y}{m:02}{d:02}*.nc'.format(
            i=instrument.upper(), s=satellite.upper(), y=year, m=month, d=day)
    return glob.glob(join(path, files))

def regression_params_from_xml(xmlfile):
    """ Read regression parameters for UTH calculation from xml file. Returns
    slopes and intercepts.
    
    Parameters:
        xmlfile (str): full path to xml file containing regression parameters
    """
    
    params = xml.load(xmlfile)
    slopes = params[0]
    intercepts = params[1]
    
    return slopes, intercepts

def get_grid_lons_lats(lon_limits, lat_limits, resolution):
    """ Returns vectors of longitudes and latitudes.
    
    Parameters:
        lon_limits (list): lower and upper limit for longitues in deg (e.g. [-180, 180])
        lat_limits (list): lower and upper limit for latitudes in deg (e.g. [-30, 30])
        resolution (numeric): grid resolution in degrees 
    """
    
    lons = np.arange(lon_limits[0], lon_limits[1] + resolution, resolution)
    lats = np.arange(lat_limits[0], lat_limits[1] + resolution, resolution)
    
    return lons, lats