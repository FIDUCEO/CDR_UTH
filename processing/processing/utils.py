#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 09:20:34 2018

@author: Theresa Lang
"""
import numpy as np
import pandas as pd
from os.path import join
from typhon.arts import xml
import glob
import datetime

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

def getFileInfo(filename):
    """ Extracts information about instrument and time from FCDR filename.
    
    Parameters:
        filename (str): filename of FCDR
    """
    info = filename.split('_')
    
    file_info = {}
    file_info['instrument'] = ''.join([i for i in info[3].upper() if not i.isdigit()])
    if file_info['instrument'] == 'SSMT':
        file_info['instrument'] = 'SSMT2'
    file_info['satellite'] = info[4].upper()
    start_time_str = info[5]
    end_time_str = info[6]
    file_info['start_time'] = datetime.datetime.strptime(start_time_str, '%Y%m%d%H%M%S')
    file_info['end_time'] = datetime.datetime.strptime(end_time_str, '%Y%m%d%H%M%S')
    
    return file_info

def collectFCDRData(FCDRs, u_types, uth_channel=3):
    """ Collects Data from several FCDRs (e.g. one day of data) in one dictionary.
    
    Parameters:
        FCDRs (list): List of FCDR objects
        u_types (list): List of uncertainty types that should be extracted 
            (e.g. ['structured', 'independent', 'common'])
        uth_channel (int): channel number used for uth-calculation, 
            brightness temperature will only be extracted for this channel
    """
    branches = ['ascending', 'descending']
    files = []
    latitude = {b: np.array([], dtype=np.float16) for b in branches}
    longitude = {b: np.array([], dtype=np.float16) for b in branches}
    latitude_unfiltered = {b: np.array([], dtype=np.float16) for b in branches}
    longitude_unfiltered = {b: np.array([], dtype=np.float16) for b in branches}
    brightness_temp = {b: np.array([]) for b in branches}
    brightness_temp_unfiltered = {b: np.array([]) for b in branches} 
    uth = {b: np.array([]) for b in branches}
    u_Tb = {t: {b: np.array([]) for b in branches} for t in u_types}
    u_Tb_unfiltered = {t: {b: np.array([]) for b in branches} for t in u_types}
    u_uth = {t: {b: np.array([]) for b in branches} for t in u_types}
    scanline = {b: np.array([]) for b in branches}
    scanline_max = {b: 0 for b in branches}
    scanline_unfiltered = {b: np.array([]) for b in branches}
    scanline_max_unfiltered = {b: 0 for b in branches}
    second_of_day = {b: np.array([]) for b in branches} 
    acquisition_time = {b: np.array([]) for b in branches} 
#    latitude_diff = {b: np.array([], dtype=np.float16) for b in branches}
#    longitude_diff = {b: np.array([], dtype=np.float16) for b in branches}
    
    node_mask = {}
    for f in FCDRs:
        files.append(f.file)
        # get all masks
        if not hasattr(f, 'node_mask'):
            f.generateNodeMask()
        if not hasattr(f, 'total_mask'):
            f.generateTotalMask()
        if not hasattr(f, 'cloud_mask'):
            f.generateCloudMask()
        if not hasattr(f, 'quality_and_issue_mask'):
            f.generateQualityAndIssueMask();
        
        node_mask['ascending'] = f.node_mask
        node_mask['descending'] = ~f.node_mask
        total_mask = f.total_mask
        quality_and_issue_mask = f.quality_and_issue_mask
#        diff_mask = np.logical_and(np.logical_not(quality_and_issue_mask), total_mask)

        # distinguish between ascending and descending node and apply all masks
        for b in branches: 
            if f.latitude[np.logical_and(~total_mask, node_mask[b])].size:
                latitude[b] = np.append(latitude[b], f.latitude[np.logical_and(~total_mask, node_mask[b])])
                longitude[b] = np.append(longitude[b], f.longitude[np.logical_and(~total_mask, node_mask[b])])
                latitude_unfiltered[b] = np.append(latitude_unfiltered[b], f.latitude[np.logical_and(~quality_and_issue_mask, node_mask[b])])
                longitude_unfiltered[b] = np.append(longitude_unfiltered[b], f.longitude[np.logical_and(~quality_and_issue_mask, node_mask[b])])
                brightness_temp[b] = np.append(brightness_temp[b], f.brightness_temp[uth_channel][np.logical_and(~total_mask, node_mask[b])])
                brightness_temp_unfiltered[b] = np.append(brightness_temp_unfiltered[b], f.brightness_temp[uth_channel][np.logical_and(~quality_and_issue_mask, node_mask[b])])
                uth[b] = np.append(uth[b], f.uth[np.logical_and(~total_mask, node_mask[b])])
                second_of_day[b] = np.append(second_of_day[b], f.second_of_day[np.logical_and(~total_mask, node_mask[b])])
                acquisition_time[b] = np.append(acquisition_time[b], f.acquisition_time[np.logical_and(~total_mask, node_mask[b])])
#                latitude_diff[b] = np.append(latitude_diff[b], f.latitude[np.logical_and(node_mask[b], diff_mask)])
#                longitude_diff[b] = np.append(longitude_diff[b], f.longitude[np.logical_and(node_mask[b], diff_mask)])
                
                # scanline information (needed for structured uncertainties)
                scanline[b] = np.append(scanline[b], scanline_max[b] + f.scanline[np.logical_and(~total_mask, node_mask[b])])
                scanline_unfiltered[b] = np.append(scanline_unfiltered[b], scanline_max_unfiltered[b] + f.scanline[np.logical_and(~quality_and_issue_mask, node_mask[b])])
                scanline_max[b] = np.max(scanline[b])
                scanline_max_unfiltered[b] = np.max(scanline_unfiltered[b])
                # get 3 types of uncertainties
                for t in u_types:
                    u_Tb[t][b] = np.append(u_Tb[t][b], f.u_Tb[t][uth_channel][np.logical_and(~total_mask, node_mask[b])])
                    u_Tb_unfiltered[t][b] = np.append(u_Tb_unfiltered[t][b], f.u_Tb[t][uth_channel][np.logical_and(~quality_and_issue_mask, node_mask[b])])
                    u_uth[t][b] = np.append(u_uth[t][b], f.u_uth[t][np.logical_and(~total_mask, node_mask[b])])
            
    
    collected_data = {}
    collected_data['latitude'] = latitude
    collected_data['longitude'] = longitude
    collected_data['brightness_temp'] = brightness_temp
    collected_data['uth'] = uth
    collected_data['second_of_day'] = second_of_day
    collected_data['acquisition_time'] = acquisition_time
    collected_data['scanline'] = scanline    
    for t in u_types:
        collected_data['u_Tb_{}'.format(t)] = u_Tb[t]
        collected_data['u_uth_{}'.format(t)] = u_uth[t]
    
    collected_data_unfiltered = {}
    collected_data_unfiltered['latitude'] = latitude_unfiltered
    collected_data_unfiltered['longitude'] = longitude_unfiltered
    collected_data_unfiltered['brightness_temp'] = brightness_temp_unfiltered 
    collected_data_unfiltered['scanline'] = scanline_unfiltered
    for t in u_types:
        collected_data_unfiltered['u_Tb_{}'.format(t)] = u_Tb_unfiltered[t]
    
#    collected_data_diff = {}
#    collected_data_diff['latitude'] = latitude_diff
#    collected_data_diff['longitude'] = longitude_diff
    
    return collected_data, collected_data_unfiltered, files

def binData(data, lat_bins, lon_bins, lat_centers, lon_centers):
    """ Bins data in the dataframe data (which contains the attributes 
    'latitude' and 'longitude') to latitude/ longitude bins. 
    
    Parameters:
        data (Dataframe): Pandas dataframe containing the attributes 'latitude'
            and 'longitude'
        lat_bins: latitude bin edges
        lon_bins: longitude bin edges
        lat_centers: latitude bin centers
        lon_centers: longitude bin centers
    """
    lon2ind = lambda lon: lon_centers.index(lon)
    lat2ind = lambda lat: lat_centers.index(lat)
    
    data['lat_bin'] = pd.cut(data.latitude, lat_bins, labels=lat_centers).apply(lat2ind)
    data['lon_bin'] = pd.cut(data.longitude, lon_bins, labels=lon_centers).apply(lon2ind)
    data = data.drop(['latitude', 'longitude'], 1)

    return data

def flattenDict(datadict, branch):
    """ Returns a dictionary containing data for one branch of the original 
    dictionary only.
    
    Parameters:
        datadict (dict): nested dictionary, every sub-dictionary has branch as
            one of its keys
        branch (string): key of sub-dictionaries that you want to extract
    """
    flattend_dict = {}
    for k, v in datadict.items():
        flattend_dict[k] = v[branch]
    
    return flattend_dict

def getLatLonBins(lat_boundaries, lon_boundaries, resolution):
    """ Takes upper and lower boundaries for longitudes and latitudes and grid
    resolution as input and returns bin edges and centers for latitudes and longitudes.
    
    Parameters:
        lat_boundaries (list): upper and lower boundaries for latitudes
        lon_boundaries (list): upper and lower boundaries for longitudes
        resolution (numeric): Lat/Lon-grid resolution in degree
    """
    lat_bins = np.arange(lat_boundaries[0] - 0.5 * resolution, lat_boundaries[-1] + 0.5 * resolution + 1., resolution)
    lat_centers = [(a + b) / 2 for a, b in zip(lat_bins, lat_bins[1:])]
    lon_bins = np.arange(lon_boundaries[0] - 0.5 * resolution, lon_boundaries[-1] + 0.5 * resolution + 1, resolution)
    lon_centers = [(a + b) / 2 for a, b in zip(lon_bins, lon_bins[1:])]
    
    return lat_bins, lon_bins, lat_centers, lon_centers

def getSecondOfDay(acquisition_time):
    """ Transforms pixel acquisition time to time in seconds of day.
    
        acquisition_time (array): vector containing (e.g. scanline) acquisition times
    """
    start_time = datetime.datetime.utcfromtimestamp(acquisition_time[0])
    midnight = datetime.datetime(start_time.year, start_time.month, start_time.day, 0, 0, tzinfo=datetime.timezone.utc).timestamp()            
    second_of_day = acquisition_time - midnight
    second_of_day[np.where(second_of_day >= 86400)] = second_of_day[np.where(second_of_day >= 86400)] - 86400
    
    return second_of_day

def getFirstNotEmptyCDR(CDRs):
    """ Returns the first CDR which is not empty from a list of CDRs.
    Parameters:
        CDRs: List of CDR objects
    """
    if (np.any([not c.is_empty for c in CDRs])):
        return next(c for c in CDRs if not c.is_empty)
    else:
        return None

def getLastNotEmptyCDR(CDRs):
    """ Returns the last CDR which is not empty from a list of CDRs.
    Parameters:
        CDRs: List of CDR objects
    """
    if (np.any([not c.is_empty for c in CDRs])):
        return next(c for c in CDRs[::-1] if not c.is_empty)
    else:
        return None

def getCDRQualityMask(observation_counts, overpass_counts):
    """ Creates a quality bitmask for the monthly mean CDR based on observation
    and overpass counts.
    
    Bit 2 is set to 1 (=use with caution), if there are less than 6 overpasses 
    in a grid cell there are less than 150 observations in a grid cell
    
    Parameters:
        observation_counts: number of pixels that were used to retrieve UTH
            in every grid cell in this month
        overpass_counts: number of satellite overpasses in every grid cell
            in one month
    """
    bitmask = np.zeros(observation_counts.shape)
    
    bitmask[overpass_counts < 6] = bitmask[overpass_counts < 6] + 2
    bitmask[observation_counts < 150] = bitmask[observation_counts < 150] + 2
    
    return bitmask

def getCorrMatrix(scanlines, correlation_vector):
    """ Calculates the correlation matrix for structured uncertainties from the
    scanline information of all pixels in one grid cell.
    
    Parameters:
        scanlines (array of integers): Scanlines of all pixels in the
            grid cell
        correlation_vector (array): Correlation vector containing correlations
            dependent on scanline differences (first entry: correlation for 
            pixels in same scanline, second entry: correlation for pixels 1 
            scanline apart, etc.)
    """
    # maximum difference of scanlines that still leads to a correlation
    max_diff = len(correlation_vector) - 1
    # when the scanline difference is larger than max_diff, the correlation is 0
    correlation_vector = np.append(correlation_vector, 0)
    # calculate scanline differences between all pixels in grid cell
    scanlines = np.reshape(scanlines, (len(scanlines), 1))
    scanline_diffs = np.abs(scanlines - scanlines.T)
    # differences larger than max_diff are set to max_diff + 1
    scanline_diffs[scanline_diffs > max_diff] = max_diff + 1
    # scanline differences are used as indices for correlation_vector to get 
    # the correlations as a function of scanline difference
    corr_matrix = np.take(correlation_vector, scanline_diffs) 

    return corr_matrix
    
    