#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 16:52:03 2018

@author: u300740
"""
import numpy as np
import os.path
from netCDF4 import Dataset
from pyresample import kd_tree, geometry
import pandas as pd
from scipy.sparse import csr_matrix, csc_matrix, diags, bmat, block_diag
import matplotlib.pyplot as plt
import datetime

class FCDR:
    def __init__(self, latitudes=None, longitudes=None, viewing_angles=None,
                 brightness_temp=None, channels=None, scanline=None, quality_mask=None, 
                 u_Tb=None, u_uth=None, quality_issue=None, start_time=None, 
                 end_time=None, acquisition_time=None,instrument=None, 
                 satellite=None, file=None, uth=None):
        self.file = file
        self.latitudes = latitudes
        self.longitudes = longitudes
        self.viewing_angles = viewing_angles
        self.brightness_temp = brightness_temp
        self.channels = channels
        self.scanline = scanline
        self.quality_mask = quality_mask
        self.u_Tb = u_Tb
        self.u_uth = u_uth
        self.quality_issue = quality_issue
        self.start_time = start_time
        self.end_time = end_time
        self.acquisition_time = acquisition_time
        self.instrument = instrument
        self.satellite = satellite
        self.uth = uth
    
    @classmethod    
    def from_netcdf(cls, fcdr_path, filename, viewing_angles=None):
        """ Read FCDR from NetCDF file.
        
        Parameters:
            fcdr_path (str): path to directory with NetCDF file
            filename (str): name of NetCDF file 
            instrument (str): instrument, e.g. Amsub or MHS
            viewing_angles (list): list of viewing angles/ pixels that should 
                be used (e.g. [44, 45] for innermost angles only)
            satellite (str): satellite, e.g. Metopb
        """
        ret = cls()
        # read data from NetCDF file
        f = Dataset(os.path.join(fcdr_path, filename))
        
        # get information from file name
        file = os.path.basename(filename)
        ret.file = file 
        file_info = file.split('_')
        
        instrument = file_info[3].upper()
        ret.instrument = ''.join([i for i in instrument if not i.isdigit()])
        ret.satellite = file_info[4].upper()

        start_time_str = file_info[5]
        end_time_str = file_info[6]
        start_time = datetime.datetime.strptime(start_time_str, '%Y%m%d%H%M%S')
        end_time = datetime.datetime.strptime(end_time_str, '%Y%m%d%H%M%S')
        ret.start_time = start_time
        ret.end_time = end_time
       
        # get longitudes and latitudes

        longitudes = f.variables['longitude'][:]
        
        # if no viewing angles are specified, use all
        if not viewing_angles:
            viewing_angles = np.arange(0, longitudes.shape[1])
        # else, use those that are specified
        else:
            longitudes = longitudes[:, viewing_angles]
        ret.viewing_angles = viewing_angles
        numangles = len(viewing_angles)
        latitudes = f.variables['latitude'][:, viewing_angles]
        
        # some variables are different in HIRS FCDRs and MW FCDRs:
        if ret.instrument == 'HIRS':
            u_types = ['independent', 'structured']
            acquisition_time = f.variables['time'][:].data
            ret.quality_mask = np.tile(f.variables['quality_scanline_bitmask'], (numangles, 1)).T
        else:
            u_types = ['independent', 'structured', 'common']
            acquisition_time = f.variables['Time'][:]
            ret.quality_mask = f.variables['quality_pixel_bitmask'][:, viewing_angles]
        
        # longitudes and latitudes can be either masked arrays or not
        scale_factor = {}
        scale_factor['HIRS'] = 0.001 # for some reason latitudes and longitudes have to be scaled in HIRS FCDRs
        scale_factor['AMSUB'] = 1
        scale_factor['MHS'] = 1
        if isinstance(longitudes, np.ma.MaskedArray):
            ret.longitudes = longitudes.data * scale_factor[ret.instrument]
        else:
            ret.longitudes = longitudes * scale_factor[ret.instrument]
        if isinstance(latitudes, np.ma.MaskedArray):
            ret.latitudes = latitudes.data * scale_factor[ret.instrument]
        else:
            ret.latitudes = latitudes * scale_factor[ret.instrument]

        # make an array containing the scan line number for every pixel
        scanline = np.arange(0, ret.latitudes.shape[0])
        ret.scanline = np.tile(scanline, (numangles, 1)).T
        
        # make an array containing acquisition time for every pixel
        midnight = int(datetime.datetime(start_time.year, start_time.month, start_time.day, 0, 0).strftime('%s'))            
        second_of_day = acquisition_time - midnight
        second_of_day = np.tile(second_of_day, (numangles, 1)).T
        if start_time.day < end_time.day:
            second_of_day[np.where(second_of_day >= 86400)] = second_of_day[np.where(second_of_day >= 86400)] - 86400
        ret.second_of_day = second_of_day
        ret.acquisition_time = np.tile(acquisition_time, (numangles, 1)).T

        # get brightness temperature, uncertainties and quality issue mask 
        brightness_temp = {}
        uncertainty = {t: {} for t in u_types}
        quality_issue = {}
        
        channels = np.arange(1, f.dimensions['channel'].size + 1)
        ret.channels = channels
        if ret.instrument == 'MHS':
            channel_names = ['1', '2', '3', '4', '5']
        elif ret.instrument == 'AMSUB':
            channel_names = ['16', '17', '18', '19', '20']
        
        if ret.instrument == 'MHS' or ret.instrument == 'AMSUB':
            for channel, channel_name in zip(channels, channel_names):
                brightness_temp[channel] = f.variables['Ch{}_BT'.format(channel_name)][:, viewing_angles]
                uncertainty['common'][channel] = f.variables['u_common_Ch{}_BT'.format(channel_name)][:, viewing_angles]
                uncertainty['structured'][channel] = f.variables['u_structured_Ch{}_BT'.format(channel_name)][:, viewing_angles]
                uncertainty['independent'][channel] = f.variables['u_independent_Ch{}_BT'.format(channel_name)][:, viewing_angles].data
                quality_issue[channel] = f.variables['quality_issue_pixel_Ch{}_bitmask'.format(channel_name)][:, viewing_angles]
        else:
            for channel in channels:
                brightness_temp[channel] = f.variables['bt'][channel - 1, :, viewing_angles].data
                for t in u_types:
                    uncertainty[t][channel] = f.variables['u_{}'.format(t)][:, viewing_angles, channel - 1].data

        ret.brightness_temp = brightness_temp
        ret.u_Tb = {}
        for t in u_types:
            ret.u_Tb[t] = uncertainty[t]

        if instrument.upper() == 'MHS' or instrument.upper() == 'AMSUB':
            ret.quality_issue = quality_issue
        
        return ret
    
    def calc_uth(self, slope_params, intercept_params):
        """ Calculate Upper Tropospheric Humidity (UTH) from brightness 
        temperatures.
        
        Parameters:
            slope_params (np.array): slope regression parameters for all 
                                     viewing angles (from innermost to outermost)
            intercept_params (np.array): intersept regression parameters for
                                         all viewing angles
            uth_channel (numeric): Tb channel used for UTH calculation (default: 3)
        """
        if self.instrument == 'HIRS':
            u_types = ['independent', 'structured']
            uth_channel = 12
        else:
            u_types = ['independent', 'common', 'structured']
            uth_channel = 3

        self.uth_channel = uth_channel
        viewing_angles = self.viewing_angles
        brightness_temp = self.brightness_temp[uth_channel].data
        u_Tb = {}
        for t in u_types:
            u_Tb[t] = self.u_Tb[t][uth_channel].data
        numel = brightness_temp.shape[0]
        slope_params = np.hstack((slope_params, slope_params[::-1]))[viewing_angles]
        intercept_params = np.hstack((intercept_params, intercept_params[::-1]))[viewing_angles]
        
        
        # calculate UTH (in %)
        uth = np.exp(brightness_temp * slope_params + intercept_params) * 100
        
        # calculate accompanying uncertainties:
        u_uth = {}
        for t in u_types:
            u_uth[t] = np.sqrt((np.tile(slope_params, (numel, 1)) * uth * u_Tb[t]) ** 2)

        self.uth = uth
        self.u_uth = u_uth
        
    def generate_cloud_mask(self):
        #TODO: thresholds are only valid for AMSU-B so far
        """ Find pixels that are contaminated with clouds.
        Reference for microwave sensors: 
            Buehler et al. (2007): A cloud filtering method for microwave 
            upper tropospheric humidity measurements
        """
        viewing_angles = self.viewing_angles
        instrument = self.instrument
        
        if instrument == 'AMSUB' or instrument == 'MHS':
            # brightness temperature thresholds for AMSU-B channel 18 (here: channel 3)
            Tb18_thresholds = [233.3, 233.9, 234.4, 234.9, 235.2, 235.5, 235.8, 236.1, 236.4,
                               236.6, 236.7, 237.0, 237.2, 237.4, 237.6, 237.8, 238.0, 238.2,
                               238.3, 238.5, 238.6, 238.7, 238.8, 239.0, 239.1, 239.2, 239.2,
                               239.3, 239.4, 239.5, 239.6, 239.6, 239.8, 239.8, 239.9, 239.9,
                               239.9, 239.9, 240.1, 240.1, 240.1, 240.1, 240.1, 240.1, 240.1]
            
            # clouds are there, if either the Tb of channel 18 is smaller than a threshold, or...
            Tb18_mask = self.brightness_temp[3] < np.hstack((Tb18_thresholds, Tb18_thresholds[::-1]))[viewing_angles]
            # ... the Tb in channel 19 is smaller than the Tb in channel 18
            deltaTb19_mask = (self.brightness_temp[4] - self.brightness_temp[3]) < 0.
            # combine both criteria to one cloud mask
            cloud_mask = np.array(np.logical_or(Tb18_mask, deltaTb19_mask))
            
        elif instrument == 'HIRS':
            Tb12_threshold = 240
            cloud_mask = self.brightness_temp[12] < Tb12_threshold
        
        self.cloud_mask = cloud_mask
    
    def generate_quality_and_issue_mask(self):
        """ Combines quality mask and quality issue mask to one mask, only for
        UTH channel.
        """
        instrument = self.instrument
        if instrument == 'HIRS':
            uth_channel = 12
        else:
            uth_channel = 3
        quality_and_issue_mask = {}
        
        if instrument == 'AMSUB' or instrument == 'MHS':
            # total mask is combination of cloud mask, general quality mask
            # and quality issue mask of the UTH channel
            quality_and_issue_mask = np.logical_or(
                    np.logical_or(self.quality_mask & 1, self.longitudes < -180.), 
                    self.quality_issue[uth_channel] >= 4)
        # HIRS FCDRs do not have a quality issue bitmask per channel
        elif instrument == 'HIRS':
            quality_and_issue_mask = np.logical_or(self.quality_mask & 1, self.longitudes < -180.)
                
        self.quality_and_issue_mask = quality_and_issue_mask
  
    def generate_total_mask(self):
        """ Combines quality mask and cloud mask to one total mask for the UTH
        channel.
        """
        total_mask = {}
        
        if ~hasattr(self, 'cloud_mask'):
            self.generate_cloud_mask()
        if ~hasattr(self, 'quality_and_issue_mask'):
            self.generate_quality_and_issue_mask()
            # total mask is combination of cloud mask, general quality mask
            # and quality issue mask of the UTH channel
        total_mask = np.logical_or(self.cloud_mask, self.quality_and_issue_mask)
        
        self.total_mask = total_mask
        
    def generate_node_mask(self):
        """ Generates a mask indicating where the satellite was in ascending/
        descending node. True indicates ascending, False indicates descending.
        """
        is_ascending = np.diff(self.latitudes, axis=0) >= 0 
        is_ascending = np.append(is_ascending, [is_ascending[-1,:]], axis=0)
        
        self.node_mask = is_ascending