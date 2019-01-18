#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 16:52:03 2018

@author: u300740
"""
import numpy as np
import os.path
from netCDF4 import Dataset
import pandas as pd
from scipy.sparse import csr_matrix, csc_matrix, diags, bmat, block_diag
import matplotlib.pyplot as plt
import datetime
import processing.utils as utils

class FCDR:
    """ An FCDR object contains the data from one FCDR file (one satellite
    orbit), which is needed for the creation of a UTH CDR. """
    def __init__(self, latitude=None, longitude=None, viewing_angles=None,
                 brightness_temp=None, channels=None, scanline=None, quality_mask=None, 
                 u_Tb=None, u_uth=None, quality_issue=None, start_time=None, 
                 end_time=None, acquisition_time=None,instrument=None, 
                 satellite=None, file=None, uth=None):
        self.file = file
        self.latitude = latitude
        self.longitude = longitude
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
    def fromNetcdf(cls, fcdr_path, filename, viewing_angles):
        """ Create FCDR from a FIDUCEO FCDR in NetCDF format.
        
        Parameters:
            fcdr_path (str): path to directory with NetCDF file
            filename (str): name of NetCDF file 
            viewing_angles (list): list of viewing angles/ pixels that should 
                be used (e.g. [44, 45] for innermost angles only)
        """
        ret = cls()
        # read data from NetCDF file
        f = Dataset(os.path.join(fcdr_path, filename))
        
        # get information from file name
        filename = os.path.basename(filename)
        file_info = utils.getFileInfo(filename)
        # number of pixels/ viewing angles
        numangles = len(viewing_angles)
        # get longitudes and latitudes
        longitude = f.variables['longitude'][:][:, viewing_angles]
        latitude = f.variables['latitude'][:, viewing_angles]
        
        # some variables are different in HIRS FCDRs and MW FCDRs:
        u_types = ['independent', 'structured', 'common']
        if file_info['instrument'] == 'HIRS':
            time_attr = 'time'
        else:
            time_attr = 'Time'
        if isinstance(f.variables[time_attr][:], np.ma.MaskedArray):
            acquisition_time = f.variables[time_attr][:].data
        else:
            acquisition_time = f.variables[time_attr][:]
            
        ret.quality_mask = f.variables['quality_pixel_bitmask'][:, viewing_angles]
        # for some reason latitudes and longitudes have to be scaled in HIRS FCDRs
        #TODO: check whether this will change with newer HIRS FCDR version!
        scale_factor = {}
        scale_factor['HIRS'] = 0.001 
        scale_factor['AMSUB'] = 1
        scale_factor['MHS'] = 1
        scale_factor['SSMT2'] = 1
        
        # longitudes and latitudes can be either masked arrays or not
        if isinstance(longitude, np.ma.MaskedArray):
            ret.longitude = longitude.data * scale_factor[file_info['instrument']]
        else:
            ret.longitude = longitude * scale_factor[file_info['instrument']]
        if isinstance(latitude, np.ma.MaskedArray):
            ret.latitude = latitude.data * scale_factor[file_info['instrument']]
        else:
            ret.latitude = latitude * scale_factor[file_info['instrument']]

        # make an array containing the scan line number for every pixel
        scanline = np.arange(0, ret.latitude.shape[0])
        scanline_array = np.tile(scanline, (numangles, 1)).T
        
        # make an array containing acquisition time for every pixel
        second_of_day = utils.getSecondOfDay(acquisition_time)
        second_of_day = np.tile(second_of_day, (numangles, 1)).T
        acquisition_time = np.tile(acquisition_time, (numangles, 1)).T

        # get brightness temperature, uncertainties and quality issue mask 
        brightness_temp = {}
        uncertainty = {t: {} for t in u_types}
        quality_issue = {}
        corr_vector = {}
        
        channels = np.arange(1, f.dimensions['channel'].size + 1)
        
        # channel names are different for MHS and AMSUB
        if file_info['instrument'] == 'MHS':
            channel_names = ['1', '2', '3', '4', '5']
        elif file_info['instrument'] == 'AMSUB':
            channel_names = ['16', '17', '18', '19', '20']
        elif file_info['instrument'] == 'SSMT2':
            channel_names = ['4', '5', '2', '1', '3']
        
        # extract brightness temperatures, uncertainties and correlation vectors
        if file_info['instrument'] in ['MHS', 'AMSUB', 'SSMT2']:
            for channel, channel_name in zip(channels, channel_names):
                brightness_temp[channel] = f.variables['Ch{}_BT'.format(channel_name)][:, viewing_angles].filled(np.nan)
                uncertainty['common'][channel] = f.variables['u_common_Ch{}_BT'.format(channel_name)][:, viewing_angles].filled(np.nan)
                uncertainty['structured'][channel] = f.variables['u_structured_Ch{}_BT'.format(channel_name)][:, viewing_angles].filled(np.nan)
                uncertainty['independent'][channel] = f.variables['u_independent_Ch{}_BT'.format(channel_name)][:, viewing_angles].filled(np.nan)
                quality_issue[channel] = f.variables['quality_issue_pixel_Ch{}_bitmask'.format(channel_name)][:, viewing_angles]
                #corr_vector[channel] = f.variables['cross_line_correlation_coefficients'][:, channel-1]
                corr_vector[channel] = []
        else:
            for channel in channels:
                brightness_temp[channel] = f.variables['bt'][channel - 1, :, viewing_angles]
                quality_issue[channel] = np.tile(f.variables['quality_channel_bitmask'][:, channel-1].data, (numangles, 1)).T
                corr_vector[channel] = f.variables['cross_line_correlation_coefficients'][:, channel-1]
                for t in u_types:
                    uncertainty[t][channel] = f.variables['u_{}'.format(t)][:, viewing_angles, channel - 1]
        
        # return variables
        ret.brightness_temp = brightness_temp
        
        ret.u_Tb = {}
        for t in u_types:
            ret.u_Tb[t] = uncertainty[t]
            
        ret.quality_issue = quality_issue
        
        ret.file = filename 
        ret.instrument = ''.join([i for i in file_info['instrument']])
        ret.satellite = file_info['satellite']
        ret.start_time = file_info['start_time']
        ret.end_time = file_info['end_time']
        ret.viewing_angles = viewing_angles
        ret.scanline = scanline_array
        ret.channels = channels
        ret.second_of_day = second_of_day
        ret.acquisition_time = acquisition_time
        ret.corr_vector = corr_vector
        
        return ret
    
    def calcUTH(self, slope_params, intercept_params):
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
            uth_channel = 12
        else:            
            uth_channel = 3
        u_types = ['independent', 'common', 'structured']
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
        
    def generateCloudMask(self):
        """ Find pixels that are contaminated with clouds.
        Reference for microwave sensors: 
            Buehler et al. (2007): A cloud filtering method for microwave 
            upper tropospheric humidity measurements
            
        For HIRS a simple threshold is used for the ch12 brightness temperature
        so far.
            
        """
        viewing_angles = self.viewing_angles
        instrument = self.instrument
        
        # brightness temperature thresholds for AMSU-B channel 18 (here: channel 3)
        Tb18_thresholds = [233.3, 233.9, 234.4, 234.9, 235.2, 235.5, 235.8, 236.1, 236.4,
                           236.6, 236.7, 237.0, 237.2, 237.4, 237.6, 237.8, 238.0, 238.2,
                           238.3, 238.5, 238.6, 238.7, 238.8, 239.0, 239.1, 239.2, 239.2,
                           239.3, 239.4, 239.5, 239.6, 239.6, 239.8, 239.8, 239.9, 239.9,
                           239.9, 239.9, 240.1, 240.1, 240.1, 240.1, 240.1, 240.1, 240.1]
        
        # SSMT2 has fewer viewing angles. Use values from nearest angles of AMSU-B
        if instrument == 'SSMT2':
            nearest_amsu_views = [8, 10, 13, 16, 19, 21, 24, 27, 29, 32, 35, 38, 40, 43]
            Tb18_thresholds = [Tb18_thresholds[i] for i in nearest_amsu_views]
        
        if instrument in ['MHS', 'AMSUB', 'SSMT2']:
            # clouds are there, if either the Tb of channel 18 is smaller than a threshold, or...
            Tb18_mask = self.brightness_temp[3] < np.hstack((Tb18_thresholds, Tb18_thresholds[::-1]))[viewing_angles]
            # ... the Tb in channel 19 is smaller than the Tb in channel 18
            deltaTb19_mask = (self.brightness_temp[4] - self.brightness_temp[3]) < 0.
            # combine both criteria to one cloud mask
            cloud_mask = np.array(np.logical_or(Tb18_mask, deltaTb19_mask))
            
        elif instrument == 'HIRS':
            #Tb8_threshold = 235
            #delta_Tb8_Tb12_threshold = 25
            #cloud_mask = np.logical_or(
            #self.brightness_temp[8] <= Tb8_threshold,
            #(self.brightness_temp[8] - self.brightness_temp[12]) <= delta_Tb8_Tb12_threshold)
            Tb12_threshold = 240
            cloud_mask = self.brightness_temp[12] < Tb12_threshold
        
        self.cloud_mask = cloud_mask
    
    def generateQualityAndIssueMask(self):
        """ Combines quality mask and quality issue mask to one mask, only for
        UTH channel.
        """
        instrument = self.instrument
        if instrument == 'HIRS':
            uth_channel = 12
        else:
            uth_channel = 3
        quality_and_issue_mask = {}
        # total mask is combination of cloud mask, general quality mask,
        # quality issue mask of the UTH channel and all other pixels that
        # are masked in the FCDR file (nanmask)
        nanmask = np.isnan(self.brightness_temp[uth_channel])
        if instrument in ['AMSUB', 'MHS', 'SSMT2']:
            quality_and_issue_mask = np.logical_or(
                    np.logical_or(self.quality_mask & 1, self.longitude < -180.), 
                    np.logical_or(self.quality_issue[uth_channel] >= 4, nanmask))
        elif instrument == 'HIRS':
            quality_mask = np.logical_or(self.quality_mask & 1, self.longitude < -180.)
            # do not use pixel that have the following flags set in the channel
            # specific mask:
            # do_not_use (bit 1), uncertainty_suspicious (bit 2),
            # self_emission_fails (bit 3), calibration_impossible (bit 4)
            issue_mask = np.logical_or(
                    self.quality_issue[uth_channel] & 1, 
                    np.logical_and(self.quality_issue[uth_channel] >= 2, self.quality_issue[uth_channel] < 16))
            quality_and_issue_mask = np.logical_or(quality_mask, issue_mask)

                
        self.quality_and_issue_mask = quality_and_issue_mask
  
    def generateTotalMask(self):
        """ Combines quality mask and cloud mask to one total mask for the UTH
        channel.
        """
        total_mask = {}
        
        if ~hasattr(self, 'cloud_mask'):
            self.generateCloudMask()
        if ~hasattr(self, 'quality_and_issue_mask'):
            self.generateQualityAndIssueMask()
            # total mask is combination of cloud mask, general quality mask
            # and quality issue mask of the UTH channel
        total_mask = np.logical_or(self.cloud_mask, self.quality_and_issue_mask)
        
        self.total_mask = total_mask
        
    def generateNodeMask(self):
        """ Generates a mask indicating where the satellite was in ascending/
        descending node. True indicates ascending, False indicates descending.
        """
        # latitude increases --> ascending
        is_ascending = np.diff(self.latitude, axis=0) >= 0 
        is_ascending = np.append(is_ascending, [is_ascending[-1,:]], axis=0)
        
        self.node_mask = is_ascending