#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 10:04:57 2018

@author: Theresa Lang
"""
import numpy as np
import pandas as pd
from os.path import join, exists
from os import makedirs
from scipy.sparse import csr_matrix, csc_matrix, diags, bmat, block_diag
import matplotlib.pyplot as plt 
import time
from datetime import datetime
from calendar import monthrange
from fiduceo.cdr.writer.cdr_writer import CDRWriter
from netCDF4 import Dataset
import processing.utils as utils

class CDR:
    def __init__(self, source=None, lat=None, lon=None, BT=None,
                 brightness_temp_std=None, uth_std=None, 
                 uth=None, u_Tb=None, u_uth=None,
                 instrument=None, satellite=None,
                 time_coverage_start=None, time_coverage_end=None,
                 time_ranges=None,
                 observation_count=None, observation_count_all=None,
                 overpass_count=None):
        self.source = source
        self.instrument = instrument
        self.satellite = satellite
        self.time_coverage_start = time_coverage_start
        self.time_coverage_end = time_coverage_end
        self.time_ranges = time_ranges
        self.lat = lat
        self.lon = lon
        self.BT = BT
        self.uth = uth
        self.u_Tb = u_Tb
        self.u_uth = u_uth
        self.observation_count = observation_count
        self.observation_count_all = observation_count_all
        self.overpass_count = overpass_count

    @classmethod
    def GriddedCDRFromFCDRs(cls, FCDRs,
                            lat_boundaries=[-90, 90], lon_boundaries =[-179, 180], 
                            resolution=1.):
        """ Creates one CDR from from a list of FCDRs.
        The CDR contains mean brightness temperatures and UTH binned to a 
        lat/lon grid as well as 3 uncertainty classes for every grid cell. 
        Uncertainty correlations between grid cells are NOT propagated!
        
        Parameters:
            FCDRs (list): list fo FCDR objects
            lat_boundaries (list): lower and upper boundary for latitude-grid
            lon_boundaries (list): lower and upper boundary for longitude-grid
            resolution (float): resolution of lat/lon grid [degree] (default: 1.0)
        """
        #t1 = time.clock()
        
        instrument = FCDRs[0].instrument
        u_types = ['independent', 'common', 'structured']
        if instrument == 'HIRS':
            uth_channel = 12
            corr_vector = FCDRs[0].corr_vector[uth_channel][FCDRs[0].corr_vector[uth_channel] > 0]
        elif instrument in ['MHS', 'AMSUB', 'SSMT2']:
            uth_channel = 3            
            #corr_vector = [1.0, 0.8465, 0.5134, 0.2231, 0.0695, 0.0155, 0.0025]
            corr_vector = FCDRs[0].corr_vector[uth_channel]
    
        branches = ['ascending', 'descending']
        
        # latitude and longitude bins and their centers
        lat_bins, lon_bins, lat_centers, lon_centers = utils.getLatLonBins(
                lat_boundaries, lon_boundaries, resolution)
        
        # initialization:
        # dictionaries for grouped quantities
        group_variances_Tb = dict.fromkeys(u_types)
        group_variances_uth = dict.fromkeys(u_types)
        # dictionaries for gridded fields that will be returned
        Tb_gridded = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        Tb_unfiltered_gridded = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        UTH_gridded = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        u_Tb_gridded = {t: {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches} for t in u_types}
        u_Tb_unfiltered_gridded = {t: {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches} for t in u_types}
        u_uth_gridded = {t: {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches} for t in u_types}
        count = {b: np.zeros((len(lat_centers), len(lon_centers))) for b in branches}
        count_all = {b: np.zeros((len(lat_centers), len(lon_centers))) for b in branches}
        count_overpasses = {b: np.zeros((len(lat_centers), len(lon_centers))) for b in branches}
        second_of_day_min = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        second_of_day_max = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        time_ranges = dict.fromkeys(branches)
        start_time = dict.fromkeys(branches)
        end_time = dict.fromkeys(branches)
        time_coverage_start = np.nan
        time_coverage_end = np.nan
        is_empty = False
 
        # collect data from all FCDRs
        collected_data, collected_data_unfiltered, collected_files = utils.collectFCDRData(
                FCDRs, u_types, uth_channel=uth_channel)

        for b in branches:
            collected_data['longitude'][b][collected_data['longitude'][b] < -180 + 0.5 * resolution] = 180.
            collected_data_unfiltered['longitude'][b][collected_data_unfiltered['longitude'][b] < -180 + 0.5 * resolution] = 180.
            # combine all values of this branch in a pandas dataframe
            data = pd.DataFrame(utils.flattenDict(collected_data, b))
            data_unfiltered = pd.DataFrame(utils.flattenDict(collected_data_unfiltered, b))

            # throw away data outside the specified new grid
            data = data[np.logical_and(
                    (data['latitude'] <= lat_bins[-1]),
                    (data['latitude'] > lat_bins[0])
                    )].reset_index()
            data = data[np.logical_and(
                    (data['longitude'] <= lon_bins[-1]),
                    (data['longitude'] > lon_bins[0])
                    )].reset_index()
            data_unfiltered = data_unfiltered[np.logical_and(
                    (data_unfiltered['latitude'] <= lat_bins[-1]),
                    (data_unfiltered['latitude'] > lat_bins[0])
                    )].reset_index()
            data_unfiltered = data_unfiltered[np.logical_and(
                    (data_unfiltered['longitude'] <= lon_bins[-1]),
                    (data_unfiltered['longitude'] > lon_bins[0])
                    )].reset_index()

#            fig, ax = plt.subplots()
#            ax.set_title(b)
#            ax.scatter(data.longitude, data.latitude, c=data.brightness_temp, s=1)
#            ax.set_xlim(lon_bins[0], lon_bins[-1])
#            ax.set_ylim(lat_bins[0], lat_bins[-1])

            if len(data.index) <= 1 and len(data_unfiltered.index) <=1:
                # no overpass in the selected grid
                print('dataframe is empty!')
                is_empty = True
            else:
                # bin data to latitude and longitude bins
                data = utils.binData(
                        data, lat_bins, lon_bins, lat_centers, lon_centers)
                data_unfiltered = utils.binData(
                        data_unfiltered, lat_bins, lon_bins, lat_centers, lon_centers)
                                
                # group data by latitude and longitude bins
                data_grouped = data.groupby([data.lat_bin, data.lon_bin], sort=False)
                data_unfiltered_grouped = data_unfiltered.groupby([data_unfiltered.lat_bin, data_unfiltered.lon_bin], sort=False)
                
                # get time of first and last data point going into this CDR
                start_time[b] = datetime.fromtimestamp(np.min(data.acquisition_time))
                end_time[b] = datetime.fromtimestamp(np.max(data.acquisition_time))
                
                # go through all groups and propagate uncertainties
                for name, group in data_grouped:
                    #group_ind = groups.get(group).values
                    lat_ind = name[0]
                    lon_ind = name[1]
                    group_size = len(group)
                    
                    # put averages of current group in the right cell on the grid 
                    Tb_gridded[b][lat_ind, lon_ind] = group.brightness_temp.mean()
                    UTH_gridded[b][lat_ind, lon_ind] = group.uth.mean()

                    # put count 
                    count[b][lat_ind, lon_ind] = group_size
                    count_overpasses[b][lat_ind, lon_ind] += 1
                    second_of_day_group = np.array(group.second_of_day)
                    second_of_day_min[b][lat_ind, lon_ind] = np.min(second_of_day_group)
                    second_of_day_max[b][lat_ind, lon_ind] = np.max(second_of_day_group)

                    # Get structured, independent and common uncertainties of this group
                    for t in u_types:
                        group_variances_Tb[t] = np.array(group['u_Tb_{}'.format(t)])
                        group_variances_uth[t] = np.array(group['u_uth_{}'.format(t)])
                    
                    # Create a covariance matrix for structured uncertainties:
                    # scanlines of data points in this group
                    scanlines = np.array(group.scanline, dtype=np.int)
                    # construct correlation matrix from scanline differences
                    corr = utils.getCorrMatrix(scanlines, corr_vector)
                    # calculate covariance matrix from correlation matrix
                    S_struct_Tb = np.multiply(
                            corr, 
                            np.outer(
                                    group_variances_Tb['structured'], 
                                    group_variances_Tb['structured']
                                    )
                            )
                    S_struct_uth = np.multiply(
                            corr, 
                            np.outer(
                                    group_variances_uth['structured'],
                                    group_variances_uth['structured']
                                    )
                            )
                    
                    # Calculate uncertainties (standard deviations) for the current group: 
                    # independent uncertainty for this grid cell 
                    # (Use Law of the Propagation of Uncertainties for independent uncertainties only)
                    u_Tb_gridded['independent'][b][lat_ind, lon_ind] = np.sqrt(
                            np.sum(group_variances_Tb['independent'] ** 2)
                            ) / group_size
                    u_uth_gridded['independent'][b][lat_ind, lon_ind] = np.sqrt(
                            np.sum(group_variances_uth['independent'] ** 2)
                            ) / group_size
                    # common uncertainty for this grid cell 
                    # (fully correlated uncertainties --> uncertainty of average is average of uncertainties)
                    # Only MW FCDRs contain common uncertainties
                    u_Tb_gridded['common'][b][lat_ind, lon_ind] = np.mean(group_variances_Tb['common'])
                    u_uth_gridded['common'][b][lat_ind, lon_ind] = np.mean(group_variances_uth['common'])
                    # structured uncertainty for this grid cell
                    u_Tb_gridded['structured'][b][lat_ind, lon_ind] = np.sqrt(np.sum(S_struct_Tb)) / group_size
                    u_uth_gridded['structured'][b][lat_ind, lon_ind] = np.sqrt(np.sum(S_struct_uth)) / group_size
                
                # go through data that is not cloud filtered:
                for name, group in data_unfiltered_grouped:
                    lat_ind = name[0]
                    lon_ind = name[1]
                    group_size = len(group)
                    count_all[b][lat_ind, lon_ind] = group_size
                    
                    Tb_unfiltered_gridded[b][lat_ind, lon_ind] = group.brightness_temp.mean()
                    for t in u_types:
                        group_variances_Tb[t] = np.array(group['u_Tb_{}'.format(t)])
                    
                    scanlines = np.array(group.scanline, dtype=np.int)
                    corr = utils.getCorrMatrix(scanlines, corr_vector)
                    S_struct_Tb_unfiltered = np.multiply(
                            corr, 
                            np.outer(
                                    group_variances_Tb['structured'], 
                                    group_variances_Tb['structured'])
                            )
                    u_Tb_unfiltered_gridded['independent'][b][lat_ind, lon_ind] = np.sqrt(np.sum(group_variances_Tb['independent'] ** 2)) / group_size
                    u_Tb_unfiltered_gridded['common'][b][lat_ind, lon_ind] = np.mean(group_variances_Tb['common'])
                    u_Tb_unfiltered_gridded['structured'][b][lat_ind, lon_ind] = np.sqrt(np.sum(S_struct_Tb_unfiltered)) / group_size
                    
            # construct time_ranges matrix
            time_ranges[b] = np.stack((second_of_day_min[b], second_of_day_max[b]))
            
        if not is_empty:
            time_coverage_start  = min(start_time.values())
            time_coverage_end = max(end_time.values())
        
        # return variables
        ret = cls()
        # lat, lon
        ret.lat = np.array(lat_centers)
        ret.lon = np.array(lon_centers)
        ret.lat_bins = lat_bins
        ret.lon_bins = lon_bins
        ret.geospatial_lat_resolution = resolution
        ret.geospatial_lon_resolution = resolution
        # BT, UTH and observation counts
        ret.BT = Tb_gridded
        ret.BT_full = Tb_unfiltered_gridded
        ret.uth = UTH_gridded
        ret.observation_count = count
        ret.observation_count_all = count_all
        ret.overpass_count = count_overpasses
        # time information
        ret.time_coverage_start = time_coverage_start
        ret.time_coverage_end = time_coverage_end
        ret.time_ranges = time_ranges
        # FCDR files included
        ret.source = collected_files
        ret.is_empty = is_empty
        # instrument and satellite
        ret.instrument = FCDRs[0].instrument
        ret.satellite = FCDRs[0].satellite
        # uncertainties
        ret.u_Tb = u_Tb_gridded
        ret.u_Tb_full = u_Tb_unfiltered_gridded
        ret.u_uth = u_uth_gridded
        
        return ret

    @classmethod
    def AveragedCDRFromCDRs(cls, CDRs):
        """ Combines several CDRs to one averaged CDR by calculating an ordinary 
        average of brightness temperature and UTH for every grid cell. 
        Uncertainties are propagated using the Law of the Propagation of 
        Uncertainty for the case of an ordenary average.
        
        Parameters:
            CDRs (list of CDRs): List of CDR objects (e.g. 30 daily CDRs that
                 should be combined to one monthly CDR)
            uncertainties (boolean): False, if uncertainties should not be 
                propagated
        """
        ret = cls()
        
        u_types = ['independent', 'structured', 'common']

        branches = ['descending', 'ascending']
        num_timesteps = len(CDRs)
        num_lons = len(CDRs[0].lon)
        num_lats = len(CDRs[0].lat)
        Tb_mean = dict.fromkeys(branches)
        Tb_full_mean = dict.fromkeys(branches)
        UTH_mean = dict.fromkeys(branches)
        u_Tb = {t: dict.fromkeys(branches) for t in u_types}
        u_Tb_full = {t: dict.fromkeys(branches) for t in u_types}
        u_uth = {t: dict.fromkeys(branches) for t in u_types}
        Tb_std = {b: np.ones((num_lats, num_lons)) * np.nan for b in branches}
        UTH_std = {b: np.ones((num_lats, num_lons)) * np.nan for b in branches}
        count = dict.fromkeys(branches)
        count_all = dict.fromkeys(branches)
        count_overpasses = dict.fromkeys(branches)
        quality_bitmask = dict.fromkeys(branches)
        second_of_day_min = dict.fromkeys(branches)
        second_of_day_max = dict.fromkeys(branches)
        time_ranges = dict.fromkeys(branches)
        files = []
        is_empty = False
        
        first_not_empty = utils.getFirstNotEmptyCDR(CDRs)
        last_not_empty = utils.getLastNotEmptyCDR(CDRs)
        
        if first_not_empty is not None:
            time_coverage_start = first_not_empty.time_coverage_start
            time_coverage_end = last_not_empty.time_coverage_end
            time_coverage_duration = pd.Timedelta(time_coverage_end - time_coverage_start).isoformat()
    
            for b in branches:
                for i in range(num_timesteps):
                    files.extend(CDRs[i].source)
                notnan_count = np.sum([~np.isnan(CDRs[i].BT[b]) for i in range(num_timesteps)], axis=0)
                Tb_mean[b] = np.around(np.nanmean([CDRs[i].BT[b] for i in range(num_timesteps)], axis=0), 2)
                Tb_full_mean[b] = np.around(np.nanmean([CDRs[i].BT_full[b] for i in range(num_timesteps)], axis=0), 2)
                UTH_mean[b] = np.around(np.nanmean([CDRs[i].uth[b] for i in range(num_timesteps)], axis=0) , 2)

                Tb_std[b] = np.nanstd([CDRs[i].BT[b] for i in range(num_timesteps)], axis=0, ddof=1)
                UTH_std[b] = np.nanstd([CDRs[i].uth[b] for i in range(num_timesteps)], axis=0, ddof=1)
                # add counts
                count[b] = np.nansum([CDRs[i].observation_count[b] for i in range(num_timesteps)], axis=0)
                count[b][count[b] == 0] = -32767
                count_all[b] = np.nansum([CDRs[i].observation_count_all[b] for i in range(num_timesteps)], axis=0)
                count_all[b][count_all[b] == 0] = -32767
                count_overpasses[b] = np.nansum([CDRs[i].overpass_count[b] for i in range(num_timesteps)], axis=0)
                # setting new quality flags based on counts - not possible yet
                # quality_bitmask[b] = utils.getCDRQualityMask(count[b], count_overpasses[b])                
                
                # combine time ranges
                second_of_day_min[b] = np.nanmin([CDRs[i].time_ranges[b][0] for i in range(num_timesteps)], axis=0)
                second_of_day_min[b][np.isnan(second_of_day_min[b])] = 4294967295
                second_of_day_max[b] = np.nanmax([CDRs[i].time_ranges[b][1] for i in range(num_timesteps)], axis=0)
                second_of_day_max[b][np.isnan(second_of_day_max[b])] = 4294967295
                time_ranges[b] = np.stack((second_of_day_min[b], second_of_day_max[b]))
                
                # combine uncertainties
                notnan_count = np.sum([~np.isnan(CDRs[i].u_Tb['independent'][b]) for i in range(num_timesteps)], axis=0)
                notnan_count_full = np.sum([~np.isnan(CDRs[i].u_Tb_full['independent'][b]) for i in range(num_timesteps)], axis=0)
                for t in ['independent', 'structured']:
                    u_Tb[t][b] = np.sqrt(np.nansum([CDRs[i].u_Tb[t][b] ** 2 for i in range(num_timesteps)], axis=0)) / notnan_count
                    u_Tb_full[t][b] = np.sqrt(np.nansum([CDRs[i].u_Tb_full[t][b] ** 2 for i in range(num_timesteps)], axis=0)) / notnan_count_full
                    u_uth[t][b] = np.sqrt(np.nansum([CDRs[i].u_uth[t][b] ** 2 for i in range(num_timesteps)], axis=0)) / notnan_count

                u_Tb['common'][b] = np.nanmean([CDRs[i].u_Tb['common'][b] for i in range(num_timesteps)], axis=0)
                u_Tb_full['common'][b] = np.nanmean([CDRs[i].u_Tb_full['common'][b] for i in range(num_timesteps)], axis=0)
                u_uth['common'][b] = np.nanmean([CDRs[i].u_uth['common'][b] for i in range(num_timesteps)], axis=0)
                # round uncertainties to 4 digits:
                numdigits = 4
                for t in u_types:
                    u_Tb[t][b] = np.around(u_Tb[t][b], numdigits)
                    u_Tb_full[t][b] = np.around(u_Tb_full[t][b], numdigits)
                    u_uth[t][b] = np.around(u_uth[t][b], numdigits)
            
            # fill quality bitmask with fill values
            quality_bitmask = np.zeros((num_lons, num_lats)) 
            # return data
            ret.lon = CDRs[0].lon
            ret.lat = CDRs[0].lat
            ret.geospatial_lat_resolution = CDRs[0].geospatial_lat_resolution
            ret.geospatial_lon_resolution = CDRs[0].geospatial_lon_resolution
            ret.lat_bnds = np.stack((CDRs[0].lat_bins[0:-1], CDRs[0].lat_bins[1:])).T
            ret.lon_bnds = np.stack((CDRs[0].lon_bins[0:-1], CDRs[0].lon_bins[1:])).T
            ret.BT = Tb_mean
            ret.BT_full = Tb_full_mean
            ret.uth = UTH_mean
            ret.BT_inhomogeneity = Tb_std
            ret.uth_inhomogeneity = UTH_std
            ret.observation_count = count
            ret.observation_count_all = count_all
            ret.overpass_count = count_overpasses
            ret.quality_pixel_bitmask = quality_bitmask
            ret.time_coverage_start = time_coverage_start 
            ret.time_coverage_end = time_coverage_end
            ret.time_coverage_duration = time_coverage_duration
            ret.time_ranges = time_ranges
            ret.instrument = CDRs[0].instrument
            ret.satellite = CDRs[0].satellite
            ret.u_BT = u_Tb
            ret.u_BT_full = u_Tb_full
            ret.u_uth = u_uth
            ret.source = files
            ret.is_empty = is_empty
    
        else:
            print('No data in this month.')
            ret.is_empty = True
        
        return ret
    
    def toNetCDF(self, CDR_path, version, comment_on_version, overwrite=False):
        """ Saves the attributes of an averaged CDR created by 
        AveragedCDRFromCDRs to a NetCDF file using the CDRWriter by Tom Block.
        
        Parameters:
            CDR_path (str): path to directory where CDRs are saved
            version (str): Version of this CDR (for CDR filename)
            comment_on_version (str): Comment on this version of the CDR.
        """
           
        u_types = ['structured', 'common', 'independent']    
        branches = ['ascend', 'descend']
        u_quantities = ['uth', 'BT', 'BT_full']
        DATE_PATTERN = "%Y%m%d%H%M%S"
        numlons = len(self.lon)
        numlats = len(self.lat)
        start_time = self.time_coverage_start
        end_time = self.time_coverage_end
        start_time_for_filename = datetime(
                start_time.year, 
                start_time.month, 
                1, 0, 0, 0
                )
        end_time_for_filename = datetime(
                start_time.year, 
                start_time.month, 
                monthrange(start_time.year, start_time.month)[1], 
                23, 59, 59
                )
        CDR_filename = CDRWriter.create_file_name_CDR(
                'UTH', 
                sensor=self.instrument, 
                platform=self.satellite, 
                start=start_time_for_filename, 
                end=end_time_for_filename, 
                type='L3', 
                version=version
                )
        # create directory for CDRs if it does not exist
        if not exists(CDR_path):
            makedirs(CDR_path)
        
        ds = CDRWriter.createTemplate('UTH', numlons, numlats)
        simple_vars = ['lat', 'lon', 'lat_bnds', 'lon_bnds']
        branch_vars = ['observation_count', 'observation_count_all',
                       'overpass_count',
                       'uth', 'BT', 'BT_full',
                       'uth_inhomogeneity', 'BT_inhomogeneity',
                       'time_ranges']
        simple_attrs = ['source', 'time_coverage_duration',
                        'geospatial_lon_resolution', 'geospatial_lat_resolution']
        for var in simple_vars:
            ds.variables[var].data = getattr(self, var)
        
        for attr in simple_attrs:
            ds.attrs[attr] = getattr(self, attr)
            
        for var in branch_vars:
            for b in branches:
                ds.variables['{}_{}'.format(var, b)].data =\
                getattr(self, '{}'.format(var))['{}ing'.format(b)]
        
        for q in u_quantities:
            for t in u_types:
                for b in branches:
                    ds.variables['u_{}_{}_{}'.format(t, q, b)].data =\
                    getattr(self, 'u_{}'.format(q))[t]['{}ing'.format(b)]

        ds.attrs['institution'] = 'University of Hamburg'
        ds.attrs['title'] = 'Upper Tropospheric Humidity (UTH)'
        ds.attrs['time_coverage_start'] = start_time.strftime(DATE_PATTERN)
        ds.attrs['time_coverage_end'] = end_time.strftime(DATE_PATTERN)
        ds.attrs['time_coverage_resolution'] = pd.Timedelta(1, unit='M').isoformat()
        ds.attrs['geospatial_lat_units'] = 'deg'
        ds.attrs['geospatial_lon_units'] = 'deg'
        ds.attrs['history'] = 'Created on {}.'.format(
                datetime.now().strftime('%Y-%m-%d %H:%M')
                )
        ds.attrs['references'] = ''
        ds.attrs['comment'] = comment_on_version
        ds.attrs['auxiliary_data'] = ''
        ds.attrs['configuration'] = ''
        CDRWriter.write(ds, join(CDR_path, CDR_filename), overwrite=overwrite)
        
        print('CDR {} for {} on {} saved as {}'.format(
                start_time.strftime('%Y/%m'),
                self.instrument,
                self.satellite,
                CDR_filename
                ))
    
    @classmethod
    def fromNetCDF(cls, filename):
        """ Creates a CDR-Object from a NetCDF file.
        
        Parameters:
            filename: path and name of NetCDF file
        """
        ret = cls()
        f = Dataset(filename)
        DATE_PATTERN = "%Y%m%d%H%M%S"
        branches = ['ascend', 'descend']
        u_types = ['structured', 'common', 'independent']
        u_quantities = ['uth', 'BT', 'BT_full']
        simple_vars = ['lat', 'lon', 'lat_bnds', 'lon_bnds']
        branch_vars = ['observation_count', 'observation_count_all',
                       'overpass_count',
                       'uth', 'BT', 'BT_full',
                       'uth_inhomogeneity', 'BT_inhomogeneity',
                       'time_ranges']
        simple_attrs = ['source', 'time_coverage_duration',
                        'geospatial_lon_resolution', 'geospatial_lat_resolution',
                        'geospatial_lat_units', 'geospatial_lon_units']
        
        for var in simple_vars:
            setattr(ret, var, f.variables[var][:])
        
        for attr in simple_attrs:
            setattr(ret, attr, getattr(f, attr))
            
        for var in branch_vars:
            h = {}
            for b in branches:
                h['{}ing'.format(b)] = f.variables['{}_{}'.format(var, b)][:]
            setattr(ret, var, h)
        
        for q in u_quantities:
            h = {}
            for t in u_types:
                h[t] = {}
                for b in branches:
                    h[t]['{}ing'.format(b)] = f.variables['u_{}_{}_{}'.format(t, q, b)][:]
            setattr(ret, 'u_{}'.format(q), h)
        
        ret.time_coverage_start = datetime.strptime(f.time_coverage_start, DATE_PATTERN)
        ret.time_coverage_end = datetime.strptime(f.time_coverage_end, DATE_PATTERN)
        
        return ret
                
