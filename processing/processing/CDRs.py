#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 10:04:57 2018

@author: Theresa Lang
"""
import numpy as np
import pandas as pd
from os.path import join
from scipy.sparse import csr_matrix, csc_matrix, diags, bmat, block_diag
import matplotlib.pyplot as plt 
import time
from datetime import datetime
from fiduceo.cdr.writer.cdr_writer import CDRWriter
from netCDF4 import Dataset
import processing.utils as utils

class CDR:
    def __init__(self, source=None, lat=None, lon=None, BT=None,
                 brightness_temp_std=None, uth_std=None, 
                 uth=None, u_Tb=None, u_uth=None,
                 instrument=None, satellite=None,
                 time_coverage_start=None, time_coverage_end=None,
                 observation_count=None, observation_count_all=None):
        source = source
        instrument = instrument
        satellite = satellite
        time_coverage_start = time_coverage_start
        time_coverage_end = time_coverage_end
        lat = lat
        lon = lon
        BT = BT
        uth = uth
        u_Tb = u_Tb
        u_uth = u_uth
        observation_count = observation_count
        observation_count_all = observation_count_all

    @classmethod
    def GriddedCDRFromFCDRs(cls, FCDRs,
                            lat_boundaries=[-90, 90], lon_boundaries =[-179, 180], 
                            resolution=1.):
        """ Creates one CDR from from a list of FCDRs.
        The CDR contains mean brightness temperatures and UTH binned to a 
        lat/lon grid as well as 3 uncertainty classes for every grid cell. 
        Uncertainty correlations between grid cells are NOT propagated!!!
        
        Parameters:
            FCDRs (list): list fo FCDR objects
            lat_boundaries (list): lower and upper boundary for latitude-grid
            lon_boundaries (list): lower and upper boundary for longitude-grid
            resolution (float): resolution of lat/lon grid [degree] (default: 1.0)
        """
        t1 = time.clock()
        
        instrument = FCDRs[0].instrument
        if instrument == 'HIRS':
            uth_channel = 12
            u_types = ['independent', 'structured'] # this will change in newer FCDR format
        elif instrument == 'MHS' or instrument == 'AMSUB':
            uth_channel = 3
            u_types = ['independent', 'common', 'structured']
    
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
        UTH_gridded = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        u_Tb_gridded = {t: {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches} for t in u_types}
        u_uth_gridded = {t: {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches} for t in u_types}
        count = {b: np.zeros((len(lat_centers), len(lon_centers))) for b in branches}
        count_all = {b: np.zeros((len(lat_centers), len(lon_centers))) for b in branches}
        count_overpasses = {b: np.zeros((len(lat_centers), len(lon_centers))) for b in branches}
        second_of_day_min = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        second_of_day_max = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        start_time = dict.fromkeys(branches)
        end_time = dict.fromkeys(branches)
 
        # collect data from all FCDRs
        collected_data, collected_data_diff, collected_files = utils.collectFCDRData(
                FCDRs, u_types, uth_channel=uth_channel)

        for b in branches:
            print(b)
#            longitudes[b][longitudes[b] == -180] = 180.
#            longitudes_all[b][longitudes_all[b] == -180] = 180.
            #TODO: stimmt das so???
            collected_data['longitude'][b][collected_data['longitude'][b] < -180 + 0.5 * resolution] = 180.
            collected_data_diff['longitude'][b][collected_data_diff['longitude'][b] < -180 + 0.5 * resolution] = 180.
            # combine all values of this branch in a pandas dataframe
            data = pd.DataFrame(utils.flattenDict(collected_data, b))
            data_diff = pd.DataFrame(utils.flattenDict(collected_data_diff, b))
            
            # throw away data outside the specified new grid
            data = data[np.logical_and((data['latitude'] <= lat_bins[-1]),(data['latitude'] > lat_bins[0]))].reset_index()
            data = data[np.logical_and((data['longitude'] <= lon_bins[-1]),(data['longitude'] > lon_bins[0]))].reset_index()
            data_diff = data_diff[np.logical_and((data_diff['latitude'] <= lat_bins[-1]),(data_diff['latitude'] > lat_bins[0]))].reset_index()
            data_diff = data_diff[np.logical_and((data_diff['longitude'] <= lon_bins[-1]),(data_diff['longitude'] > lon_bins[0]))].reset_index()
#            fig, ax = plt.subplots()
#            ax.set_title(b)
#            ax.scatter(data.longitude, data.latitude, c=data.brightness_temp, s=1)
#            ax.set_xlim(lon_bins[0], lon_bins[-1])
#            ax.set_ylim(lat_bins[0], lat_bins[-1])

            if data.empty and data_diff.empty:
                # no overpass in the selected grid
                print('dataframe is empty!')
            else:
                # bin data to latitude and longitude bins
                data = utils.binData(
                        data, lat_bins, lon_bins, lat_centers, lon_centers)
                data_diff = utils.binData(
                        data_diff, lat_bins, lon_bins, lat_centers, lon_centers)
                                
                # group data by latitude and longitude bins
                data_grouped = data.groupby([data.lat_bin, data.lon_bin], sort=False)
                data_diff_grouped = data_diff.groupby([data_diff.lat_bin, data_diff.lon_bin], sort=False)
                
                # get time of first and last data point going into this CDR
                start_time[b] = datetime.fromtimestamp(np.min(data.acquisition_time))
                end_time[b] = datetime.fromtimestamp(np.max(data.acquisition_time))
                
                # go through all groups and propagate uncertainties
                for name, group in data_grouped:
                    #group_ind = groups.get(group).values
                    lat_ind = name[0]#group[0]
                    lon_ind = name[1]#group[1]
                    group_size = len(group)#len(group_ind)
                    
                    # put averages of current group in the right cell on the grid 
                    Tb_gridded[b][lat_ind, lon_ind] = group.brightness_temp.mean()#Tb_groupmean[(lat_ind, lon_ind)]
                    UTH_gridded[b][lat_ind, lon_ind] = group.uth.mean()#UTH_groupmean[(lat_ind, lon_ind)]

                    # put count 
                    count[b][lat_ind, lon_ind] = group_size
                    count_all[b][lat_ind, lon_ind] = group_size
                    count_overpasses[b][lat_ind, lon_ind] += 1
                    second_of_day_group = np.array(group.second_of_day)#second_of_day_grouped[(lat_ind, lon_ind)]
                    second_of_day_min[b][lat_ind, lon_ind] = np.min(second_of_day_group)
                    second_of_day_max[b][lat_ind, lon_ind] = np.max(second_of_day_group)

                    # Get structured, independent and common uncertainties of this group
                    for t in u_types:
                        group_variances_Tb[t] = np.array(group['u_Tb_{}'.format(t)])
                        group_variances_uth[t] = np.array(group['u_uth_{}'.format(t)])
                    
                    # Create a covariance matrix for structured uncertainties:
                    # scanlines of data points in this group
                    scanlines = np.array(group.scanline)
                    # construct correlation matrix from scanline differences
                    scanlines_h = np.reshape(scanlines, (len(scanlines), 1))
                    corr = np.maximum(np.zeros((group_size, group_size)), 1 - np.abs(scanlines_h - scanlines_h.T) / 8)
                    # calculate covariance matrix from correlation matrix
                    S_struct_Tb = np.multiply(corr, np.outer(group_variances_Tb['structured'], group_variances_Tb['structured']))
                    S_struct_uth = np.multiply(corr, np.outer(group_variances_uth['structured'], group_variances_uth['structured']))
                    
                    # Calculate uncertainties (standard deviations) for the current group: 
                    # independent uncertainty for this grid cell (Use Law of the Propagation of Uncertainties for independent uncertainties only)
                    u_Tb_gridded['independent'][b][lat_ind, lon_ind] = np.sqrt(np.sum(group_variances_Tb['independent'] ** 2)) / group_size
                    u_uth_gridded['independent'][b][lat_ind, lon_ind] = np.sqrt(np.sum(group_variances_uth['independent'] ** 2)) / group_size
                    # common uncertainty for this grid cell (fully correlated uncertainties --> uncertainty of average is average of uncertainties)
                    # Only MW FCDRs contain common uncertainties
                    if instrument == 'MHS' or instrument == 'AMSUB':
                        u_Tb_gridded['common'][b][lat_ind, lon_ind] = np.mean(group_variances_Tb['common'])
                        u_uth_gridded['common'][b][lat_ind, lon_ind] = np.mean(group_variances_uth['common'])
                    # structured uncertainty for this grid cell
                    u_Tb_gridded['structured'][b][lat_ind, lon_ind] = np.sqrt(np.sum(S_struct_Tb)) / group_size
                    u_uth_gridded['structured'][b][lat_ind, lon_ind] = np.sqrt(np.sum(S_struct_uth)) / group_size
                
                # go through all groups that would contain additional data without cloud filtering:
                for name, group in data_diff_grouped:
                    lat_ind = name[0]
                    lon_ind = name[1]
                    group_size = len(group)
                    count_all[b][lat_ind, lon_ind] += group_size

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
        ret.uth = UTH_gridded
        ret.observation_count = count
        ret.observation_count_all = count_all
        ret.observation_count_overpasses = count_overpasses
        # time information
        ret.time_coverage_start = min(start_time.values())
        ret.time_coverage_end = max(end_time.values())
        ret.second_of_day_min = second_of_day_min
        ret.second_of_day_max = second_of_day_max
        # FCDR files included
        ret.source = collected_files
        # instrument and satellite
        ret.instrument = FCDRs[0].instrument
        ret.satellite = FCDRs[0].satellite
        # uncertainties
        ret.u_Tb = u_Tb_gridded
        ret.u_uth = u_uth_gridded
        t2 = time.clock()
        print(t2 - t1)
        return ret

    @classmethod
    def AveragedCDRFromCDRs(cls, CDRs):
        """ Combines several CDRs to one CDR by calculating an ordinary average
        of brightness temperature and UTH for every grid cell. Uncertainties 
        are propagated using the Law of the Propagation of Uncertainty for the
        case of an ordenary average.
        
        Parameters:
            CDRs (list of CDRs): List of CDR objects (e.g. 30 daily CDRs that
                 should be combined to one monthly CDR)
            uncertainties (boolean): False, if uncertainties should not be 
                propagated
        """
        ret = cls()
        
        instrument = CDRs[0].instrument
        if instrument == 'HIRS':
            u_types = ['independent', 'structured']
        elif instrument == 'MHS' or instrument == 'AMSUB':
            u_types = ['independent', 'structured', 'common']
        
        branches = ['descending', 'ascending']
        num_timesteps = len(CDRs)
        Tb_mean = dict.fromkeys(branches)
        UTH_mean = dict.fromkeys(branches)
        u_Tb = {t: dict.fromkeys(branches) for t in u_types}
        u_uth = {t: dict.fromkeys(branches) for t in u_types}
        Tb_std = dict.fromkeys(branches)
        UTH_std = dict.fromkeys(branches)
        count = dict.fromkeys(branches)
        count_all = dict.fromkeys(branches)
        count_overpasses = dict.fromkeys(branches)
        second_of_day_min = dict.fromkeys(branches)
        second_of_day_max = dict.fromkeys(branches)
        time_ranges = dict.fromkeys(branches)
        files = []

        time_coverage_start = CDRs[0].time_coverage_start
        time_coverage_end = CDRs[-1].time_coverage_end
        time_coverage_duration = pd.Timedelta(time_coverage_end - time_coverage_start).isoformat()

        for b in branches:
            for i in range(num_timesteps):
                files.extend(CDRs[i].source)
            notnan_count = np.sum([~np.isnan(CDRs[i].BT[b]) for i in range(num_timesteps)], axis=0)
            Tb_mean[b] = np.nanmean([CDRs[i].BT[b] for i in range(num_timesteps)], axis=0)
            UTH_mean[b] = np.nanmean([CDRs[i].uth[b] for i in range(num_timesteps)], axis=0) 
            Tb_std[b] = np.nanstd([CDRs[i].BT[b] for i in range(num_timesteps)], axis=0)
            UTH_std[b] = np.nanstd([CDRs[i].uth[b] for i in range(num_timesteps)], axis=0)
            count[b] = np.nansum([CDRs[i].observation_count[b] for i in range(num_timesteps)], axis=0)
            count[b][count[b] == 0] = -32767
            count_all[b] = np.nansum([CDRs[i].observation_count_all[b] for i in range(num_timesteps)], axis=0)
            count_all[b][count_all[b] == 0] = -32767
            count_overpasses[b] = np.nansum([CDRs[i].observation_count_overpasses[b] for i in range(num_timesteps)], axis=0)
            
#            a_time_coverage_start[b] = CDRs[0].a_time_coverage_start[b]
#            a_time_coverage_end[b] = CDRs[-1].a_time_coverage_end[b]
#            a_time_coverage_duration[b] = pd.Timedelta(a_time_coverage_end[b] - a_time_coverage_start[b]).isoformat()
            
            second_of_day_min[b] = np.nanmin([CDRs[i].second_of_day_min[b] for i in range(num_timesteps)], axis=0)
            second_of_day_min[b][np.isnan(second_of_day_min[b])] = 4294967295
            second_of_day_max[b] = np.nanmax([CDRs[i].second_of_day_max[b] for i in range(num_timesteps)], axis=0)
            second_of_day_max[b][np.isnan(second_of_day_max[b])] = 4294967295
            
            time_ranges[b] = np.stack((second_of_day_min[b], second_of_day_max[b]))
            
            notnan_count = np.sum([~np.isnan(CDRs[i].u_Tb['independent'][b]) for i in range(num_timesteps)], axis=0)
            for t in ['independent', 'structured']:
                u_Tb[t][b] = np.sqrt(np.nansum([CDRs[i].u_Tb[t][b] ** 2 for i in range(num_timesteps)], axis=0)) / notnan_count
                u_uth[t][b] = np.sqrt(np.nansum([CDRs[i].u_uth[t][b] ** 2 for i in range(num_timesteps)], axis=0)) / notnan_count
                
            if instrument == 'MHS' or instrument == 'AMSUB':
                u_Tb['common'][b] = np.nanmean([CDRs[i].u_Tb['common'][b] for i in range(num_timesteps)], axis=0)
                u_uth['common'][b] = np.nanmean([CDRs[i].u_uth['common'][b] for i in range(num_timesteps)], axis=0)
        
        ret.lon = CDRs[0].lon
        ret.lat = CDRs[0].lat
        ret.geospatial_lat_resolution = CDRs[0].geospatial_lat_resolution
        ret.geospatial_lon_resolution = CDRs[0].geospatial_lon_resolution
        ret.lat_bnds = np.stack((CDRs[0].lat_bins[0:-1], CDRs[0].lat_bins[1:])).T
        ret.lon_bnds = np.stack((CDRs[0].lon_bins[0:-1], CDRs[0].lon_bins[1:])).T
        ret.BT = Tb_mean
        ret.uth = UTH_mean
        ret.BT_inhomogeneity = Tb_std
        ret.uth_inhomogeneity = UTH_std
        ret.observation_count = count
        ret.observation_count_all = count_all
        ret.observation_count_overpasses = count_overpasses
        ret.time_coverage_start = time_coverage_start 
        ret.time_coverage_end = time_coverage_end
        ret.time_coverage_duration = time_coverage_duration
        ret.time_ranges = time_ranges
        ret.instrument = CDRs[0].instrument
        ret.satellite = CDRs[0].satellite
        ret.u_BT = u_Tb
        ret.u_uth = u_uth
        ret.source = files

        return ret
    
    def toNetCDF(self, CDR_path, comment_on_version):
        """ Saves the attributes of an averaged CDR created by 
        AveragedCDRFromCDRs to a NetCDF file using the CDRWriter by Tom Block.
        
        Parameters:
            comment_on_version (str): Comment on this version of the CDR.
        """
        instrument = self.instrument
        if instrument == 'HIRS':
            u_types = ['structured', 'independent']
        elif instrument == 'MHS' or instrument == 'AMSUB':
            u_types = ['structured', 'common', 'independent']
            
        branches = ['ascend', 'descend']
        u_quantities = ['uth', 'BT']
        DATE_PATTERN = "%Y%m%d%H%M%S"
        numlons = len(self.lon)
        numlats = len(self.lat)
        start_time = self.time_coverage_start
        end_time = self.time_coverage_end
        CDR_filename = CDRWriter.create_file_name_CDR('UTH', sensor=self.instrument, platform=self.satellite, start=start_time, end=end_time, type='L3', version='1.0')
        ds = CDRWriter.createTemplate('UTH', numlons, numlats)
        simple_vars = ['lat', 'lon', 'lat_bnds', 'lon_bnds']
        branch_vars = ['observation_count', 'observation_count_all',
                       'uth', 'BT',
                       'uth_inhomogeneity', 'BT_inhomogeneity']
        simple_attrs = ['source', 'time_coverage_duration',
                        'geospatial_lon_resolution', 'geospatial_lat_resolution']
        
        for var in simple_vars:
            ds.variables[var].data = getattr(self, var)
        
        for attr in simple_attrs:
            ds.attrs[attr] = getattr(self, attr)
            
        for var in branch_vars:
            for b in branches:
                ds.variables['{}_{}'.format(var, b)].data = getattr(self, '{}'.format(var))['{}ing'.format(b)]
        
        for q in u_quantities:
            for t in u_types:
                for b in branches:
                    ds.variables['u_{}_{}_{}'.format(t, q, b)].data = getattr(self, 'u_{}'.format(q))[t]['{}ing'.format(b)]
        
        #im moment noch falsch benannt: time_ranges_ascending, soll später time_ranges_ascend heißen:
        ds.variables['time_ranges_ascending'].data = getattr(self, 'time_ranges')['ascending']
        ds.variables['time_ranges_descending'].data = getattr(self, 'time_ranges')['descending']

        ds.attrs['institution'] = 'University of Hamburg'
        ds.attrs['title'] = 'Upper Tropospheric Humidity (UTH)'
        ds.attrs['time_coverage_start'] = start_time.strftime(DATE_PATTERN)
        ds.attrs['time_coverage_end'] = end_time.strftime(DATE_PATTERN)
        ds.attrs['time_coverage_resolution'] = pd.Timedelta(1, unit='M').isoformat()
        ds.attrs['geospatial_lat_units'] = 'deg'
        ds.attrs['geospatial_lon_units'] = 'deg'
        ds.attrs['history'] = 'Created on {}.'.format(datetime.now().strftime('%Y-%m-%d %H:%M'))
        ds.attrs['references'] = ''
        ds.attrs['comment'] = ''
        ds.attrs['auxiliary_data'] = ''
        ds.attrs['configuration'] = ''
        
        CDRWriter.write(ds, join(CDR_path, CDR_filename))
        return ds
    
    @classmethod
    def fromNetCDF(cls, filename):
        """ Creates a CDR-Object from a NetCDF file.
        
        Parameters:
            filename: path and name of NetCDF file
        """
        ret = cls()
        f = Dataset(filename)
        print(f.source)
        DATE_PATTERN = "%Y%m%d%H%M%S"
        branches = ['ascend', 'descend']
        u_types = ['structured', 'common', 'independent']
        u_quantities = ['uth', 'BT']
        simple_vars = ['lat', 'lon', 'lat_bnds', 'lon_bnds']
        branch_vars = ['observation_count', 'observation_count_all',
                       'uth', 'BT',
                       'uth_inhomogeneity', 'BT_inhomogeneity']
        simple_attrs = ['source', 'time_coverage_duration',
                        'geospatial_lon_resolution', 'geospatial_lat_resolution',
                        'geospatial_lat_units', 'geospatial_lon_units']
        
        for var in simple_vars:
            setattr(ret, var, f.variables[var][:])
        
        for attr in simple_attrs:
            setattr(ret, attr, getattr(f, attr))
            
        for var in branch_vars:
            print(var)
            h = {}
            for b in branches:
                h['{}ing'.format(b)] = f.variables['{}_{}'.format(var, b)][:]
            setattr(ret, var, h)
        
        for q in u_quantities:
            h = {}
            for t in u_types:
                h[t] = {}
                for b in branches:
                    h[t]['{}ing'.format(b)] = f.variables['u_{}_{}_{}'.format(t, q, b)]
            setattr(ret, 'u_{}'.format(q), h)
        
        ret.time_ranges = {}
        print(f.variables['time_ranges_ascending'][0])
        ret.time_ranges['ascending'] = f.variables['time_ranges_ascending']
        ret.time_ranges['descending'] = f.variables['time_ranges_descending']
        ret.time_coverage_start = datetime.strptime(f.time_coverage_start, DATE_PATTERN)
        ret.time_coverage_end = datetime.strptime(f.time_coverage_end, DATE_PATTERN)
        
        return ret
        
        
    @classmethod
    def GriddedCDRFromFCDRs_fullPropagation(cls, FCDRs, lat_boundaries=[-90, 90], 
                                           lon_boundaries =[-179, 180], 
                                           channel=3, resolution=1., 
                                           uncertainties=True):
        """ Creates one brightness temperature CDR from from a list of FCDRs 
        and does the full uncertainty propagation. Only needed for examplary
        uncertainty propagation for structured uncertainties. 
        
        WARNING: for big grids this consumes too much memory!!!
        
        The CDR contains mean brightness temperatures binned to a lat/lon grid
        and uncertainty covariance matrices.
        
        Parameters:
            FCDRs (list): list fo FCDR objects
            lat_boundaries (list): lower and upper boundary for latitude-grid
            lon_boundaries (list): lower and upper boundary for longitude-grid
            channel (int): brightness temperature channel that should be used
                (default: 3)
            resolution (float): resolution of lat/lon grid [degree] (default: 1.0)
        """
        ret = cls()
        # latitude and longitude bins and their centers
        lat_bins = np.arange(lat_boundaries[0] - 0.5 * resolution, lat_boundaries[-1] + 0.5 * resolution + 1., resolution)
        lat_centers = [(a + b) / 2 for a, b in zip(lat_bins, lat_bins[1:])]
        lon_bins = np.arange(lon_boundaries[0] - 0.5 * resolution, lon_boundaries[-1] + 0.5 * resolution + 1, resolution)
        lon_centers = [(a + b) / 2 for a, b in zip(lon_bins, lon_bins[1:])]
        numcells = len(lat_centers) * len(lon_centers)
        
        start_time = min([FCDRs[i].start_time for i in range(len(FCDRs))])
        end_time = max([FCDRs[i].end_time for i in range(len(FCDRs))])
        # initialization
        branches = ['ascending', 'descending']
        node_mask = {}
        latitudes = {b: np.array([]) for b in branches}
        longitudes = {b: np.array([]) for b in branches}
        brightness_temp = {b: np.array([]) for b in branches}
        data_dict = dict.fromkeys(branches)
        data = dict.fromkeys(branches)
        data_grouped = dict.fromkeys(branches)
        Tb_groupsum = dict.fromkeys(branches)
        Tb_gridded = {b: np.ones((len(lat_centers), len(lon_centers))) * np.nan for b in branches}
        count = {b: np.zeros((len(lat_centers), len(lon_centers))) for b in branches}

        u_ind = {b: np.array([]) for b in branches}
        u_com = {b: np.array([]) for b in branches}
        u_struct = {b: np.array([]) for b in branches} 
        scanline = {b: np.array([]) for b in branches}
        s_max = {b: 0 for b in branches} 

        if uncertainties:
            Sm_ind = dict.fromkeys(branches)
            Sm_com = dict.fromkeys(branches)
            Sm_struct = dict.fromkeys(branches)
        
        
        # functions for conversion from longitudes and latitudes to grid indices
        lon2ind = lambda lon: lon_centers.index(lon)
        lat2ind = lambda lat: lat_centers.index(lat)

        # collect data from all FCDRs
        for f in FCDRs:
            # get all masks
            if not hasattr(f, 'node_mask'):
                f.generate_node_mask()
            
            node_mask['ascending'] = f.node_mask
            node_mask['descending'] = ~f.node_mask
            total_mask = f.total_mask[channel]

            # distinguish between ascending and descending node and apply all masks
            for b in branches:
                latitudes[b] = np.append(latitudes[b], f.latitudes[np.logical_and(~total_mask, node_mask[b])])
                longitudes[b] = np.append(longitudes[b], f.longitudes[np.logical_and(~total_mask, node_mask[b])])
                brightness_temp[b] = np.append(brightness_temp[b], f.brightness_temp[channel][np.logical_and(~total_mask, node_mask[b])])
                u_ind[b] = np.append(u_ind[b], f.u_Tb['independent'][channel][np.logical_and(~total_mask, node_mask[b])])
                u_com[b] = np.append(u_com[b], f.u_Tb['common'][channel][np.logical_and(~total_mask, node_mask[b])])
                u_struct[b] = np.append(u_struct[b], f.u_Tb['structured'][channel][np.logical_and(~total_mask, node_mask[b])])
                scanline[b] = np.append(scanline[b], s_max[b] + f.scanline[np.logical_and(~total_mask, node_mask[b])])
                s_max[b] = np.max(scanline[b])
            
        for b in branches:
            longitudes[b][longitudes[b] < -180 + 0.5 * resolution] = 180.
            # combine all values in a pandas dataframe
            
            if uncertainties:
                data_dict[b] = {'latitude': latitudes[b], 'longitude': longitudes[b],\
                                'Tb': brightness_temp[b], 'U_ind': u_ind[b],\
                                'U_com': u_com[b], 'U_struct': u_struct[b], 
                                'Var_tot': u_com[b] ** 2 + u_struct[b] ** 2 + u_ind[b] ** 2,
                                'scanline': scanline[b]}
            else:
                data_dict[b] = {'latitude': latitudes[b], 'longitude': longitudes[b],\
                                'Tb': brightness_temp[b]}
                
            data[b] = pd.DataFrame(data_dict[b])
            # bin data to latitude and longitude bins
            data[b] = data[b][np.logical_and((data[b]['latitude'] <= lat_bins[-1]),(data[b]['latitude'] > lat_bins[0]))].reset_index()
            data[b] = data[b][np.logical_and((data[b]['longitude'] <= lon_bins[-1]),(data[b]['longitude'] > lon_bins[0]))].reset_index()
            
            if data[b].empty:
                print('dataframe is empty!')
                Tb_gridded[b] = np.ones((len(lon_centers), len(lat_centers))) * np.nan
                Sm_struct[b] = csr_matrix((numcells, numcells))
                Sm_ind[b] = csr_matrix((numcells, numcells))
                Sm_com[b] = csr_matrix((numcells, numcells))
            else:
#                fig, ax = plt.subplots()
#                ax.scatter(data[b].longitude, data[b].latitude, c=data[b].Tb)
#                ax.set_title(b)
#                ax.set_xlim(lon_bins[0], lon_bins[-1])
#                ax.set_ylim(lat_bins[0], lat_bins[-1])
#                fig, ax = plt.subplots(2,2)
#                ax[0, 0].scatter(data[b].longitude, data[b].latitude, c=data[b].U_struct)
#                ax[0, 0].set_xlim(lon_bins[0], lon_bins[-1])
#                ax[0, 0].set_ylim(lat_bins[0], lat_bins[-1])
                data[b]['lat_bin'] = pd.cut(data[b].latitude, lat_bins, labels=lat_centers)
                data[b]['lon_bin'] = pd.cut(data[b].longitude, lon_bins, labels=lon_centers)
                data[b] = data[b].drop(['latitude', 'longitude'], 1)

                # group data by latitude and longitude bins
                data_grouped[b] = data[b].groupby([data[b].lat_bin, data[b].lon_bin])
                if uncertainties:
#                    data[b]['Tb_weighted'] = data[b].Tb / data[b].Var_tot
#                    data[b]['norm'] = 1. / data[b].Var_tot
                    # sum of weighted brightness temperature and normalization factor in each bin
#                    Tb_weighted_groupsum[b] = data_grouped[b].Tb_weighted.sum()
#                    norm_groupsum[b] = data_grouped[b].norm.sum()
                    # number of data points
                    numel = len(data[b].Tb)
                    # preallocate arrays to collect rows, columns and values for 
                    # auxiliary matrix E...
                    E_columns = np.zeros(numel, dtype=np.int)
                    E_rows = np.zeros(numel, dtype=np.int)
                    E_values = np.zeros(numel, dtype=np.float16)
                    # ... and covariance matrix with independent uncertainties
                    S_ind_values = np.zeros(numel, dtype=np.float16)
                
                Tb_groupsum[b] = data_grouped[b].Tb.sum()
                
                row = 0
                col = 0

                print('Tb_gridded, S_ind')
                groups = data_grouped[b].groups
                for lat in lat_centers:
                    for lon in lon_centers:
                        group_ind = groups.get((lat, lon))
                        if group_ind is not None:
                            group_size = len(group_ind)
                            
                            if uncertainties:
                                # values and their positions for auxiliary array E
                                E_rows[col:col+group_size] = [row for i in range(group_size)]
                                E_columns[col:col+group_size] = group_ind
                                #E_values[col:col+group_size] = [1. / group_size for i in range(group_size)]
                                E_values[col:col+group_size] = 1. / group_size
                                S_ind_values[col:col+group_size] = data[b].U_ind[group_ind] ** 2

                            Tb_gridded[b][lat2ind(lat), lon2ind(lon)] = Tb_groupsum[b][lat, lon] / group_size
                            
                            count[b][lat2ind(lat), lon2ind(lon)] += group_size
                            #Tb_gridded[b][lat2ind(lat), lon2ind(lon)] = Tb_groupmean[b][lat, lon]
                        else:
                            group_size = 0
                        
                        row += 1
                        col += group_size
                
                
                if uncertainties:
                    # compose covariance matrix with structured uncertainties
                    print('S_struct')
                    
                    scanlines = np.array(data[b].scanline)
                    # construct correlation matrix from scanline differences
                    scanlines_h = np.reshape(scanlines, (len(scanlines), 1))
                    corr = np.maximum(np.zeros((len(scanlines), len(scanlines))), 1 - np.abs(scanlines_h - scanlines_h.T) / 8)

                    # calculate covariance matrix from correlation matrix
                    S_struct = csr_matrix(np.multiply(corr, np.outer(data[b].U_struct, data[b].U_struct)))
                    
                    S_ind = diags(S_ind_values)
                    #S_struct = csr_matrix((S_struct_values, (S_struct_rows, S_struct_cols)), shape=(numel, numel))
    
                    
                    # create auxiliary matrix E
                    print('E')
                    E = csr_matrix((E_values, (E_rows, E_columns)), shape=(numcells, numel))
                    print(E.todense())
                    
#                    ax[1, 0].imshow(E.todense())
#                    ax[1, 1].imshow(S_struct.todense())
#                    fig.suptitle(b)
    
                    print('matrix multiplication')
                    # Calculate covariance matrix with independent uncertainties for all grid cells 
                    Sm_ind[b] = E.dot(S_ind.dot(E.T))
                    # Calculate covariance matrix with common uncertainties for all grid cell means
                    a1 = csr_matrix(E * data[b].U_ind[:, np.newaxis])
                    a2 = csc_matrix(data[b].U_ind[np.newaxis, :] * E.T)
                    Sm_com[b] = a1.dot(a2)
                    # Calculate covariance matrix with structured uncertainties for all grid cell means
                    as1 = E.dot(S_struct)
                    Sm_struct[b] = as1.dot(E.T)
#                    ax[0, 1].imshow(Sm_struct[b].todense())

            
        # return 
        ret.lat_grid = np.array(lat_centers)
        ret.lon_grid = np.array(lon_centers)
        ret.brightness_temp_mean = Tb_gridded
        ret.count = count
        # time information
        ret.start_time = start_time
        ret.end_time = end_time
        
        if uncertainties:

            ret.S_independent = Sm_ind
            ret.S_common = Sm_com
            ret.S_structured = Sm_struct
        
        return ret 

    @classmethod
    def AveragedCDRFromCDRs_fullPropagation(cls, CDRs, uncertainties=True):
        """ Creates a new CDR by taking the average of several CDRs. The full
        uncertainty propagation is performed for structured uncertainties.
        Mittelwerte werden OHNE Wichtung der eingehenden Werte mit ihren 
        Unsicherheiten gebildet!
        """

        ret = cls()
        
        branches = ['descending', 'ascending']
        Tb_mean = dict.fromkeys(branches)
        x_mean = dict.fromkeys(branches)
        
        if uncertainties:
            S_struct_mean = dict.fromkeys(branches)
            S_com_mean = dict.fromkeys(branches)
            S_ind_mean = dict.fromkeys(branches)
            q_struct_mean = dict.fromkeys(branches)

        for b in branches:
            Tb_mean[b] = CDRs[0].brightness_temp_mean[b]
            x_mean[b] = Tb_mean[b].reshape(Tb_mean[b].size)
            
            if uncertainties:
                S_struct_mean[b] = CDRs[0].S_structured[b]
                S_com_mean[b] = CDRs[0].S_common[b]
                S_ind_mean[b] = CDRs[0].S_independent[b]
        
        numcells = Tb_mean['ascending'].size
        grid_shape = Tb_mean['ascending'].shape
        count = {b: np.zeros(grid_shape) for b in branches}
        
        for c in CDRs[1:]:
            for b in branches:
                Tb_new = c.brightness_temp_mean[b]
                x_new = Tb_new.reshape((Tb_new.size))
                count[b] += c.count[b]
                
                if uncertainties:
                    # create new covariance matrices by appending the one
                    # from this timestep to the "mean" one
                    S_struct_new = c.S_structured[b]
                    S_struct_total = block_diag((S_struct_mean[b], S_struct_new))
                    S_com_new = c.S_common[b]
                    S_com_total = block_diag((S_com_mean[b], S_com_new))
                    S_ind_new = c.S_independent[b]
                    S_ind_total = block_diag((S_ind_mean[b], S_ind_new))
                    
                    if np.any(S_struct_total.todense() < 0):
                        print(S_struct_total)
                                                                 
                    # Create auxiliary matrix E:
                    
                    # extract some useful vectors
#                    sigma_mean = S_struct_mean[b].diagonal() + S_com_mean[b].diagonal() + S_ind_mean[b].diagonal()
#                    sigma_new = S_struct_new.diagonal() + S_com_new.diagonal() + S_ind_new.diagonal()
                    
                    # create diagonals for left part of E:
                    H1 = np.zeros(numcells)
                    H1[np.logical_and(np.logical_not(np.isnan(x_mean[b])), np.isnan(x_new))] = 1.
                    x_mean_new_ind = np.logical_and(np.logical_not(np.isnan(x_new)), np.logical_not(np.isnan(x_mean[b])))
#                    H1[x_mean_new_ind] = 1. / sigma_mean[x_mean_new_ind] / (1. / sigma_mean[x_mean_new_ind] + 1. / sigma_new[x_mean_new_ind])
                    H1[x_mean_new_ind] = 1. / 2
                    # create diagonals for right part of E
                    H2 = np.zeros(numcells)
                    H2[np.logical_and(np.logical_not(np.isnan(x_new)), np.isnan(x_mean[b]))] = 1.
#                    H2[x_mean_new_ind] = 1. / sigma_new[x_mean_new_ind] / (1. / sigma_mean[x_mean_new_ind] + 1. / sigma_new[x_mean_new_ind])
                    H2[x_mean_new_ind] = 1. / 2
                    
                    # Add both parts together 
                    E1 = diags(H1)
                    E2 = diags(H2)
                    E = bmat([[E1, E2]])
                    
                    # calculate new covariance matrices with the law of propagation 
                    # of uncertainties
                    S_struct_mean[b] = E.dot(S_struct_total.dot(E.T))
                    #if np.any(S_struct_mean[b].todense() < 0):
                    print(S_struct_mean[b])
                    norm = np.sqrt(S_struct_mean[b].diagonal())
                    norm_matrix = np.outer(norm, norm)
                    q_struct_mean[b] = S_struct_mean[b] / norm_matrix
                    S_com_mean[b] = E.dot(S_com_total.dot(E.T))
                    S_ind_mean[b] = E.dot(S_ind_total.dot(E.T))
                    
                    # look where brightness temperature is not nan                
                    x_mean_ind = np.logical_not(np.isnan(x_mean[b]))
                    x_new_ind = np.logical_not(np.isnan(x_new))
    
                    x_mean_new = np.ones(numcells) * np.nan
                    # where x_new is nan and x_mean is not, the new x_mean will be
                    # the same as the old one
                    x_mean_new[np.logical_and(x_mean_ind, np.logical_not(x_new_ind))]\
                    = x_mean[b][np.logical_and(x_mean_ind, np.logical_not(x_new_ind))]
                    # where x_mean is nan and x_new is not, x_new will be the new
                    # x_mean
                    x_mean_new[np.logical_and(x_new_ind, np.logical_not(x_mean_ind))]\
                    = x_new[np.logical_and(x_new_ind, np.logical_not(x_mean_ind))]
                    # where x_new and x_mean are not nan, the new mean is a weighted
                    # average of the two 
#                    x_mean_new[x_mean_new_ind] =\
#                    (x_mean[b][x_mean_new_ind] / sigma_mean[x_mean_new_ind]\
#                     + x_new[x_mean_new_ind] / sigma_new[x_mean_new_ind])\
#                     / (1. / sigma_mean[x_mean_new_ind] + 1. / sigma_new[x_mean_new_ind])
                    x_mean_new[x_mean_new_ind] = (x_mean[b][x_mean_new_ind] + x_new[x_mean_new_ind]) / 2
                    
                    x_mean[b] = x_mean_new
                    
                else:
                    x_mean[b] = np.nanmean(np.vstack((x_mean[b], x_new)), axis=0)
                    
        
        # reshape x to get brightness temperatures in an array shaped as the grid
        for b in branches:
            Tb_mean[b] = x_mean[b].reshape(grid_shape)
        
        ret.brightness_temp_mean = Tb_mean
        ret.count = count
        ret.lat_grid = CDRs[0].lat_grid
        ret.lon_grid = CDRs[0].lon_grid
        
        if uncertainties:
            ret.S_struct = S_struct_mean
            ret.q_struct = q_struct_mean
            ret.S_com = S_com_mean
            ret.S_ind = S_ind_mean

        return ret       
       
#%% OLD FCDR methods:
#def regrid_uth(self, resolution, radius_of_influence, epsilon, lat_limits=[-90, 90], lon_limits=[-180, 180], viewing_angles=-1):
#    """ First try to regrid UTH
#    """
#    if not hasattr(self, 'uth'):
#        print('UTH has to be calculated before it can be regridded.')
#        
#    if type(viewing_angles) == int:
#        viewing_angles = self.viewing_angles
#    
#    if not hasattr(self, 'total_mask'):
#        self.generate_total_mask()
#        
#    
#    mask = self.total_mask[self.uth_channel][:, viewing_angles]
#    latitudes = self.latitudes[:, viewing_angles]
#    longitudes = self.longitudes[:, viewing_angles]
#    # devide into ascending and descending branch and apply mask
#    is_ascending = np.diff(latitudes, axis=0) >= 0 
#    is_ascending = np.append(is_ascending, [is_ascending[-1,:]], axis=0)
#    latitudes_ascending = latitudes[np.logical_and(~mask, is_ascending)]
#    latitudes_descending = latitudes[np.logical_and(~mask, ~is_ascending)]
#    longitudes_ascending = longitudes[np.logical_and(~mask, is_ascending)]
#    longitudes_descending = longitudes[np.logical_and(~mask, ~is_ascending)]
#          
#    uth_ascending = self.uth[:, viewing_angles][np.logical_and(~mask, is_ascending)]
#    uth_descending = self.uth[:, viewing_angles][np.logical_and(~mask, ~is_ascending)]
#    
#    num_lons = int(len(np.arange(lon_limits[0], lon_limits[1])) / resolution) + 1 
#    num_lats = int(len(np.arange(lat_limits[0], lat_limits[1])) / resolution) + 1
#    # - 0.5*resolution because you have to specifiy lower left CORNER and 
#    # upper right CORNER of the lower left and upper right pixel, respectively
#    area_def = geometry.AreaDefinition('areaD', 'Earth', 'areaD',
#                                       {'lat_0': '0', 'lon_0': '0',
#                                        'proj': 'latlon'},
#                                        num_lons, num_lats,
#                                        [lon_limits[0] - 0.5*resolution, lat_limits[0] - 0.5*resolution,
#                                         lon_limits[1] + 0.5*resolution, lat_limits[1] + 0.5*resolution])
#        
#    swath_def_ascending = geometry.SwathDefinition(longitudes_ascending, latitudes_ascending)
#    swath_def_descending = geometry.SwathDefinition(longitudes_descending, latitudes_descending)
#    uth_gridded_ascending = kd_tree.resample_nearest(
#            swath_def_ascending, uth_ascending, area_def, radius_of_influence=radius_of_influence,\
#            epsilon=epsilon, fill_value=np.NaN)
#    uth_gridded_descending = kd_tree.resample_nearest(
#            swath_def_descending, uth_descending, area_def, radius_of_influence=radius_of_influence,\
#            epsilon=epsilon, fill_value=np.NaN)
##    uth_gridded, stddev, count = kd_tree.resample_gauss(
##            swath_def, uth, area_def, radius_of_influence=radius_of_influence,\
##            sigmas=12000, fill_value=np.NaN, reduce_data=False, with_uncert=True)
#    
#
##        self.grid_lons = area_def.get_lonlats()[0][0]
##        self.grid_lats = area_def.get_lonlats()[1][:,0]
##        self.uth_gridded = {}
##        self.uth_gridded['ascending'] = uth_gridded_ascending
##        self.uth_gridded['descending'] = uth_gridded_descending
#    return uth_gridded_ascending, uth_gridded_descending
    
#%% OLD CDR methods:
#@classmethod
#def from_FCDRs(cls, FCDRs, resolution=1., lat_limits=[-90, 90], 
#               lon_limits=[-180, 180], radius_of_influence=15000, 
#               epsilon=0.1):
#    """ Creates a CDR object from a list of FCDR objects. The uth of the 
#    FCDRs is binned to a lat/lon grid and the mean is calculated.
#    
#    Parameter:
#        FCDRs (list): list of FCDR objects which have an attribute uth
#        resolution (numeric): resolution of the lat/lon grid in degrees
#        lat_limits (list): lower and upper limit for latitudes in the 
#            lat/lon grid in deg (e.g. [-180, 180])
#        lon_limits (list): lower and upper limit for longitues in the 
#            lat/lon grid in deg (e.g. [-180, 180])
#        radius_of_influence (numeric): radius of influence for nearest-
#            neighbor method in m
#        epsilon (numeric): allowed uncertainty for nearest neighbor method
#    """
#    ret = cls()
#    ret.start_year = FCDRs[0].year
#    ret.end_year = FCDRs[-1].year
#    ret.start_month = FCDRs[0].month
#    ret.end_month = FCDRs[-1].month
#    ret.start_day = FCDRs[0].day
#    ret.end_day = FCDRs[-1].day
#    #ret.fcdr_files = [f.file for f in FCDRs]
#    
#    branches = ['ascending', 'descending']
#    lon_grid, lat_grid = utils.get_grid_lons_lats(
#            lon_limits, lat_limits, resolution)
#    
#    uth_gridded = {b: np.zeros((len(FCDRs), len(lat_grid), len(lon_grid))) for b in branches}
#    not_nan = {b: np.zeros((len(FCDRs), len(lat_grid), len(lon_grid))) for b in branches}
#    uth_mean = {b: np.zeros((len(lat_grid), len(lon_grid))) for b in branches}
#    uth_stddev = {b: np.zeros((len(lat_grid), len(lon_grid))) for b in branches} 
#    uth_count = {b: np.zeros((len(lat_grid), len(lon_grid))) for b in branches}
#    
#    for f, i in zip(FCDRs, range(len(FCDRs))):    
#        uth_gridded['ascending'][i], uth_gridded['descending'][i] = f.regrid_uth(
#                resolution=resolution, radius_of_influence=radius_of_influence,\
#                epsilon=epsilon, lat_limits=lat_limits, lon_limits=lon_limits,\
#                viewing_angles=-1)
#        for b in branches:
#            not_nan[b][i] = ~np.isnan(uth_gridded[b][i])
#    
#    for b in branches:    
#        uth_mean[b] = np.nanmean(uth_gridded[b], axis=0)
#        uth_stddev[b] = np.nanstd(uth_gridded[b], axis=0) 
#        uth_count[b] = np.sum(not_nan[b], axis=0)
#    
#    ret.lat_grid = lat_grid
#    ret.lon_grid = lon_grid
#    ret.uth = uth_mean
#    ret.uth_stddev = uth_stddev
#    ret.uth_count = uth_count
#    
#    return ret
        
        