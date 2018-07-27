#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 09:24:52 2018

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

version_comment = 'Very first trial version.'

fcdr_path = '/scratch/uni/u237/users/tlang/CDR/FCDRs'
#fcdr_path = '/scratch/uni/u237/user_data/ihans/FCDR/easy/v2_0fv1_1_4/'
#fcdr_path = '/scratch/uni/u237/user_data/ihans/FCDR/HIRS'
regr_params_path = '/scratch/uni/u237/users/tlang/UTH/radiative_transfer/uth_definition'
cdr_path = ''

start_date = datetime.date(2015, 1, 1)
end_date = datetime.date(2015, 1, 1)
instrument = 'MHS'
satellite = 'Noaa18'
if instrument == 'HIRS':
    viewing_angles = [i for i in range(16, 42)] # or 19-38?
else:
    viewing_angles= [i for i in range(29, 61)]
        
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
    ds = monthly_mean.toNetCDF(cdr_path, comment_on_version=version_comment)
#    # save with pickle
#    outfile = 'CDR_{}_{}_{}.pkl'.format(monthly_mean.time_coverage_start.strftime('%Y-%m'), monthly_mean.instrument, monthly_mean.satellite)
#    with open(outfile, 'wb') as output:
#        pickle.dump(monthly_mean, output, pickle.HIGHEST_PROTOCOL)
    
    last_day_of_month = datetime.date(date.year, date.month, monthrange(date.year, date.month)[1])

#%% Plot overview of important CDR variables
u_types = ['independent', 'structured', 'common']
b = 'ascending'
fig, ax = plt.subplots(5, 1, figsize=(10, 10))
fig.suptitle(monthly_mean.time_coverage_start.strftime('%Y-%m'))
lon, lat = np.meshgrid(monthly_mean.lon, monthly_mean.lat)
im = plots.plotTb(monthly_mean, b, ax=ax[0], cmap=plt.get_cmap('temperature', 20))
ax[0].set_title('Tb')
fig.colorbar(im, ax=ax[0])
im = plots.plotUTH(monthly_mean, b, ax=ax[1], cmap=plt.get_cmap('speed', 20))
ax[1].set_title('UTH')
fig.colorbar(im, ax=ax[1])
im = plots.plotCDRQuantity(monthly_mean.uth_inhomogeneity[b], monthly_mean.lat, monthly_mean.lon, ax=ax[2], cmap=plt.get_cmap('speed', 20))
ax[2].set_title('UTH standard deviation of daily averages')
fig.colorbar(im, ax=ax[2])
im = plots.plotCount(monthly_mean, node=b, ax=ax[3], cmap=plt.get_cmap('density', 20), vmin=0)
ax[3].set_title('count')
fig.colorbar(im, ax=ax[3])
im = plots.plotCDRQuantity(monthly_mean.observation_count_all[b] - monthly_mean.observation_count[b], monthly_mean.lat, monthly_mean.lon, ax=ax[4], cmap=plt.get_cmap('density', 20))
ax[4].set_title('count all - count ')
fig.colorbar(im, ax=ax[4])
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

fig, ax = plt.subplots(len(u_types), 1, figsize=(10, 10))
fig.suptitle(monthly_mean.time_coverage_start.strftime('%Y-%m'))
for i, t in enumerate(u_types):
    im = plots.plotTbUncertainty(monthly_mean, t, node=b, ax=ax[i])
    ax[i].set_title('{} uncertainties for Tb'.format(t))
    fig.colorbar(im, ax=ax[i])
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

fig, ax = plt.subplots(len(u_types), 1, figsize=(10, 10))
fig.suptitle(monthly_mean.time_coverage_start.strftime('%Y-%m'))
for i, t in enumerate(u_types):
    #im = ax[i].pcolormesh(lon, lat, monthly_mean.u_uth[t][b], cmap='Reds')
    im = plots.plotUTHUncertainty(monthly_mean, t, node=b, ax=ax[i])
    ax[i].set_title('{} uncertainties for UTH'.format(t))
    fig.colorbar(im, ax=ax[i])
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
#%% Plot time ranges
fig, ax = plt.subplots(2, 2, figsize=(18, 8))
im = plots.plotCDRQuantity(monthly_mean.time_ranges['ascending'][0] / 60 / 60, monthly_mean.lat, monthly_mean.lon, ax=ax[0, 0], cmap=plt.get_cmap('density', 24), vmin=0, vmax=24)
ax[0, 0].set_title('ascending min')
im = plots.plotCDRQuantity(monthly_mean.time_ranges['ascending'][1] / 60 / 60, monthly_mean.lat, monthly_mean.lon, ax=ax[1, 0], cmap=plt.get_cmap('density', 24), vmin=0, vmax=24)
ax[1, 0].set_title('ascending max')
im = plots.plotCDRQuantity(monthly_mean.time_ranges['descending'][0] / 60 / 60, monthly_mean.lat, monthly_mean.lon, ax=ax[0, 1], cmap=plt.get_cmap('density', 24), vmin=0, vmax=24)
ax[0, 1].set_title('descending min')
im = plots.plotCDRQuantity(monthly_mean.time_ranges['descending'][1] / 60 / 60, monthly_mean.lat, monthly_mean.lon, ax=ax[1, 1], cmap=plt.get_cmap('density', 24), vmin=0, vmax=24)
ax[1, 1].set_title('descending max')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

  

