#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 11:49:19 2018

@author: u300740

Read in UTH CDR and make some plots.
"""

import matplotlib.pyplot as plt
from processing import CDRs
from processing import plots
import numpy as np
#%% test: load NetCDF file 
monthly_mean = None
monthly_mean = CDRs.CDR.fromNetCDF('CDR_files/MHS/NOAA18/2015/FIDUCEO_CDR_UTH_MHS_NOAA18_20150101000000_20150131235959_L3_v1.0_fv1.1.5.nc')
#%% Plot overview of important CDR variables
u_types = ['independent', 'structured', 'common']
b = 'ascending'
xticklocs = np.arange(-180, 181, 60)
xticklabels = [str(i) + '°' for i in xticklocs]
yticklocs = np.arange(-30, 31, 15)
yticklabels = [str(i) + '°' for i in yticklocs]

fig, ax = plt.subplots(5, 1, figsize=(10, 10))
fig.suptitle(monthly_mean.time_coverage_start.strftime('%Y-%m'))
lon, lat = np.meshgrid(monthly_mean.lon, monthly_mean.lat)
im = plots.plotTb(monthly_mean, b, ax=ax[0], cmap=plt.get_cmap('temperature', 20))
ax[0].set_title('Brightness temperature')
ax[0].set_xlabel('Longitude')
ax[0].set_ylabel('Latitude')
ax[0].set_xticks(xticklocs)
ax[0].set_xticklabels(xticklabels)
ax[0].set_yticks(yticklocs)
ax[0].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[0])
im = plots.plotUTH(monthly_mean, b, ax=ax[1], cmap=plt.get_cmap('speed', 20), vmin=0, vmax=100)
ax[1].set_title('UTH')
ax[1].set_xlabel('Longitude')
ax[1].set_ylabel('Latitude')
ax[1].set_xticks(xticklocs)
ax[1].set_xticklabels(xticklabels)
ax[1].set_yticks(yticklocs)
ax[1].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[1])
im = plots.plotCDRQuantity(monthly_mean.uth_inhomogeneity[b], monthly_mean.lat, monthly_mean.lon, ax=ax[2], cmap=plt.get_cmap('speed', 20))
ax[2].set_title('UTH standard deviation of daily averages')
ax[2].set_xlabel('Longitude')
ax[2].set_ylabel('Latitude')
ax[2].set_xticks(xticklocs)
ax[2].set_xticklabels(xticklabels)
ax[2].set_yticks(yticklocs)
ax[2].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[2])
im = plots.plotCount(monthly_mean, node=b, ax=ax[3], cmap=plt.get_cmap('density', 20), vmin=0)
ax[3].set_title('Observation count')
ax[3].set_xlabel('Longitude')
ax[3].set_ylabel('Latitude')
ax[3].set_xticks(xticklocs)
ax[3].set_xticklabels(xticklabels)
ax[3].set_yticks(yticklocs)
ax[3].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[3])
im = plots.plotCDRQuantity(monthly_mean.observation_count_all[b] - monthly_mean.observation_count[b], monthly_mean.lat, monthly_mean.lon, ax=ax[4], cmap=plt.get_cmap('density', 20))
ax[4].set_title('count all - count ')
ax[4].set_xlabel('Longitude')
ax[4].set_ylabel('Latitude')
ax[4].set_xticks(xticklocs)
ax[4].set_xticklabels(xticklabels)
ax[4].set_yticks(yticklocs)
ax[4].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[4])
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

#%%
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

#%% plot small overview
b = 'ascending'
xticklocs = np.arange(-180, 181, 60)
xticklabels = [str(i) + '°' for i in xticklocs]
yticklocs = np.arange(-30, 31, 15)
yticklabels = [str(i) + '°' for i in yticklocs]

plt.rcParams.update({'font.size': 13}) 
fig, ax = plt.subplots(3, 1, figsize=(10, 10))
fig.suptitle(monthly_mean.time_coverage_start.strftime('%Y-%m'))
lon, lat = np.meshgrid(monthly_mean.lon, monthly_mean.lat)
im = plots.plotTb(monthly_mean, b, ax=ax[0], cmap=plt.get_cmap('temperature', 20))
ax[0].set_title('Brightness temperature')
ax[0].set_xlabel('Longitude')
ax[0].set_ylabel('Latitude')
ax[0].set_xticks(xticklocs)
ax[0].set_xticklabels(xticklabels)
ax[0].set_yticks(yticklocs)
ax[0].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[0])
im = plots.plotUTH(monthly_mean, b, ax=ax[1], cmap=plt.get_cmap('speed', 20), vmin=0, vmax=100)
ax[1].set_title('UTH')
ax[1].set_xlabel('Longitude')
ax[1].set_ylabel('Latitude')
ax[1].set_xticks(xticklocs)
ax[1].set_xticklabels(xticklabels)
ax[1].set_yticks(yticklocs)
ax[1].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[1])
im = plots.plotCount(monthly_mean, node=b, ax=ax[2], cmap=plt.get_cmap('density', 20), vmin=0)
ax[2].set_title('Observation count')
ax[2].set_xlabel('Longitude')
ax[2].set_ylabel('Latitude')
ax[2].set_xticks(xticklocs)
ax[2].set_xticklabels(xticklabels)
ax[2].set_yticks(yticklocs)
ax[2].set_yticklabels(yticklabels)
fig.colorbar(im, ax=ax[2])

plt.tight_layout(rect=[0, 0.03, 1, 0.95])