#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 12:03:00 2018

@author: u300740
"""

import numpy as np
from os.path import join
from datetime import timedelta, date
import matplotlib.pyplot as plt
import glob
from processing.CDRs import CDR
from processing.FCDRs import FCDR

#fcdr_path = '/scratch/uni/u237/user_data/ihans/FCDR/easy/v2_0fv1_1_4/'
fcdr_path = '/scratch/uni/u237/user_data/tlang/FCDRs/HIRS/'

# time period, instrument, satellite and viewing angles
start_date = date(1992, 1, 1)
end_date = date(1992, 1, 31)
instrument = 'HIRS' #'HIRS'
satellite = 'noaa11' #'METOPA'
viewing_angles = {}
viewing_angles['AMSUB'] = [i for i in range(33, 58)]
viewing_angles['MHS'] = [i for i in range(33, 58)]
viewing_angles['HIRS'] = [i for i in range(17, 41)]
lat_boundaries=[-30, 30] #[10, 24]
lon_boundaries=[100, 160] #[0, 14]

# resolution of lat/lon grid:
resolution = 1.0

filenames = []
fcdrs = []
gridded_days = []
dates = []
d = start_date
delta = timedelta(days=1)

while d <= end_date:
    print(d.strftime("%Y-%m-%d"))
    path = join(fcdr_path, satellite.lower(), str(d.year), '{m:02}'.format(m=d.month), '{d:02}'.format(d=d.day))
    filenames_new = glob.glob(join(path, '*0.7.nc'))
    fcdrs = []
    for file in filenames_new:
        #print(file)
        f = FCDR.fromNetcdf(
                 fcdr_path, file, viewing_angles[instrument])
        f.generateTotalMask()
        f.generateCloudMask()
        fcdrs.append(f)

    # gridded daily CDRs
    gridded_days.append(CDR.GriddedCDRFromFCDRs_fullPropagation(fcdrs, lat_boundaries=lat_boundaries, 
                                                lon_boundaries=lon_boundaries, resolution=resolution,
                                                uncertainties=True))
    
    d += delta


 #%%
#for i, c in enumerate(gridded_days):
#    for b in ['ascending', 'descending']:
#        if np.any(~np.isnan(c.brightness_temp_mean[b])):
#            print(i)
#            print(b)
#            print('Tb:')
#            print(c.brightness_temp_mean[b])
#            print('count:')
#            print(c.count[b])
#            print('S_struct:')
#            print(c.S_structured[b].todense())
#            print('S_common:')
#            print(c.S_common[b].todense())

#%%
print('monthly mean')
monthly_mean = CDR.AveragedCDRFromCDRs_fullPropagation(gridded_days)

#%% get correlation structure for center grid cell
numcells = len(monthly_mean.lat_grid) * len(monthly_mean.lon_grid)
center = int(np.floor(numcells/2))
branches = ['ascending', 'descending']
center_corrs = {}
for b in branches:
    center_corrs[b] = np.flipud(monthly_mean.q_struct[b][center, :].ravel().reshape((len(monthly_mean.lat_grid), len(monthly_mean.lon_grid))))
    #corr_grid[b] = np.zeros((len(monthly_mean.lon_grid), len(monthly_mean.lat_grid)))

e_folding = 1 / np.e
plt.rcParams.update({'font.size': 16}) 
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
for i, b in enumerate(branches):
    ax[i].set_title(b)
    im = ax[i].imshow(center_corrs[b], cmap=plt.get_cmap('density', 10), vmin=0, vmax=1)
    #ax[i].contour(center_corrs[b], levels=[e_folding], colors='r')
    ax[i].set_yticks([])
    ax[i].set_xticks([])
    ax[i].set_xticks(np.arange(-.5, center_corrs[b].shape[0], 1), minor=True)
    ax[i].set_yticks(np.arange(-.5, center_corrs[b].shape[0], 1), minor=True)
    ax[i].grid(which='minor', color='k', linestyle='-', linewidth=0.5)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label('correlation')

    
#%%
    
#data = gridded_orbits[1]
#lon, lat = np.meshgrid(data.lon_grid, data.lat_grid)
#
#Tb = monthly_mean.brightness_temp_mean
#fig, ax = plt.subplots()
#ax.pcolormesh(lon, lat, Tb['ascending'])
figurename = 'S_struct_{i}_{s}_{m}_{a}angles_lon_{lon1}_{lon2}_lat_{lat1}_{lat2}.png'.format(
        i=instrument, s=satellite, m=start_date.strftime('%Y%m'), a=len(viewing_angles),\
        lon1=lon_boundaries[0], lon2=lon_boundaries[-1], lat1=lat_boundaries[0], lat2=lat_boundaries[-1])
gridnames = [(lat, lon) for lat in monthly_mean.lat_grid for lon in monthly_mean.lon_grid] 
title = 'gridcells: {}\ninstrument: {}, satellite: {}, month: {},\nnumber of viewing angles: {}'.format(gridnames, instrument, satellite, start_date.strftime('%m/%Y'), len(viewing_angles))

plt.rcParams.update({'font.size': 13})
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
fig.suptitle(title)
im = ax[0].imshow(monthly_mean.q_struct['ascending'], cmap=plt.get_cmap('density', 10), vmin=0, vmax=1)
ax[0].set_xticks([])
#ax[0].set_yticks([])
ax[0].set_xticks(np.arange(-.5, monthly_mean.q_struct['ascending'].shape[0], 1), minor=True);
ax[0].set_yticks(np.arange(-.5, monthly_mean.q_struct['ascending'].shape[0], 1), minor=True);
#ax[0].grid(which='minor', color='k', linestyle='-', linewidth=0.5)
ax[0].set_title('S$_{struct}$ ascending')

im = ax[1].imshow(monthly_mean.q_struct['descending'], cmap=plt.get_cmap('density', 10), vmin=0, vmax=1)
ax[1].set_title('S$_{struct}$ descending')
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[1].set_xticks(np.arange(-.5, monthly_mean.q_struct['ascending'].shape[0], 1), minor=True);
ax[1].set_yticks(np.arange(-.5, monthly_mean.q_struct['ascending'].shape[0], 1), minor=True);
#ax[1].grid(which='minor', color='k', linestyle='-', linewidth=0.5)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
#plt.savefig(join('plots/cov_mat_examples', figurename))


fig, ax = plt.subplots(1, 2, figsize=(8, 6))
fig.suptitle(title)
im = ax[0].imshow(monthly_mean.S_com['ascending'].todense(), vmin=0, vmax=0.5)
ax[0].set_title('S$_{com}$ ascending')
ax[0].set_axis_off()
im = ax[1].imshow(monthly_mean.S_com['descending'].todense(), vmin=0, vmax=0.5)
ax[1].set_title('S$_{com}$ descending')
ax[1].set_axis_off()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

#%% plot correlation vectors for HIRS and microwave sensors:
plt.rcParams.update({'font.size': 15})
corr_vector_mw = [1.0, 0.8465, 0.5134, 0.2231, 0.0695, 0.0155, 0.0025]
x_mw = np.arange(7)

#corr_vector_hirs =
fig, ax = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
ax[0].plot(x_mw, corr_vector_mw, 'x', markersize=5, markeredgewidth=2)
ax[0].set_xlabel('Scanline difference')
ax[0].set_ylabel('Correlation')
ax[0].set_ylim(-0.02, 1.02)
ax[0].set_title('MW')

if instrument == 'HIRS':
    corr_vector_hirs = f.corr_vector[:,12]
    x_hirs = np.arange(len(corr_vector_hirs))
    ax[1].plot(x_hirs, corr_vector_hirs, 'x', markersize=1, markeredgewidth=1)
    ax[1].set_xlim(0, 230)
    ax[1].set_xlabel('Scanline difference')
    ax[1].set_title('HIRS')

plt.tight_layout()



   
