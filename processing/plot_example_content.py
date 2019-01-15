#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 09:26:58 2018

@author: Theresa Lang

Script to plot example content of UTH CDR
"""

import matplotlib.pyplot as plt
from processing import CDRs
from processing import plots
import numpy as np
#%% load NetCDF file 
example_CDR = CDRs.CDR.fromNetCDF('/scratch/uni/u237/users/tlang/CDR/CDR_UTH/CDR_files/fromFCDRv4_0/MHS/NOAA18/2012/FIDUCEO_CDR_UTH_MHS_NOAA18_20120701000000_20120731235959_L3_v1_0_fv2.0.0.nc')

# specify branch to plot
b = 'ascending'

#%%%%%%%%%%%%%%% UTH %%%%%%%%%%%%%%%%
fig, ax = plt.subplots(5, 1, figsize=(10, 10))

im = plots.plotCDRQuantity(example_CDR.uth[b], example_CDR.lat, example_CDR.lon, ax=ax[0], cmap=plt.get_cmap('speed', 14), vmin=0, vmax=100)
ax[0].set_title('uth_{}'.format(b))
cb = fig.colorbar(im, ax=ax[0])
cb.set_label('[%]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.uth_inhomogeneity[b], example_CDR.lat, example_CDR.lon, ax=ax[1], cmap=plt.get_cmap('speed', 9), vmin=0, vmax=45)
ax[1].set_title('uth_inhomogeneity_{}'.format(b))
cb = fig.colorbar(im, ax=ax[1])
cb.set_label('[%]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.u_uth['independent'][b], example_CDR.lat, example_CDR.lon, ax=ax[2], cmap=plt.get_cmap('Reds', 8), vmin=0, vmax=1.6)
ax[2].set_title('u_uth_independent_{}'.format(b))
cb = fig.colorbar(im, ax=ax[2])
cb.set_label('[%]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.u_uth['structured'][b], example_CDR.lat, example_CDR.lon, ax=ax[3], cmap=plt.get_cmap('Reds', 8), vmin=0, vmax=1.6)
ax[3].set_title('u_uth_structured_{}'.format(b))
cb = fig.colorbar(im, ax=ax[3])
cb.set_label('[%]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.u_uth['common'][b], example_CDR.lat, example_CDR.lon, ax=ax[4], cmap=plt.get_cmap('Reds', 8), vmin=0, vmax=1.6)
ax[4].set_title('u_uth_common_{}'.format(b))
cb = fig.colorbar(im, ax=ax[4])
cb.set_label('[%]', rotation=0, labelpad=20)

plt.tight_layout(pad=2.5)

#%%%%%%%%%%%%%%%% BT %%%%%%%%%%%%%%%%%%%
fig, ax = plt.subplots(5, 1, figsize=(10, 10))

im = plots.plotCDRQuantity(example_CDR.BT[b], example_CDR.lat, example_CDR.lon, ax=ax[0], cmap=plt.get_cmap('temperature', 14), vmin=240, vmax=275)
ax[0].set_title('BT_{}'.format(b))
cb = fig.colorbar(im, ax=ax[0])
cb.set_label('[K]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.BT_inhomogeneity[b], example_CDR.lat, example_CDR.lon, ax=ax[1], cmap=plt.get_cmap('temperature', 9), vmin=0, vmax=18)
ax[1].set_title('BT_inhomogeneity_{}'.format(b))
cb = fig.colorbar(im, ax=ax[1])
cb.set_label('[K]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.u_BT['independent'][b], example_CDR.lat, example_CDR.lon, ax=ax[2], cmap=plt.get_cmap('Reds', 8), vmin=0, vmax=0.3)
ax[2].set_title('u_BT_independent_{}'.format(b))
cb = fig.colorbar(im, ax=ax[2])
cb.set_label('[K]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.u_BT['structured'][b], example_CDR.lat, example_CDR.lon, ax=ax[3], cmap=plt.get_cmap('Reds', 6), vmin=0, vmax=0.3)
ax[3].set_title('u_BT_structured_{}'.format(b))
cb = fig.colorbar(im, ax=ax[3])
cb.set_label('[K]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.u_BT['common'][b], example_CDR.lat, example_CDR.lon, ax=ax[4], cmap=plt.get_cmap('Reds', 12), vmin=0, vmax=0.3)
ax[4].set_title('u_BT_common_{}'.format(b))
cb = fig.colorbar(im, ax=ax[4])
cb.set_label('[K]', rotation=0, labelpad=20)

plt.tight_layout(pad=2.5)

#%%%%%%%%%%%%%% Counts %%%%%%%%%%%%%%%%%%%%

fig, ax = plt.subplots(4, 1, figsize=(10, 10))

im = plots.plotCDRQuantity(example_CDR.overpass_count[b], example_CDR.lat, example_CDR.lon, ax=ax[0], cmap=plt.get_cmap('density', 6), vmin=0, vmax=6)
ax[0].set_title('overpass_count_{}'.format(b))
cb = fig.colorbar(im, ax=ax[0])
cb.set_label('[-]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.observation_count[b], example_CDR.lat, example_CDR.lon, ax=ax[1], cmap=plt.get_cmap('density', 8), vmin=0, vmax=160)
ax[1].set_title('observation_count_{}'.format(b))
cb = fig.colorbar(im, ax=ax[1])
cb.set_label('[-]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.observation_count_all[b], example_CDR.lat, example_CDR.lon, ax=ax[2], cmap=plt.get_cmap('density', 8), vmin=0, vmax=160)
ax[2].set_title('observation_count_all_{}'.format(b))
cb = fig.colorbar(im, ax=ax[2])
cb.set_label('[-]', rotation=0, labelpad=20)

im = plots.plotCDRQuantity(example_CDR.observation_count_all[b] - example_CDR.observation_count[b], example_CDR.lat, example_CDR.lon, ax=ax[3], cmap=plt.get_cmap('density', 9), vmin=0, vmax=90)
ax[3].set_title('observation_count_all_{} - observation_count_{}'.format(b, b))
cb = fig.colorbar(im, ax=ax[3])
cb.set_label('[-]', rotation=0, labelpad=20)


plt.tight_layout(pad=2.5)

#%%%%%%%%%%%%%% Time Ranges %%%%%%%%%%%%%%%%

fig, ax = plt.subplots(2, 1, figsize=(10, 5))

im = plots.plotCDRQuantity(example_CDR.time_ranges['ascending'][0] / 60 / 60, example_CDR.lat, example_CDR.lon, ax=ax[0], cmap=plt.get_cmap('density', 24))
ax[0].set_title('time_ranges_ascending min')

im = plots.plotCDRQuantity(example_CDR.time_ranges['descending'][0] / 60 / 60, example_CDR.lat, example_CDR.lon, ax=ax[1], cmap=plt.get_cmap('density', 24))
ax[1].set_title('time_ranges_descending min')

#fig.subplots_adjust(right=0.6)
cbar_ax = fig.add_axes([0.11, 0.08, 0.88, 0.03])
cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cb.set_label('[hours of day]', labelpad=-3)


plt.tight_layout(pad=4.5, w_pad=0.3)