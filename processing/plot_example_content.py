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
example_CDR = CDRs.CDR.fromNetCDF('/scratch/uni/u237/users/tlang/CDR/CDR_UTH/CDR_files/fromFCDRv4_1/v1_2/AMSUB/NOAA17/2005/FIDUCEO_CDR_UTH_AMSUB_NOAA17_20050701000000_20050731235959_L3_v1_2_fv2.0.0.nc')
#('/scratch/uni/u237/users/tlang/CDR/CDR_UTH/CDR_files/fromFCDRv4_1/SSMT2/F15/2004/FIDUCEO_CDR_UTH_SSMT2_F15_20041101000000_20041130235959_L3_v1_1_fv2.0.0.nc')

# specify branch to plot
b = 'ascending'

#%%%%%%%%%%%%%%% UTH %%%%%%%%%%%%%%%%
plt.rcParams.update({'font.size': 15})
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

#%%%%%%%% UTH Ascending + Descending %%%%%%%%%%%%
######### Combined uncertainty #############
fig, ax = plt.subplots(2, 1, figsize=(15, 5))

u_ind = np.sqrt(example_CDR.u_uth['independent']['ascending'] ** 2 + example_CDR.u_uth['independent']['descending'] ** 2)
u_struct = (example_CDR.u_uth['structured']['ascending'] + example_CDR.u_uth['structured']['descending']) / 2
u_com = (example_CDR.u_uth['common']['ascending'] + example_CDR.u_uth['common']['descending']) / 2
u_tot = np.sqrt(u_ind ** 2 + u_struct ** 2 + u_com ** 2)

im = plots.plotCDRQuantity(
        np.mean(np.dstack((example_CDR.uth['ascending'], example_CDR.uth['descending'])), axis=2),
        example_CDR.lat, example_CDR.lon, ax=ax[0], cmap=plt.get_cmap('speed', 10), vmin=0, vmax=100)
cb = fig.colorbar(im, ax=ax[0])
cb.set_label('UTH [%]', labelpad=20)


im = plots.plotCDRQuantity(u_tot, example_CDR.lat, example_CDR.lon, ax=ax[1], cmap=plt.get_cmap('Reds', 10), vmin=0, vmax=2.5)
cb = fig.colorbar(im, ax=ax[1])
cb.set_label('UTH uncertainty [%]', labelpad=20)


#%%%%%%%% BT Ascending + Descending %%%%%%%%%%%%
######### Combined uncertainty #############
fig, ax = plt.subplots(2, 1, figsize=(15, 5))

u_ind = np.sqrt(example_CDR.u_BT['independent']['ascending'] ** 2 + example_CDR.u_BT['independent']['descending'] ** 2)
u_struct = (example_CDR.u_BT['structured']['ascending'] + example_CDR.u_BT['structured']['descending']) / 2
u_com = (example_CDR.u_BT['common']['ascending'] + example_CDR.u_BT['common']['descending']) / 2
u_tot = np.sqrt(u_ind ** 2 + u_struct ** 2 + u_com ** 2)

im = plots.plotCDRQuantity(np.mean(np.dstack((example_CDR.BT['ascending'], example_CDR.BT['descending'])), axis=2),
                           example_CDR.lat, example_CDR.lon, ax=ax[0], cmap=plt.get_cmap('temperature', 14), vmin=240, vmax=275)
cb = fig.colorbar(im, ax=ax[0])
cb.set_label('BT [K]', labelpad=20)

im = plots.plotCDRQuantity(u_tot, example_CDR.lat, example_CDR.lon, ax=ax[1], cmap=plt.get_cmap('Reds', 10), vmin=0.3, vmax=0.55)
cb = fig.colorbar(im, ax=ax[1])
cb.set_label('BT uncertainty [K]', labelpad=20)