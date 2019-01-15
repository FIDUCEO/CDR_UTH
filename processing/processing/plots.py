#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 10:19:45 2018

@author: u300740
"""
import numpy as np
import matplotlib.pyplot as plt
import typhon
from mpl_toolkits.basemap import Basemap

def plotUTH(CDR, node, ax=None, **kwargs):
    """ plot UTH field of a CDR object.
    parameters:
        CDR: CDR-object
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    
    if ax is None:
        ax = plt.gca()
    
    default_kwargs = {}    
    if 'cmap' not in kwargs:
        default_kwargs = {
                'cmap': plt.get_cmap('speed', 20)}
    
    lon, lat = np.meshgrid(CDR.lon, CDR.lat)
        
    return ax.pcolormesh(lon, lat, CDR.uth[node], **default_kwargs, **kwargs)

def plotTb(CDR, node, ax=None, **kwargs):
    """ plot brightness temperature field of a CDR object.
    parameters:
        CDR: CDR-object
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    
    if ax is None:
        ax = plt.gca()
    
    default_kwargs = {}    
    if 'cmap' not in kwargs:
        default_kwargs = {
                'cmap': plt.get_cmap('temperature', 20)}
    
    lon, lat = np.meshgrid(CDR.lon, CDR.lat)
        
    return ax.pcolormesh(lon, lat, CDR.BT[node], **default_kwargs, **kwargs)

def plotCount(CDR, node, ax=None, **kwargs):
    """ plot count field of a CDR object.
    parameters:
        CDR: CDR-object
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    
    if ax is None:
        ax = plt.gca()
    
    default_kwargs = {}    
    if 'cmap' not in kwargs:
        default_kwargs = {
                'cmap': plt.get_cmap('density', 20)}
    
    lon, lat = np.meshgrid(CDR.lon, CDR.lat)
        
    return ax.pcolormesh(lon, lat, CDR.observation_count[node], **default_kwargs, **kwargs)

def plotCountAll(CDR, node, ax=None, **kwargs):
    """ plot count field of a CDR object.
    parameters:
        CDR: CDR-object
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    
    if ax is None:
        ax = plt.gca()
    
    default_kwargs = {}    
    if 'cmap' not in kwargs:
        default_kwargs = {
                'cmap': plt.get_cmap('density', 20)}
    
    lon, lat = np.meshgrid(CDR.lon, CDR.lat)
        
    return ax.pcolormesh(lon, lat, CDR.observation_count_all[node], **default_kwargs, **kwargs)

def plotOverpassCount(CDR, node, ax=None, **kwargs):
    """ plot overpass count field of a CDR object.
    parameters:
        CDR: CDR-object
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    
    if ax is None:
        ax = plt.gca()
    
    default_kwargs = {}    
    if 'cmap' not in kwargs:
        default_kwargs = {
                'cmap': plt.get_cmap('density', 20)}
    
    lon, lat = np.meshgrid(CDR.lon, CDR.lat)
        
    return ax.pcolormesh(lon, lat, CDR.overpass_count[node], **default_kwargs, **kwargs)

def plotTbUncertainty(CDR, uncertainty_type, node, ax=None, **kwargs):
    """ plot brightness temperature uncertainty field of a CDR object.
    parameters:
        CDR: CDR-object
        uncertainty_type (string): type of uncertainty ('independent', 'common' or 'structured')
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    if ax is None:
        ax = plt.gca()
    
    default_kwargs = {}    
    if 'cmap' not in kwargs:
        default_kwargs = {
                'cmap': 'Reds'}
    
    lon, lat = np.meshgrid(CDR.lon, CDR.lat)
        
    return ax.pcolormesh(lon, lat, CDR.u_BT[uncertainty_type][node], **default_kwargs, **kwargs)

def plotUTHUncertainty(CDR, uncertainty_type, node, ax=None, **kwargs):
    """ plot UTH uncertainty field of a CDR object.
    parameters:
        CDR: CDR-object
        uncertainty_type (string): type of uncertainty ('independent', 'common' or 'structured')
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    if ax is None:
        ax = plt.gca()
    
    default_kwargs = {}    
    if 'cmap' not in kwargs:
        default_kwargs = {
                'cmap': 'Reds'}
    
    lon, lat = np.meshgrid(CDR.lon, CDR.lat)
        
    return ax.pcolormesh(lon, lat, CDR.u_uth[uncertainty_type][node], **default_kwargs, **kwargs)

def plotCDRQuantity(CDR_quantity, latitudes, longitudes, ax=None, **kwargs):
    """ plot any given quantity of a CDR object.
    parameters:
        CDR_quantity (2DArray): Any 2D field that should be plotted
        latitudes (Array): Latitudes of CDR grid
        longitudes (Array): Longitudes of CDR grid
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    if ax is None:
        ax = plt.gca()
    
    lon, lat = np.meshgrid(longitudes, latitudes)
    ma = Basemap(projection='mill',lon_0=0, llcrnrlat=-30, urcrnrlat=30,
            llcrnrlon=-180, urcrnrlon=180, ax=ax)
    # plot coastlines, draw label meridians and parallels.
    ma.drawcoastlines()
    labels_meridians=[0,0,0,1]
    labels_parallels = [1,0,0,0]    
    ma.drawparallels(np.arange(ma.latmin, ma.latmax, 15), labels=labels_parallels, fontsize=11)
    ma.drawmeridians(np.arange(ma.lonmin,ma.lonmax,60),labels=labels_meridians, fontsize=11)
        
    return ma.pcolormesh(lon, lat, CDR_quantity, latlon=True, **kwargs)

def plotOverview(CDR, node, **kwargs):
    """ Plots an overview of the CDR-object containing brightness temperature,
    UTH and count for the specified node.
    parameters:
        CDR: CDR-object
        node (string): 'ascending' or 'descending' (satellite branch)
        ax(AxesSubplot, optional): Axes to plot in
        **kwargs: Additional keyword arguments
    """
    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    im = plotTb(CDR, node, ax=ax[0])
    cb = fig.colorbar(im, ax=ax[0])
    im = plotUTH(CDR, node, ax=ax[1])
    cb = fig.colorbar(im, ax=ax[1])
    im = plotCount(CDR, node, ax=ax[2])
    cb = fig.colorbar(im, ax=ax[2])
    for i in range(3):
        ax[i].set_xlabel('Longitudes (°)')
        ax[i].set_ylabel('Latitudes (°)')
    fig.suptitle(CDR.time_coverage_start.strftime('%Y-%m'))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    
    
    
    
    



    
        
    