#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:19:31 2020

@author: mallen
"""


import rasterio as rio
import rasterio.mask
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors
import os
from affine import Affine
from datetime import datetime

# =============================================================================
# bring in and prep data
# =============================================================================
# start time
print("start: " + str(datetime.now().time()))

### parameters
# parameters are read in from a .csv in the parameter_files directory and populated below
parameters = pd.read_csv("E:/ssse_ecostress/parameter_files/nyc_view_seen_01142021.csv")

# model settings
city = parameters[parameters['param'] == 'city']['setting'].item()
svs = parameters[parameters['param'] == 'svs']['setting'].item()
sns = parameters[parameters['param'] == 'sns']['setting'].item()
buildings_only = parameters[parameters['param'] == 'buildings_only']['setting'].item()
max_height = parameters[parameters['param'] == 'max_height']['setting'].item()
plot = parameters[parameters['param'] == 'plot']['setting'].item()

# image parameters
px = int(parameters[parameters['param'] == 'px']['setting'])
pixelsize = int(parameters[parameters['param'] == 'pixelsize']['setting'])
thresh = int(parameters[parameters['param'] == 'thresh']['setting'])

# filepaths/working directory
cover = parameters[parameters['param'] == 'cover']['setting'].item()
dsm = parameters[parameters['param'] == 'dsm']['setting'].item()
metadata = parameters[parameters['param'] == 'metadata']['setting'].item()
basepath = parameters[parameters['param'] == 'basepath']['setting'].item()

# set wd
os.chdir(basepath)

# set filename/plotting parameters
if svs == 'view' and sns == "seen":
    tag0 = 'untagged'
    tag1 = 'seen wall'
    fnappend = '_sensor'
elif svs == 'view' and sns == 'not seen':
    tag0 = 'untagged'
    tag1 = 'obscured'
    fnappend = '_sensor_notseen'
else:
    tag0 = 'sunlit'
    tag1 = 'shaded'
    fnappend = '_shade'

# filter for buildings (if buildings_only == yes)
if buildings_only == "yes":
    lc = rio.open("../data/" + cover)
    lcr = lc.read()[0,:,:]
    building_id = int(parameters[parameters['param'] == 'buildings_ID']['setting'])
    lcr[lcr != building_id] = 0
    lcr[lcr == building_id] = 1

# bring in the metadata
meta = pd.read_csv(metadata)

### note: the azimuth/solar zenith angles in ECOSTRESS metadata files are N=0 and horizon = 90
# we adjust them below so the math is easier to interpret
# adjust azimuth angles
adj = []
for i in range(len(meta)):
    if meta['solar_azimuth'][i] > 180:
        adj.append(meta['solar_azimuth'][i] - 180)
    else:
        adj.append(meta['solar_azimuth'][i] + 180)
meta['solar_azimuth_adj'] = pd.Series(adj)
adj = []
for i in range(len(meta)):
    if sns == 'seen':
        if meta['view_azimuth'][i] > 180:
            adj.append(meta['view_azimuth'][i])
        else:
            adj.append(meta['view_azimuth'][i])
    else:
        if meta['view_azimuth'][i] > 180:
            adj.append(meta['view_azimuth'][i] - 180)
        else:
            adj.append(meta['view_azimuth'][i] + 180)
meta['view_azimuth_adj'] = pd.Series(adj)

# adjust solar zenith
meta['solar_zenith_adj'] = 90 - meta['solar_zenith'] 

# remove nighttime for solar
if svs == 'solar':
    meta = meta[meta['solar_zenith_adj'] > 10].reset_index(drop = True)

# =============================================================================
# enter main loop
# =============================================================================
### meta loop
for f in range(len(meta)):
    # bring in dsm (need fresh copy every time)
    la = rio.open("../data/" + dsm) # changed to 5m version
    lar = la.read()[0,:,:]
    lameta = la.meta
    latransform = la.transform
    
    # filter
    lar = np.where(np.isnan(lar), 0, lar)
    lar[lar == -999] = 0
    lar[lar < 0] = 0
    lar[lar > 336] = 0 # tallest building in LA is 3pxm tall
    # filter for building
    larm = lar * lcr # if we only want buildings
    #larm = lar.copy() # if we want shade from all objects
    
    # get metadata for output
    metaout = la.meta.copy()
    metaout.update({'dtype': 'float64'})
    
    ### sensor/sun parameters
    va = rio.open("../data/" + city + "geo/" + "ECO1BGEO.001_Geolocation_" + svs + "_azimuth_" + meta['filename'][f][12:] + "_utm_clip.tif").read()[0,:,:]
    vz = rio.open("../data/" + city + "geo/" + "ECO1BGEO.001_Geolocation_" + svs + "_zenith_" + meta['filename'][f][12:] + "_utm_clip.tif").read()[0,:,:]
    # mask 
    va[va < -1000] = np.nan
    vz[vz < -1000] = np.nan
    
    # adjust azimuth
    if svs == 'solar':
        if meta[svs + "_azimuth"][f] > 180:
            va = va - 180
        else:
            va = va + 180
    elif svs == 'view' and sns == 'notseen':
        if meta[svs + "_azimuth"][f] > 180:
            va = va - 180
        else:
            va = va + 180
    
    # set filename
    outfn = meta['filename'][f] + fnappend
    print("file: " + outfn)
    # convert zenith to elevation angle
    vz_e = 90 - vz
    
    # initialize matrix
    shade_matrix = np.zeros([lar.shape[0], lar.shape[1]])
    
    # 1) calculate direction of shade
    # find quadrant
    va_quad = int(np.ceil(meta[svs+"_azimuth_adj"][f]/90))
    # find proportion v/h in quadrant
    v = -(meta[svs+"_azimuth_adj"][f] - 90*va_quad)/90 # flip sign so positive when moving up 
    h = 1-v # find proportion in opposite direction (sum to 1)
    # change signs and flip arrays to reflect quadrants
    if va_quad == 1:
        shade_matrix = np.fliplr(shade_matrix)
        lar = np.fliplr(lar)
        va = np.fliplr(va)
        vz_e = np.fliplr(vz_e)
        v1 = v
        h1 = h
    elif va_quad == 2: 
        shade_matrix = np.flipud(np.fliplr(shade_matrix))
        lar = np.flipud(np.fliplr(lar))
        va = np.flipud(np.fliplr(va))
        vz_e = np.flipud(np.fliplr(vz_e))
        v1 = h
        h1 = v
    elif va_quad == 3:
        shade_matrix = np.flipud(shade_matrix)
        lar = np.flipud(lar)
        va = np.flipud(va)
        vz_e = np.flipud(vz_e)
        v1 = v
        h1 = h
    elif va_quad == 4:
        shade_matrix = shade_matrix
        lar = lar
        v1 = h
        h1 = v
         
    ### enter loop
    vat = []
    vzt = []
    print("enter loop: " + str(datetime.now().time()))
    for i in range(10, lar.shape[0]-10):
        for j in range(10, lar.shape[1]-10):
            z = lar[i,j] # get the height of a pixel
            # calculate shade for surrounding pixels as a function of
            # - va, vz, z
            # if height is > 1m, keep going
            if z < thresh:
                continue
            else:
                # update va and vz
                vaij = np.nanmean(va[int(np.floor(i/px)):int(np.ceil(i/px)), int(np.floor(j/px)):int(np.ceil(j/px))])
                vz_eij = np.nanmean(vz_e[int(np.floor(i/px)):int(np.ceil(i/px)), int(np.floor(j/px)):int(np.ceil(j/px))])
                
                # fill in the edges with average values
                if np.isnan(vaij) == True:
                    vaij = meta[svs+"_azimuth_adj"][f]
                    if svs == 'view':
                        vz_eij = 90 - meta[svs+"_zenith"][f]
                    else:
                        vz_eij = meta[svs+"_zenith_adj"][f]
                    #print("edge?")
                    
                #vat.append(vaij)
                #vzt.append(vz_eij)
                # 2) calculate length of shade
                l = z/np.tan(vz_eij/360*(2*np.pi))
                # divide by pixelsize to get length in pixels,
                # this represents the MAX shade length in pixels (i.e. the length for pixels
                # parallel to the sun/sensor LOS/Azimuth)
                lpx = l/pixelsize
                
                # 3) use v and h to calculate shade plane extending from the origin pixel
                v = v1 # update with correct orientation
                h = h1
                lpx_v = lpx*v
                lpx_h = lpx*h
                
                # set up sub-matricies
                sm_shade = np.zeros([int(np.ceil(lpx_v)), int(np.ceil(lpx_h))])
                sm_lar = lar[i:i+int(np.ceil(lpx_v)), j:j+int(np.ceil(lpx_h))].copy()
                sm_dist = sm_shade.copy()
                sm_linefinder = sm_shade.copy()
                
                for ii in range(sm_shade.shape[0]):
                    for jj in range(sm_shade.shape[1]):
                        # distance for each ring around the pixel
                        sm_dist[ii, jj] = max([ii, jj])
                        # compare ratio of v/h to ii/jj, this finds the correct line
                        # line represents a comparion between points that are parallel with the sun
                        # LOS and the pixel grid
                        # this part is fuzzy and likely needs improvement
                        sm_linefinder[ii, jj] = np.absolute(v/h - (ii+1)/(jj+1)) # add one to get away from zero, this is fine
                
                # calculate shade plane height, subtract building height
                # if negative, shade plane is attenuated
                sm_shade = (l-(sm_dist*pixelsize))*np.tan(vz_eij/360*(2*np.pi)) - sm_lar
                
                # 4) filter shadeplane matrix using optimal shade "path"
                # find the optimal "path" for shade
                # copy the matricies of interest
                sm_linefound = np.zeros([int(np.ceil(lpx_v)), int(np.ceil(lpx_h))])
                # loop to tag the minimum value (i.e. the optimal path) in each "ring"
                for w in range(0, int(np.max(sm_dist)+1)):
                    # reup linefinder
                    sm_linefinderw = sm_linefinder.copy()
                    sm_distw = sm_dist.copy()
                    # define target ring
                    sm_distw[sm_distw != w] = np.nan
                    sm_distw[sm_distw == w] = 1
                    # filter for ring
                    sm_linefinderw = sm_linefinderw * sm_distw
                    # find lowest value
                    wi, wj = np.where(sm_linefinderw == np.nanmin(sm_linefinderw))
                    sm_linefound[wi, wj] = 1
                # multiply shade plane vs optimal line
                sm_shade_line = sm_shade * sm_linefound
                
                # 5) attenuate the optimized shade plane
                # convert to binary
                sm_shade_line[sm_shade_line > 0] = 1 
                sm_shade_line[sm_shade_line <= 0] = 0 
                # find first location where shade plane is negative
                vi, vj = np.where(sm_shade_line < 0)
                if (len(vi) == 0) and (len(vj) == 0):
                    sm_shade_line = sm_shade_line
                else:
                    sm_shade_line = sm_shade_line
                
                # 6) apply the shadeline to the shade_matrix
                shade_matrix[i:i+sm_shade_line.shape[0],j:j+sm_shade_line.shape[1]] = shade_matrix[i:i+sm_shade_line.shape[0],j:j+sm_shade_line.shape[1]] + sm_shade_line
    
    # 7) flip back around
    if va_quad == 1:
        shade_matrix_filt = np.fliplr(shade_matrix.copy())
        larf = np.fliplr(lar.copy())
    elif va_quad == 2:
        shade_matrix_filt = np.flipud(np.fliplr(shade_matrix.copy()))
        larf = np.flipud(np.fliplr(lar.copy()))
    elif va_quad == 3:
        shade_matrix_filt = np.flipud(shade_matrix.copy())
        larf = np.flipud(lar.copy())
    elif va_quad == 4:
        shade_matrix_filt = shade_matrix.copy()
        larf = lar.copy()
    
    shade_matrix_filt[shade_matrix_filt > 1] = 1
    print("exit loop: " + str(datetime.now().time()))
    
# =============================================================================
# plot
# =============================================================================
    if plot == 'no':
        continue
    else:
        # create colormap
        colors = ['0.9','dodgerblue']
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors, 2)
        
        # plot
        fig, [ax1, ax0] = plt.subplots(2, 1)
        
        p0 = ax0.imshow(shade_matrix_filt[5660:5780, 5580:5740], vmin = 0, vmax = 1, cmap = cmap)
        cb0 = plt.colorbar(p0, ax = ax0, ticks = [0.25, 0.75])
        #cb0.ax.set_yticklabels(["Sunlit", "Shaded"]) # for sun
        cb0.ax.set_yticklabels([tag0, tag1])
        p1 = ax1.imshow(larf[5660:5780, 5580:5740], vmin = 0, vmax = 45) ##########
        ax1.text(0.03, 
             1.03, 
             "elev: " + str(round(vz_eij, 1)) + ", azim: " + str(round(vaij, 1)) + ", time: " + str(meta['hourpst'][f]) + ":" + str(meta['minute'][f]),
             va = "bottom",
             ha = "left", 
             c = '0',
             transform = ax1.transAxes)
        plt.colorbar(p1, ax = ax1, label = "Height Above Ground, m")
        
        plt.tight_layout()
        plt.savefig("../plots/" + outfn + ".png", dpi = 500)
    
# =============================================================================
# output
# =============================================================================
    with rio.open("../data/" + outfn + ".tif", 'w', **metaout) as dst:
        dst.write(shade_matrix_filt[None,:,:])