#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:19:31 2020

@author: mallen
"""


import rasterio as rio  
import rasterio.mask
import rasterio.crs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime

# =============================================================================
# bring in and prep data
# =============================================================================
# start time
print("start: " + str(datetime.now().time()))

### parameters
# parameters are read in from a .csv in the parameter_files directory and populated below
parametersfn = input("input path to parameter file, e.g. /path/to/file.csv: \n")
parameters = pd.read_csv(parametersfn)

# model settings
city = parameters[parameters['param'] == 'city']['setting'].item()
svs = parameters[parameters['param'] == 'svs']['setting'].item()
sns = parameters[parameters['param'] == 'sns']['setting'].item()
buildings_only = parameters[parameters['param'] == 'buildings_only']['setting'].item()
max_height = int(parameters[parameters['param'] == 'max_height']['setting'].item())
plot = parameters[parameters['param'] == 'plot']['setting'].item()
perpixel_angles = parameters[parameters['param'] == 'perpixel_angles']['setting'].item()

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
    lc = rio.open(cover)
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
    la = rio.open(dsm) # changed to 5m version
    lar = la.read()[0,:,:]
    lameta = la.meta
    latransform = la.transform
    
    # filter
    lar = np.where(np.isnan(lar), 0, lar)
    lar[lar == -999] = 0
    lar[lar < 0] = 0
    lar[lar > max_height] = max_height
    # filter for building
    if buildings_only == "yes":
        lar = lar * lcr # if we only want buildings
    else:
        lar = lar.copy() # if we want shade from all objects
    
    # get metadata for output
    metaout = la.meta.copy()
    metaout.update({'dtype': 'float64'})
    
    ### sensor/sun parameters
    if perpixel_angles == 'yes':
        va = rio.open("../data/" + city + "geo/" + "ECO1BGEO.001_Geolocation_" + svs + "_azimuth_" + meta['filename'][f][12:] + "_utm_clip.tif").read()[0,:,:]
        vz = rio.open("../data/" + city + "geo/" + "ECO1BGEO.001_Geolocation_" + svs + "_zenith_" + meta['filename'][f][12:] + "_utm_clip.tif").read()[0,:,:]
        # mask 
        va[va < -1000] = np.nan
        vz[vz < -1000] = np.nan
    else:
        if svs == 'view':
            vz = meta['view_zenith'][f]
            va = meta['view_azimuth_adj'][f]
        elif svs == 'solar':
            vz = meta['solar_zenith'][f] # note, not using adjusted because it makes it easier
            va = meta['solar_azimuth_adj'][f]
    
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
    # convert zenith to elevation angle (note, because we don't use adjusted, both are actually just the zenith angle)
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
        if perpixel_angles == 'yes':
            va = np.fliplr(va)
            vz_e = np.fliplr(vz_e)
        v1 = v
        h1 = h
    elif va_quad == 2: 
        shade_matrix = np.flipud(np.fliplr(shade_matrix))
        lar = np.flipud(np.fliplr(lar))
        if perpixel_angles == 'yes':
            va = np.flipud(np.fliplr(va))
            vz_e = np.flipud(np.fliplr(vz_e))
        v1 = h
        h1 = v
    elif va_quad == 3:
        shade_matrix = np.flipud(shade_matrix)
        lar = np.flipud(lar)
        if perpixel_angles == 'yes':
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
    print("enter loop: " + str(datetime.now().time()))
    for i in range(20, lar.shape[0]-20):
        for j in range(20, lar.shape[1]-20):
            z = lar[i,j] # get the height of a pixel
            # calculate shade for surrounding pixels as a function of
            # - va, vz, z
            # if height is > threshold height (default 2m), keep going
            if z < thresh:
                continue
            else:
                # update va and vz if perpixel angles is yes
                if perpixel_angles == 'yes':
                    vaij = np.nanmedian(va[int(np.floor(i/px)):int(np.ceil(i/px)), int(np.floor(j/px)):int(np.ceil(j/px))]) # changed to median feb 19, 2021
                    vz_eij = np.nanmedian(vz_e[int(np.floor(i/px)):int(np.ceil(i/px)), int(np.floor(j/px)):int(np.ceil(j/px))])
                    
                    # fill in the edges with average values
                    if np.isnan(vaij) == True:
                        vaij = meta[svs+"_azimuth_adj"][f]
                        if svs == 'view':
                            vz_eij = 90 - meta[svs+"_zenith"][f]
                        else:
                            vz_eij = meta[svs+"_zenith_adj"][f]
                        #print("edge?")
                else: 
                    # otherwise, use global values
                    vaij = va
                    vz_eij = vz_e
                    
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
                # find first location where shade plane is negative (on margin)
                vi, vj = np.where(sm_shade_line < -0.01)
                if (len(vi) == 0) and (len(vj) == 0):
                    sm_shade_line = sm_shade_line
                else:
                    sm_shade_line = sm_shade_line[:min(vi), :min(vj)]
                # convert to binary
                sm_shade_line[sm_shade_line > 0] = 1 
                sm_shade_line[sm_shade_line <= 0] = 0 
                
                # 6a) if seen = yes, collapse shadeline onto target pixel
                if (svs == 'view') and (sns == 'seen'):
                    if sm_shade_line.size > 0:
                        sm_shade_line = sum(sum(sm_shade_line)).reshape([1,1]) # sum over both axes
                    else:
                        sm_shade_line = np.array(0).reshape([1,1])
                    
                
                # 6b) apply the shadeline to the shade_matrix
                # only apply if matrix is larger than 0 (i.e. not truncated at source pixel)
                if len(sm_shade_line) > 0:
                    shade_matrix[i:i+sm_shade_line.shape[0],j:j+sm_shade_line.shape[1]] = shade_matrix[i:i+sm_shade_line.shape[0],j:j+sm_shade_line.shape[1]] + sm_shade_line
                else:
                    continue
                
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
    if (svs == 'view') and (sns == 'seen'):
        shade_matrix_filt = shade_matrix_filt.copy()
    else:
        shade_matrix_filt[shade_matrix_filt > 1] = 1
    print("exit loop: " + str(datetime.now().time()))

    # 8) output
    with rio.open("./example_data/shadeseen_output/" + outfn + ".tif", 'w', **metaout) as dst:
        dst.write(shade_matrix_filt[None,:,:])

# =============================================================================
# plot
# =============================================================================
    if plot == 'no':
        continue
    else:
        # plot
        fig, [ax1, ax0] = plt.subplots(2, 1)
        
        p0 = ax0.imshow(shade_matrix_filt)
        cb0 = plt.colorbar(p0, ax = ax0)
        p1 = ax1.imshow(larf, vmin = 0, vmax = 45)
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
