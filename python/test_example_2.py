#!/usr/bin/python
# -*- coding: utf8 -*-
'''
Testing for drought aggregation routines.

Example 2 - 3D test to identify dummy coherent regions ######################################## 

A convoluted example involving interlocking shapes with time.

Nb relative paths - assumes code is run from [repo]/drought/python directory 
where [repo] is wherever is was cloned into.

'''
__author__="Ben Lloyd-Hughes"
__date__ ="$Mar 25, 2013 2:27:03 AM$"


#external libraries
from mayavi import mlab
import matplotlib.colors as mpc
from netCDF4 import Dataset
import numpy as np
import time
import os

#internal libraries
import structure
import plotting

#select some test data (for grid only - dummy data will be inserted below)
#e.g. land fraction from the HiGEM climate model.
infile = '../data/higem_land_qrparm.mask_frac_npsea.nc'

#Step 1 remap data to a predefined geodesic grid #####################################################
#e.g http://kiwi.atmos.colostate.edu/BUGS/geodesic/data/index.geodesic.html
#and http://kiwi.atmos.colostate.edu/BUGS/geodesic/interpolate.html
gridfile = '../data/C02562.global.nc' #predefined geodesic grid 
datafile = '..'+infile.split('.')[2]+'.'+infile.split('.')[3]+'.'+'C02562.nc'

#nb external call to cdo - 'climate data operators'
#https://code.zmaw.de/projects/cdo
cmd = 'cdo remapcon,'+gridfile +' '+ infile+' '+datafile
os.system(cmd)



#Step 2 - dummy data###########################
#
#read grid
nc = Dataset(datafile)
latv = nc.variables['lat_vertices'][:,:] + np.pi/2 # shift range to (0,pi)
lonv = nc.variables['lon_vertices'][:,:]
latc = nc.variables['lat'][:] + np.pi/2 # shift range to (0,pi)
lonc = nc.variables['lon'][:]
land = nc.variables['lsm'][:,:,:]
land = land.reshape((1,land.shape[2]))
nc.close()

#dummy five years worth of data initialized to zero
indat = np.zeros((60,2562))

#dummy some structure into the array ##############################

#start with a patch at theta=0, phi=0 of radius 1 (i.e.  7 cells on the north, prime meridian)
#then move it east one cell at each timestep then drop straight
#back again to the surface at south pole

#find neighbours (for lookup)
nn = structure.nind(datafile,1)
#nn = structure.nind(datafile,2)
#nn = structure.nind(datafile,3)

ts = np.arange(0,indat.shape[0])
#first loop up
cnt = 0
for t in ts:
    if 0.0+4*t < 180:
        wpt = structure.phitheta2cell(0, 0.0+4*t, latc, lonc)        
        nni = nn[wpt]
        for ni in nni:
            indat[t,ni] = 1.0
        cnt  = cnt+1

#then back down
wpt = structure.phitheta2cell(0, 180, latc, lonc)
for t in np.arange(0,cnt):
    nni = nn[wpt]
    for ni in nni:
        indat[t,ni] = 1.0


#ring around the equator
for lon in np.arange(0,360,1):
    wpt = structure.phitheta2cell(lon, 90, latc, lonc)
    nni = nn[wpt]
    for ni in nni:
        indat[5:10,ni] = 1.0


#simple vertical stack
wpt = structure.phitheta2cell(30, 60, latc, lonc)
nni = nn[wpt]
for ni in nni:
    indat[0:30,ni] = 1.0


#simple patch (close but not touching the ring)
wpt = structure.phitheta2cell(0, 90, latc, lonc)
nni = nn[wpt] #nb nn3
for ni in nni:
    indat[0:3,ni] = 1.0




#########################################################################

#Detect structures

#set thresh
thresh = 0.5
morethan = True
blobs = structure.aggregate(indat, nn, thresh, morethan=morethan, verbose=True)#should give sizes [1210,  617,  210,   21]



#create colormap (and normalisation funtion for lookups on blob size)
#cmap = mpc.LinearSegmentedColormap.from_list(name='mycmap', colors=['m','r','y','w','g','c','b'])
#cnorm = mpc.Normalize(vmin=np.min(blobs[1]), vmax=np.max(blobs[1]))
#
#simple lookup works better if we're just plotting the first n-largest blobs
cmap = mpc.ListedColormap(['b','r','y','w','g','c','k'])
cmap2 = mpc.ListedColormap(['#909090','#bbbbbb'])
n = 5

reload(plotting)
ib = np.arange(1,n)
#ib = [3]
Rd = plotting.blobplotMayvi(blobs, latv, lonv, cmap, ib, R=1)
#Rm = plotting.blobplotMayvi(landblobs, latvl, lonvl, cmap2, [1,2], opacity=1, R=1) #add in land below the blobs
#Rm2 = plotting.blobplotMayvi(landblobs, latvl, lonvl, cmap2, [1,2], opacity=0.1, R=Rd) #add in land above the blobs
mlab.show()




