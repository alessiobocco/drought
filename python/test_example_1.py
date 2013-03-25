#!/usr/bin/python
# -*- coding: utf8 -*-
'''
Testing for drought aggregation routines.

Example 1 - 2D test to identify cohenerent regions of land ######################################## 

Nb relative paths - assumes code is run from [repo]/drought/python directory 
where [repo] is wherever is was cloned into.

'''
__author__="Ben Lloyd-Hughes"
__date__ ="$Mar 25, 2013 12:47:56 AM$"


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

#select some test data
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



#Step 2 read data and locate coherent regions ########################################################
nc = Dataset(datafile)
latv = nc.variables['lat_vertices'][:,:] + np.pi/2 # shift range to (0,pi)
lonv = nc.variables['lon_vertices'][:,:]
indat = nc.variables['lsm'][:,:,:]
indat = indat.reshape((1,indat.shape[2])) #drop unused dimension on lsm (format for aggregate() is 2D with time as shape[0])
nc.close()

#build map of nearest neighbours
nn = structure.nind(datafile,1) #with thresh=0.5 this gives blobs of size [406, 127, 89, 70, 39, 3, 3, 3, 2, 2, 1, 1]
#nn = structure.nind(datafile,10) #enormous radius gives one big blob

#set thresh (land are coded as 1 and sea is 0 in 'lsm' variable.
thresh = 0.5
morethan = True
blobs = structure.aggregate(indat, nn, thresh, morethan=morethan, verbose=True)


#Step 3 plot the results ############################################################################
#create colormap (and normalisation funtion for lookups on blob size)
#cmap = mpc.LinearSegmentedColormap.from_list(name='mycmap', colors=['m','r','y','w','g','c','b'])
#cnorm = mpc.Normalize(vmin=np.min(blobs[1]), vmax=np.max(blobs[1]))
#
#simple lookup works better if we're just plotting the first n-largest blobs
cmap = mpc.ListedColormap(['b','r','y','w','g','c','k'])
n = 7
reload(plotting)
ib = np.arange(0,n)
s = plotting.blobplotMayvi(blobs, latv, lonv, cmap, ib, opacity=1)
mlab.show()




