'''
Methods to quantify drought structure.

'''
__author__="Ben Lloyd-Hughes"
__date__ ="$Aug 21, 2012 11:00:56 AM$"


import netCDF4
import numpy as np

def nind(ncfile,r):
    '''
    Nearest neighbours within radius r-cells 
    on a geodesic grid (supplied as a netcdf)

    nb grids need to specified in radians lats (-pi/2,pi/2) lons (0, 2pi)
    (nb note very crude fuzzy matching - can tolerate some grids being a bit off)

    '''
    
    #read grid
    nc = netCDF4.Dataset(ncfile)
    latc = nc.variables['lat'][:]
    lonc = nc.variables['lon'][:]
    nc.close()
    
    #sanity check the grid
    #separation of lons (lonc assumed to run from 0 to 2pi)
    if (min(lonc) < -0.001) | (max(lonc) > 2.001*np.pi):
        raise Exception("Longitude out of range (0,2pi)")
    
    #separation of lats (latc is assumed to run from -pi/2 to pi/2)
    if (min(latc) < -1.001*np.pi/2) | (max(latc) > 1.001*np.pi/2):
        raise Exception("Latitude out of range (-pi/2,pi/2)")
    
    #transform from polar to cartesian coords
    R = 1
    xc = np.zeros(len(latc))
    yc = np.zeros(len(latc))
    zc = np.zeros(len(latc))
    for i in np.arange(0,len(latc)):
        xc[i] = R * np.sin(latc[i]+ np.pi/2 ) * np.cos(lonc[i]) #nb shift lat into range(0,pi)
        yc[i] = R * np.sin(latc[i]+ np.pi/2 ) * np.sin(lonc[i])
        zc[i] = R * np.cos(latc[i]+ np.pi/2 )
        
    #loop through cells and find neighbours
    neighbours = []
    for i in np.arange(0,len(xc)):
        
        #Euclidean distance
        dx = xc - xc[i]
        dy = yc - yc[i]
        dz = zc - zc[i]
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        
        #number of neighbours is 6r + 6(r-1) + ...
        if (i == 0) or (i == 1):
            nv = 5 #nverts can be 5 if at pole
        else:
            nv = 6
        npr = 0
        for j in np.arange(0,r)+1:
            npr = npr + nv*j
        
        ni = np.argsort(dist)[0:npr+1] #nb start from 0 to collect 'self'
            
        neighbours.append(ni)
        
    return(neighbours)




def aggregate (indat, nn, thresh, tstep=1, morethan=False, verbose=False):
    '''
    A  method to identify 'coherent' groups in data based upon neigbours being 
    below the threshold 'thresh'.
    
    Nb stores the output in single array of sequentially numbered 'blobs'
    
    Keyword 'morethan' switches this to above thresh.

    indat are a 2D np array with time in the shape[0] dimension as output from 
    cdo geodesic regridding. See test_aggregation.py for example usage.

    nn is the index of neighbours for the indat grid as output from nind()
    this is a list of length indat.shape[1]. Each element of nn contains the 
    indices of the neighbours to point 'i' in indat.

    thresh is just a number.

    tstep controls the degree of temporal agregation. This defaults to immediate 
    neighbours in time.

     Returns a tuple of (blobs, sizes) where blobs is an array the same size as 
    indat with non zero elements where the blob exists. The non zero elements 
    are equal to the rank of the blob in descending order of size. ie blob id
    1 is the largest.
    '''
    nt = indat.shape[0] #how many time steps?
    npt = indat.shape[1] #how many data points?

    blobs = np.zeros(indat.shape,dtype='int') #storage for output 
    cnt = 1 #counter to keep track of blob ids

    for i in np.arange(0,nt):
      for j in np.arange(0,npt):

        if verbose: print str(i) + ' ' + str(j) + ' ' + str(cnt)

        #Does this point qualify as a blob?
        test = indat[i,j] < thresh
        if morethan:
            test = indat[i,j] > thresh

        if not test:
            continue
        else:
            blobi = np.zeros(indat.shape,dtype='int') #temporary storage for this blob
            blobi[i,j] = cnt #assign an id to this blob

            #what do the neighbours look like? Also looking backward and forward in time.
            for k in np.arange(0, len(nn[j])):
               for it in np.arange(-1*tstep, tstep+1, 1):
                   nindex = nn[j][k]
                   nit = np.clip(i+it, 0, nt - 1) #need to be a bit careful with the time range at the ends of the indat array
                   if not morethan:
                       testn = indat[nit,nindex] < thresh #Does this neighbouring point qualify as a blob?
                   else:
                       testn = indat[nit,nindex] > thresh
                   if not testn:
                      continue  #jump out quickly if nothing to do
                   else:
                      blobi[nit,nindex] = cnt #count this point as a blob
                      #does this blob intersect with any other previously stored in blobs?
                      #what do the neighbours in previous blobs look like? Also looking backward and forward in time.
                      testb = blobs[nit,nindex] > 0 #Does the intersect with any others?
                      if not testb:
                          continue
                      else: #if so, copy the previous blob into the current blob
                        cntb = blobs[nit,nindex]
                        wb = np.where(blobs == cntb)
                        blobi[wb] = cnt
            
            #store this blob in the output blob array
            wb = np.where(blobi > 0)
            blobs[wb] = cnt

            #increment blob id counter
            cnt = cnt + 1
    
    #tidy up the output - order by size and give blobs sequential ids
    #return(blobs)
    ids = np.unique(blobs)
    ids = ids[np.where(ids > 0)]
    bsizes = []
    for id in ids:
        bsizes.append(len(np.where(blobs == id)[1]))
    bsi = np.argsort(bsizes)[::-1] #indices of sorted values the [::-1] reverses the order
    blobso = np.zeros(indat.shape,dtype='int') #storage for output
    for i in np.arange(0, len(bsi)): #loop by size (largest first)
        wb = np.where(blobs == ids[bsi[i]])
        for j in np.arange(0, len(wb[0])):
            blobso[wb[0][j], wb[1][j]] = i+1 #assign sequential id (starting at 1) 
    
    #all done
    return((blobso, np.take(bsizes,bsi)))





def wrapphi(x):
    '''
    Wraps longitude into range 0,360 (wrap is a saw tooth wave)
    '''
    if (np.min(x) < 0.0):
            raise Exception("x must be positive")
    else:
            x = x % 360
    
    return(x)

def wraptheta(y):
    '''
    Wraps latitude into range 0,180 (wrap is a triangluar wave)
    [a bit clunky looking but it seems to work...]
    '''
    y = (-1) * (np.abs(((y) % 360) - 180) - 180);
    return(y)


#test wrapping routines ########################
#xt = np.arange(0,4*360)
#xw = []
#for x in xt:
#    xw.append(wrapphi(x))
#
#plt.plot(xt, xw)
#
#yt = np.arange(0,2*360)
#yw = []
#for y in yt:
#    yw.append(wraptheta(y))
#
#plt.plot(yt, yw)


def phitheta2cell (x, y, latc, lonc):
    '''
    Method to find grid cell number for a give (lat,lon) position.
    x (phi) is lon which runs 0 to 360 (west to east)
    y (theta) is lat which runs 0 to 180 (north to south)
    '''
    
    #range check x and wrap where necessary
    x = wrapphi(x)
    
    #range check y and wrap where necessary
    y = wraptheta(y)
    
    #sanity check input grid
    #lonc assumed to run from 0 to 2pi (nb note very crude fuzzy matching - can tolerate some grids being a bit off)
    if (np.min(lonc) < -0.001) | (np.max(lonc) > 2.001*np.pi):
        raise Exception("Longitude out of range (0,2pi)")
    
    #latc is assumed to run from 0 to pi
    if (np.min(latc) < -0.001) | (np.max(latc) > 1.001*np.pi):
        raise Exception("Latitude out of range (0, pi)")
    
    #convert x and y to radians
    xr = (x/360.0) * 2*np.pi
    yr = (y/360.0) * 2*np.pi
    
    #dists
    dx2 = (lonc - xr)**2
    dy2 = (latc - yr)**2
    d = np.sqrt(dx2+dy2)
    wm = np.where(d == np.min(d)) #nearest
    
    return(wm[0])
    
    
    
    
    

            

            

    
