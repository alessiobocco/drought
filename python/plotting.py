'''
Methods for plotting drought structure data.
'''
__author__="Ben Lloyd-Hughes"
__date__ ="$Aug 22, 2012 2:43:03 PM$"
import numpy as np


#Transform from angular to cartesian coords:
#(for sets of verts - see ang2cartc if converting grid cell centres only))
def ang2cartv(latv,lonv,R=1):
    
  xv = latv.copy()
  yv = latv.copy()
  zv = latv.copy()
  for i in np.arange(0,latv.shape[0]):
    for j in np.arange(0,latv.shape[1]):
      xv[i,j] = R * np.sin(latv[i,j]) * np.cos(lonv[i,j])
      yv[i,j] = R * np.sin(latv[i,j]) * np.sin(lonv[i,j])
      zv[i,j] = R * np.cos(latv[i,j])
  
  return((xv,yv,zv))

#Transform from angular to cartesian coords:
#(grid cell centres only)
def ang2cartc(latc,lonc,R=1):
  
  x = []
  y = []
  z = []
  for i in np.arange(0,latc.shape[0]):
    x.append(R * np.sin(latc[i]) * np.cos(lonc[i]))
    y.append(R * np.sin(latc[i]) * np.sin(lonc[i]))
    z.append(R * np.cos(latc[i]))
    
  return((x,y,z))
  



#Workout the triangles that are needed to cover the voxel
def trimesh(latvnn, lonvnn, R=1, zscl=1):
    
    pvox = []
    
    for i in np.arange(len(latvnn)): #loop each set of verts
        
        #bottom points to cartesian coords
        bpts = ang2cartc(latvnn[i,:],lonvnn[i,:],R=R) #nb using ang2cartc here since we've only got one line of verts
        bx = bpts[0]
        by = bpts[1]
        bz = bpts[2]
        #scale 'height' between grids to be roughly the same as the spatial width
        dz = np.max(np.abs(bx - bx[0]))*zscl
        #upper pts to Cartesian coords nb adjusted radius to account for 'height'
        tpts = ang2cartc(latvnn[i,:],lonvnn[i,:],R=R+dz)
        tx = tpts[0]
        ty = tpts[1]
        tz = tpts[2]
        
        #next stitch vertices into 3 long vectors
        vx = np.zeros((len(bx)*2)) #storage nb needs to be twice as big to store top and bottom of voxel
        vx[0:len(bx)] = bx
        vx[len(bx):] = tx 
        vy = np.zeros((len(bx)*2))
        vy[0:len(bx)] = by
        vy[len(bx):] = ty 
        vz = np.zeros((len(bx)*2))
        vz[0:len(bx)] = bz
        vz[len(bx):] = tz 
         
        #construct list of triplets that define the triangle that cover the polygon 
        triangles = []
        #lowerside
        for i in np.arange(2, len(bx)):
            triangles.append((0,i-1,i))
        
        #upperside
        for i in np.arange(2+len(bx), len(tx)+len(bx)):
            triangles.append((len(bx),i-1,i))
            
        #sides
        for i in np.arange(0, len(bx)):
            triangles.append((i,i+1,i+len(bx)))
            triangles.append((i,len(bx)+i-1,i+len(bx)))
        
        pvox.append((vx,vy,vz,triangles,R,dz)) #nb return R and dz so these can be fed back in repeated calls to build up vertical structure
    
    
    return(pvox)


def flattenTrimesh(tris, cols=False):
    '''
    Trimesh returns a list of triangles (20 per voxel)
    This can be slow to plot. This method flattens that list so
    that the whole lot can be plotted at once.
    
    If a list of colours is supplied (one element per voxel) then the
    list of triangles in the results tuple will have a corresponding list of colors.
    Nb input colours are expected to be of class  <type 'numpy.ndarray'> [r,g,b,a] whereas mayavi
    expects colours a (r,g,b) tuple. WebGL expects alpha  - so need to include this format too.
    '''
    nvox = len(tris)        #how many voxels
    nvert = len(tris[0][0]) #how many vertices per voxel (normally 12 i.e. hexagonal)
    ntri = len(tris[0][3])  #how many triangles to cover the voxel (20 for hexagonal voxel)
    
    #output storage
    x = np.empty(nvox*nvert)
    x[:] = np.nan
    y = np.empty(nvox*nvert)
    y[:] = np.nan
    z = np.empty(nvox*nvert)
    z[:] = np.nan
    triOut = [(np.nan,np.nan,np.nan)]*nvox*ntri
    if np.any(cols):
        colsOut = [(np.nan,np.nan,np.nan)]*nvox*ntri
        colsOutgl = [[np.nan,np.nan,np.nan,np.nan]]*nvox*ntri
    else:
        colsOut = False
        colsOutgl = False
    
    #loop
    cnt = 0
    cntTri= 0
    for ti, tri in enumerate(tris):
        nv = len(tri[0])
        x[0+cnt:nv+cnt] = tri[0][:] #nb nv can be different to nvert - but only if it is less (ditto nt below)
        y[0+cnt:nv+cnt] = tri[1][:]
        z[0+cnt:nv+cnt] = tri[2][:]
        
        nt = len(tri[3])
        triOutT = [(np.nan,np.nan,np.nan)]*nt
        if np.any(cols):
            colTi = cols[ti]
            colT = (colTi[0], colTi[1], colTi[2])
            colOutT = [colT]*nt
            colOutTgl = [colTi]*nt
        
        for i,t, in enumerate(tri[3]):
            triOutT[i] = (t[0]+cnt, t[1]+cnt, t[2]+cnt)
        triOut[0+cntTri:nt+cntTri] = triOutT
        if np.any(cols):
            colsOut[0+cntTri:nt+cntTri] = colOutT
            colsOutgl[0+cntTri:nt+cntTri] = colOutTgl
        
        cnt = cnt+nv 
        cntTri = cntTri + nt
    
    return((x,y,z,triOut, colsOut, colsOutgl))






def blobplotMayvi(blobs, latv, lonv, cmap, ib=[1], opacity=0.3, R=1.0):
    '''
    Renders an interactive MayaVi scene which contains the subset of blobs defined by the list ib.
    
    Default is to just plot the first blob.
    
    Colours are defined by a mpc.ListedColormap
    
    '''
    from mayavi import mlab

    for j,i in enumerate(ib):
        
        #extract data for this blob
        wdata = np.where(blobs[0] == i)
        
        #assign a color to this blob according to id number
        bcol = cmap(j)
        
        #inner loop over the time steps dts (need to adjust R by dz each time)
        dts = np.unique(wdata[0])
        Rin = R
        for t in dts:
            
            wtdata = np.where(wdata[0] == t)
            widata = np.take(wdata[1], wtdata)
            latvnn = np.take(latv, widata, axis=0)[0] #extract verts for this blob for this time
            lonvnn = np.take(lonv, widata, axis=0)[0]
            
            res = trimesh(latvnn, lonvnn, R)
            R = res[0][4] + res[0][5]    #adjust radius to next level (nb not perfect since this changes with lat - maybe refine?)
            resFlat = flattenTrimesh(res)
        
            #assing a color to this blob according to id number
            bcol = cmap(j)
            x = resFlat[0]
            y = resFlat[1]
            z = resFlat[2]
            triangles = resFlat[3]
            s = mlab.triangular_mesh(x, y, z*(-1), triangles, color=bcol[0:3], opacity=opacity) #nb need to flip z
        R = Rin
    mlab.show()

    return(R)
        
