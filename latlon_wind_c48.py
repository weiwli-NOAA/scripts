import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os


#datadir = '/work2/noaa/da/cmartin/CI/GDASApp/data/lowres/gdas.20210323/12/atmos/RESTART/'
datadir = '/work/noaa/da/weiwli/GDASApp/build/fv3-jedi/test/Data/'
griddir = '/work/noaa/da/weiwli/c48/gfs_land_c48/grid'

# one tile for now
#tile=1

#fname = os.path.join(datadir, f"20210323.180000.fv_srf_wnd.res.tile{tile}.nc")
#oroname = os.path.join(griddir, f"C48_oro_data.tile{tile}.nc")

# read data and grid
#ncdf = nc.Dataset(fname)
#ncgf = nc.Dataset(oroname)

#lat = ncgf.variables['geolat'][:]
#lon = ncgf.variables['geolon'][:]
#lon[np.where(lon>180)] = lon[np.where(lon>180)]-360

#t2m = ncdf.variables['t2m'][0,:,:]

# set up map
#ax = plt.axes(projection=ccrs.PlateCarree())

# define grid
#gnomonic = ccrs.Gnomonic(central_latitude=lat[384,384], central_longitude=lon[384,384])

# plot data
#ax.pcolormesh(lon, lat, t127)
##ax.pcolormesh(lon, lat, t127, transform=gnomonic)
#ax.coastlines(zorder=10)
#ax.set_global()
#plt.show()

# loop and try all six tiles
lon3d=np.empty((6,25,48))
lat3d=np.empty((6,25,48))
u3d=np.empty((6,25,48))
v3d=np.empty((6,25,48))

for tile in range(1,1): 
       
    #fname = os.path.join(datadir, f"20210323.180000.fv_srf_wnd.res.tile{tile}.nc")
    fname = os.path.join(datadir, f"gfs.bkg.lonlat.20201215_000000z.nc4")

    oroname = os.path.join(griddir, f"C48_oro_data.tile{tile}.nc")


    
    #read data and grid
    ncdf = nc.Dataset(fname)
    ncgf = nc.Dataset(oroname)
    lat = ncgf.variables['geolat'][:]
    lon = ncgf.variables['geolon'][:]

    lon[np.where(lon>180)] = lon[np.where(lon>180)]-360
    u = ncdf.variables['u_srf'][0,:,:]
    v = ncdf.variables['v_srf'][0,:,:]
    print(tile)
    print(u)

    #put these 2D arrays in the 3D arrays
    lat3d[tile-1,...] = lat
    lon3d[tile-1,...] = lon
    
    u3d[tile-1,...] = u
    v3d[tile-1,...] = v


#get min/max for plotting
#minval = np.nanmin(arr2d)
#maxval = np.nanmax(arr2d)

#set up map
ax = plt.axes(projection =ccrs.PlateCarree())

# loop and plot data
for itile in range(0,0):
    print(itile)    
    
    ax.quiver(lon3d[itile], lat3d[itile], u3d[itile], v3d[itile], transform = ccrs.PlateCarree(), angles = "xy", color= "red", regrid_shape = 20)

ax.coastlines()
ax.set_title("Surface Wind Vector C48")
plt.show()
