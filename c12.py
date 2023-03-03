import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os


datadir = '/work/noaa/da/weiwli/c12/'
griddir = '/work/noaa/da/weiwli/c12/oro/'

# one tile for now
tile=1

fname = os.path.join(datadir, f"20201215.000000.fv_core.res.tile{tile}.nc")
oroname = os.path.join(griddir, f"C12_oro_data.tile{tile}.nc")

# read data and grid
ncdf = nc.Dataset(fname)
ncgf = nc.Dataset(oroname)

lat = ncgf.variables['geolat'][:]
lon = ncgf.variables['geolon'][:]
lon[np.where(lon>180)] = lon[np.where(lon>180)]-360

t127 = ncdf.variables['T'][0,-1,:,:]

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

# loop and try all six and get the min/max values
arr3d=np.empty((6,12,12))
lat3d=np.empty((6,12,12))
lon3d=np.empty((6,12,12))
for tile in range(1,7):
    fname = os.path.join(datadir, f"20201215.000000.fv_core.res.tile{tile}.nc")
    oroname = os.path.join(griddir, f"C12_oro_data.tile{tile}.nc")
    # read data and grid
    ncdf = nc.Dataset(fname)
    ncgf = nc.Dataset(oroname)
    lat = ncgf.variables['geolat'][:]
    lon = ncgf.variables['geolon'][:]
    lon[np.where(lon>180)] = lon[np.where(lon>180)]-360
    t127 = ncdf.variables['T'][0,-1,:,:]
    # put these 2D arrays in the 3D arrays
    arr3d[tile-1,...] = t127
    lat3d[tile-1,...] = lat
    lon3d[tile-1,...] = lon
# get min/max for plotting
minval = np.nanmin(arr3d)
maxval = np.nanmax(arr3d)

# set up map
ax = plt.axes(projection=ccrs.PlateCarree())

# loop and plot data
cmaps=['jet', 'brg', 'terrain', 'gnuplot', 'coolwarm', 'RdGy']
for itile in range(0,6):
    ax.pcolormesh(lon3d[itile], lat3d[itile], arr3d[itile], vmin=minval, vmax=maxval)
    #ax.pcolormesh(lon3d[itile], lat3d[itile], arr3d[itile], vmin=minval, vmax=maxval, cmap=cmaps[itile])
ax.coastlines()
plt.show()
