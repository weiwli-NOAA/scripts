import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os


datadir = '/work/noaa/da/weiwli/c48/gfs_land_c48/mem001/'
griddir = '/work/noaa/da/weiwli/c48/gfs_land_c48/grid/'

# one tile for now
tile=1

fname = os.path.join(datadir, f"20191215.180000.sfc_data.tile{tile}.nc")
oroname = os.path.join(griddir, f"C48_oro_data.tile{tile}.nc")

# read data and grid
ncdf = nc.Dataset(fname)
ncgf = nc.Dataset(oroname)

lat = ncgf.variables['geolat'][:]
lon = ncgf.variables['geolon'][:]
lon[np.where(lon>180)] = lon[np.where(lon>180)]-360

t2m = ncdf.variables['t2m'][0,:,:]

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
arr3d=np.empty((6,48,48))
lat3d=np.empty((6,48,48))
lon3d=np.empty((6,48,48))

tilenum = 6
tile = tilenum
#for tile in range(1,7,2):
print(tile)
fname = os.path.join(datadir, f"20191215.180000.sfc_data.tile{tile}.nc")
oroname = os.path.join(griddir, f"C48_oro_data.tile{tile}.nc")
    # read data and grid
ncdf = nc.Dataset(fname)
ncgf = nc.Dataset(oroname)
lat = ncgf.variables['geolat'][0:48,0:48]
lon = ncgf.variables['geolon'][0:48,0:48]
lon[np.where(lon>180)] = lon[np.where(lon>180)]-360
t2m  = ncdf.variables['t2m'][0,:,:]
    # put these 2D arrays in the 3D arrays
arr3d[tile-1,...] = t2m
lat3d[tile-1,...] = lat
lon3d[tile-1,...] = lon
# get min/max for plotting
minval = np.nanmin(arr3d)
maxval = np.nanmax(arr3d)

# set up map
ax = plt.axes(projection=ccrs.SouthPolarStereo())
#ax = plt.axes(projection=ccrs.SouthPolarStereo())

ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())

# loop and plot data
cmaps=['jet', 'brg', 'terrain', 'gnuplot', 'coolwarm', 'RdGy']

itilenum = 5
itile = itilenum
print(itile)
ax.pcolormesh(lon3d[itile], lat3d[itile], arr3d[itile], vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), cmap='terrain')

#ax.pcolormesh(lon3d[itile], lat3d[itile], arr3d[itile], vmin=minval, vmax=maxval, cmap=cmaps[itile])
ax.coastlines()
ax.gridlines()
ax.set_title("C48 South Pole T_2m")

plt.show()
