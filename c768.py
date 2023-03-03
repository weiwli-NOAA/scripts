import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os


datadir = '/work/noaa/da/weiwli/c768/restart/'
griddir = '/work/noaa/da/weiwli/c768/oro/'

# one tile for now
tile=1

fname = os.path.join(datadir, f"20210801.000000.sfc_data.tile{tile}.nc")
oroname = os.path.join(griddir, f"C768_oro_data.tile{tile}.nc")

# read data and grid
ncdf = nc.Dataset(fname)
ncgf = nc.Dataset(oroname)

lat = ncgf.variables['geolat'][:]
lon = ncgf.variables['geolon'][:]
lon[np.where(lon>180)] = lon[np.where(lon>180)]-360

q2m = ncdf.variables['q2m'][0,:,:]

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
#arr3d=np.empty((6,768,768))

#lat3d=np.empty((6,768,768))
#lon3d=np.empty((6,768,768))

tilenum = 3
tile = tilenum
fname = os.path.join(datadir, f"20210801.000000.sfc_data.tile{tile}.nc")
oroname = os.path.join(griddir, f"C768_oro_data.tile{tile}.nc")
# read data and grid
ncdf = nc.Dataset(fname)
ncgf = nc.Dataset(oroname)
lat = ncgf.variables['geolat'][:]
lon = ncgf.variables['geolon'][:]
lon[np.where(lon>180)] = lon[np.where(lon>180)]-360
t2m = ncdf.variables['t2m'][0,:,:]
# put these 2D arrays in the 3D arrays
#arr3d[tile-1,...] = q2m
arr2d=t2m
print(lon)
print(lat)
print(arr2d)
#lat3d[tile-1,...] = lat
#lon3d[tile-1,...] = lon
# get min/max for plotting
minval = np.nanmin(arr2d)
maxval = np.nanmax(arr2d)

#set up map
ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-180, 180, -90, 0], ccrs.PlateCarree())

#ax.plot('lon', 'lat', 'arr2d')


#ax = plt.subplot(1, 1, 1, projection = ccrs.NorthPolarStereo())
#ax.set_extent([-180, 180, -90, -60])
#ransform=ccrs.PlateCarree()
#ax = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(0, 90))

# loop and plot data

#cmaps=['terrain', 'gnuplot', 'coolwarm', 'RdGy']
#for itile in range(0,6):
   #print(itile)
#    ax.pcolormesh(lon3d[itile], lat3d[itile], arr3d[itile], vmin=minval, vmax=maxval)
plt.pcolormesh(lon, lat, arr2d, vmin=minval, vmax=maxval, transform = ccrs.PlateCarree())

ax.coastlines()
ax.gridlines()
ax.set_title("South Polar Stereographic Projection")
plt.show()
