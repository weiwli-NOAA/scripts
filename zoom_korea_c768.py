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
arr3d=np.empty((6,768,768))
lat3d=np.empty((6,768,768))
lon3d=np.empty((6,768,768))

tilenum = [2,3,4]


for tile in tilenum:
#for tile in range(1,7,2):
    print(tile)
    fname = os.path.join(datadir, f"20210801.000000.sfc_data.tile{tile}.nc")
    oroname = os.path.join(griddir, f"C768_oro_data.tile{tile}.nc")
    # read data and grid
    ncdf = nc.Dataset(fname)
    ncgf = nc.Dataset(oroname)
    lat = ncgf.variables['geolat'][0:768,0:768]
    lon = ncgf.variables['geolon'][0:768,0:768]
    lon[np.where(lon>180)] = lon[np.where(lon>180)]-360
    q2m  = ncdf.variables['q2m'][0,:,:]
    # put these 2D arrays in the 3D arrays
    arr3d[tile-1,...] = q2m
    lat3d[tile-1,...] = lat
    lon3d[tile-1,...] = lon
# get min/max for plotting
minval = np.nanmin(arr3d)
maxval = np.nanmax(arr3d)

# set up map
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([120, 135, 30, 45])

# loop and plot data
cmaps=['jet', 'brg', 'terrain', 'gnuplot', 'coolwarm', 'RdGy']


itilenum = [1,2,3]

for itile in itilenum:
    print(itile)
    #ax.pcolormesh(lon3d[itile], lat3d[itile], arr3d[itile], vmin=minval, vmax=maxval)
    ax.pcolormesh(lon3d[itile], lat3d[itile], arr3d[itile], vmin=minval, vmax=maxval, cmap=cmaps[itile])
    ax.coastlines()
plt.show()
