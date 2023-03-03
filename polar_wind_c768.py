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

fname = os.path.join(datadir, f"20210801.000000.fv_srf_wnd.res.tile{tile}.nc")
oroname = os.path.join(griddir, f"C768_oro_data.tile{tile}.nc")

# read data and grid
ncdf = nc.Dataset(fname)
ncgf = nc.Dataset(oroname)

lat = ncgf.variables['geolat'][:]
lon = ncgf.variables['geolon'][:]
lon[np.where(lon>180)] = lon[np.where(lon>180)]-360

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
lon3d=np.empty((6,768,768))
lat3d=np.empty((6,768,768))
u3d=np.empty((6,768,768))
v3d=np.empty((6,768,768))

for tile in range(1,7): 
    fname = os.path.join(datadir, f"20210801.000000.fv_srf_wnd.res.tile{tile}.nc")
    oroname = os.path.join(griddir, f"C768_oro_data.tile{tile}.nc")

    
    #read data and grid
    ncdf = nc.Dataset(fname)
    ncgf = nc.Dataset(oroname)
    lat = ncgf.variables['geolat'][:]
    lon = ncgf.variables['geolon'][:]


    lon[np.where(lon>180)] = lon[np.where(lon>180)]-360
    u = ncdf.variables['u_srf'][0,:,:]
    v = ncdf.variables['v_srf'][0,:,:]

    #put these 2D arrays in the 3D arrays
    lat3d[tile-1,...] = lat
    lon3d[tile-1,...] = lon
    
    u3d[tile-1,...] = u
    v3d[tile-1,...] = v


#get min/max for plotting
#minval = np.nanmin(arr2d)
#maxval = np.nanmax(arr2d)

#set up map
ax = plt.axes(projection =ccrs.SouthPolarStereo())
# loop and plot data
for itile in range(0,6):
    print(itile)    
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    
    ax.quiver(lon3d[itile], lat3d[itile], u3d[itile], v3d[itile], transform = ccrs.PlateCarree(), angles = "xy", color="red", regrid_shape=20)

ax.coastlines()
ax.set_title("Surface Wind Vector South Pole C768")
plt.show()
