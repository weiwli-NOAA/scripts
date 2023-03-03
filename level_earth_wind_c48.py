import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import xarray as xr
import metpy
from metpy.units import units

datadir = '/work2/noaa/da/cmartin/CI/GDASApp/data/lowres/gdas.20210323/12/atmos/RESTART/'
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
lon3d=np.empty((6,48,48))
lat3d=np.empty((6,48,48))
u3d=np.empty((6,48,48))
v3d=np.empty((6,48,48))

for tile in range(1,7): 
    
    fname = os.path.join(datadir, f"20210323.180000.fv_core.res.tile{tile}.nc")
    oroname = os.path.join(griddir, f"C48_oro_data.tile{tile}.nc")
   
#read grid
    #ncdf = nc.Dataset(fname)    
    ncgf = nc.Dataset(oroname)
    
    lat = ncgf.variables['lat'][:]
    lon = ncgf.variables['lon'][:]

    geolat = ncgf.variables['geolat'][:]
    geolon = ncgf.variables['geolon'][:]

    lon[np.where(lon>180)] = lon[np.where(lon>180)]-360
    geolon[np.where(geolon>180)] = geolon[np.where(geolon>180)]-360


#read data file
    windds = xr.open_dataset(fname)

    #ugrid = windds['u_srf'].assign_coords({'ylat', 'xlon'})
    #vgrid = windds['v_srf'].assign_coords({'ylat', 'xlon'})

   # yaxis_1 = windds['yaxis_1']
   # xaxis_1 = windds['xaxis_1']
    
    
   

    #ugrid = windds['u_srf'].assign_coords({'latitude': lat, 'longitude': lon})
    #vgrid = windds['v_srf'].assign_coords({'latitude': lat, 'longitude': lon})
    #print(ugrid) 
    #windds['ylat'] = lat[:]
    #windds['xlon'] = lon[:]
    
#    ugrid = windds['u_srf'].metpy.assign_crs(grid_mapping_name = 'gnomonic_cubic')
#    vgrid = windds['v_srf'].metpy.assign_crs(grid_mapping_name = 'gnomonic_cubic')


    ugrid = windds['ua'].metpy.assign_crs(grid_mapping_name='latitude_longitude',earth_radius=6371229.0)
    vgrid = windds['va'].metpy.assign_crs(grid_mapping_name='latitude_longitude',earth_radius=6371229.0)
    #windds = windds.metpy.assign_crs
    
    ugrid = ugrid[0,99,:,:]
    vgrid = vgrid[0,99,:,:]

    #print(ugrid)
    #windds = windds.metpy.assign_crs(grid_mapping_name='')
    x = windds['xaxis_1'].values
    y = windds['yaxis_2'].values 
    #print(y)
    #print(vgrid)
    

    #vgrid = windds['v_srf'].assign_coords({'time': 'time','y': 'yaxis_1', 'x': 'xaxix_1'})
    #windds = windds.metpy.assign_crs(grid_mapping_name='x_y')


    data_crs = ugrid.metpy.cartopy_crs

    #parse full dataset
    #windds = windds.metpy.parse_cf()
    #windds = windds.metpy.assign_crs(cf_attributes=cfparams)
    #windds =windds.metpy. assign_y_x(self, force=False, tolerance=None)
    #windds = windds.metpy.assign_y_x()   

    #ugrid = windds.metpy.parse_cf('u_srf', coordinates={'y': 'yaxis_1', 'x': 'xaxis_1'}) 
    #vgrid = windds.metpy.parse_cf('v_srf', coordinates={'y': 'yaxis_1', 'x': 'xaxis_1'}) 


  
    #grab wind components
  
  

   # print(ugrid)


    #grab Cartopy CRS from metadata for plotting wind vectors
    #data_crs = ugrid.metpy.cartopy_crs
    #data_crs = windds['u_srf'].metpy.cartopy_crs()

    #x=ncdf.variables['xaxis_1'][:]
#    y=ncdf.variables['yaxis_1'][:]

#    ugrid = ncdf.variables['u_srf'][0,:,:]
#    print(ugrid)

#    vgrid = ncdf.variables['v_srf'][0,:,:]

    print(tile)

    #convert grid-relative u and v into earth relative u and v
  
#   data_crs =ugrid.metpy.cartopy_crs  

    xx, yy = np.meshgrid(x, y)
    #print(xx.shape)
    #print(yy.shape)

    ugrid = ugrid[:,:]
    vgrid = vgrid[:,:]

    #print(ugrid)
    #print(ugrid.values.shape)
    #print(vgrid.values.shape)

    uearth, vearth = ccrs.PlateCarree().transform_vectors(data_crs, xx, yy, ugrid.values, vgrid.values)
    print(ugrid)

    print(uearth)
   
  
    #put these 2D arrays in the 3D arrays
    lat3d[tile-1,...] = geolat
    lon3d[tile-1,...] = geolon
    
    u3d[tile-1,...] = uearth
    v3d[tile-1,...] = vearth
   
 
#get min/max for plotting
#minval = np.nanmin(arr2d)
#maxval = np.nanmax(arr2d)

#set up map
ax = plt.axes(projection =ccrs.PlateCarree())

# loop and plot data
for itile in range(0,6):
    print(itile)    
    
    ax.quiver(lon3d[itile], lat3d[itile], u3d[itile], v3d[itile], transform = ccrs.PlateCarree(), angles = "xy", color= "red", regrid_shape = 20)

ax.coastlines()
ax.set_title("Surface Wind Vector C48")
plt.show()
