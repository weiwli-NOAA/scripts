Cartopy can not correctly plot a tile when the tile contains both positive and negative longitude values. Tile 4 is one example. 
More work is needed to split Tile 4 into 2 tiles.
Grid-relative winds need to be converted into earth-related, in order to plot with PlateCarree projection.
