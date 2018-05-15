################################################################
# Modern concepts and methods inMacroecology and Biogeography #
#                                                              #
# Day 3: Raster layers, data extraction and predictors of      #
#        species richness                                      #
#                                                              #


###### Content

# 1. the latitudinal gradient model
# 2. environmental (raster) layers 



######################################################
# 1. the latitudinal gradient model
rm(list=ls()) # remove everything 

library(plyr)
library(raster)
library(rgdal)
library(rgeos)

if (file.exists("D:/rastertemp")){rasterOptions(tmpdir="D:/rastertemp/")} # raster saves large temporary files; You can define where 

# read in RData files from hard drive
load("data/species_num.RData")

grid <- readOGR("data/30min_grid_select50%.shp", integer64="allow.loss")

grid@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs") # remember to include coordinate system

# load raster library

# Import global digital elevation model at 30 arc seconds resolution
elev <- raster("data/mn30_grd")
par(mfrow=c(1,1))
plot(elev)

# The raster layer is larger than it needs to be for Neotropical palm data. We can crop it based on the palm gridcell shapefile 
extent(grid)
elev <- crop(elev, extent(grid)) # crop raster layers to smaller extent
plot(elev)



# Import climate data from Chelsa at 30 arc seconds resolution
bio1 <- raster("data/CHELSA_bio10_1.tif")
bio1 <- crop(bio1, extent(grid)) # crop raster layers to smaller extent
bio12 <- raster("data/CHELSA_bio10_12.tif")
bio12 <- crop(bio12, extent(grid)) # crop raster layers to smaller extent

# you can do calculations with the raster layers
bio1 <- bio1/10 # Careful! this is computer intense and takes a long time

par(mfrow=c(1,2))
plot(bio1)
plot(bio12)


#####################################################
### extract environmental information for the palm
### grid polygons
par(mfrow=c(1,1))
plot(bio1, ylim=c(-6,-2), xlim=c(-82,-76))
plot(grid, add=TRUE)

### Temperature
bio1_values <- extract(bio1,grid[c(1:10),])
bio1_values
bio1_mean <- extract(bio1,grid[c(1:10),], fun=mean)
bio1_mean

head(grid@data)


# How long would it take for all gridcells?
nrow(grid)
nrow(grid)/10*system.time(bio1_mean <- extract(bio1,grid[c(1:10),], fun=mean))/60
# 20 min?

### faster work-around

# aggregate raster layer to the resolution of the palm grid
bio1_30min <- aggregate(bio1, fact=60, fun=mean, na.rm=TRUE)
dim(bio1)
extent(bio1)
dim(bio1_30min)
extent(bio1_30min)


# look at it
par(mfrow=c(1,2))

plot(bio1, ylim=c(-6,-2), xlim=c(-82,-76))
plot(grid, add=TRUE)

plot(bio1_30min, ylim=c(-6,-2), xlim=c(-82,-76))
plot(grid, add=TRUE)

grid_centroids <- gCentroid(grid, byid=TRUE)
plot(grid_centroids, add=TRUE)

# extract values based on spatial points shapefile
temp_mean <- extract(bio1_30min,grid_centroids)
#temp_mean <- temp_mean/10
temp_mean[1:10]

### Precipitation
bio12_30min <- aggregate(bio12, fact=60, fun=mean, na.rm=TRUE)
prec_mean <- extract(bio12_30min,grid_centroids)
prec_mean[1:10]

### Elevation range
elev_30min_min <- aggregate(elev, fact=60, fun=min, na.rm=TRUE)
elev_30min_max <- aggregate(elev, fact=60, fun=max, na.rm=TRUE)
elev_min <- extract(elev_30min_min,grid_centroids)
elev_max <- extract(elev_30min_max,grid_centroids)
elev_range <- elev_max - elev_min
elev_range[1:10]
hist(elev_range)



elev_range[1:10]
prec_mean[1:10]
temp_mean[1:10]

#head(grid@data)
grid@data$ID[1:10]


# join environmental data and species numbers


#model

#predictions plotting





