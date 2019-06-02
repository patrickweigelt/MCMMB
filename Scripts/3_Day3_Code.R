# Day 3

# load packages

library(plyr)
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(BIEN)
library(ROCR)
library(performance)

# load data

species <- read.csv("Data/palms_species_per_gridcell.csv")

species$species <- paste(species$GENUS,species$EPITHET, sep = " ") # make new column combining 'Genus' and 'Epithet'

species <- species[which(species$species %in% c("Geonoma pauciflora","Geonoma weberbaueri","Phytelephas tenuicaulis","Mauritia flexuosa")),]

unique(species$species)

species$species[which(species$species == "Geonoma weberbaueri")] <- "Geonoma undata"

# load shape files

grid <- readOGR("Data/30min_grid_select50%.shp", integer64="allow.loss")
grid@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs") # assign projection

species <- species[,c("grid_id","species")]

species_grid <- as.data.frame.matrix(table(species))
species_grid <- data.frame(grid_id = as.integer(row.names(species_grid)), species_grid)

names(grid@data)[1] <- "grid_id"

# join occurrences with spatial data
grid@data <- join(grid@data, species_grid, by="grid_id", type="left")

grid@data[is.na(grid@data)] <- 0 # convert NAs to zeroes 
grid@data[1:5,1:5] # inspect data frame

load("Data/BIEN_data.RData")

#write.csv(spec_occ, "Data/BIEN_data.csv",row.names=FALSE)

# exercise

# remove NAs from longitude/latitude

spec_occ <- spec_occ[which(!(is.na(spec_occ$latitude)|is.na(spec_occ$longitude))),] 
spec_occ <- spec_occ[which(!is.na(spec_occ$latitude) & !is.na(spec_occ$longitude)),]
head(spec_occ)

# Part 2
table(spec_occ$scrubbed_species_binomial)

#

par(mfrow=c(2,2))

plot(grid, col=ifelse(grid@data$Geonoma.pauciflora==1,"#ef8a62","#7fbf7b"), border=FALSE, main = "Geonoma pauciflora")
#species 1
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Geonoma pauciflora")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Geonoma pauciflora")], pch=19,
       col="lightgrey")

#species 2
plot(grid, col=ifelse(grid@data$Geonoma.undata==1,"#ef8a62","#7fbf7b"), border=FALSE, main = "Geonoma undata")
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Geonoma undata")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Geonoma undata")],  pch=19,
       col="lightgrey")

#species 3
plot(grid, col=ifelse(grid@data$Phytelephas.tenuicaulis==1,"#ef8a62","#7fbf7b"), border=FALSE, main = "Phytelephas tenuicaulis")
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Phytelephas tenuicaulis")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Phytelephas tenuicaulis")],  pch=19,
       col="lightgrey")

#species 4
plot(grid, col=ifelse(grid@data$Mauritia.flexuosa==1,"#ef8a62","#7fbf7b"), border=FALSE, main = "Mauritia flexuosa")
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Mauritia flexuosa")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Mauritia flexuosa")],  pch=19,
       col="lightgrey")

# create background

grid_centroids <- gCentroid(grid, byid=TRUE)
grid_centroids <- SpatialPointsDataFrame(grid_centroids, data = data.frame(grid_centroids@coords, present = 0))
names(grid_centroids@data) <- c("longitude","latitude","present")
par(mfrow=c(1,1))
plot(grid_centroids)

# load raster files

elev <- raster("Data/mn30_grd")
bio1 <- raster("Data/CHELSA_bio10_01.tif") # temperature*10
bio12 <- raster("Data/CHELSA_bio10_12.tif")

# create spatial points data frame
americas <- readOGR("Data/americas.shp")
americas@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs") 

# geonoma pauciflora

geo_pau <- SpatialPointsDataFrame(coords = spec_occ[which(spec_occ$scrubbed_species_binomial == "Geonoma pauciflora"),c("longitude","latitude")], data = data.frame(spec_occ[which(spec_occ$scrubbed_species_binomial == "Geonoma pauciflora"),c("longitude","latitude")], present = 1))
geo_pau@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs") 

# pseudo-absences
set.seed(100)
geo_pau_abs <- rbind(grid_centroids[round(runif(nrow(geo_pau), min = 1, max = nrow(grid_centroids))),],geo_pau)

par(mfrow=c(2,2))

plot(geo_pau_abs, col=ifelse(geo_pau_abs@data$present==1,"#ef8a62","#7fbf7b"))
plot(americas,add=TRUE)

geo_und <- SpatialPointsDataFrame(coords = spec_occ[which(spec_occ$scrubbed_species_binomial == "Geonoma undata"),c("longitude","latitude")], data = data.frame(spec_occ[which(spec_occ$scrubbed_species_binomial == "Geonoma undata"),c("longitude","latitude")], present = 1))

geo_und@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs") 

# add pseudo-absences
set.seed(100)
geo_und_abs <- rbind(grid_centroids[round(runif(nrow(geo_und), min = 1, max = nrow(grid_centroids))),],geo_und)

phy_ten <- SpatialPointsDataFrame(coords = spec_occ[which(spec_occ$scrubbed_species_binomial == "Phytelephas tenuicaulis"),c("longitude","latitude")], data = data.frame(spec_occ[which(spec_occ$scrubbed_species_binomial == "Phytelephas tenuicaulis"),c("longitude","latitude")], present = 1))

phy_ten@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs") 

# add pseudo-absences
set.seed(100)
phy_ten_abs <- rbind(grid_centroids[round(runif(nrow(phy_ten), min = 1, max = nrow(grid_centroids))),],phy_ten)

mau_fle <- SpatialPointsDataFrame(coords = spec_occ[which(spec_occ$scrubbed_species_binomial == "Mauritia flexuosa"),c("longitude","latitude")], data = data.frame(spec_occ[which(spec_occ$scrubbed_species_binomial == "Mauritia flexuosa"),c("longitude","latitude")], present = 1))

mau_fle@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs") 

# add pseudo-absences
set.seed(100)
mau_fle_abs <- rbind(grid_centroids[round(runif(nrow(mau_fle), min = 1, max = nrow(grid_centroids))),],mau_fle)

# add in env. layers
geo_pau_abs@data$elev <- extract(elev,geo_pau_abs)
geo_pau_abs@data$bio1 <- extract(bio1,geo_pau_abs)
geo_pau_abs@data$bio12 <- extract(bio12,geo_pau_abs)

geo_und_abs@data$elev <- extract(elev,geo_und_abs)
geo_und_abs@data$bio1 <- extract(bio1,geo_und_abs)
geo_und_abs@data$bio12 <- extract(bio12,geo_und_abs)

phy_ten_abs@data$elev <- extract(elev,phy_ten_abs)
phy_ten_abs@data$bio1 <- extract(bio1,phy_ten_abs)
phy_ten_abs@data$bio12 <- extract(bio12,phy_ten_abs)

mau_fle_abs@data$elev <- extract(elev,mau_fle_abs)
mau_fle_abs@data$bio1 <- extract(bio1,mau_fle_abs)
mau_fle_abs@data$bio12 <- extract(bio12,mau_fle_abs)

# fit sdm

logreg_geo_pau <- glm(present ~ elev + bio1 + bio12, data = geo_pau_abs@data, family = "binomial")
summary(logreg_geo_pau)

logreg_geo_und <- glm(present ~ elev + bio1 + bio12, data = geo_und_abs@data, family = "binomial")

logreg_phy_ten <- glm(present ~ elev + bio1 + bio12, data = phy_ten_abs@data, family = "binomial")

logreg_mau_fle <- glm(present ~ elev + bio1 + bio12, data = mau_fle_abs@data, family = "binomial")

# model evaluation part I

p <- predict(logreg_mau_fle, newdata=mau_fle_abs@data, type="response")
pr <- prediction(p, mau_fle_abs@data$present)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

#
new_extent <- extent(-97, -35, -30, 17) # crop the raster layers to smaller extent

elev <- crop(elev, new_extent) # crop raster layers to smaller extent
elev <- aggregate(elev, fact=4, fun=mean) # aggregate data

bio1 <- crop(bio1, new_extent)
bio1 <- aggregate(bio1, fact=4, fun=mean)
bio12 <- crop(bio12, new_extent)
bio12 <- aggregate(bio12, fact=4, fun=mean)


# predictions

new_data <- data.frame(elev = getValues(elev), bio1 = getValues(bio1), bio12 = getValues(bio12))

# predict presences
geo_pau_predict <- predict(logreg_geo_pau, newdata = new_data, type = "response")
geo_und_predict <- predict(logreg_geo_und, newdata = new_data, type = "response")
phy_ten_predict <- predict(logreg_phy_ten, newdata = new_data, type = "response")
mau_fle_predict <- predict(logreg_mau_fle, newdata = new_data, type = "response")

par(mfrow=c(1,1),mar=c(2,2,2,2))

# create raster layers with SDM predictions
geo_pau_predict_rst <- raster(new_extent, nrows=nrow(elev), ncols=ncol(elev), vals = geo_pau_predict)
geo_und_predict_rst <- raster(new_extent, nrows=nrow(elev), ncols=ncol(elev), vals = geo_und_predict)
phy_ten_predict_rst <- raster(new_extent, nrows=nrow(elev), ncols=ncol(elev), vals = phy_ten_predict)
mau_fle_predict_rst <- raster(new_extent, nrows=nrow(elev), ncols=ncol(elev), vals = mau_fle_predict)

# plot pretty stuff

par(mfrow=c(4,2), mar=c(3,1,1,1))

plot(grid, col=ifelse(grid@data$Geonoma.pauciflora==1,"#7fbf7b","white"), border=FALSE, main = "Geonoma pauciflora")
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Geonoma pauciflora")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Geonoma pauciflora")], pch=3,cex=0.3)
plot(americas,add=TRUE)


plot(geo_pau_predict_rst, main="SDM predictions")
plot(americas,add=TRUE)

plot(grid, col=ifelse(grid@data$Geonoma.undata==1,"#7fbf7b","white"),border=FALSE, main = "Geonoma undata")
plot(americas,add=TRUE)
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Geonoma undata")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Geonoma undata")], pch=3,cex=0.3)

plot(geo_und_predict_rst,main="SDM predictions")
plot(americas,add=TRUE)


plot(grid, col=ifelse(grid@data$Phytelephas.tenuicaulis==1,"#7fbf7b","white"), border=FALSE, main = "Phytelephas tenuicaulis")
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Phytelephas tenuicaulis")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Phytelephas tenuicaulis")], pch=3,cex=0.3)
plot(americas,add=TRUE)

plot(phy_ten_predict_rst,main="SDM predictions")
plot(americas,add=TRUE)

plot(grid, col=ifelse(grid@data$Mauritia.flexuosa==1,"#7fbf7b","white"), border=FALSE, main = "Mauritia flexuosa")
points(spec_occ$longitude[which(spec_occ$scrubbed_species_binomial == "Mauritia flexuosa")],spec_occ$latitude[which(spec_occ$scrubbed_species_binomial == "Mauritia flexuosa")], pch=3,cex=0.3)
plot(americas,add=TRUE)

plot(mau_fle_predict_rst,main="SDM predictions")
plot(americas,add=TRUE)