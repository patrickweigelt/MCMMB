##############
# Day 5 ######
##############

elev <- raster("Data/mn30_grd")
bio1 <- raster("Data/CHELSA_bio10_01.tif") # temperature*10
bio12 <- raster("Data/CHELSA_bio10_12.tif") # precipitation

grid <- readOGR("Data/30min_grid_select50%.shp", integer64="allow.loss")

grid@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs")

americas <- readOGR("Data/americas.shp")
americas@proj4string <- CRS("+proj=longlat +ellps=WGS84 +no_defs")

extent(grid)

elev <- crop(elev, extent(grid)) # crop raster layers to smaller extent
bio1 <- crop(bio1, extent(grid))
bio12 <- crop(bio12, extent(grid))

# test aggregation
bio1_values <- extract(bio1,grid[c(1:10),])

mean(bio1_values[[1]][1:100])

# calculate mean for each grid cell
bio1_mean <- sapply(bio1_values, mean)

bio1_mean <- extract(bio1,grid[c(1:10),], fun=mean) 

#  how slow am I today?
nrow(grid)/10*system.time(bio1_mean <- extract(bio1,grid[c(1:10),], fun=mean))/60

# lets do this faster
bio1_30min <- aggregate(bio1, fact=60, fun=mean, na.rm=TRUE)
dim(bio1)
extent(bio1)

grid_centroids <- gCentroid(grid, byid=TRUE)

#extract raster data 

temp_mean <- extract(bio1_30min,grid_centroids)
temp_mean<-temp_mean/10

bio12_30min <- aggregate(bio12, fact=60, fun=mean, na.rm=TRUE)
prec_mean <- extract(bio12_30min,grid_centroids)

elev_30min_min <- aggregate(elev, fact=60, fun=min, na.rm=TRUE) # minimum elevation
elev_30min_max <- aggregate(elev, fact=60, fun=max, na.rm=TRUE) # maximum elevation
elev_min <- extract(elev_30min_min,grid_centroids)

elev_max <- extract(elev_30min_max,grid_centroids)
elev_range <- elev_max - elev_min # calculate range

grid@data$temp <- temp_mean
grid@data$prec <- prec_mean
grid@data$elev <- elev_range

# check collinearity

cor(grid@data[,c("temp","prec","elev")])

par(mfrow=c(1,3))
hist(grid@data$temp, main="", xlab="MAT (C)")
hist(grid@data$prec, main="", xlab="MAP (mm/yr)")
hist(grid@data$elev, main="", xlab="Elevation (m)")

load("Data/species_num.RData")

# exercise
colnames(species_num)
colnames(grid@data)

names(grid@data)[1] <- "grid_id" # rename ID column in gridcell shapefile

species_num <- join(species_num, grid@data, by="grid_id")
head(species_num)

# fit models

lm_elev<-lm(spec_num~elev, data=species_num)
lm_precip<-lm(spec_num~prec, data=species_num)
lm_temp<-lm(spec_num~temp, data=species_num)

# compare models
AIC(lm_elev, lm_precip, lm_temp)

#
plot(species_num$spec_num ~ species_num$elev, xlab="Elevation a.s.l. (m)",  ylab="Species number")
abline(lm_elev, col="orange")

plot(species_num$spec_num ~ species_num$temp, xlab="Elevation a.s.l. (m)",  ylab="Species number")
abline(lm_temp, col="orange")

plot(species_num$spec_num ~ species_num$prec, xlab="Elevation a.s.l. (m)",  ylab="Species number")
abline(lm_precip, col="purple")

model_full <- lm(spec_num ~ temp + prec + elev,data=species_num)
summary(model_full)

AIC(model_full, lm_elev, lm_precip, lm_temp)


crPlots(model_full, layout=c(1, 3), ask=FALSE, main="Partial residual plots", grid=FALSE, ylab="Partial residuals (species richness)")

species_num$pred <- predict(model_full)
species_num$resi <- residuals(model_full)

grid@data <- join(grid@data, 
                  species_num[,c("grid_id","spec_num","pred","resi")], type="left", by="grid_id")


par(mfrow=c(1,3), mar=c(1,1,2,1))
#Create color scheme

my.class.fr<-classIntervals(grid@data$spec_num, n=10, style="equal", dataPrecision=0) # bin data into n quantiles
my.pal<-matlab.like(10)
my.col.fr<-findColours(my.class.fr,my.pal) # ramp colors based on classInts

plot(grid, col=my.col.fr, border=NA, main="Palm species richness")
plot(americas, add=TRUE)

legend("bottomleft", # position
       legend = names(attr(my.col.fr, "table")),
       title = "Species number",
       fill = attr(my.col.fr, "palette"),
       cex = 1,
       bty = "n") # no box around it


my.class.fr<-classIntervals(grid@data$pred, n=10, style="equal", dataPrecision=0) # bin data into n quantiles


my.pal<-matlab.like(10)
my.col.fr<-findColours(my.class.fr,my.pal) # ramp colors based on classInts

plot(grid, col=my.col.fr, border=NA, main="Model predictions")
plot(americas, add=TRUE)

legend("bottomleft", # position
       legend = names(attr(my.col.fr, "table")),
       title = "Species number",
       fill = attr(my.col.fr, "palette"),
       cex = 1,
       bty = "n") # no box around it



my.class.fr<-classIntervals(grid@data$resi, n=10, style="equal", dataPrecision=0)

my.pal<-green2red(10)
my.col.fr<-findColours(my.class.fr,my.pal) # ramp colors based on classInts

plot(grid, col=my.col.fr, border=NA, main="Model residuals")
plot(americas, add=TRUE)

legend("bottomleft", # position
       legend = names(attr(my.col.fr, "table")),
       title = "Species number",
       fill = attr(my.col.fr, "palette"),
       cex = 1,
       bty = "n") # no box around it