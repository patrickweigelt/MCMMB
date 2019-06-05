# Day 7 code

library(plyr) # basic data manipulation

library(rgdal) #  Geospatial fun
library(rgeos) #  Geospatial fun
library(maptools)# Geospatial fun

library(Taxonstand) # Taxonomic standardisation of plant species names
library(ape) # Analyses of phylogenetics and evolution
library(picante) # Phylogenies and ecology
library(pez) # Phylogenies and ecology

library(psych) # basic stats stuff
library(classInt) # for plotting nice colors
library(colorRamps) # 

palmtrees <- read.nexus("Data/Phylogeny_ConservativeCon_Checklist.nex")

palmtree <- palmtrees[[1]]
str(palmtree)

head(palmtree$tip.label)

length(palmtree$tip.label)

plot(palmtree,type="fan",cex=0.3, edge.color="gray70",tip.color="#ef8a62")

# get palm distribution data

species <- read.csv("Data/palms_species_per_gridcell.csv",sep=",", stringsAsFactors = FALSE)
head(species)

# exercise

#1
species$species<-paste(species$GENUS,species$EPITHET, sep="_")

species$species[1:5]

#2 
length(unique(species$species))

#3 %in%

length(unique(species$species[which(!species$species %in% palmtree$tip.label)]))

a<-data.frame(unique(species$species))
colnames(a)[1]<-"species"
a$Americas<-"present"
  
b<-data.frame(palmtree$tip.label) # larger data frame
colnames(b)[1]<-"species"

c<-plyr::join(b,a,by="species",type="left")
550-summary(as.factor(c$Americas))[1]

# taxonomic standardization

specmissing <- unique(species$species[which(!species$species %in% palmtree$tip.label)])

specmissing <- gsub("_"," ",specmissing) # # replace underscore with a space
 
taxstand <- TPL(specmissing, diffchar = 2, max.distance = 1)

#

for (i in 1:nrow(taxstand)){
  species$GENUS[which(species$GENUS == taxstand$Genus[i] & species$EPITHET == taxstand$Species[i])] <- taxstand$New.Genus[i]
  species$EPITHET[which(species$GENUS == taxstand$Genus[i] & species$EPITHET == taxstand$Species[i])] <- taxstand$New.Species[i]
}

species$species <- paste(species$GENUS,species$EPITHET, sep="_")

specmissing <- unique(species$species[which(!species$species %in% palmtree$tip.label)])

species <- species[which(!species$species %in% specmissing),]


# 

americas <- readOGR("Data/americas.shp")

grid <- readOGR("Data/30min_grid_select50%.shp", integer64="allow.loss")

names(grid@data)[1] <- "grid_id"

grid_centroids <- gCentroid(grid, byid=TRUE) # make spatial points dataframe from grid

grid_centroids <- SpatialPointsDataFrame(grid_centroids, data = data.frame(grid_centroids@coords))
names(grid_centroids@data) <- c("longitude","latitude")
par(mfrow=c(1,1))
plot(grid_centroids)
plot(americas, add=TRUE)

grid@data <- cbind(grid@data,grid_centroids@data)
head(grid@data)

range(grid@data$longitude)
range(grid@data$latitude)

longitude <- seq(-116,-34, by=2)
latitude  <- seq(-35,35, by=2)

new_ID_matrix <- matrix(c(1:(length(longitude)*length(latitude))), nrow=length(longitude), ncol=length(latitude) , dimnames = list(longitude, latitude), byrow=TRUE)
grid@data$new_ID <- NA

for(i in 1:nrow(grid@data)){
  grid@data$new_ID[i] <-   new_ID_matrix[which(longitude < grid@data$longitude[i] & (longitude + 2) > grid@data$longitude[i]), which(latitude < grid@data$latitude[i] & (latitude + 2) > grid@data$latitude[i])]
}

for(i in 1:nrow(grid@data)){
  grid@data$new_ID[i] <-   new_ID_matrix[which(longitude < grid@data$longitude[i] & (longitude + 2) > grid@data$longitude[i]), which(latitude < grid@data$latitude[i] & (latitude + 2) > grid@data$latitude[i])]
}

grid_2degrees <- unionSpatialPolygons(grid, grid@data$new_ID)

uniqueIDs <- unique(grid@data$new_ID)
grid_2degrees <- SpatialPolygonsDataFrame(grid_2degrees, data=data.frame(new_ID=uniqueIDs, row.names = uniqueIDs))
head(grid_2degrees@data)


species <- plyr::join(species, grid@data, by="grid_id", type="inner")

# 
palmtree_pruned <- drop.tip(palmtree,palmtree$tip.label[which(!palmtree$tip.label %in% unique(species$species))])

length(palmtree_pruned$tip.label)

length(unique(species$species))

write.nexus(palmtree_pruned, file="Data/palmtree_pruned.nex")

# load pruned phylogeny

# build spp x site matricx

species <- unique(species[,c("new_ID","species")])
species_grid <- as.data.frame.matrix(table(species))
species_grid<-as.matrix(species_grid)

write.csv(species_grid, file = "Data/palms_specsxsites_phylo.csv")

# Faith's PD

pd_palms <- pd(species_grid, palmtree_pruned)

pd_model <- lm(PD ~ SR, data = pd_palms)
pd_palms$residuals <- residuals(pd_model)

par(mfrow=c(1,2))
plot(pd_palms$SR, pd_palms$PD, main="Species richness vs Faith's PD", xlab="Species number",  ylab="Faith's PD", cex.main=0.8,cex.axis=0.8)
abline(pd_model)
plot(pd_palms$SR, pd_palms$residuals, main="Species richness vs Faith's PD residuals", xlab="Species number",  ylab="Faith's PD residuals", cex.main=0.8,cex.axis=0.8)

c.data <- comparative.comm(palmtree_pruned , species_grid)

pd.null<-generic.null(c.data, c(.pd,.mpd,.mntd),null.model = "richness",comp.fun=.ses,permute=100)

colnames(pd.null) <- c("Faith_PD","Corrected_FaithPD","MPD", "MNTD")
rownames(pd.null) <- sites(c.data)
pd.null<-data.frame(pd.null)
pd.null$new_ID<-as.numeric(rownames(pd.null))

pd.out    <- pd.null[,c(21, 1,13,3,15,4,16)]  # select columns for final data set
pd.out [is.na(pd.out )] <- 0

par(mfrow=c(1,2))
plot(pd_palms$SR, pd_palms$residuals, main="Species richness vs Faith's PD residuals", xlab="Species number",  ylab="Faith's PD residuals", cex.main=0.8,cex.axis=0.8)
plot(pd.out$Faith_PD.observed, pd.out$Faith_PD.SES,main="Faith's PD vs standardised Faith's PD",cex.main=0.8,xlab="observed Faith's PD",ylab="standardised Faith's PD",cex.axis=0.8)


# exercise 2

#1 Create new data frame

pd_palms$newID<-as.numeric(rownames(pd_palms))

palmdiversity<-plyr::join(pd_palms, pd.out,by="new_ID")

#2 plot correlations
pairs.panels(palmdiversity[,c(3:10)], method="pearson",density=F,ellipses=F, hist.col = "white")

# corrplot is a nice option 
#3 join 'palmdiversity' with grid_2degrees

grid_2degrees@data <- plyr::join(grid_2degrees@data, palmdiversity, by = "new_ID", type="left")
head(grid_2degrees@data)

#4 make pretty maps of observed and standardized pd metrics

