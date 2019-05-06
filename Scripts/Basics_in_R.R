################################################################
# Modern concepts and methods inMacroecology and Biogeography #
#                                                              #
# Day 1: Basics in R and the latitudinal gradient              #
#                                                              #
# www.r-project.org/                                           #


### General remark

# Set the format your comupter uses to display numbers to 
# English format for exchanging data with others

# decimal symbol: "."
# Digit grouping symbol: ","
# List seperator: ","


###### Content

# 1. getting started
# 2. data structures and indexing
# 3. loading data into R



### 1. Getting started

# use the console to do simple calculations
2 + 2

2 * 2

8 / 2

6 - 2

sqrt(16) # sqrt() is a function

16^(1/2)

2^2

1/2 ^ (-2)

log10(10000)


# assign a value to an object

a <- 4
a

b <- sqrt(16)
b

number_3 <- 2 * 2
number_3

text_a <- "Hello"
text_a

text_a <- "Good bye" # overwrite the object
text_a

# see the objects in memory; In R Studio the objects are summed up at the top-right corner
ls()


# getting help

# ?function.name
# ??keyword

?sqrt
??squareroot

?ls




### 2. data structures and indexing

## values
mode(a)
str(a)

mode(text_a)
str(text_a)


## vectors (a sequence of values of the same mode)

vector_1 <- c(1,5,6,3.4,10,11)
vector_1

vector_1[3]
vector_1[3] <- 5
vector_1[3]

vector_2 <- c(1:6)
vector_2

length(vector_2)
mode(vector_2)

mean(vector_2)
sd(vector_2)

vector_3 <- rep(c(1,2),3)
vector_3

vector_4 <- c("blue","red","green","yellow","grey","black")
mode(vector_4)

vector_4[4]


# you can do calculations with numeric vectors just like with values
vector_2
vector_2 * 2
vector_2 ^ 2

vector_1 + vector_2
vector_5 <- vector_1 + vector_2


## arrays and matrices (tables of values of the same mode with n dimensions (matrix = 2 dimensions))
?matrix
matrix_1 <- matrix(NA, nrow = 6, ncol = 6)
matrix_1

matrix_1 <- matrix(1:(6^2), nrow = 6, ncol = 6)
matrix_1

matrix_2 <- matrix(1:(6^2), nrow = 6, ncol = 6, byrow = TRUE)
matrix_2

matrix_2[3,5]
matrix_2[3,]
matrix_2[,5]

# transpose matrix
t(matrix_1)
matrix_2


# you can do calculations with numeric matrices just like with values and vectors
matrix_1 + 1
matrix_1 ^ 2
mean(matrix_1)
colMeans(matrix_1)

matrix_1 + matrix_2

matrix_1 * vector_2 # vector is too short and is recycled!


## data.frames (tables of values of one mode per column)
df_1 <- data.frame(ID = vector_2, number_1 = vector_1, colour = vector_4)
df_1

str(df_1) # colour is a factor; makes sense if certain values are repeated a lot of times
df_1[,3]

df_1 <- data.frame(ID = vector_2, number_1 = vector_1, colour = vector_4, stringsAsFactors = FALSE)
str(df_1) # colour is a character vector
df_1[,3]
df_1$colour

plot(df_1$ID,df_1$number_1, col=df_1$colour)



######################################################
## data.frame exercise

# Imagine a field study with 5 plots (1,2,3,...) in each of 3 vegetation types (forest, shrubland, grassland)
# create a dataframe with 1 row per plot and 2 columns for vegetation_type, plot_ID
# Every second day of a month a student visits one of the plots. Add a third variable called day starting from 1 and ending at 29 using the function seq()



######################################################
######################################################
### 3. loading data into R
ls()
rm(list=ls()) # remove everything 
setwd("D:/Teaching/Macroecology/Practicals/Day1/data") # adjust to your needs
getwd()

## libraries / packages
?install.packages
?update.packages

#install.packages("plyr")
library(plyr)
search() # Gives a list of attached packages and R objects
help(package="plyr")


######################################################
## reading in a data.frame
?read.table

# Data from  Kreft, H., Sommer, J.H. & Barthlott, W. (2006)
#            The significance of geographic range size for spatial diversity
#            patterns in Neotropical palms. Ecography, 29, 21-30.


species <- read.csv("palms_species_per_gridcell.csv",sep=",")
View(species)
str(species)
summary(species)
head(species)
tail(species)

nrow(species)
ncol(species)
dim(species)

species$GENUS
unique(species$GENUS)
class(species$GENUS)
levels(species$GENUS)

colnames(species)
rownames(species)


## sorting
order(species$grid_id)
species <- species[order(species$grid_id),]
View(species)


## Adding a column for the full species name
species$species <- paste(species$GENUS,species$EPITHET, sep = " ")


## subsetting the dataframe
species_unique <- unique(species[,c("spp_Id","GENUS","EPITHET","species")])


# count species per gridcell
?table
species_num <- data.frame(table(species$grid_id))
head(species_num)
names(species_num) <- c("grid_id","spec_num")
head(species_num)

str(species_num)
species_num$grid_id <- as.numeric(as.character(species_num$grid_id))
head(species_num)

str(species_num)



##########################################################
## loading a GIS shapefile into R
library(rgdal)
library(rgeos)
library(maptools)

### Exercise
# Open the two shapefiles provided into QGIS and have a look at them


# read a polygon shapefile
grid <- readOGR("30min_grid_select50%.shp", integer64="allow.loss")
plot(grid)


### calculate latitude and longitude from the grids polygons
# convert the polygons to points (mass centroids)
grid_centroids <- gCentroid(grid, byid=TRUE)

plot(grid_centroids)

# extract the coordinates of the points
coordinates <- data.frame(grid_centroids@coords)
head(coordinates)

# access the attribute table
head(grid@data)

# add grid_id to coordinates table
coordinates <- data.frame(grid_id = grid@data$ID, Long = coordinates$x, Lat = coordinates$y)



##########################################################
## joins species numbers and coordinates

species_num <- join(species_num, coordinates, by="grid_id", type="inner")

View(species_num)


## join species numbers to shapefile
names(grid@data)[1] <- "grid_id"
grid@data <- join(grid@data, species_num, by="grid_id", type="left")


## Make a species richness map

#Create color scheme      
library(classInt)
library(colorRamps)

my.class.fr<-classIntervals(grid@data$spec_num,n=10,style="equal") # bin data into n quantiles
my.pal<-matlab.like(10)
my.col.fr<-findColours(my.class.fr,my.pal) # ramp colors based on classInts

plot(grid, col=my.col.fr)






