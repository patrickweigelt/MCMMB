#Day 6 code: island biogeography!

# packages
library(maptools)
library(maps)
library(maptools)
library(classInt)
data(worldMapEnv)
library(psych)
library(car)


# load island data set

islands <- read.csv("Data/islands.csv")

islands <- islands[order(islands$spec_num),] # arrange islands by species richness (low to high)

# choose pretty colors
my.class.fr.S<-classIntervals(islands$spec_num,n=9, style="fixed", fixedBreaks=c(0,25,50,100,200,400,800,1600,3200,6400)) # bin data into n quantiles
my.pal<-c("midnightblue","dodgerblue4","dodgerblue2","chartreuse4","gold","orangered1","darkred")
my.col.fr.S<-findColours(my.class.fr.S,my.pal) # ramp colors based on classInts
my.col.fr.S[islands$spec_num < 1] <- "#666666" # assign grey for zero species

if(.Platform$OS.type == "windows"){ 
  windows(width=20*0.3937008, height=12*0.3937008)
}else{
  X11(width=20*0.3937008, height=12*0.3937008,type="cairo")
}

map('world', fill = FALSE)
title(main = "Angiosperm species richness")
points(islands$longitude, islands$latitude, pch=ifelse(islands$spec_num >= 1,1,4), cex = (log10(islands$spec_num+1)/2.6)^2+1, col = my.col.fr.S, lwd=2)

  dev.copy2pdf(file = "Data/islands.pdf")
  
  
#
  
pairs.panels(islands[,c(4:7,11:14)], density=F,ellipses=F, hist.col = "white")

# SARs

SAR0 <- lm(spec_num ~ area, data=islands)
summary(SAR0)

plot(spec_num~area,data=islands)
regLine(SAR0, col="darkgreen")

# non-linear
SAR1 <- nls(spec_num ~ c + area^z, data = islands, 
            start = list(c = 0, z = 1))

# make some predictions
newdata <- data.frame(area=seq(min(islands$area),max(islands$area),length.out = 1000))


predicted_richness_SAR1 <- predict(SAR1, newdata = newdata)

plot(islands$area, islands$spec_num, xlab = "Island area (km)", ylab="Species number", main="Species-area relationship")
points(newdata$area,predicted_richness_SAR1, type="l", col="darkred", lwd=2)

# Exercise 1

par(mfrow=c(1,2))
# 1. look at the histograms of species richness and area
hist(islands$area)
hist(islands$spec_num)

# data aren't normal - what to do?

islands$log_area <- log10(islands$area + 1)
islands$log_spec_num <- log10(islands$spec_num + 1)

par(mfrow=c(1,2))
hist(islands$log_area)
hist(islands$log_spec_num)

SAR2 <- lm(log_spec_num ~ log_area, data = islands)
summary(SAR2)

# plot linear model 2

plot(islands$log_area, islands$log_spec_num, xlab = "log10 island area (log10 km)", ylab="log10 species number", main="log-log Species-area relationship")
regLine(SAR2)

newdata_SAR3 <- data.frame(log_area=seq(min(islands$log_area),max(islands$log_area),length.out = 1000))
predicted_richness_SAR3 <- predict(SAR3, newdata = newdata_SAR3, type="response")
points(newdata_SAR3$log_area,log10(predicted_richness_SAR3+1), type="l", col="purple", lwd=2)

#
plot(islands$area, islands$spec_num, xlab = "Island area (km?)", ylab="Species number", main="Species-area relationship")

# plot the simple linear model
regLine(SAR0, col="darkgreen")

# plot the non-linear power law model
points(newdata$area,predicted_richness_SAR1, type="l", col="darkred", lwd=2)

# plot the log-log space model in untransformed space
newdata_SAR2 <- data.frame(log_area=log10(seq(min(islands$area),max(islands$area),length.out = 100)+1))
predicted_richness_SAR2 <- predict(SAR2, newdata = newdata_SAR2)
points((10^newdata_SAR2$log_area)-1,(10^predicted_richness_SAR2)-1, type="l", col="purple", lwd=2)
points(newdata$area,(10^predicted_richness_SAR2)-1, type="l", col="darkgreen", lwd=2, lty=2)

# plot the Poisson model in untransformed space
predicted_richness_SAR3 <- predict(SAR3, newdata = newdata_SAR3, type="response")
points((10^newdata_SAR3$log_area)-1,predicted_richness_SAR3, type="l", col="darkblue", lwd=2)

AIC(SAR0, SAR1, SAR2, SAR3)

# Exercise 2

# data transformation

islands$log_area<-log10(islands$area)
islands$log_spec_num<-log10(islands$spec_num+1)

islands$log_max_elev <- log10(islands$max_elev + 1)
islands$log_mean_prec <- log10(islands$mean_prec)
islands$log_mean_temp <- log10(islands$mean_temp+5)

# 

islands$geology<-as.character(islands$geology)
islands$geology[which(islands$geology == "floor" | islands$geology == "volcanic")] <- "oceanic"
islands$geology[which(islands$geology %in% c("floor" , "volcanic"))] <- "oceanic"

islands$geology[which(islands$geology == "floor")] <- "oceanic"
islands$geology[which(islands$geology == "volcanic")] <- "oceanic"

summary(as.factor(islands$geology))

#

new_islands<-islands[,c(6,8,15:19)] # create new database with only explanatory variable of interest
new_islands<-na.omit(new_islands) # remove observations with NAs

# full model


Multi.mod <- lm(log_spec_num ~ log_area + dist + log_max_elev+ log_mean_temp + log_mean_prec + geology, data= new_islands)

summary(Multi.mod)

crPlots(Multi.mod)

