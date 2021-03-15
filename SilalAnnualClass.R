library(raster)
library(rgdal)
library(sp)
library(randomForest)
library(RStoolbox)
library(e1071)
library(reshape2)
library(cluster)
library(terra)
library(rpart)
library(plyr)
library(caTools)
library(gridExtra)
library(fields)
library(lattice)
library(dplyr)
library(purrr) 
library(magick)
library(caret)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(ggsn)
setwd("C:\\UserData\\wellingm\\OneDrive - Australian National University/DEAfrica Outputs")
set.seed(1)

#load the 2019 geomedian
SilalAnn <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\Annual_Geomeds\\Silalatshani2019sMAD.tif")

summary(SilalAnn)
s2bands <- list("green", "red", "blue", "nir", "swir_1", "swir_2", "SMAD", "NDVI", "NDWI", "BSI")
SilalAnn <- setNames(SilalAnn, s2bands)
SilalAnnplot <- plotRGB(SilalAnn, r=2, g=1, b=3, stretch="hist")
SilalAnn


#load training polygons
SilalTrain <- readOGR("C:\\UserData\\wellingm\\OneDrive - Australian National University\\QGIS\\Silalatshani.shp")
summary(SilalTrain)
plot(SilalTrain)
SilalTrain$id <- as.factor(SilalTrain$id)
SilalTrain$Class <- as.factor(SilalTrain$id)
classes <- c("Irrigated","Other")
classesdf <- data.frame(classnum = c(1,2), classnames=classes)
levels(SilalTrain$Class) <- c("Irrigated","Other")

#' Divide data into training and validation
Silal_smp_size <- floor(0.8*nrow(SilalTrain))
Silal_train_ind <- sample(seq_len(nrow(SilalTrain)), size=Silal_smp_size)
SilalTrain_tr <- SilalTrain[Silal_train_ind,]
SilalTrain_val <- SilalTrain[-Silal_train_ind,]

#' RStoolbox superclass
#superclassSilalAnn <- superClass(SilalAnn, SilalTrain, valData = NULL, responseCol = "id",
                              #nSamples = 100000, polygonBasedCV = TRUE, trainPartition = 0.8,
                              #model = "rf", tuneLength = 2, kfold = 5, minDist = 1,
                             # mode = "classification", predict = FALSE, predType = "raw",
                             # filename = NULL, verbose=TRUE, overwrite = TRUE)
#superclassSilalAnn


#' Generate point samples on training data
pt_silaltrain <- spsample(SilalTrain_tr, 80000, type='random')
pt_silaltrain$Class <- over(pt_silaltrain, SilalTrain_tr)$Class
pt_silaltrain <- vect(pt_silaltrain)

#' Extract spectral vals on training data
xy_train <- as.matrix(geom(pt_silaltrain)[,c('x','y')])
df_silaltrain <- extract(SilalAnn, xy_train)
head(df_silaltrain)
data_silaltrain <- data.frame(Class=pt_silaltrain$Class, df_silaltrain)
summary(data_silaltrain)
data_silaltrain <- na.omit(data_silaltrain)

#' Generate point samples on validation data
pt_silalval <- spsample(SilalTrain_val, 20000, type='random')
pt_silalval$Class <- over(pt_silalval, SilalTrain_val)$Class
pt_silalval <- vect(pt_silalval)

#' Extract spectral vals on validation data
xy_val <- as.matrix(geom(pt_silalval)[,c('x','y')])
df_silalval <- extract(SilalAnn, xy_val)
head(df_silalval)
data_silalval <- data.frame(Class=pt_silalval$Class, df_silalval)
summary(data_silalval)
data_silalval <- na.omit(data_silalval)

#' Create training and validation data
SilalTrainingAnn = data_silaltrain
SilalValidAnn = data_silalval

#' RF classifier
Silal_Ann_Random_Forest_Classification <- randomForest(x=SilalTrainingAnn[-1], y=as.factor(SilalTrainingAnn$Class), ntree=500, importance=TRUE)
plot(Silal_Ann_Random_Forest_Classification)
Silal_Ann_Random_Forest_Classification

#' Variable importance
SilalVIP <- varImpPlot(Silal_Ann_Random_Forest_Classification, type=1)
silalImp <- as.data.frame(SilalVIP)
silalImp$varnames <- rownames(silalImp)
rownames(silalImp) <- NULL
silalImp
SilalVIPgg <- ggplot(silalImp, aes(x=reorder(varnames, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseAccuracy))+
  xlab("Variable")+ylab("Mean Decrease Accuracy")+coord_flip()

#' Prediction stats and confusion matrix
SilalPredictionAnn <- predict(Silal_Ann_Random_Forest_Classification, newdata=SilalValidAnn)
conMatSilalAnn <- confusionMatrix(SilalPredictionAnn, as.factor(SilalValidAnn$Class))
conMatSilalAnn

#' Look at the map
SilalMapAnn <- predict(SilalAnn, Silal_Ann_Random_Forest_Classification, filename="RF_Silal.img", type="response", 
                    index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
plot(SilalMapAnn)

SilalProb <- predict(SilalAnn, Silal_Ann_Random_Forest_Classification, filename="RF_SilalProb.img", type="prob", 
                  index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
SilalMapGGProb <- as.data.frame(SilalProb, xy=TRUE) %>%
  mutate(Probability=RF_SilalProb)

SilalMapAnnRat <- ratify(SilalMapAnn)
rat <- levels(SilalMapAnnRat)[[1]]
rat$legend <- classesdf$classnames
levels(SilalMapAnnRat) <- rat
SilalLevelMap <- levelplot(SilalMapAnnRat, scales=list(x=list(at=NULL), y=list(at=NULL)))
SilalLevelMap
SilalMapGG <- as.data.frame(SilalMapAnn, xy=TRUE) %>%
  mutate(Classification=factor(RF_Silal))
SilalMapGG
bicol <- c("grey44", "orange3")


#'Image Comparison
SilalAnn
SilalRGB <- ggRGB(SilalAnn, r=2, g=1, b=3, stretch = "hist")+ggsn::scalebar(x.min=2823740,
  x.max= 2829530,y.min=-2599410, y.max=-2592230,dist=2,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230, location = "bottomleft", symbol=1)+
  theme_void()
SilalRGB
grid.arrange(SilalRGB,SilalLevelMap, ncol=2)

#' Overlaid transparent classified
SilalMapGGclass <- ggRGB(SilalAnn, r=2, g=1, b=3, stretch = "hist") +geom_raster(data=SilalMapGG, aes(x=x, y=y, 
  fill=Classification), alpha=0.5)+
  scale_fill_manual(values=viridis::viridis(2), breaks=1:2,labels=c("Irrigated","Other"))+
  ggsn::scalebar(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230,dist=2,transform=FALSE, 
  dist_unit="km",model="WGS84")+north(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230,
  location = "bottomleft", symbol=1)+theme_void()
SilalMapGGclass

SilalMapGGprob <- ggplot() +geom_raster(data=SilalMapGGProb , aes(x=x, y=y, 
  fill=Probability))+
  ggsn::scalebar(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230,dist=2,transform=FALSE, 
                 dist_unit="km",model="WGS84")+north(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230,
   location = "bottomleft", symbol=1)+theme_void()+coord_fixed()+theme(legend.position = c(0.85,0.24),
   legend.background = element_rect(fill="white", colour="white",size=5))
SilalMapGGprob

grid.arrange(KnMapGGprob, SilalMapGGprob, ncol=2)

SilalSMAD <- ggplot()+geom_raster(data=SilalAnn,aes(x=x,y=y,fill=SMAD))+ggsn::scalebar(x.min=2823740,x.max= 2829530,
 y.min=-2599410, y.max=-2592230,dist=2,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230, 
  location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.86,0.26),
  legend.background = element_rect(fill="white", colour="white", size=5))
SilalSMAD
SilalNDVI <- ggplot()+geom_raster(data=SilalAnn,aes(x=x,y=y,fill=NDVI))+ggsn::scalebar(x.min=2823740,x.max= 2829530,
  y.min=-2599410, y.max=-2592230,dist=2,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230, 
        location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.86,0.26),
    legend.background = element_rect(fill="white", colour="white", size=5))
SilalNDVI
grid.arrange(SilalRGB, SilalNDVI, SilalSMAD, ncol=3)
grid.arrange(SilalSMAD, SilalMapGGclass, ncol=2)
grid.arrange(KnMapGGprob, SilalMapGGprob, ncol=2)

#' Load WaPOR
Africa <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\RS_FebIssue\\WaPOR\\L1_LCC_19.tif")
Africa

#' Crop WaPOR to area
aoim <- extent(2823740,2829530,-2599410,-2592230)
aoi <- extent(29.19,29.39,-20.89,-20.69)
SilalWapCrop <- crop(Africa, aoi)
SilalWapPrj <- projectRaster(SilalWapCrop, crs=crs(SilalAnn))
SilalWapPrj <- crop(SilalWapPrj,aoim)
SilalWapPrj

#' Allocate classes to WaPOR LCC product
SilalWapdf <- as.data.frame(SilalWapPrj, xy=TRUE)
summary(SilalWapdf)
SilalWapdf$Land_Cover_Class <- cut(SilalWapdf$layer,
                                   breaks=c(0, 15,25,35,41.48,42.51,43.5,55,65,75,80.5,85,95,111.5,112.5,114.5,
                                            115.5,116.5,121.5,122.5,123.5,124.5,125.5,150,225,Inf),
                                   labels=c("n.a","Shrubland","Grassland","Cropland, rainfed", "Cropland, irrigated", "Cropland, fallow",
                                            "Built-up", "Bare/Sparse vegetation", "Snow/ice", "water bodies", "Temp. water bodies",
                                            "Shrub cover, flooded", "Tree-cover, evergreen", "Tree-cover, evergreen broadleaved",
                                            "Tree-cover, closed deciduous", "Tree-cover, closed mixed", "Tree cover, closed unknown",
                                            "Tree-cover, open evergreen", "Tree cover, open evergreen broadleaf", "Tree cover, open deciduous",
                                            "Tree-cover, open deciduous broadleaf", "Tree-cover, open mixed",  "Tree cover, open unknown",
                                            "Sea water", "No data"))


mycols <- c("grey44", "orange3", "lawngreen", "olivedrab4", "orangered3", "gold",
            "gray10", "gray95", "azure", "royalblue2", "royalblue4",
            "paleturquoise", "seagreen1", "seagreen2",
            "seagreen3", "seagreen4", "seagreen",
            "palegreen", "palegreen1", "palegreen2",
            "palegreen3", "palegreen4",  "greenyellow",
            "navy", "gray100")
SilaWapMap <- ggplot()+geom_raster(data=SilalWapdf, aes(x=x, y=y, fill=Land_Cover_Class))
SilaWapMap

#' Make WaPOR classes binary and plot
SilalWapdf$Classification <- cut(SilalWapdf$layer,
                                 breaks=c(0, 15,25,35,41.48,42.51,43.5,55,65,75,80.5,85,95,111.5,112.5,114.5,
                                          115.5,116.5,121.5,122.5,123.5,124.5,125.5,150,225,Inf),
                                 labels=c("Other","Other","Other","Other", "Cropland, irrigated", "Other",
                                          "Other", "Other", "Other", "Other", "Other",
                                          "Other", "Other", "Other",
                                          "Other", "Other", "Other",
                                          "Other", "Other", "Other",
                                          "Other", "Other",  "Other",
                                          "Other", "Other"))
bicol <- c("purple","yellow")
SilalWapdf$Classification <- factor(SilalWapdf$Classification , levels=c("Cropland, irrigated","Other"))
SilalGMRastMapBi <- ggRGB(SilalAnn, r=2, g=1, b=3, stretch = "hist")+
  geom_raster(data=SilalWapdf, aes(x=x, y=y, fill=Classification), alpha=0.5)+
  ggsn::scalebar(x.min=2823740, x.max= 2829530,y.min=-2599410, y.max=-2592230,dist=2,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2823740,x.max= 2829530,y.min=-2599410, y.max=-2592230, location = "bottomleft", symbol=1)+
  theme_void()+
  scale_fill_manual(values=bicol)
SilalGMRastMapBi

grid.arrange(SilalMapGGclass, SilalGMRastMapBi, ncol=2)

#' Check training data confusion for  other product

#' get xy from valid here
Validxy <- data.frame(SilalValidAnn$x, SilalValidAnn$y)
names(Validxy) <- c("x", "y")

SilWaPORconf <- extract(SilalWapPrj, xy_val, cellnumbers=TRUE)
head(SilWaPORconf)
SilWaPORconf <- data.frame(SilWaPORconf)
summary(SilWaPORconf)
SilWaPORconf <- na.omit(SilWaPORconf)


SilWaPORconf$WaPOR_Class <- cut(SilWaPORconf$layer,
                                      breaks=c(0, 15,25,35,41.48,42.51,43.5,55,65,75,80.5,85,95,111.5,112.5,114.5,
                                               115.5,116.5,121.5,122.5,123.5,124.5,125.5,150,225,Inf),
                                      labels=c("Other","Other","Other","Other", "Cropland, irrigated", "Other",
                                               "Other", "Other", "Other", "Other", "Other",
                                               "Other", "Other", "Other",
                                               "Other", "Other", "Other",
                                               "Other", "Other", "Other",
                                               "Other", "Other",  "Other",
                                               "Other", "Other"))
SilWaPORconf$WaPOR_Class <- factor(SilWaPORconf$WaPOR_Class, levels=c("Cropland, irrigated","Other"))
levels(SilWaPORconf$WaPOR_Class)[levels(SilWaPORconf$WaPOR_Class)=="Cropland, irrigated"]<- "Irrigated"


SilalAnnWapTab <- table(SilWaPORconf$WaPOR_Class,SilalValidAnn$Class)
SilalAnnWapTab
#' WAPOR IS VERT COL, TEST DAT IS HORIZ ROW
conmatSilWap <- confusionMatrix(SilalAnnWapTab)
conmatSilWap

