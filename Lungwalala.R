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
library(ggmap)
library(caTools)
library(fields)
library(dplyr)
library(purrr) 
library(magick)
library(caret)
library(viridis)
library(dplyr)
library(lattice)
library(mapproj)
library(ggsn)
setwd("C:\\UserData\\wellingm\\OneDrive - Australian National University/DEAfrica Outputs")
set.seed(1)

#load the 2019 geomedian
LungAnn <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\Annual_Geomeds\\Lungwalala2019sMAD.tif")
LungAnn
summary(LungAnn)
s2bands <- list("green", "red", "blue", "nir", "swir_1", "swir_2", "SMAD", "NDVI", "NDWI", "BSI")
LungAnn <- setNames(LungAnn, s2bands)
LungAnnplot <- plotRGB(LungAnn, r=2, g=1, b=3, stretch="hist")
summary(LungAnn)

#load training polygons
LungTrain <- readOGR("C:\\UserData\\wellingm\\OneDrive - Australian National University\\QGIS\\Lungwalala.shp")
summary(LungTrain)
plot(LungTrain)
LungTrain$Class <- as.factor(LungTrain$id)
classes <- c("Irrigated","Other")
classesdf <- data.frame(classnum = c(1,2), classnames=classes)
levels(LungTrain$Class) <- c("Irrigated","Other")

summary(LungTrain@data)
ggplot()+geom_polygon(data=LungTrain, aes(x=long, y=lat, group=id))

#' Divide data into training and validation
Lung_smp_size <- floor(0.8*nrow(LungTrain))
Lung_train_ind <- sample(seq_len(nrow(LungTrain)), size=Lung_smp_size)
LungTrain_tr <- LungTrain[Lung_train_ind,]
LungTrain_val <- LungTrain[-Lung_train_ind,]

#' RStoolbox superclass
#superclassLungAnn <- superClass(LungAnn, LungTrain, valData = NULL, responseCol = "id",
 #                               nSamples = 100000, polygonBasedCV = TRUE, trainPartition = 0.8,
  #                              model = "rf", tuneLength = 2, kfold = 5, minDist = 1,
   #                             mode = "classification", predict = FALSE, predType = "raw",
    #                            filename = NULL, verbose=TRUE, overwrite = TRUE)
#superclassLungAnn


#' Generate point samples on training data
pt_Lungtrain <- spsample(LungTrain_tr, 80000, type='random')
pt_Lungtrain$Class <- over(pt_Lungtrain, LungTrain_tr)$Class
pt_Lungtrain <- vect(pt_Lungtrain)

#' Extract spectral vals on training data
xy_train <- as.matrix(geom(pt_Lungtrain)[,c('x','y')])
df_Lungtrain <- extract(LungAnn, xy_train)
head(df_Lungtrain)
data_Lungtrain <- data.frame(Class=pt_Lungtrain$Class, df_Lungtrain)
summary(data_Lungtrain)
data_Lungtrain <- na.omit(data_Lungtrain)

#' Generate point samples on validation data
pt_Lungval <- spsample(LungTrain_val, 20000, type='random')
pt_Lungval$Class <- over(pt_Lungval, LungTrain_val)$Class
pt_Lungval <- vect(pt_Lungval)

#' Extract spectral vals on validation data
xy_val <- as.matrix(geom(pt_Lungval)[,c('x','y')])
df_Lungval <- extract(LungAnn, xy_val)
head(df_Lungval)
data_Lungval <- data.frame(Class=pt_Lungval$Class, df_Lungval)
summary(data_Lungval)
data_Lungval <- na.omit(data_Lungval)

#' Create training and validation data
LungTrainingAnn = data_Lungtrain
LungValidAnn = data_Lungval

#' RF classifier
Lung_Ann_Random_Forest_Classification <- randomForest(x=LungTrainingAnn[-1], y=as.factor(LungTrainingAnn$Class), ntree=500, importance=TRUE)
plot(Lung_Ann_Random_Forest_Classification)
Lung_Ann_Random_Forest_Classification

#' Variable importance
LungVIP <- varImpPlot(Lung_Ann_Random_Forest_Classification, type=1)
LungImp <- as.data.frame(LungVIP)
LungImp$varnames <- rownames(LungImp)
rownames(LungImp) <- NULL
LungImp
LungVIPgg <- ggplot(LungImp, aes(x=reorder(varnames, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseAccuracy))+
  xlab("Variable")+ylab("Mean Decrease Accuracy")+coord_flip()
LungVIPgg

#' Prediction stats and confusion matrix
LungPredictionAnn <- predict(Lung_Ann_Random_Forest_Classification, newdata=LungValidAnn)
conMatLungAnn <- confusionMatrix(LungPredictionAnn, as.factor(LungValidAnn$Class))
conMatLungAnn

#' Look at the map
LungMapAnn <- predict(LungAnn, Lung_Ann_Random_Forest_Classification, filename="RF_Lung.img", type="response", 
                      index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
LungProb <- predict(LungAnn, Lung_Ann_Random_Forest_Classification, filename="RF_LungProb.img", type="prob", 
                    index=1, na.rm=TRUE, progress="window", overwrite=TRUE)

plot(LungMapAnn)
summary(LungMapAnn)
LungMapAnnRat <- ratify(LungMapAnn)
rat <- levels(LungMapAnnRat)[[1]]
rat$legend <- classesdf$classnames
levels(LungMapAnnRat) <- rat
LungLevelMap <- levelplot(LungMapAnnRat,scales=list(x=list(at=NULL), y=list(at=NULL)))
LungLevelMap
LungMapAnn
LungMapGG <- as.data.frame(LungMapAnn, xy=TRUE) %>%
  mutate(Classification=factor(RF_Lung))
LungMapGG


#'Image Comparison
LungAnn
LungRGB <- ggRGB(LungAnn, r=2, g=1, b=3, stretch = "hist")+ggsn::scalebar(x.min=2657330,
  x.max= 2661200,y.min=-2254540, y.max=-2249670,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2657330,x.max= 2661200,y.min=-2254540, y.max=-2249670, location = "bottomleft", symbol=1)+
  theme_void()
LungRGB
grid.arrange(LungRGB,LungLevelMap, ncol=2)

LungMapGGProb <- as.data.frame(LungProb, xy=TRUE) %>%
  mutate(Probability=RF_LungProb)

LungMapGG <- as.data.frame(LungMapAnn, xy=TRUE) %>%
  mutate(Classification=factor(RF_Lung))

LungMapGGprob <- ggplot()+geom_raster(data=LungMapGGProb, aes(x=x, y=y, fill=Probability))+
  ggsn::scalebar(x.min=2657330,x.max= 2661200,y.min=-2254540, y.max=-2249670,dist=1,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2657330,x.max= 2661200,y.min=-2254540, y.max=-2249670, location = "bottomleft", symbol=1)+
  theme_void()+coord_fixed()+theme(legend.position = c(0.85,0.26),
                                   legend.background = element_rect(fill="white", size=5, colour="white"))
LungMapGGprob
grid.arrange(NabuMapGGprob, LungMapGGprob, ncol=2)

#' Overlaid transparent classified

LungMapGGclass <- ggRGB(LungAnn, r=2, g=1, b=3, stretch = "hist")+geom_raster(data=LungMapGG, aes(x=x, y=y, fill=Classification), alpha=0.5)+
  scale_fill_manual(values=viridis::viridis(2), breaks=1:2,labels=c("Irrigated","Other"))+
  theme_void()+ggsn::scalebar(x.min=2657330,
  x.max= 2661200,y.min=-2254540, y.max=-2249670,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2657330,x.max= 2661200,y.min=-2254540, y.max=-2249670, location = "bottomleft", symbol=1)
LungMapGGclass

LungSMAD <-  ggplot()+geom_raster(data=LungAnn,aes(x=x,y=y,fill=SMAD))+ggsn::scalebar(x.min=2657330,
  x.max= 2661200,y.min=-2254540, y.max=-2249670,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2657330,x.max= 2661200,y.min=-2254540, y.max=-2249670, 
  location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.86,0.26),
  legend.background = element_rect(fill="white", colour="white", size=5))
LungSMAD
LungNDVI <- ggplot()+geom_raster(data=LungAnn,aes(x=x,y=y,fill=NDVI))+ggsn::scalebar(x.min=2657330,
  x.max= 2661200,y.min=-2254540, y.max=-2249670,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2657330,x.max= 2661200,y.min=-2254540, y.max=-2249670, 
        location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.86,0.26),
 legend.background = element_rect(fill="white", colour="white", size=5))
LungNDVI
grid.arrange(LungRGB, LungNDVI, LungSMAD, ncol=3)
grid.arrange(LungRGB, LungNDVI, LungSMAD, LungMapGGclass, ncol=4)
grid.arrange(LungRGB,LungMapGGclass, ncol=2)
grid.arrange(KnLevelMap, LungLevelMap, NabuLevelMap, LungLevelMap, ncol=2, nrow=2)
grid.arrange(KnMapGGclass, LungMapGGclass, ncol=2)
grid.arrange(NabuMapGGclass, LungMapGGclass, ncol=2)
grid.arrange(NabuLevelMap, LungLevelMap, ncol=2)

grid.arrange(KnVIPgg, SilalVIPgg, ncol=2)
grid.arrange(NabuVIPgg, LungVIPgg, ncol=2)
grid.arrange(SilalQtrVIPgg, NabuQtrVIPgg, ncol=2)
grid.arrange(LungQtrVIPgg, ncol=2)

#' Load WaPOR
Africa <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\RS_FebIssue\\WaPOR\\L1_LCC_19.tif")
Africa

#' Crop WaPOR to area
LungAnn
Lungaoim <- extent(2657330,2661200,-2254540,-2249670)
Lungaoi <- extent(26,29,-18.9,-16)
LungWapCrop <- crop(Africa, Lungaoi)
LungWapPrj <- projectRaster(LungWapCrop, crs=crs(LungAnn))
LungWapPrj <- crop(LungWapPrj,Lungaoim)
LungWapPrj

#' Allocate classes to WaPOR LCC product
LungWapdf <- as.data.frame(LungWapPrj, xy=TRUE)
summary(LungWapdf)
LungWapdf$Land_Cover_Class <- cut(LungWapdf$layer,
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
LungWapMap <- ggplot()+geom_raster(data=LungWapdf, aes(x=x, y=y, fill=Land_Cover_Class))+coord_fixed()
LungWapMap

#' Make WaPOR classes binary and plot
LungWapdf$Classification <- cut(LungWapdf$layer,
                                breaks=c(0, 15,25,35,41.48,42.51,43.5,55,65,75,80.5,85,95,111.5,112.5,114.5,
                                         115.5,116.5,121.5,122.5,123.5,124.5,125.5,150,225,Inf),
                                labels=c("Other","Other","Other","Other", "Cropland, irrigated", "Other",
                                         "Other", "Other", "Other", "Other", "Other",
                                         "Other", "Other", "Other",
                                         "Other", "Other", "Other",
                                         "Other", "Other", "Other",
                                         "Other", "Other",  "Other",
                                         "Other", "Other"))
bicol <- c( "purple","yellow")
LungWapdf$Classification <- factor(LungWapdf$Classification , levels=c("Cropland, irrigated","Other"))
LungGMRastMapBi <- ggRGB(LungAnn, r=2, g=1, b=3, stretch = "hist")+
  geom_raster(data=LungWapdf, aes(x=x, y=y, fill=Classification), alpha=0.5)+
  ggsn::scalebar(x.min=2657330,
                 x.max= 2661200,y.min=-2254540, y.max=-2249670,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2657330,
        x.max= 2661200,y.min=-2254540, y.max=-2249670, location = "bottomleft", symbol=1)+
  theme_void()+
  scale_fill_manual(values=bicol)
LungGMRastMapBi
grid.arrange(LungMapGGclass, LungGMRastMapBi, ncol=2)

LungAnn
#' Check training data confusion for  other product
LungValidxy <- data.frame(LungValidAnn$x, LungValidAnn$y)
names(LungValidxy) <- c("x", "y")

LungWaPORconf <- extract(LungWapPrj, xy_val, cellnumbers=TRUE)
head(LungWaPORconf)
LungWaPORconf <- data.frame(LungWaPORconf)
summary(LungWaPORconf)
LungWaPORconf <- na.omit(LungWaPORconf)


LungWaPORconf$WaPOR_Class <- cut(LungWaPORconf$layer,
                                 breaks=c(0, 15,25,35,41.48,42.51,43.5,55,65,75,80.5,85,95,111.5,112.5,114.5,
                                          115.5,116.5,121.5,122.5,123.5,124.5,125.5,150,225,Inf),
                                 labels=c("Other","Other","Other","Other", "Cropland, irrigated", "Other",
                                          "Other", "Other", "Other", "Other", "Other",
                                          "Other", "Other", "Other",
                                          "Other", "Other", "Other",
                                          "Other", "Other", "Other",
                                          "Other", "Other",  "Other",
                                          "Other", "Other"))
LungWaPORconf$WaPOR_Class <- factor(LungWaPORconf$WaPOR_Class, levels=c("Cropland, irrigated","Other"))
levels(LungWaPORconf$WaPOR_Class)[levels(LungWaPORconf$WaPOR_Class)=="Cropland, irrigated"]<- "Irrigated"
LungAnnWapTab <- table(LungWaPORconf$WaPOR_Class,LungValidAnn$Class)
LungAnnWapTab
#' WAPOR IS VERT COL, TEST DAT IS HORIZ ROW
conmatLungWap <- confusionMatrix(LungAnnWapTab)
conmatLungWap
