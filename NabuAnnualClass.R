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
library(caTools)
library(fields)
library(dplyr)
library(purrr) 
library(magick)
library(caret)
library(dplyr)
library(viridis)
library(gridExtra)
library(ggsn)
setwd("")
set.seed(1)
#load the 2019 geomedian
NabuAnn <- stack(".tif")

summary(NabuAnn)
s2bands <- list("green", "red", "blue", "nir", "swir_1", "swir_2", "SMAD", "NDVI", "NDWI", "BSI")
NabuAnn <- setNames(NabuAnn, s2bands)
NabuAnnplot <- plotRGB(NabuAnn, r=2, g=1, b=3, stretch="hist")
summary(NabuAnn)

#load training polygons
NabuTrain <- readOGR(".shp")
summary(NabuTrain)
plot(NabuTrain)
NabuTrain$id <- as.factor(NabuTrain$id)
NabuTrain$Class <- as.factor(NabuTrain$id)
classes <- c("Irrigated","Other")
classesdf <- data.frame(classnum = c(1,2), classnames=classes)
levels(NabuTrain$Class) <- c("Irrigated","Other")

#' Divide data into training and validation
Nabu_smp_size <- floor(0.8*nrow(NabuTrain))
Nabu_train_ind <- sample(seq_len(nrow(NabuTrain)), size=Nabu_smp_size)
NabuTrain_tr <- NabuTrain[Nabu_train_ind,]
NabuTrain_val <- NabuTrain[-Nabu_train_ind,]

#' RStoolbox superclass
#superclassNabuAnn <- superClass(NabuAnn, NabuTrain, valData = NULL, responseCol = "id",
 #                                nSamples = 100000, polygonBasedCV = TRUE, trainPartition = 0.8,
  #                               model = "rf", tuneLength = 2, kfold = 5, minDist = 1,
   #                              mode = "classification", predict = FALSE, predType = "raw",
    #                             filename = NULL, verbose=TRUE, overwrite = TRUE)
#superclassNabuAnn

#' Generate point samples on training data
pt_Nabutrain <- spsample(NabuTrain_tr, 80000, type='random')
pt_Nabutrain$Class <- over(pt_Nabutrain, NabuTrain_tr)$Class
pt_Nabutrain <- vect(pt_Nabutrain)

#' Extract spectral vals on training data
xy_train <- as.matrix(geom(pt_Nabutrain)[,c('x','y')])
df_Nabutrain <- extract(NabuAnn, xy_train)
head(df_Nabutrain)
data_Nabutrain <- data.frame(Class=pt_Nabutrain$Class, df_Nabutrain)
summary(data_Nabutrain)
data_Nabutrain <- na.omit(data_Nabutrain)

#' Generate point samples on validation data
pt_Nabuval <- spsample(NabuTrain_val, 20000, type='random')
pt_Nabuval$Class <- over(pt_Nabuval, NabuTrain_val)$Class
pt_Nabuval <- vect(pt_Nabuval)

#' Extract spectral vals on validation data
xy_val <- as.matrix(geom(pt_Nabuval)[,c('x','y')])
df_Nabuval <- extract(NabuAnn, xy_val)
head(df_Nabuval)
data_Nabuval <- data.frame(Class=pt_Nabuval$Class, df_Nabuval)
summary(data_Nabuval)
data_Nabuval <- na.omit(data_Nabuval)

#' Create training and validation data
NabuTrainingAnn = data_Nabutrain
NabuValidAnn = data_Nabuval

#' RF classifier
Nabu_Ann_Random_Forest_Classification <- randomForest(x=NabuTrainingAnn[-1], y=as.factor(NabuTrainingAnn$Class), ntree=500, importance=TRUE)
plot(Nabu_Ann_Random_Forest_Classification)
Nabu_Ann_Random_Forest_Classification

#' Variable importance
NabuVIP <- varImpPlot(Nabu_Ann_Random_Forest_Classification, type=1)
NabuImp <- as.data.frame(NabuVIP)
NabuImp$varnames <- rownames(NabuImp)
rownames(NabuImp) <- NULL
NabuImp
NabuVIPgg <- ggplot(NabuImp, aes(x=reorder(varnames, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseAccuracy))+
  xlab("Variable")+ylab("Mean Decrease Accuracy")+coord_flip()
NabuVIPgg

#' Prediction stats and confusion matrix
NabuPredictionAnn <- predict(Nabu_Ann_Random_Forest_Classification, newdata=NabuValidAnn)
conMatNabuAnn <- confusionMatrix(NabuPredictionAnn, as.factor(NabuValidAnn$Class))
conMatNabuAnn

#' Look at the map
NabuMapAnn <- predict(NabuAnn, Nabu_Ann_Random_Forest_Classification, filename="RF_nabu.img", type="response", 
                       index=1, na.rm=TRUE, progress="window", overwrite=TRUE)

NabuProb <- predict(NabuAnn, Nabu_Ann_Random_Forest_Classification, filename="RF_NabuProb.img", type="prob", 
                     index=1, na.rm=TRUE, progress="window", overwrite=TRUE)

plot(NabuMapAnn)
NabuMapAnnRat <- ratify(NabuMapAnn)
rat <- levels(NabuMapAnnRat)[[1]]
rat$legend <- classesdf$classnames
levels(NabuMapAnnRat) <- rat
NabuLevelMap <- levelplot(NabuMapAnnRat, scales=list(x=list(at=NULL), y=list(at=NULL)))
NabuLevelMap
NabuMapGG <- as.data.frame(NabuMapAnn, xy=TRUE) %>%
  mutate(Classification=factor(RF_nabu))
NabuMapGG


#' Image comparison
NabuAnn
nabuRGB <- ggRGB(NabuAnn, r=2, g=1, b=3, stretch = "hist")+ggsn::scalebar(x.min=2705790,
  x.max= 2709660,y.min=-2196710, y.max=-2191830,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2705790, x.max= 2709660,y.min=-2196710, y.max=-2191830, location = "bottomleft", symbol=1)+
  theme_void()
nabuRGB
grid.arrange(nabuRGB,nabuLevelMap, ncol=2)

#' Overlaid transparent classified
NabuMapGGclass <- ggRGB(NabuAnn, r=2, g=1, b=3, stretch = "hist")+geom_raster(data=NabuMapGG, aes(x=x, y=y, fill=Classification), alpha=0.5)+
  scale_fill_manual(values=viridis::viridis(2), breaks=1:2,labels=c("Irrigated","Other"))+
  ggsn::scalebar(x.min=2705790, x.max= 2709660,y.min=-2196710, y.max=-2191830,dist=1,
  transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2705790, x.max= 2709660,y.min=-2196710, y.max=-2191830, location = "bottomleft", symbol=1)+
  theme_void()
NabuMapGGclass

NabuMapGGProb <- as.data.frame(NabuProb, xy=TRUE) %>%
  mutate(Probability=RF_NabuProb)

NabuMapGGprob <- ggplot()+geom_raster(data=NabuMapGGProb, aes(x=x, y=y, fill=Probability))+
  ggsn::scalebar(x.min=2705790, x.max= 2709660,y.min=-2196710, y.max=-2191830,dist=1,
  transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2705790, x.max= 2709660,y.min=-2196710, y.max=-2191830, location = "bottomleft", symbol=1)+
  theme_void()+coord_fixed()+theme(legend.position = c(0.85,0.26),
  legend.background = element_rect(fill="white", colour="white", size=5))
NabuMapGGprob

nabuSMAD <- ggplot()+geom_raster(data=NabuAnn,aes(x=x,y=y,fill=SMAD))+ggsn::scalebar(x.min=2705790, x.max= 2709660,y.min=-2196710,
   y.max=-2191830,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2705790, x.max= 2709660,y.min=-2196710, y.max=-2191830, 
 location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.86,0.26),
   legend.background = element_rect(fill="white", colour="white", size=5))
nabuSMAD
nabuNDVI <- ggplot()+geom_raster(data=NabuAnn,aes(x=x,y=y,fill=NDVI))+ggsn::scalebar(x.min=2705790, x.max= 2709660,y.min=-2196710,
  y.max=-2191830,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2705790, x.max= 2709660,y.min=-2196710, y.max=-2191830, 
        location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.86,0.26),
legend.background = element_rect(fill="white", colour="white", size=5))
nabuNDVI
grid.arrange(nabuRGB, nabuNDVI, nabuSMAD, ncol=3)

#' Load WaPOR
Africa <- stack(".tif")
Africa

#' Crop WaPOR to area
Nabaoim <- extent(2705790,2709660,-2196710,-2191830)
Nabaoi <- extent(27,29,-17.9,-17.1)
NabuWapCrop <- crop(Africa, Nabaoi)
NabuWapPrj <- projectRaster(NabuWapCrop, crs=crs(NabuAnn))
NabuWapPrj <- crop(NabuWapPrj,Nabaoim)
NabuWapPrj

#' Allocate classes to WaPOR LCC product
NabuWapdf <- as.data.frame(NabuWapPrj, xy=TRUE)
summary(NabuWapdf)
NabuWapdf$Land_Cover_Class <- cut(NabuWapdf$layer,
                                   breaks=c(0, 15,25,35,41.48,42.51,43.5,55,65,75,80.5,85,95,111.5,112.5,114.5,
                                            115.5,116.5,121.5,122.5,123.5,124.5,125.5,150,225,Inf),
                                   labels=c("n.a","Shrubland","Grassland","Cropland, rainfed", "Cropland, irrigated", "Cropland, fallow",
                                            "Built-up", "Bare/Sparse vegetation", "Snow/ice", "water bodies", "Temp. water bodies",
                                            "Shrub cover, flooded", "Tree-cover, evergreen", "Tree-cover, evergreen broadleaved",
                                            "Tree-cover, closed deciduous", "Tree-cover, closed mixed", "Tree cover, closed unNabuown",
                                            "Tree-cover, open evergreen", "Tree cover, open evergreen broadleaf", "Tree cover, open deciduous",
                                            "Tree-cover, open deciduous broadleaf", "Tree-cover, open mixed",  "Tree cover, open unNabuown",
                                            "Sea water", "No data"))


mycols <- c("grey44", "orange3", "lawngreen", "olivedrab4", "orangered3", "gold",
            "gray10", "gray95", "azure", "royalblue2", "royalblue4",
            "paleturquoise", "seagreen1", "seagreen2",
            "seagreen3", "seagreen4", "seagreen",
            "palegreen", "palegreen1", "palegreen2",
            "palegreen3", "palegreen4",  "greenyellow",
            "navy", "gray100")
NabuWapMap <- ggplot()+geom_raster(data=NabuWapdf, aes(x=x, y=y, fill=layer))+coord_fixed()
NabuWapMap

#' Make WaPOR classes binary and plot
NabuWapdf$Classification <- cut(NabuWapdf$layer,
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
NabuWapdf$Classification <- factor(NabuWapdf$Classification , levels=c("Cropland, irrigated","Other"))
NabuGMRastMapBi <- ggRGB(NabuAnn, r=2, g=1, b=3, stretch = "hist")+
  geom_raster(data=NabuWapdf, aes(x=x, y=y, fill=Classification), alpha=0.5)+
  ggsn::scalebar(x.min=2705790,
                 x.max= 2709660,y.min=-2196710, y.max=-2191830,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2705790,
        x.max= 2709660,y.min=-2196710, y.max=-2191830, location = "bottomleft", symbol=1)+
  theme_void()+
  scale_fill_manual(values=bicol)

grid.arrange(NabuMapGGclass,NabuGMRastMapBi,ncol=2)

NabuAnn
#' Check training data confusion for  other product
NabuValidxy <- data.frame(NabuValidAnn$x, NabuValidAnn$y)
names(NabuValidxy) <- c("x", "y")

NabuWaPORconf <- extract(NabuWapPrj, xy_val, cellnumbers=TRUE)
head(NabuWaPORconf)
NabuWaPORconf <- data.frame(NabuWaPORconf)
summary(NabuWaPORconf)
NabuWaPORconf <- na.omit(NabuWaPORconf)


NabuWaPORconf$WaPOR_Class <- cut(NabuWaPORconf$layer,
                                breaks=c(0, 15,25,35,41.48,42.51,43.5,55,65,75,80.5,85,95,111.5,112.5,114.5,
                                         115.5,116.5,121.5,122.5,123.5,124.5,125.5,150,225,Inf),
                                labels=c("Other","Other","Other","Other", "Cropland, irrigated", "Other",
                                         "Other", "Other", "Other", "Other", "Other",
                                         "Other", "Other", "Other",
                                         "Other", "Other", "Other",
                                         "Other", "Other", "Other",
                                         "Other", "Other",  "Other",
                                         "Other", "Other"))
NabuWaPORconf$WaPOR_Class <- factor(NabuWaPORconf$WaPOR_Class, levels=c("Cropland, irrigated","Other"))
levels(NabuWaPORconf$WaPOR_Class)[levels(NabuWaPORconf$WaPOR_Class)=="Cropland, irrigated"]<- "Irrigated"
NabuAnnWapTab <- table(NabuWaPORconf$WaPOR_Class,NabuValidAnn$Class)
NabuAnnWapTab
#' WAPOR IS VERT COL, TEST DAT IS HORIZ ROW
conmatNabuWap <- confusionMatrix(NabuAnnWapTab)
conmatNabuWap

#' Aquastat
Astat <- raster(".asc")
projection(Astat) <- "+init=EPSG:4326"
Nabuaoi <- extent(0,100,-80,0)
NabuAstat <- crop(Astat,Nabuaoi)
plot(NabuAstat)
AstatPrj <- projectRaster(NabuAstat, crs=crs(NabuAnn))
AstatPrj
summary(AstatPrj)
Nabuaoim <- extent(2705790, 2709660,-2196710, -2191830)
NabuAstatPrj <- crop(AstatPrj, Nabuaoim)
plot(NabuAstatPrj)

#' Aquastat plot
NabuAstatdf <- as.data.frame(NabuAstatPrj, xy=TRUE)
summary(NabuAstatdf)

#' Aquastat categorise for plotting
NabuAstatdf$gmia_pct_cat <- cut(NabuAstatdf$gmia_v5_aei_pct,
                              breaks=c(-Inf,1, 5,10,20,35,50,75,100),
                              labels=c("0%","1-5%","5-10%", "10-20%", "20-35%", "35-50%",
                                       "50-75%", "75-100%"))
summary(NabuAstatdf$gmia_pct_cat)

#' Plot
Astatcols = c("yellow", "green", "orange", "purple", "aquamarine", "coral3", "darkgreen", "red")
NabuGMAstatMap <- ggRGB(NabuAnn, r=2, g=1, b=3, stretch = "hist")+ggsn::scalebar(x.min=2705790,
    x.max= 2709660,y.min=-2196710, y.max=-2191830,dist=1,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=2705790,
        x.max= 2709660,y.min=-2196710, y.max=-2191830, location = "bottomleft", symbol=1)+
  theme_void()+
  geom_raster(data=NabuAstatdf, aes(x=x, y=y, fill=gmia_pct_cat), alpha=0.5)+
  scale_fill_manual(values=Astatcols)+
  labs(fill= "Area Equipped for Irrigation (%)")
NabuGMAstatMap

#' AQUASTAT confusion matrix
NabuAstatPrj
NabuAstatconf <- extract(NabuAstatPrj, xy_val)
summary(NabuAstatconf)
NabuAstatconf <- data.frame(class=pt_Nabutrain$Class, NabuAstatconf)
summary(NabuAstatconf)
NabuAstatconf <- na.omit(NabuAstatconf)
NabuAstatconf$AstatClassification <- cut(NabuAstatconf$NabuAstatconf,
                                       breaks=c(-Inf,0.1, 5,10,20,35,50,75,100),
                                       labels=c("Other","Irrigated","Irrigated", "Irrigated", "Irrigated", "Irrigated",
                                                "Irrigated", "Irrigated"))
NabuAstatconf$class <- as.factor(NabuAstatconf$class)
NabuAstatconf$AstatClassification <- as.factor(NabuAstatconf$AstatClassification)
NabuAstatconf$AstatClassification <- factor(NabuAstatconf$AstatClassification, levels=c("Irrigated","Other"))
summary(NabuAstatconf)
NabuAstatTab <- table(NabuAstatconf$class, NabuAstatconf$AstatClassification)
NabuAstatTab
NabuAstatconMat <- confusionMatrix(NabuAstatTab)
NabuAstatconMat
