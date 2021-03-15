library(raster)
library(rgdal)
library(sp)
library(randomForest)
library(RStoolbox)
library(e1071)
library(reshape2)
library(cluster)
library(terra)
library(stringr)
library(rpart)
library(caTools)
library(fields)
library(dplyr)
library(maptools)
library(purrr) 
library(magick)
library(viridis)
library(broom)
library(caret)
library(rasterVis)
library(ggplot2)
library(gridExtra)
library(ggsn)
library(plyr)
library(sf)
library(dplyr)
library(rfinterval)
library(spatialEco)
setwd("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs")
set.seed(1)

#load the 2017 geomedian
KnGM <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs/ls8gm_tmad/Kununurra_2017_gm.tif")
ls8bands <- list("blue", "green", "red", "nir", "swir1", "swir2", "NDVI", "NDWI", "BSI", "SMAD","EMAD", "bcdev")
KnGM <- setNames(KnGM, ls8bands)
summary(KnGM)
KnPlot <- plotRGB(KnGM, r=3, g=2, b=1, stretch="hist")

#load training polygons
KnTrain <- readOGR("C:\\UserData\\wellingm\\OneDrive - Australian National University\\QGIS/Kununurra_3577.shp")
summary(KnTrain)
plot(KnTrain)
KnTrain$Class <- as.factor(KnTrain$id)
classes <- c("Irrigated","Other")
classesdf <- data.frame(classnum = c(1,2), classnames=classes)
levels(KnTrain$Class) <- c("Irrigated","Other")
KnTrain
KnTrainF <- tidy(KnTrain, region="Class")
KnTrainF

#' Divide data into training and validation
Kn_smp_size <- floor(0.8*nrow(KnTrain))
Kn_train_ind <- sample(seq_len(nrow(KnTrain)), size=Kn_smp_size)
KnTrain_tr <- KnTrain[Kn_train_ind,]
KnTrain_val <- KnTrain[-Kn_train_ind,]


KnTrainDemo <-ggRGB(KnGM, r=3, g=2, b=1, stretch = "lin")+
  geom_polygon(data=KnTrainF, aes(x=long, y=lat, group=group, fill=id), alpha=0.7)+
  labs(fill="Class")+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()
KnTrainDemo

KnTrain_trF <- tidy(KnTrain_tr, region="Class")

KnTrainDemo_tr <-ggRGB(KnGM, r=3, g=2, b=1, stretch = "lin")+
  geom_polygon(data=KnTrain_trF, aes(x=long, y=lat, group=group, fill=id), alpha=0.7)+
  labs(fill="Class")+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()
KnTrainDemo_tr

KnTrain_valF <- tidy(KnTrain_val, region="Class")

KnTrainDemo_val <-ggRGB(KnGM, r=3, g=2, b=1, stretch = "lin")+
  geom_polygon(data=KnTrain_valF, aes(x=long, y=lat, group=group, fill=id), alpha=0.7)+
  labs(fill="Class")+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()
KnTrainDemo_val

#' RStoolbox superclass NOT USED IN FINAL
#superclassKn <- superClass(KnGM, KnTrain, valData = NULL, responseCol = "id",
#                              nSamples = 100000, polygonBasedCV = TRUE, trainPartition = 0.8,
#                              model = "rf", tuneLength = 2, kfold = 5, minDist = 1,
#                              mode = "classification", predict = FALSE, predType = "raw",
#                              filename = NULL, verbose=TRUE, overwrite = TRUE)
#superclassKn

#superpredict <- predict(superclassKn, KnGM)
#superclassKn
#plot(superpredict)
#superclassimp <- varImp(superclassKn$model)
#plot(superclassimp, type = 0)
#varImpPlot(superclassKn$model, type=1)

#' Generate point samples on training data
pt_kntrain <- spsample(KnTrain_tr, 80000, type='random')
pt_kntrain$Class <- over(pt_kntrain, KnTrain_tr)$Class
pt_kntrain <- vect(pt_kntrain)

#' Extract spectral vals on training data
xy_train <- as.matrix(geom(pt_kntrain)[,c('x','y')])
df_kntrain <- extract(KnGM, xy_train)
head(df_kntrain)
data_kntrain <- data.frame(Class=pt_kntrain$Class, df_kntrain)
summary(data_kntrain)
data_kntrain <- na.omit(data_kntrain)

#' Generate point samples on validation data
pt_knval <- spsample(KnTrain_val, 20000, type='random')
pt_knval$Class <- over(pt_knval, KnTrain_val)$Class
pt_knval <- vect(pt_knval)

#' Extract spectral vals on validation data
xy_val <- as.matrix(geom(pt_knval)[,c('x','y')])
df_knval <- extract(KnGM, xy_val)
head(df_knval)
data_knval <- data.frame(Class=pt_knval$Class, df_knval)
summary(data_knval)
data_knval <- na.omit(data_knval)

#' Create training and validation data
knTraining = data_kntrain
knValid = data_knval

#' RF classifier
Random_Forest_Classification <- randomForest(x=knTraining[-1], y=as.factor(knTraining$Class), ntree=500, 
                                             importance=TRUE, probability=TRUE)
plot(Random_Forest_Classification)
Random_Forest_Classification
partialPlot(Random_Forest_Classification,knTraining,"SMAD")

#' Variable importance
KnVIP <- varImpPlot(Random_Forest_Classification, type=1)
KnImp <- as.data.frame(KnVIP)
KnImp$varnames <- rownames(KnImp)
rownames(KnImp) <- NULL
KnImp
KnVIPgg <- ggplot(KnImp, aes(x=reorder(varnames, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseAccuracy))+
  xlab("Variable")+ylab("Mean Decrease Accuracy")+coord_flip()
KnVIPgg

#' Prediction stats and confusion matrix
knPrediction <- predict(Random_Forest_Classification, newdata=knValid)
KnconMat <- confusionMatrix(knPrediction, as.factor(knValid$Class))
KnconMat
prodacc <- 6823/(6823+70)
prodacc
useracc <- 6823/(6823+68)
useracc


#' Look at the map
KnMap <- predict(KnGM, Random_Forest_Classification, filename="RF_Kn.img", type="response", 
                                  index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
KnProb <- predict(KnGM, Random_Forest_Classification, filename="RF_KnProb.img", type="prob", 
                  index=1, na.rm=TRUE, progress="window", overwrite=TRUE)

plot(KnMap)
plot(KnProb)
KnMapRat <- ratify(KnMap)
rat <- levels(KnMapRat)[[1]]
rat$legend <- classesdf$classnames
levels(KnMapRat) <- rat
KnMapRat
KnLevelMap <- levelplot(KnMapRat, scales=list(x=list(at=NULL), y=list(at=NULL)))
KnLevelMap

KnMapGG <- as.data.frame(KnMap, xy=TRUE) %>%
  mutate(Classification=factor(RF_Kn))

KnMapGGclass <- ggplot()+geom_raster(data=KnMapGG, aes(x=x, y=y, fill=Classification))+
  scale_fill_manual(values=viridis::viridis(2), breaks=1:2,labels=c("Irrigated","Other"))+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
  transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
  location = "bottomleft", symbol=1)+
  theme_void()+coord_fixed()
KnMapGGclass

KnMapGGProb <- as.data.frame(KnProb, xy=TRUE) %>%
  mutate(Probability=RF_KnProb)
KnMapGGprob <- ggplot()+geom_raster(data=KnMapGGProb, aes(x=x, y=y, fill=Probability))+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()+coord_fixed()
KnMapGGprob

KnMapGGprobtrans <- ggRGB(KnGM, r=3, g=2, b=1, stretch = "lin")+geom_raster(data=KnMapGGProb, aes(x=x, y=y, fill=Probability), alpha=0.5)+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()
KnMapGGprobtrans

KnMapGGclasstrans <- ggRGB(KnGM, r=3, g=2, b=1, stretch = "lin")+geom_raster(data=KnMapGG, aes(x=x, y=y, fill=Classification),alpha=0.5)+
  scale_fill_manual(values=viridis::viridis(2), breaks=1:2,labels=c("Irrigated","Other"))+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()
KnMapGGclasstrans
grid.arrange(KnMapGGclasstrans, KnGMAstatMap, ncol=2)

KnMapGGprobdf <- as.data.frame(KnProb, xy=TRUE)%>%
  mutate(Probability=RF_KnProb)

KnMapGGprob <- ggplot()+geom_raster(data=KnMapGGprobdf, aes(x=x, y=y, fill=Probability))+
  ggsn::scalebar(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,
                 transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700, x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()+coord_fixed()+theme(legend.position = c(0.87,0.24),
  legend.background = element_rect(fill="white", size=5, colour="white"))
KnMapGGprob

#'Image Comparison
KnRGB <- ggRGB(KnGM, r=3, g=2, b=1, stretch = "lin")+ggsn::scalebar(x.min=-360700,
  x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700,
        x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+
  theme_void()
KnRGB
grid.arrange(KnRGB,KnLevelMap, ncol=2)

KnSMAD <- ggplot()+geom_raster(data=KnGM,aes(x=x,y=y,fill=SMAD))+ggsn::scalebar(x.min=-360700,
  x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700,
        x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.88,0.28),
                                    legend.background = element_rect(fill="white", size=5, colour="white"))
KnSMAD
KnNDVI <- ggplot()+geom_raster(data=KnGM,aes(x=x,y=y,fill=NDVI))+ggsn::scalebar(x.min=-360700,
     x.max= -338525,y.min=-1673300, y.max=-1650750,dist=5,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750, 
        location = "bottomleft", symbol=1)+coord_fixed()+theme_void()+theme(legend.position = c(0.88,0.28),
            legend.background = element_rect(fill="white", size=5, colour="white"))
KnNDVI
grid.arrange(KnRGB,KnNDVI,KnSMAD,ncol=3)
grid.arrange(KnRGB,KnMapGGclass, ncol=2)

#' Aquastat
Astat <- raster("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\RS_FebIssue\\WaPOR\\gmia_v5_aei_pct.asc")
projection(Astat) <- "+init=EPSG:4326"
Knaoi <- extent(100,170,-50,0)
KnAstat <- crop(Astat,Knaoi)
plot(KnAstat)
plot(Astat)
AstatPrj <- projectRaster(KnAstat, crs=crs(KnGM))
AstatPrj

summary(AstatPrj)
Knaoim <- extent(-357000,-338525,-1673300, -1650750)
KnAstatPrj <- crop(AstatPrj, Knaoim)
plot(KnAstatPrj)

#' Aquastat plot
KnAstatdf <- as.data.frame(KnAstatPrj, xy=TRUE)
summary(KnAstatdf)

#' Aquastat categorise for plotting
KnGMLarge <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\KnDS\\KnGmLARGE.tif")
KnGMLarge
KnAstatdf$gmia_pct_cat <- cut(KnAstatdf$gmia_v5_aei_pct,
                                 breaks=c(-Inf,1, 5,10,20,35,50,75,100),
                                 labels=c("0%","1-5%","5-10%", "10-20%", "20-35%", "35-50%",
                                          "50-75%", "75-100%"))
summary(KnAstatdf$gmia_pct_cat)

#' Plot
Astatcols = c("yellow", "green", "orange", "purple", "aquamarine", "coral3", "darkgreen", "red")
KnGMAstatMap <- ggRGB(KnGM, r=1, g=2, b=3, stretch = "hist")+ggsn::scalebar(x.min=-360700,x.max= -338525,
                  y.min=-1673300, y.max=-1650750,dist=5,transform=FALSE, dist_unit="km",model="WGS84")+
  north(x.min=-360700,x.max= -338525,y.min=-1673300, y.max=-1650750, location = "bottomleft", symbol=1)+
  theme_void()+
  geom_raster(data=KnAstatdf, aes(x=x, y=y, fill=gmia_pct_cat), alpha=0.5)+
  scale_fill_manual(values=Astatcols)+
  labs(fill= str_wrap("Area Equipped for Irrigation (%)",20))
KnGMAstatMap



#' AQUASTAT confusion matrix
KnAstatPrj
KnAstatconf <- extract(KnAstatPrj, xy_val)
summary(KnAstatconf)
KnAstatconf <- data.frame(class=pt_knval$Class, KnAstatconf)
summary(KnAstatconf)
KnAstatconf <- na.omit(KnAstatconf)
KnAstatconf$AstatClassification <- cut(KnAstatconf$KnAstatconf,
                                        breaks=c(-Inf,0.1, 5,10,20,35,50,75,100),
                                        labels=c("Other","Irrigated","Irrigated", "Irrigated", "Irrigated", "Irrigated",
                                                 "Irrigated", "Irrigated"))
KnAstatconf$class <- as.factor(KnAstatconf$class)
KnAstatconf$AstatClassification <- as.factor(KnAstatconf$AstatClassification)
KnAstatconf$AstatClassification <- factor(KnAstatconf$AstatClassification, levels=c("Irrigated","Other"))
summary(KnAstatconf)
KnAstatTab <- table(KnAstatconf$class, KnAstatconf$AstatClassification)
KnAstatTab
KnAstatconMat <- confusionMatrix(KnAstatTab)
KnAstatconMat

