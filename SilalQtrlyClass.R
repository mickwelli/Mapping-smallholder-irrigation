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
setwd("C:\\UserData\\wellingm\\OneDrive - Australian National University/DEAfrica Outputs")
set.seed(1)

#load the 2019 geomedians
SilalQ1 <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\Qtry_Geomeds\\Silalatshani\\Silal2019sMADQ1.tif")
SilalQ2 <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\Qtry_Geomeds\\Silalatshani\\Silal2019sMADQ2.tif")
SilalQ3 <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\Qtry_Geomeds\\Silalatshani\\Silal2019sMADQ3.tif")
SilalQ4 <- stack("C:\\UserData\\wellingm\\OneDrive - Australian National University\\DEAfrica Outputs\\Qtry_Geomeds\\Silalatshani\\Silal2019sMADQ4.tif")
summary(SilalQ1)
s2bands <- list("green", "red", "blue", "nir", "swir_1", "swir_2", "SMAD", "NDVI", "NDWI", "BSI")
SilalQ1 <- setNames(SilalQ1, s2bands)
SilalQ1plot <- plotRGB(SilalQ1, r=2, g=1, b=3, stretch="hist")
summary(SilalQ1)
SilalQ2 <- setNames(SilalQ2, s2bands)
SilalQ3 <- setNames(SilalQ3, s2bands)
SilalQ4 <- setNames(SilalQ4, s2bands)

#StackQtrs
SilalQtrs <- stack(SilalQ1,SilalQ2,SilalQ3,SilalQ4)

#load training polygons
SilalTrain <- readOGR("C:\\UserData\\wellingm\\OneDrive - Australian National University\\QGIS\\Silalatshani.shp")
summary(SilalTrain)
plot(SilalTrain)
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
#superclassSilal <- superClass(SilalQtrs, SilalTrain, valData = NULL, responseCol = "id",
 #                          nSamples = 100000, polygonBasedCV = TRUE, trainPartition = 0.8,
  #                         model = "rf", tuneLength = 2, kfold = 5, minDist = 1,
   #                        mode = "classification", predict = FALSE, predType = "raw",
    #                       filename = NULL, verbose=TRUE, overwrite = TRUE)
#superclassSilal

#superpredictSilal <- predict(superclassSilal, SilalQtrs)
#plot(superpredictSilal)
#superclassimp <- varImp(superclassSilal$model)
#plot(superclassimp, type = 0)
#varImpPlot(superclassSilal$model, type=1)

#' Generate point samples on training data
pt_silaltrain <- spsample(SilalTrain_tr, 80000, type='random')
pt_silaltrain$Class <- over(pt_silaltrain, SilalTrain_tr)$Class
pt_silaltrain <- vect(pt_silaltrain)

#' Extract spectral vals on training data
xy_train <- as.matrix(geom(pt_silaltrain)[,c('x','y')])
df_silaltrain <- extract(SilalQtrs, xy_train)
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
df_silalval <- extract(SilalQtrs, xy_val)
head(df_silalval)
data_silalval <- data.frame(Class=pt_silalval$Class, df_silalval)
summary(data_silalval)
data_silalval <- na.omit(data_silalval)

#' Create training and validation data
SilalTraining = data_silaltrain
SilalValid = data_silalval

#' RF classifier
Silal_Random_Forest_Classification <- randomForest(x=SilalTraining[-1], y=as.factor(SilalTraining$Class), ntree=500, importance=TRUE)
plot(Silal_Random_Forest_Classification)
Silal_Random_Forest_Classification

#' Variable importance
SilalQtrVIP <- varImpPlot(Silal_Random_Forest_Classification, type=1)
SilalQtrImp <- as.data.frame(SilalQtrVIP)
SilalQtrImp$varnames <- rownames(SilalQtrImp)
rownames(SilalQtrImp) <- NULL
SilalQtrImp
SilalQtrVIPgg <- ggplot(SilalQtrImp, aes(x=reorder(varnames, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseAccuracy))+
  xlab("Variable")+ylab("Mean Decrease Accuracy")+coord_flip()
SilalQtrVIPgg

#' Prediction stats and confusion matrix
SilalPrediction <- predict(Silal_Random_Forest_Classification, newdata=SilalValid)
SilalQtrconMat <- confusionMatrix(SilalPrediction, as.factor(SilalValid$Class))
SilalQtrconMat

#' Look at the map
SilalMap <- predict(SilalQtrs, Silal_Random_Forest_Classification, filename="RF_Silal.img", type="response", 
                 index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
plot(SilalMap)

