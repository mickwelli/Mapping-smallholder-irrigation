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
library(ggplot2)
library(rasterVis)
setwd("")
set.seed(1)

#load the 2019 geomedians
NabuQ1 <- stack("")
NabuQ2 <- stack("")
NabuQ3 <- stack(".tif")
NabuQ4 <- stack("tif")
summary(NabuQ1)
s2bands <- list("green", "red", "blue", "nir", "swir_1", "swir_2", "SMAD", "NDVI", "NDWI", "BSI")
NabuQ1 <- setNames(NabuQ1, s2bands)
NabuQ1plot <- plotRGB(NabuQ1, r=2, g=1, b=3, stretch="hist",axes=TRUE)
ggRGB(NabuQ1, r=7, g=NULL, b=NULL)
summary(NabuQ1)
NabuQ2 <- setNames(NabuQ2, s2bands)
NabuQ3 <- setNames(NabuQ3, s2bands)
NabuQ4 <- setNames(NabuQ4, s2bands)

#StackQtrs
NabuQtrs <- stack(NabuQ1,NabuQ2,NabuQ3,NabuQ4)

#load training polygons
NabuTrain <- readOGR(".shp")
summary(NabuTrain)
plot(NabuTrain)
NabuTrain$Class <- as.factor(NabuTrain$id)
classes <- c("Irrigated","Other")
classesdf <- data.frame(classnum = c(1,2), classnames=classes)
levels(NabuTrain$Class) <- c("Irrigated","Other")

#' Divide data into training and validation
Nabu_smp_size <- floor(0.8*nrow(NabuTrain))
Nabu_train_ind <- sample(seq_len(nrow(NabuTrain)), size=Nabu_smp_size)
NabuTrain_tr <- NabuTrain[Nabu_train_ind,]
NabuTrain_val <- NabuTrain[-Nabu_train_ind,]

#' RStoolbox superclass NORMALLY SKIP THIS SECTION (it is an alternative rf modelling method)
#superclassNabu <- superClass(NabuQtrs, NabuTrain, valData = NULL, responseCol = "id",
 #                            nSamples = 100000, polygonBasedCV = TRUE, trainPartition = 0.8,
  #                           model = "rf", tuneLength = 2, kfold = 5, minDist = 1,
   #                          mode = "classification", predict = FALSE, predType = "raw",
    #                         filename = NULL, verbose=TRUE, overwrite = TRUE)
#superclassNabu

#superpredictNabu <- predict(superclassNabu, NabuQtrs)
#plot(superpredictNabu)
#superclassimp <- varImp(superclassNabu$model)
#plot(superclassimp, type = 0)
#varImpPlot(superclassNabu$model, type=1)

#' Generate point samples on training data
pt_Nabutrain <- spsample(NabuTrain_tr, 80000, type='random')
pt_Nabutrain$Class <- over(pt_Nabutrain, NabuTrain_tr)$Class
pt_Nabutrain <- vect(pt_Nabutrain)

#' Extract spectral vals on training data
xy_train <- as.matrix(geom(pt_Nabutrain)[,c('x','y')])
df_Nabutrain <- extract(NabuQtrs, xy_train)
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
df_Nabuval <- extract(NabuQtrs, xy_val)
head(df_Nabuval)
data_Nabuval <- data.frame(Class=pt_Nabuval$Class, df_Nabuval)
summary(data_Nabuval)
data_Nabuval <- na.omit(data_Nabuval)

#' Create training and validation data
NabuTraining = data_Nabutrain
NabuValid = data_Nabuval

#' RF classifier
Nabu_Random_Forest_Classification <- randomForest(x=NabuTraining[-1], y=as.factor(NabuTraining$Class), ntree=500, importance=TRUE)
plot(Nabu_Random_Forest_Classification)
Nabu_Random_Forest_Classification

#' Variable importance
NabuQtrVIP <- varImpPlot(Nabu_Random_Forest_Classification, type=1)
NabuQtrImp <- as.data.frame(NabuQtrVIP)
NabuQtrImp$varnames <- rownames(NabuQtrImp)
rownames(NabuQtrImp) <- NULL
NabuQtrImp
NabuQtrVIPgg <- ggplot(NabuQtrImp, aes(x=reorder(varnames, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseAccuracy))+
  xlab("Variable")+ylab("Mean Decrease Accuracy")+coord_flip()
NabuQtrVIPgg

grid.arrange(NabuQtrVIPgg, NabuQtrVIPgg, ncol=2)
grid.arrange(LungQtrVIPgg, ncol=2)
#' Prediction stats and confusion matrix
NabuPrediction <- predict(Nabu_Random_Forest_Classification, newdata=NabuValid)
NabuconMat <- confusionMatrix(NabuPrediction, as.factor(NabuValid$Class))
NabuconMat

#' Look at the map
NabuMap <- predict(NabuQtrs, Nabu_Random_Forest_Classification, filename="RF_Nabu.img", type="response", 
                   index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
plot(NabuMap)

NabuMap <- ratify(NabuMap)
rat <- levels(NabuMap)[[1]]
rat$legend <- classesdf$classnames
levels(NabuMap) <- rat
NabuLevelMap <- levelplot(NabuMap)
NabuLevelMap
