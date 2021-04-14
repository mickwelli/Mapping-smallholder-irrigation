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
setwd("")
set.seed(1)

#load the 2019 geomedians
LungQ1 <- stack(".tif")
LungQ2 <- stack(".tif")
LungQ3 <- stack(".tif")
LungQ4 <- stack(".tif")
summary(LungQ1)
s2bands <- list("green", "red", "blue", "nir", "swir_1", "swir_2", "SMAD", "NDVI", "NDWI", "BSI")
LungQ1 <- setNames(LungQ1, s2bands)
LungQ1plot <- plotRGB(LungQ1, r=2, g=1, b=3, stretch="hist")
summary(LungQ1)
LungQ2 <- setNames(LungQ2, s2bands)
LungQ3 <- setNames(LungQ3, s2bands)
LungQ4 <- setNames(LungQ4, s2bands)

#StackQtrs
LungQtrs <- stack(LungQ1,LungQ2,LungQ3,LungQ4)

#load training polygons
LungTrain <- readOGR(".shp")
summary(LungTrain)
plot(LungTrain)
LungTrain$Class <- as.factor(LungTrain$id)
classes <- c("Irrigated","Other")
classesdf <- data.frame(classnum = c(1,2), classnames=classes)
levels(LungTrain$Class) <- c("Irrigated","Other")

#' Divide data into training and validation
Lung_smp_size <- floor(0.8*nrow(LungTrain))
Lung_train_ind <- sample(seq_len(nrow(LungTrain)), size=Lung_smp_size)
LungTrain_tr <- LungTrain[Lung_train_ind,]
LungTrain_val <- LungTrain[-Lung_train_ind,]

#' RStoolbox superclass
#superclassLung <- superClass(LungQtrs, LungTrain, valData = NULL, responseCol = "id",
                             # nSamples = 100000, polygonBasedCV = TRUE, trainPartition = 0.8,
                              #model = "rf", tuneLength = 2, kfold = 5, minDist = 1,
                             # mode = "classification", predict = FALSE, predType = "raw",
                             # filename = NULL, verbose=TRUE, overwrite = TRUE)
#superclassLung

#superpredictLung <- predict(superclassLung, LungQtrs)
#plot(superpredictLung)
#superclassimp <- varImp(superclassLung$model)
#plot(superclassimp, type = 0)
#varImpPlot(superclassLung$model, type=1)

#' Generate point samples on training data
pt_Lungtrain <- spsample(LungTrain_tr, 80000, type='random')
pt_Lungtrain$Class <- over(pt_Lungtrain, LungTrain_tr)$Class
pt_Lungtrain <- vect(pt_Lungtrain)

#' Extract spectral vals on training data
xy_train <- as.matrix(geom(pt_Lungtrain)[,c('x','y')])
df_Lungtrain <- extract(LungQtrs, xy_train)
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
df_Lungval <- extract(LungQtrs, xy_val)
head(df_Lungval)
data_Lungval <- data.frame(Class=pt_Lungval$Class, df_Lungval)
summary(data_Lungval)
data_Lungval <- na.omit(data_Lungval)

#' Create training and validation data
LungTraining = data_Lungtrain
LungValid = data_Lungval

#' RF classifier
Lung_Random_Forest_Classification <- randomForest(x=LungTraining[-1], y=as.factor(LungTraining$Class), 
    ntree=500, importance=TRUE)
plot(Lung_Random_Forest_Classification)
Lung_Random_Forest_Classification

#' Variable importance
LungQtrVIP <- varImpPlot(Lung_Random_Forest_Classification, type=1)
LungQtrImp <- as.data.frame(LungQtrVIP)
LungQtrImp$varnames <- rownames(LungQtrImp)
rownames(LungQtrImp) <- NULL
LungQtrImp
LungQtrVIPgg <- ggplot(LungQtrImp, aes(x=reorder(varnames, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=MeanDecreaseAccuracy))+
  xlab("Variable")+ylab("Mean Decrease Accuracy")+coord_flip()
LungQtrVIPgg

#' Prediction stats and confusion matrix
LungPrediction <- predict(Lung_Random_Forest_Classification, newdata=LungValid)
LungconMat <- confusionMatrix(LungPrediction, as.factor(LungValid$Class))
LungconMat

#' Look at the map
LungMap <- predict(LungQtrs, Lung_Random_Forest_Classification, filename="RF_Lung.img", type="response", 
                    index=1, na.rm=TRUE, progress="window", overwrite=TRUE)
plot(LungMap)
LungMap <- ratify(LungMap)
LungMap
levels(LungMap$RF_Lung) <- c("Irrigated","Other")
library(rasterVis)
LungMap <- ratify(LungMap)
rat <- levels(LungMap)[[1]]
rat$legend <- classesdf$classnames
levels(LungMap) <- rat
LungLevelMap <- levelplot(LungMap)

