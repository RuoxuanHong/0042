library(sf)
library(tmap)
library(tmaptools)
library(tidymodels)
library(spdep)
library(RColorBrewer)
library(here)
library(reshape)
library(lattice)
library(tidyverse)
library(ggplot2)
library(raster)
library(maptools)
library(sp)
library(readr)
library(randomForest)
library(dplyr)
library(pROC)
library(rpart)
library(caret)
library(e1071)


# Load Data---------------------------------------------------------------------
AQI <- read_csv("16stations.csv", 
                col_names = TRUE, 
                locale = locale(encoding = 'Latin1'))

AQI2019 <- AQI[,c(4,22:381)]

AQIstation_matrix<-data.matrix(AQI[,5:ncol(AQI)])
AQI2019station_matrix<-data.matrix(AQI2019[,2:ncol(AQI2019)])

# Basic Anlaysis----------------------------------------------------------------
mu = mean(AQIstation_matrix)
mu

rowMeans=rowMeans(AQIstation_matrix)


sdev = sd(AQIstation_matrix)
sdev

rowsd = apply(AQIstation_matrix,1,sd)
rowsd

me=median(AQIstation_matrix)
apply(AQIstation_matrix, 1, median)


# Distribution------------------------------------------------------------------
hist(AQIstation_matrix,main="Distribution of AQI")
abline(v=mu, col="red")
abline(v=me, col="blue")


qqnorm(AQIstation_matrix)
qqline(AQIstation_matrix, col="red")



# Shanghai AQI time series Line Chart------------------------------------------------------------

plot(colMeans(AQI2019station_matrix), xlab = "Year-Month", ylab = "AQI", type="l", xaxt="n", main="AQI Time Series 2019")
axis(1, at=c(15,46,74,106,136,167,197,227,258,287,315,345),labels=paste(c("1","2", "3", "4", "5", "6", "7", "8", "9","10","11","12"),sep=""))


# spatial Elements -----------------------------------------------------------------------
station <- paste("sta",1:nrow(AQI2019),sep="")
AQI2019<-cbind(station, AQI2019)
colnames(AQI2019)[3:ncol(AQI2019)] <- as.character(c(1:360))
AQIstation_matrix<-data.matrix(AQI2019[,3:ncol(AQI2019)])


newAQI <- melt(AQI2019, id.vars=c("District","station"), measure.vars=3:ncol(AQI2019))
colnames(newAQI)[3:4] <- c("Date", "AQI")

# every station time series analysis----------------------------------------------

xyplot(AQI ~ Date | District, xlab = "Month", type = "l",
       layout = c(4, 4),
       data=newAQI,
       scales=list(x=list(at=seq(15,360,30), labels=paste(c("1","2", "3", "4", "5", "6", "7", "8", "9","10","11","12")))),
       main ="AQI in Shanghai stations 2019" )            

# every station monthly average 2019--------------------------------------------
daymean<-colMeans(AQI2019station_matrix)

Janmean <- mean(daymean[1:31])
Febmean <- mean(daymean[32:59])
Marmean <- mean(daymean[60:90])
Aprmean <- mean(daymean[91:120])
Maymean <- mean(daymean[121:151])
Junmean <- mean(daymean[152:181])
JUlmean <- mean(daymean[182:211])
Augmean <- mean(daymean[212:242])
Sepmean <- mean(daymean[243:272])
Octmean <- mean(daymean[273:300])
Novmean <- mean(daymean[301:329])
Decmean <- mean(daymean[330:360])



Month <- c("1","2","3","4","5","6","7","8","9","10","11","12")
month<- factor(Month, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))
Monthmean <- c(Janmean,Febmean,Marmean,Aprmean,Maymean,Junmean,JUlmean,Augmean,Sepmean,Octmean,Novmean,Decmean)
studentData <- data.frame(month, Monthmean)
ggplot(studentData, mapping= aes(x=month, y=Monthmean)) + 
  geom_bar(stat="identity")



# Map period mean for every district-------------------------------------------------------------

shanghaimap <- st_read(here::here("shanghaidistricts", "Shanghai.shp"))
qtm(shanghaimap)   


AQI$periodmean <- rowMeans
Stationprofile<- AQI[c(1,2,3,4,832)]

AirqualitySHdistricts <- shanghaimap%>%
  left_join(.,
            Stationprofile, 
            by = c("PAC" = "PAC"))


tmap_mode("view")
qtm(AirqualitySHdistricts, 
    fill = "periodmean",
    borders = "grey",  
    fill.palette = "Blues" )




#Spatial Weight Matrix ---------------------------------------------------------
coordsW <- AirqualitySHdistricts%>%
  st_centroid()%>%
  st_geometry()

plot(coordsW)


knn_wards <-coordsW %>%
  knearneigh(., k=4)
LWard_knn <- knn_wards %>%
  knn2nb()
plot(LWard_knn, st_geometry(coordsW), col="blue")

Lward.knn_4_weight <- LWard_nb %>%
  nb2listw(., style="C")

# Moran's I-------------------------------------------------------------------
AQImoran <- AirqualitySHdistricts %>%
  st_drop_geometry()%>%
  dplyr::select(periodmean)%>%
  pull()%>%
  moran.test(., Lward.knn_4_weight)
AQImoran

# Geary's C--------------------------------------------------
AQIGeary <- AirqualitySHdistricts %>%
  st_drop_geometry()%>%
  dplyr::select(periodmean)%>%
  pull()%>%
  geary.test(., Lward.knn_4_weight)
AQIGeary

# Getis Ord General G---------------------------------------------------------
AQIG <- AirqualitySHdistricts %>%
  st_drop_geometry()%>%
  dplyr::select(periodmean)%>%
  pull()%>%
  globalG.test(., Lward.knn_4_weight)
AQIG



# local Moran's I-------------------------------------------------------------
AQIlocalmoran <- AirqualitySHdistricts %>%
  st_drop_geometry()%>%
  dplyr::select(periodmean)%>%
  pull()%>%
  localmoran(., Lward.knn_4_weight)%>%
  as_tibble()
AQIlocalmoran


# map local moran's statistic
local <- localmoran(x = AirqualitySHdistricts$periodmean, listw = Lward.knn_4_weight)
moran.map <- cbind(AirqualitySHdistricts, local)

tm_shape(moran.map) +
  tm_fill(col = "Ii",
          style = "quantile",
          title = "local moran's I statistic") +
  tm_borders(col="grey")




# Train set and Test set------------------------------------------------------------
Qingpu <- read_csv("QPtrain.csv", 
                   col_names = TRUE)

Qingputest <- read_csv("QPtest.csv", 
                       col_names = TRUE)

# Decision Tree---------------------------------------------------------------------

# Decision Tree on train set
set.seed(100)
model_dt = train(quality.level ~ ., data = Qingpu, method = "rpart")
model_dt_1 = predict(model_dt, data = Qingpu)
table(model_dt_1, Qingpu$quality.level)

mean(model_dt_1 == Qingpu$quality.level)


# Running on Validation Set

model_dt_vs = predict(model_dt, newdata = Qingputest)
table(model_dt_vs, Qingputest$quality.level)

mean(model_dt_vs == Qingputest$quality.level)




# Random Forest----------------------------------------------------------------------

# choose mtry
set.seed(100)
accuracy=c()
i=5
for (i in 3:8){
  model3 <- randomForest(as.factor(quality.level) ~., data = Qingpu,importance = T,ntree=1000,mtry=i)
  AQpred<-predict(model3, newdata=Qingputest)
  accuracy[i-2] = mean(AQpred == Qingputest$quality.level)
}

accuracy

plot(3:8,accuracy)


# choose ntree
set.seed(100)
r=randomForest(as.factor(quality.level)~.,data=Qingpu,mytree=5,ntree=1000,importance=TRUE,do.trace=100)

set.seed(100)
rf_train = randomForest(as.factor(quality.level) ~ .,data = Qingpu,importance = T,mtry=5,ntree=1000)
rf_train 
plot(rf_train) 

# train model
set.seed(100)
rf_train = randomForest(as.factor(quality.level) ~ .,data = Qingpu,importance = T,mtry=5,ntree=400)
rf_train 
plot(rf_train)


treesize(rf_train)
hist(treesize(rf_train))



# test dataset prediction--------------------------------------------------------

set.seed(100)
AQ_pred<-predict(rf_train, newdata=Qingputest)
cm<-table(AQ_pred,Qingputest$quality.level)
cm


# Model performance evaluation--------------------------------------------------

importance(rf_train)
varImpPlot(rf_train)


# Accuracy, Precision, Recall, F1
n = sum(cm) # number of instances
nc = nrow(cm) # number of classes
diag = diag(cm) # number of correctly classified instances per class 
rowsums = apply(cm, 1, sum) # number of instances per class
colsums = apply(cm, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes

accuracy = sum(diag) / n 
accuracy

precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1) 

macroPrecision = mean(precision)
macroRecall = mean(recall)
macroF1 = mean(f1)
data.frame(macroPrecision, macroRecall, macroF1)















