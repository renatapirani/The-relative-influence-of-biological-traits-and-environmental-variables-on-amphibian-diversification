## by Fabricius Domingos  05/2024
########################################
## 7. Random Forest script for Amphibians - Anura, Gymnophiona and Caudata 
########################################

#######################
#### Load packages ####
#######################
library(VSURF)
library(randomForest)
library(ggplot2)
library(ggpubr)
library(foreach)
library(doParallel)
library(parallel)
library(RRF)
library(graphics)
library(ape)
library(dplyr)

##############################
#### Load and verify data ####
##############################

## upload the data
setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/input")
data <- read.csv("Traits_AMPHI_test.csv")
head(data)
summary(data)
str(data)

# two NA's in a few collumns. We will have to drop them first.
amphi <- na.omit(data)
summary(amphi)
str(amphi)

## log transformed the body size and mass
amphi$BodyLength_mm <- log(amphi$BodyLength_mm)
amphi$BodyMass_g <- log(amphi$BodyMass_g)

# Subset the data frame by Order
amphi_anura <- amphi[amphi$Order == "Anura", ]
print(amphi_anura)

amphi_caudata <- amphi[amphi$Order == "Caudata", ]
print(amphi_caudata)

amphi_gymno <- amphi[amphi$Order == "Gymnophiona", ]
print(amphi_gymno)

###################################### FOR ALL AMPHIBIANS ###########################

################################################
#### Variable selection using Random Forest ####
################################################

setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/Random_Forest/All_data/")

### Using all continuous variables, because RF does not handle categorical variables
# Dependent variables is lambda_ClaDS
# Running on parallel, using 6 cores
amphi.vsurf <- VSURF(x = amphi[,-(1:5)], y = amphi[,4], 
                     nfor.thres=100,
                     nfor.interp=50,
                     nfor.pred=50,
                     mtry = 100,
                     parallel = TRUE,
                     ncores = 6)

# There was only one warning message: Warning message:
# In randomForest.default(x = x, y = y, ntree = ntree.thres, mtry = mtry,  
#    :invalid mtry: reset to within valid range

# Checking the results
amphi.vsurf
summary(amphi.vsurf)
par(las=2)
plot(amphi.vsurf, var.names=T)

# quick and dirty just to have a look
dev.copy(jpeg,filename="amphi_vsurf_Bio.jpg");
dev.off ();

# Threshold for variable elimination based on Variable Importance (VI)
amphi.vsurf$varselect.thres
plot(amphi.vsurf, step='thres', imp.sd=F, var.names=T, lty=1.5)

# Variable selection for interpretation
# These are the variables which considerably explain the data, and have the smallest error.
amphi.vsurf$varselect.interp
plot(amphi.vsurf, step='interp', imp.sd=F, var.names=T)

# Variable selection for prediction
# This is selected based on the previous variables, in a stepwise selection but using successive RF models. This step selects the most discriminant variables, with smaller correlation, and that better predicts the variable of interest.
amphi.vsurf$varselect.pred
plot(amphi.vsurf, step='pred', imp.sd=F, var.names=T)

# Variable importance for table (supplementary material, probably)
amphi.vsurf$imp.varselect.thres
amphi.vsurf$imp.mean.dec
amphi.vsurf$imp.mean.dec.ind

###########################################
#### Classic RF and Partial dependence ####
###########################################

### Importance of variables
# Only using the variables selected by VSURF

# First we will run a RRF using only the data selected by VSURF
# BodyLength_mm, Longitude, BodyMass_g,TempSeasonality
set.seed(666)
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi[,c(4,6,7,8,9,10,11)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi[,c(4,6,7,8,9,10,11)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
# Mean of squared residuals: 4.542417e-06
# % Var explained: 49.37
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#                %IncMSE     IncNodePurity
#BodyLength_mm   369.0073     1.3790771
#BodyMass_g      332.2173     1.1402160
#Longitude       462.6383     1.3295000
#TempSeasonality 353.4255     0.9178485

# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi[,c(4,6,7,10,12)][,-1],
                     amphi[,c(4,6,7,10,12)][, 1],
                      cv.fold = 10,
                      step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi[,c(4,6,7,10,12)][,-1],
        amphi[,c(4,6,7,10,12)][, 1],
        cv.fold = 10,
        step = step)
}


#summary(rrfcv.100.amphi)
(error.cv <- sapply(rrfcv.100.amphi, "[[", "error.cv"))

pdf(file="CV_amphi.pdf", width = 8, height = 8)

matplot(rrfcv.100.amphi[[1]]$n.var,
        cbind(error.cv, rowMeans(error.cv)),
        type = "l",
        lwd = c(rep(0.5, ncol(error.cv)), 2),
        col = c(rep(rgb(69, 117, 180, 100, max = 255), ncol(error.cv)), rgb(215, 48, 39, max = 255)),
        lty = 1,
        bty = "n",
        las = 1,
        axes=FALSE,
        xlab = "Number of variables",
        ylab = "Cross validation error",cex.axis=0.8)
axis(1, at = seq(0,11,1))
axis(2, at= seq(16000,30000,2000))
dev.off()        

pdf(file="Importance_amphi.pdf", width = 8, height = 8)
RRF::varImpPlot(amphi.grrf, type = 1, main = NULL, n.var = 4, cex.axis=0.8)
dev.off()        


pdf(file="Partial_amphi.pdf", width = 10, height = 3)  
par(mfrow=c(1,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines

# Longitude
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,10,12)],
                 Longitude, "lambda_ClaDS",
                 xlab= "Longitude", ylab="Lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyLength_mm
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,10,12)],
                 BodyLength_mm, "lambda_ClaDS", 
                 xlab= "Body size (mm)", ylab="",
                 cex.lab = 1, main = NULL, col="royalblue2")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BIO4
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,10,12)],
                 TempSeasonality, "lambda_ClaDS",
                 xlab= "BIO4", ylab="",
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyMass_g
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,10,12)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="",  
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

dev.off()


###################################### ORDER GYMNOPHIONA ###########################

################################################
#### Variable selection using Random Forest ####
################################################

setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/Random_Forest/Order_gymno")

### Using all continuous variables, because RF does not handle categorical variables
# Dependent variables is lambda_ClaDS
# Running on parallel, using 6 cores
amphi.vsurf <- VSURF(x = amphi_gymno[,-c(1:5)], y = amphi_gymno[,4], 
                     nfor.thres=100,
                     nfor.interp=50,
                     nfor.pred=50,
                     mtry = 100,
                     parallel = TRUE,
                     ncores = 6)

# There was only one warning message: Warning message:
# In randomForest.default(x = x, y = y, ntree = ntree.thres, mtry = mtry,  
#    :invalid mtry: reset to within valid range

# Checking the results
amphi.vsurf
summary(amphi.vsurf)
par(las=2)
plot(amphi.vsurf, var.names=T)

# quick and dirty just to have a look
dev.copy(jpeg,filename="amphi_vsurf.jpg");
dev.off ();

# Threshold for variable elimination based on Variable Importance (VI)
amphi.vsurf$varselect.thres
plot(amphi.vsurf, step='thres', imp.sd=F, var.names=T, lty=1.5)

# Variable selection for interpretation
# These are the variables which considerably explain the data, and have the smallest error.
amphi.vsurf$varselect.interp
plot(amphi.vsurf, step='interp', imp.sd=F, var.names=T)

# Variable selection for prediction
# This is selected based on the previous variables, in a stepwise selection but using successive RF models. This step selects the most discriminant variables, with smaller correlation, and that better predicts the variable of interest.
amphi.vsurf$varselect.pred
plot(amphi.vsurf, step='pred', imp.sd=F, var.names=T)

# Variable importance for table (supplementary material, probably)
amphi.vsurf$imp.varselect.thres
amphi.vsurf$imp.mean.dec
amphi.vsurf$imp.mean.dec.ind

###########################################
#### Classic RF and Partial dependence ####
###########################################

### Importance of variables
# Only using the variables selected by VSURF

# First we will run a RRF using only the data selected by VSURF
# BodyMass_g   Nocturnality     Longitude    
set.seed(666)
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi_gymno[,c(4,7,8,10)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi_gymno[,c(4,7,8,10)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
#  Mean of squared residuals: 1.074272e-07
# % Var explained:  42.3
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#                %IncMSE IncNodePurity
#BodyMass_g   101.0872   0.001291949
#Nocturnality 139.4062   0.000739566
#Longitude    133.8455   0.001457334

# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi_gymno[,c(4,7,8,10)][,-1],
                     amphi_gymno[,c(4,7,8,10)][, 1],
                     cv.fold = 10,
                     step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi_gymno[,c(4,7,8,10)][,-1],
        amphi_gymno[,c(4,7,8,10)][, 1],
        cv.fold = 10,
        step = step)
}


#summary(rrfcv.100.amphi)
(error.cv <- sapply(rrfcv.100.amphi, "[[", "error.cv"))

pdf(file="CV_amphi.pdf", width = 8, height = 8)

matplot(rrfcv.100.amphi[[1]]$n.var,
        cbind(error.cv, rowMeans(error.cv)),
        type = "l",
        lwd = c(rep(0.5, ncol(error.cv)), 2),
        col = c(rep(rgb(69, 117, 180, 100, max = 255), ncol(error.cv)), rgb(215, 48, 39, max = 255)),
        lty = 1,
        bty = "n",
        las = 1,
        axes=FALSE,
        xlab = "Number of variables",
        ylab = "Cross validation error",cex.axis=0.8)
axis(1, at = seq(0,11,1))
axis(2, at= seq(16000,30000,2000))
dev.off()        

pdf(file="Importance_amphi.pdf", width = 8, height = 8)
RRF::varImpPlot(amphi.grrf, type = 1, main = NULL, n.var = 3, cex.axis=0.8)
dev.off()        

## Partial plots
pdf(file="Partial_amphi.pdf", width = 8, height = 5)  
par(mfrow=c(2,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines

# Nocturnality
RRF::partialPlot(amphi.grrf, amphi_gymno[,c(4,7,8,10)],
                 Nocturnality, "lambda_ClaDS", 
                 xlab= "Nocturnality", ylab="Lambda ClaDS",
                 cex.lab = 1, main = NULL, col="mediumpurple")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Longitude
RRF::partialPlot(amphi.grrf, amphi_gymno[,c(4,7,8,10)],
                 Longitude, "lambda_ClaDS",
                 xlab= "Longitude", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="mediumpurple")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyMass_g
RRF::partialPlot(amphi.grrf, amphi_gymno[,c(4,7,8,10)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="mediumpurple")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

dev.off()


###################################### ORDER ANURA ###########################

################################################
#### Variable selection using Random Forest ####
################################################

setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/Random_Forest/Order_anura/")

### Using all continuous variables, because RF does not handle categorical variables
# Dependent variables is lambda_ClaDS
# Running on parallel, using 6 cores
amphi.vsurf <- VSURF(x = amphi_anura[,-c(1:5)], y = amphi_anura[,4], 
                     nfor.thres=100,
                     nfor.interp=50,
                     nfor.pred=50,
                     mtry = 100,
                     parallel = TRUE,
                     ncores = 6)

# There was only one warning message: Warning message:
# In randomForest.default(x = x, y = y, ntree = ntree.thres, mtry = mtry,  
#    :invalid mtry: reset to within valid range

# Checking the results
amphi.vsurf
summary(amphi.vsurf)
par(las=2)
plot(amphi.vsurf, var.names=T)

# quick and dirty just to have a look
dev.copy(jpeg,filename="amphi_vsurf.jpg");
dev.off ();

# Threshold for variable elimination based on Variable Importance (VI)
amphi.vsurf$varselect.thres
plot(amphi.vsurf, step='thres', imp.sd=F, var.names=T, lty=1.5)

# Variable selection for interpretation
# These are the variables which considerably explain the data, and have the smallest error.
amphi.vsurf$varselect.interp
plot(amphi.vsurf, step='interp', imp.sd=F, var.names=T)

# Variable selection for prediction
# This is selected based on the previous variables, in a stepwise selection but using successive RF models. This step selects the most discriminant variables, with smaller correlation, and that better predicts the variable of interest.
amphi.vsurf$varselect.pred
plot(amphi.vsurf, step='pred', imp.sd=F, var.names=T)

# Variable importance for table (supplementary material, probably)
amphi.vsurf$imp.varselect.thres
amphi.vsurf$imp.mean.dec
amphi.vsurf$imp.mean.dec.ind

###########################################
#### Classic RF and Partial dependence ####
###########################################

### Importance of variables
# Only using the variables selected by VSURF

# First we will run a RRF using only the data selected by VSURF
# BodyLength_mm, Longitude, BodyMass_g,TempSeasonality
set.seed(666)
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi_anura[,c(4,6,7,8,9,10,11,12)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi_anura[,c(4,6,7,8,9,10,11,12)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
#  Mean of squared residuals: 4.517875e-06
# % Var explained: 40.13
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#              %IncMSE IncNodePurity
#BodyLength_mm   263.8935     0.4372371
#BodyMass_g      272.8328     0.4924183
#Nocturnality    252.6098     0.1111768
#Verticality     372.7620     0.2081072
#Longitude       367.9320     0.5299338
#Latitude        281.2982     0.4987869
#TempSeasonality 259.1692     0.4408220


# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi_anura[,c(4,6,7,8,9,10,11,12)][,-1],
                     amphi_anura[,c(4,6,7,8,9,10,11,12)][, 1],
                     cv.fold = 10,
                     step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi_anura[,c(4,6,7,8,9,10,11,12)][,-1],
        amphi_anura[,c(4,6,7,8,9,10,11,12)][, 1],
        cv.fold = 10,
        step = step)
}


#summary(rrfcv.100.amphi)
(error.cv <- sapply(rrfcv.100.amphi, "[[", "error.cv"))

pdf(file="CV_amphi.pdf", width = 8, height = 8)

matplot(rrfcv.100.amphi[[1]]$n.var,
        cbind(error.cv, rowMeans(error.cv)),
        type = "l",
        lwd = c(rep(0.5, ncol(error.cv)), 2),
        col = c(rep(rgb(69, 117, 180, 100, max = 255), ncol(error.cv)), rgb(215, 48, 39, max = 255)),
        lty = 1,
        bty = "n",
        las = 1,
        axes=FALSE,
        xlab = "Number of variables",
        ylab = "Cross validation error",cex.axis=0.8)
axis(1, at = seq(0,11,1))
axis(2, at= seq(16000,30000,2000))
dev.off()


# Importance plot
pdf(file="Importance_amphi.pdf", width = 8, height = 8)
RRF::varImpPlot(amphi.grrf, type = 1, main = NULL, n.var = 7, cex.axis=0.8)
dev.off() 

# Partial plots
pdf(file="Partial_amphi.pdf", width = 8, height = 5)  
par(mfrow=c(2,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines

# Verticality
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,10,11,12)],
                 Verticality, "lambda_ClaDS",
                 xlab= "Verticality", ylab="Lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Longitude
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,10,11,12)],
                 Longitude, "lambda_ClaDS", 
                 xlab= "Longitude", ylab="",
                 cex.lab = 1, main = NULL, col="olivedrab2")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Latitude
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,10,11,12)],
                 Latitude, "lambda_ClaDS",
                 xlab= "Latitude", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyMass_g
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,10,11,12)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyLength_mm
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,10,11,12)],
                 BodyLength_mm, "lambda_ClaDS",
                 xlab= "Body size (mm)", ylab="Lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# TempSeasonality
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,10,11,12)],
                 TempSeasonality, "lambda_ClaDS",
                 xlab= "BIO4", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Nocturnality
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,10,11,12)],
                 Nocturnality, "lambda_ClaDS",
                 xlab= "Nocturnality", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

dev.off()


###################################### ORDER CAUDATA ###########################

################################################
#### Variable selection using Random Forest ####
################################################

setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/Random_Forest/Order_caudata/")

### Using all continuous variables, because RF does not handle categorical variables
# Dependent variables is lambda_ClaDS
# Running on parallel, using 6 cores
amphi.vsurf <- VSURF(x = amphi_caudata[,-c(1:5)], y = amphi_caudata[,4], 
                     nfor.thres=100,
                     nfor.interp=50,
                     nfor.pred=50,
                     mtry = 100,
                     parallel = TRUE,
                     ncores = 6)

# There was only one warning message: Warning message:
# In randomForest.default(x = x, y = y, ntree = ntree.thres, mtry = mtry,  
#    :invalid mtry: reset to within valid range

# Checking the results
amphi.vsurf
summary(amphi.vsurf)
par(las=2)
plot(amphi.vsurf, var.names=T)

# quick and dirty just to have a look
dev.copy(jpeg,filename="amphi_vsurf.jpg");
dev.off ();

# Threshold for variable elimination based on Variable Importance (VI)
amphi.vsurf$varselect.thres
plot(amphi.vsurf, step='thres', imp.sd=F, var.names=T, lty=1.5)

# Variable selection for interpretation
# These are the variables which considerably explain the data, and have the smallest error.
amphi.vsurf$varselect.interp
plot(amphi.vsurf, step='interp', imp.sd=F, var.names=T)

# Variable selection for prediction
# This is selected based on the previous variables, in a stepwise selection but using successive RF models. This step selects the most discriminant variables, with smaller correlation, and that better predicts the variable of interest.
amphi.vsurf$varselect.pred
plot(amphi.vsurf, step='pred', imp.sd=F, var.names=T)

# Variable importance for table (supplementary material, probably)
amphi.vsurf$imp.varselect.thres
amphi.vsurf$imp.mean.dec
amphi.vsurf$imp.mean.dec.ind

###########################################
#### Classic RF and Partial dependence ####
###########################################

### Importance of variables
# Only using the variables selected by VSURF

# First we will run a RRF using only the data selected by VSURF
set.seed(666)
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi_caudata[,c(4,7,11,12,6,10,14,9)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi_caudata[,c(4,7,11,12,6,10,14,9)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
#  Mean of squared residuals: 6.314227e-06
# % Var explained: 0.96
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#                %IncMSE IncNodePurity
# BodyLength_mm     39.30269  0.0008787689
# Nocturnality      35.13756  0.0005470041
# Latitude          27.85079  0.0008991927
# PrecipSeasonality 33.14919  0.0009860312


# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi_caudata[,c(4,7,11,12,6,10,14,9)][,-1],
                     amphi_caudata[,c(4,7,11,12,6,10,14,9)][, 1],
                     cv.fold = 10,
                     step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi_caudata[,c(4,7,11,12,6,10,14,9)][,-1],
        amphi_caudata[,c(4,7,11,12,6,10,14,9)][, 1],
        cv.fold = 10,
        step = step)
}


#summary(rrfcv.100.amphi)
(error.cv <- sapply(rrfcv.100.amphi, "[[", "error.cv"))

pdf(file="CV_amphi.pdf", width = 8, height = 8)

matplot(rrfcv.100.amphi[[1]]$n.var,
        cbind(error.cv, rowMeans(error.cv)),
        type = "l",
        lwd = c(rep(0.5, ncol(error.cv)), 2),
        col = c(rep(rgb(69, 117, 180, 100, max = 255), ncol(error.cv)), rgb(215, 48, 39, max = 255)),
        lty = 1,
        bty = "n",
        las = 1,
        axes=FALSE,
        xlab = "Number of variables",
        ylab = "Cross validation error",cex.axis=0.8)
axis(1, at = seq(0,11,1))
axis(2, at= seq(16000,30000,2000))
dev.off()        

pdf(file="Importance_amphi.pdf", width = 8, height = 8)
RRF::varImpPlot(amphi.grrf, type = 1, main = NULL, n.var = 4, cex.axis=0.8)
dev.off() 

# Partial plots
pdf(file="Partial_amphi.pdf", width = 8, height = 5)  
par(mfrow=c(2,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines


# Body mass
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,7,11,12,6,10,14,9)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="Lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Latitude
RRF::partialPlot(amphi.grrf,amphi_caudata[,c(4,7,11,12,6,10,14,9)],
                 Latitude, "lambda_ClaDS", 
                 xlab= "Latitude", ylab="",
                 cex.lab = 1, main = NULL, col="orange2")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# TempSeasonality
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,7,11,12,6,10,14,9)],
                 TempSeasonality, "lambda_ClaDS",
                 xlab= "BIO4", ylab="",  
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyLength_mm
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,7,11,12,6,10,14,9)],
                 BodyLength_mm, "lambda_ClaDS",
                 xlab= "Body size (mm)", ylab="",  
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Longitude
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,7,11,12,6,10,14,9)],
                 Longitude, "lambda_ClaDS",
                 xlab= "Longitude", ylab="Lambda ClaDS", 
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# AnnuMeanTemp
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,7,11,12,6,10,14,9)],
                 AnnuMeanTemp, "lambda_ClaDS",
                 xlab= "BIO1", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Verticality
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,7,11,12,6,10,14,9)],
                 Verticality, "lambda_ClaDS",
                 xlab= "Verticality", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

dev.off()


######end




