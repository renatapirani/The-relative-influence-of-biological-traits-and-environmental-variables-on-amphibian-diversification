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
data <- read.csv("TetrapodTraits_v2.0.1.csv") # I changed it for the figure BAMM
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
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi[,c(4,6,7,8,9,11,12,13,15)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

#% Var explained: 50.1

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi[,c(4,6,7,8,9,11,12,13,15)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
# Mean of squared residuals: 0.0003548056
# % Var explained: 50.45
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#                %IncMSE IncNodePurity
#BodyLength_mm   409.4785     1.0997795
#BodyMass_g      370.3333     1.0191138
#Nocturnality    215.8147     0.1776580
#Verticality     237.0232     0.2682718
#AnnuMeanTemp    238.2031     0.5854794
#AnnuPrecip      227.2837     0.4849366
#TempSeasonality 234.9960     0.6412888
#Elevation       215.6944     0.5106816

# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi[,c(4,6,7,8,9,11,12,13,15)][,-1],
                     amphi[,c(4,6,7,8,9,11,12,13,15)][, 1],
                      cv.fold = 10,
                      step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi[,c(4,6,7,8,9,11,12,13,15)][,-1],
        amphi[,c(4,6,7,8,9,11,12,13,15)][, 1],
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
RRF::varImpPlot(amphi.grrf, type = 1, main = NULL, n.var = 8, cex.axis=0.8)
dev.off()        


pdf(file="Partial_amphi.pdf", width = 10, height = 6)  
par(mfrow=c(2,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines

# BodyLength_mm
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 BodyLength_mm, "lambda_ClaDS", 
                 xlab= "Body size (mm)", ylab="Lambda ClaDS",
                 cex.lab = 1, main = NULL, col="royalblue2")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyMass_g
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="",  
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BIO1 - AnnuMeanTemp
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 AnnuMeanTemp, "lambda_ClaDS",
                 xlab= "BIO1", ylab="",
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Verticality
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 Verticality, "lambda_ClaDS",
                 xlab= "Verticality", ylab="",  
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BIO4 -  TempSeasonality
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 TempSeasonality, "lambda_ClaDS",
                 xlab= "BIO4", ylab="Lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BIO12 - AnnuPrep
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 AnnuPrecip, "lambda_ClaDS",
                 xlab= "BIO12", ylab="",  
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Elevation
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 Elevation, "lambda_ClaDS",
                 xlab= "Elevation", ylab="",  
                 cex.lab = 1, main = NULL, col="royalblue2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Nocturnality
RRF::partialPlot(amphi.grrf, amphi[,c(4,6,7,8,9,11,12,13,15)],
                 Nocturnality, "lambda_ClaDS",
                 xlab= "Nocturnality", ylab="",  
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
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi_gymno[,c(4,6,7,8,15)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi_gymno[,c(4,6,7,8,15)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
#  Mean of squared residuals: 2.157819e-05
# % Var explained: 28.1
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#                %IncMSE IncNodePurity
#BodyLength_mm 102.22889   0.001042046
#BodyMass_g     98.68141   0.001228991
#Nocturnality   91.94154   0.000468940
#Elevation      88.40263   0.001210897

# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi_gymno[,c(4,6,7,8,15)][,-1],
                     amphi_gymno[,c(4,6,7,8,15)][, 1],
                     cv.fold = 10,
                     step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi_gymno[,c(4,6,7,8,15)][,-1],
        amphi_gymno[,c(4,6,7,8,15)][, 1],
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


## Partial plots
pdf(file="Partial_amphi.pdf", width = 10, height = 6)  
par(mfrow=c(2,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines

# Nocturnality
RRF::partialPlot(amphi.grrf, amphi_gymno[,c(4,6,7,8,15)],
                 Nocturnality, "lambda_ClaDS", 
                 xlab= "Nocturnality", ylab="Lambda ClaDS",
                 cex.lab = 1, main = NULL, col="mediumpurple")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Elevation
RRF::partialPlot(amphi.grrf, amphi_gymno[,c(4,6,7,8,15)],
                 Elevation, "lambda_ClaDS",
                 xlab= "Elevation", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="mediumpurple")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyMass_g
RRF::partialPlot(amphi.grrf, amphi_gymno[,c(4,6,7,8,15)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="mediumpurple")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# Bodysize
RRF::partialPlot(amphi.grrf, amphi_gymno[,c(4,6,7,8,15)],
                 BodyLength_mm, "lambda_ClaDS",
                 xlab= "Body size (mm)", ylab="",  # No y-axis label
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
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
#Mean of squared residuals: 0.0003250946
#% Var explained: 31.21
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#              %IncMSE IncNodePurity
#odyLength_mm     212.2647    0.40550045
#BodyMass_g        230.0653    0.41664470
#Nocturnality      175.2833    0.09734054
#Verticality       310.8631    0.17795781
#AnnuMeanTemp      218.9060    0.34070476
#AnnuPrecip        217.8600    0.31644102
#TempSeasonality   321.7548    0.44685341
#PrecipSeasonality 207.9317    0.30842441
#Elevation         188.8207    0.30224683


# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)][,-1],
                     amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)][, 1],
                     cv.fold = 10,
                     step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)][,-1],
        amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)][, 1],
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
RRF::varImpPlot(amphi.grrf, type = 1, main = NULL, n.var = 9, cex.axis=0.8)
dev.off() 

# Partial plots
pdf(file="Partial_amphi.pdf", width = 10, height = 7)  
par(mfrow=c(3,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines

# Verticality
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 Verticality, "lambda_ClaDS",
                 xlab= "Verticality", ylab="Lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# TempSeasonality - BIO4
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 TempSeasonality, "lambda_ClaDS", 
                 xlab= "BIO4", ylab="",
                 cex.lab = 1, main = NULL, col="olivedrab2")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyLength_mm
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 BodyLength_mm, "lambda_ClaDS",
                 xlab= "Body size (mm)", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyMass_g
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# AnnuPrecip - 
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 AnnuPrecip, "lambda_ClaDS",
                 xlab= "BIO12", ylab="lambda ClaDS",  # No y-axis label
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


#AnnuMeanTemp
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 AnnuMeanTemp, "lambda_ClaDS",
                 xlab= "BIO1", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# PrecipSeasonality
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 PrecipSeasonality, "lambda_ClaDS",
                 xlab= "BIO15", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# Elevation
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 Elevation, "lambda_ClaDS",
                 xlab= "Elevation", ylab="",  
                 cex.lab = 1, main = NULL, col="olivedrab2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# Nocturnality
RRF::partialPlot(amphi.grrf, amphi_anura[,c(4,6,7,8,9,11,12,13,14,15)],
                 Nocturnality, "lambda_ClaDS",
                 xlab= "Nocturnality", ylab="lambda ClaDS",  
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
amphi.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=amphi_caudata[,c(4,6,7,9,11,12,13,14)])
amphi.RRF
plot(amphi.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(amphi.RRF$importance/(max(amphi.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
amphi.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=amphi_caudata[,c(4,6,7,9,11,12,13,14)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
amphi.grrf
#MMean of squared residuals: 0.0003903136
# % Var explained: 77.73
plot(amphi.grrf)

RRF::importance(amphi.grrf)
#                %IncMSE IncNodePurity
#BodyLength_mm     118.2643    0.12216132
#BodyMass_g        179.5174    0.20460553
#Verticality       106.2226    0.03154278
#AnnuMeanTemp      147.5114    0.24531886
#AnnuPrecip        119.6747    0.10291116
#TempSeasonality   237.1378    0.33236880
#PrecipSeasonality 120.2694    0.08213524


# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.amphi <- rrfcv(amphi_caudata[,c(4,6,7,9,11,12,13,14)][,-1],
                     amphi_caudata[,c(4,6,7,9,11,12,13,14)][, 1],
                     cv.fold = 10,
                     step = step)
summary(rrfcv.amphi)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.100.amphi <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(amphi_caudata[,c(4,6,7,9,11,12,13,14)][,-1],
        amphi_caudata[,c(4,6,7,9,11,12,13,14)][, 1],
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
RRF::varImpPlot(amphi.grrf, type = 1, main = NULL, n.var = 7, cex.axis=0.8)
dev.off() 

# Partial plots
pdf(file="Partial_amphi.pdf", width = 10, height = 6)  
par(mfrow=c(2,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines


# TempSeasonality - BIO4
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,6,7,9,11,12,13,14)],
                 TempSeasonality, "lambda_ClaDS",
                 xlab= "BIO4", ylab="Lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# Body mass
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,6,7,9,11,12,13,14)],
                 BodyMass_g, "lambda_ClaDS",
                 xlab= "Body mass (g)", ylab="",  
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# AnnuMeanTemp
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,6,7,9,11,12,13,14)],
                 AnnuMeanTemp, "lambda_ClaDS",
                 xlab= "BIO1", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# PrecipSeasonality
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,6,7,9,11,12,13,14)],
                 PrecipSeasonality, "lambda_ClaDS",
                 xlab= "BIO15", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# BodyLength_mm
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,6,7,9,11,12,13,14)],
                 BodyLength_mm, "lambda_ClaDS",
                 xlab= "Body size (mm)", ylab="lambda ClaDS",  
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# AnnuPrecip
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,6,7,9,11,12,13,14)],
                 AnnuPrecip, "lambda_ClaDS",
                 xlab= "BIO12", ylab="", 
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# Verticality
RRF::partialPlot(amphi.grrf, amphi_caudata[,c(4,6,7,9,11,12,13,14)],
                 Verticality, "lambda_ClaDS",
                 xlab= "Verticality", ylab="",  # No y-axis label
                 cex.lab = 1, main = NULL, col="orange2")  # Remove yaxt='n'
box(which = "plot", bty = "l")  # Add only left and bottom box lines

dev.off()


######end




