## by Fabricius Domingos/ Renata Pirani  05/2024
########################################
## 7. Random Forest script for Biogeographical regions
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

##############################
#### Load and verify data ####
##############################
### upload the phylogeny
setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/BAMM_ClaDS")
tree <- read.tree("Amph.final.tre")
tree <- drop.tip(tree, "Homo_sapiens")

##upload the data
setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/input")
Bio_data <- read.csv("Traits_AMPHI_total.csv")
head(Bio_data)
summary(Bio_data)
str(Bio_data)

# remove columns
Bio_data <- subset(Bio_data, select = -AssessedStatus)
Bio_data <- subset(Bio_data, select = -lambda_BAMM)

# two NA's in a few columns. We will have to drop them first.
Bio_data <- na.omit(Bio_data)
summary(Bio_data)
str(Bio_data)

## log transformed the body size and mass
Bio_data$BodyLength_mm <- log(Bio_data$BodyLength_mm)
Bio_data$BodyMass_g <- log(Bio_data$BodyMass_g)


## Change the imputed for 0 or 1 for BIOREGIONS column
Bio_data$Neotropic [Bio_data$Neotropic > 0.5] <- 1
Bio_data$Neotropic [Bio_data$Neotropic < 0.5] <- 0

Bio_data$Afrotropic [Bio_data$Afrotropic > 0.5] <- 1
Bio_data$Afrotropic [Bio_data$Afrotropic < 0.5] <- 0

Bio_data$Australasia [Bio_data$Australasia > 0.5] <- 1
Bio_data$Australasia [Bio_data$Australasia < 0.5] <- 0

Bio_data$IndoMalay [Bio_data$IndoMalay > 0.5] <- 1
Bio_data$IndoMalay [Bio_data$IndoMalay < 0.5] <- 0

Bio_data$Nearctic [Bio_data$Nearctic > 0.5] <- 1
Bio_data$Nearctic [Bio_data$Nearctic < 0.5] <- 0

Bio_data$Oceania [Bio_data$Oceania > 0.5] <- 1
Bio_data$Oceania [Bio_data$Oceania < 0.5] <- 0

Bio_data$Palearctic [Bio_data$Palearctic > 0.5] <- 1
Bio_data$Palearctic [Bio_data$Palearctic < 0.5] <- 0

## filter the data by BIOREGIONS (IUCN) 
data_Neo <- subset(Bio_data, Bio_data$Neotropic == 1)
data_Afro <- subset(Bio_data, Bio_data$Afrotropic == 1)
data_Aus <- subset(Bio_data, Bio_data$Australasia == 1)
data_Ind <- subset(Bio_data, Bio_data$IndoMalay == 1)
data_Near <- subset(Bio_data, Bio_data$Nearctic == 1)
data_Oce <- subset(Bio_data, Bio_data$Oceania == 1) ## did not use this data, only 3 sp
data_Pale <- subset(Bio_data, Bio_data$Palearctic == 1)


## remove unnecessary columns for each df
data_Neo <- data_Neo[, !(names(data_Neo) %in% c("X","Palearctic","Oceania","Afrotropic",
                                                "Nearctic","IndoMalay","Australasia",
                                                "Neotropic"))]

data_Afro <- data_Afro[, !(names(data_Afro) %in% c("X","Palearctic","Oceania","Afrotropic",
                                                "Nearctic","IndoMalay","Australasia",
                                                "Neotropic"))]

data_Aus <- data_Aus[, !(names(data_Aus) %in% c("X","Palearctic","Oceania","Afrotropic",
                                                   "Nearctic","IndoMalay","Australasia",
                                                   "Neotropic"))]

data_Ind <- data_Ind[, !(names(data_Ind) %in% c("X","Palearctic","Oceania","Afrotropic",
                                                "Nearctic","IndoMalay","Australasia",
                                                "Neotropic"))]

data_Near <- data_Near[, !(names(data_Near) %in% c("X","Palearctic","Oceania","Afrotropic",
                                                "Nearctic","IndoMalay","Australasia",
                                                "Neotropic"))]

data_Oce <- data_Oce[, !(names(data_Oce) %in% c("X","Palearctic","Oceania","Afrotropic",
                                                   "Nearctic","IndoMalay","Australasia",
                                                   "Neotropic"))]

data_Pale <- data_Pale[, !(names(data_Pale) %in% c("X","Palearctic","Oceania","Afrotropic",
                                                "Nearctic","IndoMalay","Australasia",
                                                "Neotropic"))]


# drop tips of the tree
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% data_Neo$tiplabel])

# run PIC for each group separately 
variables <- c("lambda_ClaDS", "BodyLength_mm", "BodyMass_g", "Nocturnality",
               "Verticality", "Longitude", "Latitude", "TempSeasonality", "RangeSize",
               "AnnuMeanTemp", "AnnuPrecip", "PrecipSeasonality", "Elevation")

Neo_pics <- lapply(variables, function(var) {
  named <- data_Neo[[var]]
  names(named) <- data_Neo$tiplabel
  pic(named, tree)
})
names(Neo_pics) <- variables
Neo_pics <- data.frame(Neo_pics)


################################################
#### Variable selection using Random Forest ####
################################################

setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/Random_Forest/BIO_regions/")

### Using all continuous variables, because RF does not handle categorical variables
# Dependent variables is lambda_ClaDS
# Running on parallel, using 6 cores
BIO.vsurf <- VSURF(x = Neo_pics[,-c(1:5)], y = Neo_pics[,4], 
                     nfor.thres=100,
                     nfor.interp=50,
                     nfor.pred=50,
                     mtry = 100,
                     parallel = TRUE,
                     ncores = 6)

# Checking the results
BIO.vsurf
summary(BIO.vsurf)
par(las=2)
plot(BIO.vsurf, var.names=T)


# quick and dirty just to have a look - change the name for each region
dev.copy(jpeg,filename="Near_vsurf.jpg");
dev.off ();

# Threshold for variable elimination based on Variable Importance (VI)
BIO.vsurf$varselect.thres
plot(BIO.vsurf, step='thres', imp.sd=F, var.names=T, lty=1.5)

# Variable selection for interpretation
# These are the variables which considerably explain the data, and have the smallest error.
BIO.vsurf$varselect.interp
plot(BIO.vsurf, step='interp', imp.sd=F, var.names=T)

# Variable selection for prediction
# This is selected based on the previous variables, in a stepwise selection but using successive RF models. This step selects the most discriminant variables, with smaller correlation, and that better predicts the variable of interest.
BIO.vsurf$varselect.pred
plot(BIO.vsurf, step='pred', imp.sd=F, var.names=T)

# Variable importance for table (supplementary material, probably)
BIO.vsurf$imp.varselect.thres
BIO.vsurf$imp.mean.dec
BIO.vsurf$imp.mean.dec.ind

###########################################
#### Classic RF and Partial dependence ####
###########################################

### Importance of variables
# Only using the variables selected by VSURF

# First we will run a RRF using only the data selected by VSURF
set.seed(666)
BIO.RRF <- RRF(lambda_ClaDS ~ ., flagReg=0, data=data_Near[,c(4,14,6,13,16,10,11,7,8,9)])
BIO.RRF
plot(BIO.RRF)
gamma <- 0.5
coefReg.estimated <- (1-gamma)+(gamma*(BIO.RRF$importance/(max(BIO.RRF$importance))))

# Now, run the GRRF on the coefficient and find the final importance of the predictors
BIO.grrf <- RRF(lambda_ClaDS ~ ., flagReg=1, data=data_Near[,c(4,14,6,13,16,10,11,7,8,9)], coefReg=coefReg.estimated, importance=T, ntree=10000, type=1)
BIO.grrf
# Mean of squared residuals: 0.0003893558
# % Var explained: 55.34
plot(BIO.grrf)

RRF::importance(BIO.grrf)
#              %IncMSE IncNodePurity
# BodyLength_mm 279.0638     0.7086939
# BodyMass_g    272.6015     0.6314494
# Verticality   236.6533     0.1914988
# Longitude     247.3908     0.5645824


# Cross-validation
set.seed(667)
step <-  1 - (1 / 10)

# testing the script
rrfcv.BIO <- rrfcv(data_Near[,c(4,14,6,13,16,10,11,7,8,9)][,-1],
                   data_Near[,c(4,14,6,13,16,10,11,7,8,9)][, 1],
                     cv.fold = 10,
                     step = step)
summary(rrfcv.BIO)


# 100 replicates of cross-validation - estou fazendo com 10
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
trials <- 10
rrfcv.10.BIO <- foreach(icount(trials), .packages="RRF") %do% {
  rrfcv(data_Near[,c(4,14,6,13,16,10,11,7,8,9)][,-1],
        data_Near[,c(4,14,6,13,16,10,11,7,8,9)][, 1],
        cv.fold = 10,
        step = step)
}


#summary(rrfcv.100.amphi)
(error.cv <- sapply(rrfcv.10.BIO, "[[", "error.cv"))

pdf(file="Near_CV.pdf", width = 8, height = 8)

matplot(rrfcv.10.BIO[[1]]$n.var,
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


# Importance plot - change the region names
pdf(file="Near_Importance_amphi.pdf", width = 8, height = 8)
RRF::varImpPlot(BIO.grrf, type = 1, main = NULL, n.var = 9, cex.axis=0.8)
dev.off() 

# Partial plots
pdf(file="Near_Partial_amphi.pdf", width = 8, height = 7)  
par(mfrow=c(3,4), bty='l')  # Setting layout to one row, three columns, and remove top/right box lines

# BodyLength_mm
RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 BodyLength_mm, "lambda_ClaDS", 
                 xlab= "Body size (mm)", ylab="lambda_ClaDS",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BIO1
RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 AnnuMeanTemp, "lambda_ClaDS", 
                 xlab= "BIO1", ylab="",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines


# Range size
RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 RangeSize, "lambda_ClaDS", 
                 xlab= "Range size", ylab="",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

# BodyLength_mm
RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 Latitude, "lambda_ClaDS", 
                 xlab= "Latitude", ylab="",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines




# Verticality
RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 Verticality, "lambda_ClaDS", 
                 xlab= "Verticality", ylab="lambda_ClaDS",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines


RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 Longitude, "lambda_ClaDS", 
                 xlab= "Longitude", ylab="",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines


RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 Nocturnality, "lambda_ClaDS", 
                 xlab= "Nocturnality", ylab="",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 PrecipSeasonality, "lambda_ClaDS", 
                 xlab= "BIO15", ylab="",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

RRF::partialPlot(BIO.grrf, data_Near[,c(4,14,6,13,16,10,11,7,8,9)],
                 BodyMass_g, "lambda_ClaDS", 
                 xlab= "Body mass", ylab="lambda_ClaDS",
                 cex.lab = 1, main = NULL, col="black")
box(which = "plot", bty = "l")  # Add only left and bottom box lines

dev.off()


### end





