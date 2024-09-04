## by Renata Pirani 03/2024
###############################################
## 2.  BAMM priors - Load the packages and the phylogenetic tree
###############################################

# Here I run the script for the Jetz & Pyron (2018) phylogeny

setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/BAMM/")
library(phytools)
library(ape)
library(stringr)
library(phylotools)
library(tidyverse)
library(dplyr)
library(plotrix)
library(tidyr)
library(phangorn)
library(BAMMtools)

## Load Jetz & Pyron (2018) Phylogeny
basepath <- '~/Documents/UCLA/Project/Tree_Analyses/VertLife_Supplemental_Figures_S1-12/tree/'
Amph.tree     <- paste0(basepath, "amph_shl_new_Posterior_7238.1.tre")
Amph.tree     <- read.tree(Amph.tree)

head(Amph.tree$tip.label)
is.ultrametric(Amph.tree) # FASLE but has to be TRUE
is.binary(Amph.tree) # TRUE
is.rooted(Amph.tree) # TRUE
min(Amph.tree$edge.length) # check min branch length, has to be higher than 0

## Plot the Tree in fan format
plotTree(Amph.tree,type="fan",fsize=0.1,ftype="i",lwd=c(0.5,0.5))


##  Force the tree to be ultrametric (<http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html>)
force.ultrametric<-function(Amph.tree,method=c("extend")){
  method<-method[1]
  if(method=="nnls") Amph.tree<-nnls.tree(cophenetic(Amph.tree),Amph.tree,
                                          rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(Amph.tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(Amph.tree),function(x,y) which(y==x),
               y=Amph.tree$edge[,2])
    Amph.tree$edge.length[ii]<-Amph.tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  Amph.tree
}

Amph.tree<-force.ultrametric(Amph.tree) ## default method
is.ultrametric(Amph.tree) # has to the TRUE
is.binary(Amph.tree) # has to the TRUE
is.rooted(Amph.tree) # has to the TRUE
min(Amph.tree$edge.length) # check min branch length, has to be higher than 0

## save Amph.tree to run BAMM
write.tree(Amph.tree, file = "Amph.final.tre", append = FALSE,
           digits = 7239, tree.names = FALSE)


##Setting up BAMM prior (<http://bamm-project.org/>) for Amph.tree

setBAMMpriors(read.tree("Amph.final.tre"))
#Results at the file myPriors.txt
#Prior block chosen by BAMMtools::setBAMMpriors
#expectedNumberOfShifts = 1.0
#lambdaInitPrior = 8.54206680321215
#lambdaShiftPrior = 0.00328966696212668
#muInitPrior = 8.54206680321215

##  Change these parameters in a .txt file or you can do it here and run BAMM
generateControlFile('divcontrol.txt', type = 'diversification', params = list(
  treefile = 'Amph.final.tre',
  globalSamplingFraction = '0.83',           # according to AmphibiaWeb, total of 8,729 species
  numberOfGenerations = '50000000',          # test with 50 million
  overwrite = '0',
  lambdaInitPrior = '8.54206680321215',      # parameter for the speciation rate
  lambdaShiftPrior = '0.00328966696212668',  # standard deviation of the normal distribution
  muInitPrior = '8.54206680321215',          # extinction rate
  expectedNumberOfShifts = '50'))            # since the tree has 7.238 tips, I tested with 50

### end

