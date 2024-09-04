## by Renata Pirani 03/2024
#####################################
## 3. BAMM results
###################################
setwd("/Users/renatapirani/Documents/UCLA/Project/Speciation_Rate/analyses_R/output/Pyron_Phylogeny/")

library(BAMMtools)
tree <- read.tree("Amph.final.tre")
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1)

## Assessing MCMC convergence
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

## Matrix of Bayes factors. Total number of shifts = 38
computeBayesFactors(mcmcout, expectedNumberOfShifts=50, burnin=0.1)

## Graphic showing the priors results
plotPrior(mcmcout, expectedNumberOfShifts=50)

## Discart 10% of burnin
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

## Check the effective sample sizes of the log-likelihood and the number of shift events present in each sample ( ideal is to be > 200)
library(coda)
effectiveSize(postburn$N_shifts) # = 404.7757
effectiveSize(postburn$logLik) # = 273.5901 


## Calculation of extinction and Speciation rates and for all the amphibians
allrates <- getCladeRates(edata)

## mean speciation rate for the data and estimate the 90% highest posterior density 
mean(allrates$lambda) # = 0.06066985, low speciation
mean(allrates$mu) # = 0.005595196
mean(edata$meanTipMu) # = 0.006922332
quantile(allrates$lambda, c(0.05, 0.95))


## Posterior probabilities of each rate shift count observed during simulation of the posterior.
shift_probs <- summary(edata)

## Prior distribution on the number of rate shifts
css <- credibleShiftSet(edata, expectedNumberOfShifts=50, threshold=5, set.limit = 0.95)
css$number.distinct # = 8544
css$indices[[1]] # = 1003 1004
summary(css)
#plot.credibleshiftset(css, plot.max = 9)

## Finding the single best shift configuration
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=50)
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=1.5)

## Plot Rate variation through time: color density plot
plot.new()
#par(mfrow=c(1,3))
st <- max(branching.times(tree))
plotRateThroughTime(edata, intervalCol="red", avgCol="red", start.time=st, ylim=c(0,0.15), cex.axis=1)
text(x=200, y= 0.12, label="Neotropical", font=2, cex=2.0, pos=4)


## Get the tip rates using BAMM to be able to run PGLS
meanlam <- getTipRates(edata, returnNetDiv = FALSE,
                       statistic = 'mean')$lambda.avg
meanlam

# Create a vector of values for the meanlam values from BAMM
Sp_meanlam <- data.frame(lambda_BAMM = meanlam)

# Duplicate the index column so you have a new column named tiplabel
Sp_meanlam$tiplabel <- rownames(Sp_meanlam)


## Include the lamba values as a column into the Traits_Neo table and match species names
#Traits_Neo <- read.csv("Traits_NEOTROPICS.csv", header = TRUE, row.names = TRUE)
Traits_Neo <- merge(data, Sp_meanlam, by = "tiplabel")

##save the data
write.csv(Traits_Neo, file = "BAMM_amphibian.csv", row.names = TRUE)

## end



