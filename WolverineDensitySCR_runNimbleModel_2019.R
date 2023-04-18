################################################################################
##### ----- SPATIAL DETERMINANTS OF WOLVERINE DENSITY IN SCANDINAVIA ----- ##### 
#####                 MONITORING SEASON: DEC 2018 - JUN 2019               #####
#####                  SPATIAL CAPTURE-RECAPTURE ANALYSIS                  #####
################################################################################

## Clean
rm(list = ls())
cat("\014")
gc()

## Library
library(nimble)
library(nimbleSCR) 


## ----------------------------------------------------------------------------------------------
## ------ I. LOAD NIMBLE INPUT FILES -----
## ----------------------------------------------------------------------------------------------

setwd("") #--define to read/source the required files

## Load a custom nimble function
source("dbin_LESS_Cached_MultipleCovResponse.R")

## Load the input data for each sex
# Female
load("SCR_Female_2019_1.RData")
# # Male
# load("SCR_Male_2019_1.RData")


## -----------------------------------------------------------------------------
## ------ II. MODEL SETTING AND RUNNING ------- 
## -----------------------------------------------------------------------------

### ==== 1. NIMBLE MODEL DEFINITION ====
modelCode <- nimbleCode({
  
  ##----------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##
  
  ## Priors on covariate effects on density
  for (h in 1:n.habCovs) {  
    betaHabCovs[h] ~ dnorm(0.0, 0.01)
  }#h
  
  habIntensity[1:numHabWindows] <- exp(
    
    ## Zone 1 is an "implied" intercept
      betaHabCovs[1]  * habCovs[1:numHabWindows, 1] +  #--Zone 2
      betaHabCovs[2]  * habCovs[1:numHabWindows, 2] +  #--Zone 3
      betaHabCovs[3]  * habCovs[1:numHabWindows, 3] +  #--Zone 4
      
      betaHabCovs[4]  * habCovs[1:numHabWindows, 4] +  #--Distance from the relict range in Zone 1
      
      betaHabCovs[5]  * habCovs[1:numHabWindows, 5] + #--TRI
      betaHabCovs[6]  * habCovs[1:numHabWindows, 6] + #--Snow
      betaHabCovs[7]  * habCovs[1:numHabWindows, 7] + #--Forest
      betaHabCovs[8]  * habCovs[1:numHabWindows, 8] + #--Human settlements
      betaHabCovs[9]  * habCovs[1:numHabWindows, 9] + #--Moose
      
      ## Interaction btw distance from the relic range and zones
      betaHabCovs[10] * habCovs[1:numHabWindows, 4] * habCovs[1:numHabWindows, 1] +
      betaHabCovs[11] * habCovs[1:numHabWindows, 4] * habCovs[1:numHabWindows, 2] +
      betaHabCovs[12] * habCovs[1:numHabWindows, 4] * habCovs[1:numHabWindows, 3] 
  )
  #-----------------------------------------------------------------------------
  
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for (i in 1:n.individuals) {
    
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max
    )
    
  }#i
  
  
  ##----------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##
  
  psi ~ dunif(0,1)
  probDetBefore ~ dunif(0,1)
  for (i in 1:n.individuals) {
    z[i] ~ dbern(psi)
    detResponse[i] ~ dbern(probDetBefore) #--Individual-level covariate on detection
  }#i
  
  
  ##----------------------------------------------------------------------------
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------##
  
  sigma ~ dunif(0,4)
  
  ## Priors on covariate effects on detection probability
  for (c in 1:n.covs) {
    betaCovs[c] ~ dunif(-5,5)
  }#c
  
  betaResponse ~ dunif(-5,5)
  
  ## Management region specific intercept for detection probability
  for (c in 1:n.countries) {
    p0[c] ~ dunif(0,1)
  }#c  
  
  for (i in 1:n.individuals) {
    y.alive[i, 1:nMaxDetectors] ~ dbin_LESS_Cached_MultipleCovResponse(sxy = sxy[i, 1:2],
                                                                       sigma = sigma,
                                                                       nbDetections[i],
                                                                       yDets = yDets[i, 1:nMaxDetectors],
                                                                       detector.xy =  detector.xy[1:n.detectors, 1:2],
                                                                       trials = trials[1:n.detectors],
                                                                       detectorIndex = detectorIndex[1:n.cellsSparse, 1:maxNBDets],
                                                                       nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                       ResizeFactor = ResizeFactor,
                                                                       maxNBDets = maxNBDets,
                                                                       habitatID = habitatIDDet[1:y.maxDet, 1:x.maxDet],
                                                                       indicator = z[i],
                                                                       p0[1:n.countries],
                                                                       detCountries[1:n.detectors],
                                                                       detCov = detCovs[1:n.detectors, 1:n.covs],
                                                                       betaCov = betaCovs[1:n.covs],
                                                                       BetaResponse = betaResponse,
                                                                       detResponse = detResponse[i]
    )
    
  }#i
  
  
  ##----------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  
  N <- sum(z[1:n.individuals])
  
})


### ----------------------------------------------------------------------------
nimParams <- c("N", "psi","betaHabCovs", 
               "p0", "sigma", "betaCovs","betaResponse","probDetBefore",
               "habIntensity","z", "sxy")



## -----------------------------------------------------------------------------
## ------ III. MODEL FITTING ------- WITHOUT VARIABLE SELECTION
## -----------------------------------------------------------------------------

## Create nimble model
model <- nimbleModel(code = modelCode,
                     constants = nimConstants,
                     inits = nimInits,
                     data = nimData,
                     check = FALSE,
                     calculate = FALSE) 
model$calculate()

## Compile and Run MCMC 
conf <- configureMCMC(model,
                      monitors = nimParams,
                      print = FALSE)

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(Rmcmc, project = model)
MCMC_runtime <- system.time(
  ## A very short test run!
  samples <- runMCMC(Cmcmc,
                     niter = 100, 
                     nburnin = 10, 
                     nchains = 2, 
                     samplesAsCodaMCMC = TRUE)
)



## -----------------------------------------------------------------------------
## ------ IV. MODEL FITTING ------- WITH VARIABLE SELECTION: RJMCMC
## -----------------------------------------------------------------------------

## Create nimble objects
lmIndicatorCode <- nimbleCode({
  
  ##----------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##
  
  ## Priors on covariate effects on density
  for (h in 1:n.habCovs) {
    betaHabCovs[h] ~ dunif(-5, 5)
  }#h
  
  for (v in 1:n.habCovs) {
    w[v] ~ dbern(omega)  #--Indicator variable for each coefficient
  }#v
  omega ~ dbeta(2,8)     #--Prior on inclusion probability
  # This can be considered an informative prior; dbeta(2,8) puts an average prior probability of 0.2 on retaining each covariate (with SD = 0.12)
  # alternatively, consider dunif(0,1), a non-informative prior, which assigns an average prior probability of 0.5 (with higher SD = 0.29)
  
  
  ## Zone 1 is an "implied" intercept
  habIntensity[1:numHabWindows] <- exp(
    
      betaHabCovs[1]  * w[1] * habCovs[1:numHabWindows, 1] +
      betaHabCovs[2]  * w[2] * habCovs[1:numHabWindows, 2] +
      betaHabCovs[3]  * w[3] * habCovs[1:numHabWindows, 3] +
      
      betaHabCovs[4]  * w[4] * habCovs[1:numHabWindows, 4] +
      
      betaHabCovs[5]  * w[5] * habCovs[1:numHabWindows, 5] +
      betaHabCovs[6]  * w[6] * habCovs[1:numHabWindows, 6] +
      betaHabCovs[7]  * w[7] * habCovs[1:numHabWindows, 7] +
      betaHabCovs[8]  * w[8] * habCovs[1:numHabWindows, 8] +
      betaHabCovs[9]  * w[9] * habCovs[1:numHabWindows, 9] +
      
      ## Interaction btw distance from the relic range and zones
      betaHabCovs[10] * w[10] * habCovs[1:numHabWindows, 4] * habCovs[1:numHabWindows, 1] +
      betaHabCovs[11] * w[11] * habCovs[1:numHabWindows, 4] * habCovs[1:numHabWindows, 2] +
      betaHabCovs[12] * w[12] * habCovs[1:numHabWindows, 4] * habCovs[1:numHabWindows, 3] 
  )
  
  ## Constraints: interaction terms are included only if the additive effects are included
  c1 ~ dconstraint(w[1] >= w[10] & w[4] >= w[10])
  c2 ~ dconstraint(w[2] >= w[11] & w[4] >= w[11])
  c3 ~ dconstraint(w[3] >= w[12] & w[4] >= w[12])
  
  
  #-----------------------------------------------------------------------------
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for (i in 1:n.individuals) {
    sxy[i, 1:2] ~ dbernppAC(lowerCoords = lowerHabCoords[1:numHabWindows, 1:2], 
                            upperCoords = upperHabCoords[1:numHabWindows, 1:2], 
                            logIntensities = logHabIntensity[1:numHabWindows], 
                            logSumIntensity = logSumHabIntensity,
                            habitatGrid = habitatGrid[1:y.max, 1:x.max], 
                            numGridRows = y.max, numGridCols = x.max)
  }#i
  
  
  ##----------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##
  
  psi ~ dunif(0, 1)
  probDetBefore ~ dunif(0,1)
  
  for (i in 1:n.individuals) {
    z[i] ~ dbern(psi)
    detResponse[i] ~ dbern(probDetBefore)
  }#i
  
  
  ##----------------------------------------------------------------------------
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------##
  
  sigma ~ dunif(0, 4)
  for (c in 1:n.covs) {
    betaCovs[c] ~ dunif(-5, 5)
  }#c
  
  betaResponse ~ dunif(-5, 5)
  for (c in 1:n.countries) {
    p0[c] ~ dunif(0, 1)
  }#c
  
  for (i in 1:n.individuals) {
    y.alive[i, 1:nMaxDetectors] ~ dbin_LESS_Cached_MultipleCovResponse(sxy = sxy[i, 1:2], 
                                                                       sigma = sigma, 
                                                                       nbDetections = nbDetections[i], 
                                                                       yDets = yDets[i, 1:nMaxDetectors], 
                                                                       detector.xy = detector.xy[1:n.detectors, 1:2], 
                                                                       trials = trials[1:n.detectors], 
                                                                       detectorIndex = detectorIndex[1:n.cellsSparse, 1:maxNBDets], 
                                                                       nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse], 
                                                                       ResizeFactor = ResizeFactor, maxNBDets = maxNBDets, 
                                                                       habitatID = habitatIDDet[1:y.maxDet, 1:x.maxDet], 
                                                                       indicator = z[i], p0State = p0[1:n.countries], 
                                                                       detCountries = detCountries[1:n.detectors], 
                                                                       detCov = detCovs[1:n.detectors, 1:n.covs], 
                                                                       betaCov = betaCovs[1:n.covs], 
                                                                       BetaResponse = betaResponse, 
                                                                       detResponse = detResponse[i])
  }#i
  
  
  ##----------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  
  N <- sum(z[1:n.individuals])
  
})

##------------------------------------------------------------------------------										  
## Set up the model
# Constants
lmIndicatorConstants <- nimConstants

# Initial values
lmIndicatorInits <- nimInits
lmIndicatorInits$omega <- .5
lmIndicatorInits$w <- rep(1, nimConstants$n.habCovs)
lmIndicatorInits$betaHabCovs <- runif(nimConstants$n.habCovs, -0.1, 0.1)

# Input data
nimData$c1 <- 1
nimData$c2 <- 1
nimData$c3 <- 1
lmIndicatorData <- nimData


##------------------------------------------------------------------------------										  
lmIndicatorModel <- nimbleModel(code = lmIndicatorCode,
                                constants = lmIndicatorConstants,
                                inits = lmIndicatorInits,
                                data = lmIndicatorData)

lmIndicatorConf <- configureMCMC(lmIndicatorModel)
lmIndicatorConf$addMonitors('w','omega','betaHabCovs','N')

configureRJ(lmIndicatorConf,
            targetNodes = c('betaHabCovs'),
            indicatorNodes = c('w'),
            control = list(mean = 0, scale = .2))

## Check the assigned samplers
lmIndicatorConf$printSamplers(c("w", "betaHabCovs"))

mcmcIndicatorRJ <- buildMCMC(lmIndicatorConf)

cIndicatorModel <- compileNimble(lmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, project = lmIndicatorModel)

## A very short test run!
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ,
                                        samplesAsCodaMCMC = TRUE,
                                        niter = 100,
                                        nburnin = 10))


