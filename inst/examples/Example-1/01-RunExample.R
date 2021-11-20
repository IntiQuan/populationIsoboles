#### HEADER ================================================================
#
# S001-ImplementModel.R
#
# [PURPOSE]
# Demonstrate the use of the population isoboles API
#
# [AUTHOR]
# Daniel Lill (daniel.lill@intiquan.com)
#
# [CLEANING]
rm(list = ls())
#
# [INPUT]
#
# [OUTPUT]
.outputFolder <- "Outputs/"
#
# [OTHER]

## Preliminaries ====

try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

library(populationIsoboles)
library(ggplot2)
if(!require("MASS")) {install.packages("MASS"); library(MASS)}
if(!require("dMod")) {install.packages("dMod"); library(dMod)}

# -------------------------------------------------------------------------#
# Load Model and Parameter information ----
# -------------------------------------------------------------------------#
# .. Model -----
model <- source("Resources/modelDeparsed.R")$value
print(model)
# .. Parameter information -----
parameterInfo <- source("Resources/parametersDeparsed.R")$value
print(parameterInfo$theta) # Theta = Point estimates of fixed effects
print(parameterInfo$S)     # S     = Covariance matrix for Theta
print(parameterInfo$Omega) # Omega = Covariance matrix of random effects

# .. Show the basic PKPD model for one parameter -----
pars_doses <- c(AMTx1 = 500, AMTx2 = 500)
simtime    <- seq(0,24*28,24)
sim <- model(simtime, c(parameterInfo$theta, pars_doses))
plot(sim) +  
  geom_hline(yintercept = log(10), linetype = 2, color = "blue") + 
  annotate("text", 600, 3, label = "LLOQ", color = "blue") + 
  labs(x = "Time [h]", y  = "log Parasites / mL") + 
  guides(color = FALSE)

# -------------------------------------------------------------------------#
# Define functions to simulate populations ----
# -------------------------------------------------------------------------#
# .. Parameter sampling -----

#' Sample realizations of population parameters
samplePopPars <- function(Nsubj, theta, S) {
  thetapars <- MASS::mvrnorm(1, mu = theta , Sigma = S)
  thetapars <- matrix(thetapars, nrow = Nsubj, ncol = length(theta), byrow = TRUE)
  thetapars
}

#' Sample realizations of random effects
sampleEtas <- function(Nsubj, Omega) {
  etas <- MASS::mvrnorm(Nsubj, mu = rep(0, dim(Omega)[1]) , Sigma = Omega)
  etas
}

#' Function to simulate individual parameters
#'
#' Draws a set of population parameters and adds random effects on log-scale
#'
#' @param Nsubj Numeric, number of subjects to sample
#' @param theta Population parameter vector
#' @param S     Covariance matrix for theta
#' @param Omega Covariance matrix for random effects
#'
#' @return matrix[individuals, parameters]
sampleIndividuals <- function(Nsubj, theta, S, Omega) {
  exp(log(samplePopPars(Nsubj,theta,S)) + sampleEtas(Nsubj, Omega))
}

#' Add dosing information to parameters
#'
#' Replicates parsIndiv for each dose combination
#'
#' @param AMTx1,AMTx2 Vectors of sample length with the respective doses
#' e.g. AMTx1[1] and AMTx2[1] are the first dose combination
#' @param parsIndiv output of [sampleIndividuals()]
#'
#' @return matrix with length(AMTx1) * nrow(parsIndiv) rows
expandParametersWithDosing <- function(AMTx1, AMTx2, parsIndiv) {
  doses_expanded <- matrix(rep(c(AMTx1,AMTx2), each = nrow(parsIndiv)), ncol = 2)
  doses_expanded <- `colnames<-`(doses_expanded, c("AMTx1", "AMTx2"))
  pars_expanded <- do.call(rbind, lapply(1:length(AMTx1), function(x) parsIndiv))
  pars <- cbind(doses_expanded, pars_expanded)
}


# .. Model simulation -----

#' simulate the model for one parameter set and decide if the patient is cured
#'
#' Cure criterion: log Parasite concentration BLOQ (< log10) at 28 days
#'
#' @param pars0
#'
#' @return logical(1L), TRUE = "cure", FALSE = "no cure"
simModelBinaryOutput <- function(pars0) {
  prediction <- (prd)(seq(0,24*28, 24), pars0, deriv = FALSE)[[1]]
  prediction <- prediction[nrow(prediction), "OUTPUT1", drop  = TRUE]
  cure28 <- prediction < log(10)
  cure28
}

#' Simulate model for a population
#'
#' @param pars matrix[individuals,c(dosingParameters, individualParameters)]
#'   Each row corresponds to one simulation, determined by all parameters in this row
#'   The columns of this matrix need to be all model parameters +  the dosing parameters AMTx1, AMTx2
#'
#' @return logical(nrow(individuals)) indicating if a subject is cured or not
simModelMultiplePars <- function(pars) {
  ncores <- 4
  iscured <- parallel::mclapply(
    X = seq_len(nrow(pars)), mc.cores = ncores, FUN = function(i) {
      pars <- pars[i,,drop = TRUE]
      simModelBinaryOutput(pars)
    })
  do.call(c, iscured)
}

#' Final objective function to evaluate the cure rate
#'
#' @param AMTx1,AMTx2 Vectors of sample length with the respective doses
#'   e.g. AMTx1[1] and AMTx2[1] are the first dose combination
#' @param parsIndiv Output from [sampleIndividuals()]
#'
#' @return Vector of length(AMTx1), indicating the success rate at the dose
#'   combinations given by AMTx1,AMTx2
#' @export
objectiveFunction <- function(AMTx1,AMTx2, parsIndiv) {
  pars <- expandParametersWithDosing(AMTx1, AMTx2, parsIndiv)
  simres <- simModelMultiplePars(pars)
  cureRates <- vapply(seq_along(AMTx1), function(dose_i) {
    mean(simres[(dose_i-1)*nrow(parsIndiv) + 1:nrow(parsIndiv)])},
    FUN.VALUE = 0.1)
  cureRates
}



# -------------------------------------------------------------------------#
# Illustrate parameter sampling and model simulation ----
# -------------------------------------------------------------------------#

# Parameters of individual subjects
set.seed(1)
indivPars <- do.call(sampleIndividuals, c(list(N = 10), parameterInfo))
# Add dosing information, dose drug1 400mg, drug2 200 mg
pars <- cbind(AMTx1 = 400, AMTx2 = 200, indivPars)
simModelMultiplePars(pars) # 3 out of 10 subjects are cured
objectiveFunction(AMTx1 = c(400,400), AMTx2 = c(200, 400), indivPars)
# -------------------------------------------------------------------------#
# Additional Information for isobole simulation ----
# -------------------------------------------------------------------------#
# Determine maximum doses
gridInfo <- gridInfo_default(AMTx1Max = 5000, AMTx2Max = 5000, imax = 5, offset = 0)
objectiveValue <- 0.95 # Target cure rate in each population

# -------------------------------------------------------------------------#
# Isoboles function ----
# -------------------------------------------------------------------------#
# .. Population isobole simulation function -----
# Loop over populations:
# 1. Sample Individuals of this population
# 2. Run fastIsoboles, the fast isobole algorithm
# 3. Save results into .outputFolder/Simulations
runPopulationIsoboles <- function(
  Npop, Nsubj,
  gridInfo, objectiveValue,
  objectiveFunction,
  sampleIndividuals,
  parameterInfo,
  .outputFolder
) {
  dir.create(file.path(.outputFolder, "Simulations"), recursive = TRUE)
  # Save gridInfo for later reuse
  saveRDS(gridInfo, file.path(.outputFolder, "Simulations", "001-GridInfo.rds"))

  for (k in 1:Npop) {
    # Sample parameters
    parsIndiv <- do.call(sampleIndividuals, c(list(Nsubj = Nsubj), parameterInfo))
    # Define objective function which only takes arguments (AMTx1, AMTx2) and 
    #   has access to parsIndiv via the environment
    obj <- function(AMTx1, AMTx2) {
      objectiveFunction(AMTx1,AMTx2, parsIndiv)
    }
    # Run isobole algorithm for this population
    fastIsoboles(objfun = obj, objvalue = objectiveValue,
                    gridmin = gridInfo$gridmin, gridmax = gridInfo$gridmax,
                    imin = 5, imax = 5, FLAGverbose = FALSE, # Force 5 iterations
                    k = k, .outputFolder = .outputFolder)
    PI_zip(PI_files(.outputFolder, k)) # This step is useful if many populations are run
  }
}

# .. Run the population isoboles function -----
set.seed(1234)
runPopulationIsoboles(
  Npop = 3, Nsubj = 20, # Values should be higher to appropriately sample the parameter space.
                        # Needs powerful computer though
  gridInfo          = gridInfo,
  objectiveValue    = objectiveValue,
  objectiveFunction = objectiveFunction,
  sampleIndividuals = sampleIndividuals,
  parameterInfo     = parameterInfo,
  .outputFolder     = .outputFolder
)

# -------------------------------------------------------------------------#
# Plot results ----
# -------------------------------------------------------------------------#

## Spaghetti plot of isoboles
plotSpaghetti(file.path(.outputFolder, "Simulations"))

## Spaghetti with confidence levels overlaid
plotQuantiles(file.path(.outputFolder, "Simulations"))

## Confidence level response surface
plotCLRS(file.path(.outputFolder, "Simulations"))

# -------------------------------------------------------------------------#
# Clean up output folder ----
# -------------------------------------------------------------------------#
unlink(.outputFolder, recursive = TRUE)

# Exit ----

