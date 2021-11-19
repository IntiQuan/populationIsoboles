#### HEADER ================================================================
#
# Example 2
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
library(IQRtools)

aux_version("IQRtools", minVersion = "1.4")
# -------------------------------------------------------------------------#
# Load Model and Parameter information ----
# -------------------------------------------------------------------------#
# .. Model -----
model <- IQRmodel("Resources/Model.txt")

# .. Parameter information -----
parameterGPF <- load_GPF("Resources/Parameters.xlsx")

# .. Show the basic PKPD model for one parameter -----
dosing <- IQRdosing(c(0,0), c(1,2), c(500,500))
pars <- sampleIndParamValues(parameterGPF)$typicalIndParamValues
pars <- pars[setdiff(names(pars), c("ID", "ID.POP"))]
pars <- cbind(pars, PLbase = 9.)
simtime    <- seq(0,24*28,24)
sim <- sim_IQRmodel(model, simtime = simtime,
                    parameters = pars,
                    dosingTable = dosing, FLAGoutputsOnly = TRUE)
plot(sim) +
  geom_hline(yintercept = log(10), linetype = 2, color = "blue") +
  annotate("text", 600, 3, label = "LLOQ", color = "blue") +
  labs(x = "Time [h]", y  = "log Parasites / mL") +
  guides(color = FALSE)

# -------------------------------------------------------------------------#
# Define functions to simulate populations ----
# -------------------------------------------------------------------------#

#' Simulate model and determin cure rate
#'
#' @param AMTx1,AMTx2 Vector of doses
#' @param parsIndiv output from [IQRtools::sampleIndParamValues()]
#'
#' @return
#' @export
#'
#' @examples
objectiveFunction <- function(AMTx1,AMTx2, parsIndiv) {
  pars <- parsIndiv$indParamValues
  pars <- pars[setdiff(names(pars), "ID.POP")]
  # Optimal opportunity to parallelize simulations
  ncores <- 4
  cureRates <- parallel::mclapply(
    X = seq_along(AMTx1),
    mc.cores = ncores, FUN = function(i) {
      dosing <- IQRdosing(rep(0, length(AMTx1[i])*2),
                          ADM = rep(c(1,2), each = length(AMTx1[i])),
                          AMT = c(AMTx1[i], AMTx2[i]))
      eventTable <- create_IQReventTable(dosing, pars)
      sim <- sim_IQRmodel(model,simtime = seq(0,24*28,24),
                          eventTable = eventTable,
                          FLAGoutputsOnly = TRUE)
      PLfinal <- sim[sim$TIME==24*28, "OUTPUT1"]
      mean(PLfinal < log(10)) # Parasitemia less than 10parasites/mL
    })
  unlist(cureRates)
}

# -------------------------------------------------------------------------#
# Illustrate objective Function ----
# -------------------------------------------------------------------------#
parsIndiv <- sampleIndParamValues(parameterGPF, NULL, 10, Npop = 1)
AMTx1 <- AMTx2 <- c(100,500)
objectiveFunction(AMTx1, AMTx2, parsIndiv)

# -------------------------------------------------------------------------#
# Isoboles function ----
# -------------------------------------------------------------------------#
set.seed(1234)
# Nsubj and Npop Values should be higher to appropriately
# sample the parameter space. Needs powerful computer though
Npop <- 3
Nsubj <- 20
objectiveValue <- 0.95 # Target cure rate in each population
runPopulationIsoboles(
  Npop                  = Npop,
  Nsubj                 = Nsubj,
  gridInfo              = list(AMTx1Max = 5000, AMTx2Max = 5000),
  objectiveValue        = objectiveValue,
  objectiveFunction     = objectiveFunction,
  sampleIndividuals     = sampleIndParamValues,
  argsSampleIndividuals = list(spec = parameterGPF, Nsamples = Nsubj),
  FLAGverbose           = TRUE,
  .outputFolder         = .outputFolder
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

