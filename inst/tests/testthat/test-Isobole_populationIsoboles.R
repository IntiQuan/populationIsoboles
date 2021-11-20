test_that("Vignette results are unchanged", {

  if (require("IQRtools")) {
    library(populationIsoboles)
    aux_version("IQRtools", minVersion = "1.4")
    # -------------------------------------------------------------------------#
    # Load Model and Parameter information ----
    # -------------------------------------------------------------------------#
    # .. Model -----
    model <- IQRmodel("testdata/Resources/Model.txt")

    # .. Parameter information -----
    parameterGPF <- load_GPF("testdata/Resources/Parameters.xlsx")

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
    # Isoboles function ----
    # -------------------------------------------------------------------------#
    set.seed(1234)
    # Nsubj and Npop Values should be higher to appropriately
    # sample the parameter space. Needs powerful computer though
    Npop <- 3
    Nsubj <- 20
    objectiveValue <- 0.95 # Target cure rate in each population
    .outputFolder <- tempdir()
    runPopulationIsoboles(
      Npop                  = Npop,
      Nsubj                 = Nsubj,
      gridInfo              = list(AMTx1Max = 5000, AMTx2Max = 5000),
      objectiveValue        = objectiveValue,
      objectiveFunction     = objectiveFunction,
      sampleIndividuals     = sampleIndParamValues,
      argsSampleIndividuals = list(spec = parameterGPF, Nsamples = Nsubj),
      FLAGverbose           = FALSE,
      .outputFolder         = .outputFolder
    )

    CLRS <- aggregateIsoboles(file.path(.outputFolder, "Simulations"))

    # Compare CLRS to expectation
    CLRS_expectation <- readRDS("testdata/Expectations/test-populationIsoboles-CLRS.rds")

    expect_equivalent(CLRS, CLRS_expectation)
  }
})
