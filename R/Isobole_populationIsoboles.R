
#' Run population isoboles simulation
#'
#' @param Npop Number of populations
#' @param Nsubj Number of subjects per population
#' @param objectiveValue Target value
#' @param objectiveFunction Function to calculate target value. Needs to be function(AMTx1, AMTx2, parsIndiv),
#'  where AMTx1, AMTx2 are numeric vectors and parsIndiv is the output of the supplied function `sampleIndividuals`
#' @param sampleIndividuals Function to sample individual parameters
#' @param argsSampleIndividuals List of arguments for `sampleIndividuals`
#' @param gridInfo A list(AMTx1Max, AMTx2Max,offset), specifying the searched grid dimensions
#' @param itermin,itermax Number of iterations for fastIsoboles algorithm
#' @param .outputFolder Path to store results in
#' @param FLAGverbose print messages to console?
#'
#' @return Nothing, all results are written to disc
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
#' @export
#'
#' @examples
#' \dontrun{
#' # See examples in inst/examples/Example-01/01-RunExample.R
#' # See examples in inst/examples/Example-02/01-RunExample.R
#' }
runPopulationIsoboles <- function(
  Npop, Nsubj, objectiveValue,
  objectiveFunction,
  sampleIndividuals,
  argsSampleIndividuals,
  gridInfo = list(AMTx1Max = 1000, AMTx2Max = 1000, offset = 0),
  itermin = 4, itermax = 5,
  FLAGverbose = FALSE,
  .outputFolder = "."
) {
  # Ensure output folder exists
  dir.create(file.path(.outputFolder, "Simulations"), recursive = TRUE)
  # Augment gridInfo including lower boundaries and distances.
  # Save for later use in plotting
  gridInfo <- gridInfo_default(AMTx1Max = gridInfo$AMTx1Max,
                               AMTx2Max = gridInfo$AMTx2Max,
                               imax = itermax, offset = gridInfo$offset)
  saveRDS(gridInfo, file.path(.outputFolder, "Simulations", "001-GridInfo.rds"))

  # Loop over populations
  for (k in 1:Npop) {
    if (FLAGverbose) cat("Population ", k," ------------------- \n")
    # Sample parameters
    parsIndiv <- do.call(sampleIndividuals, argsSampleIndividuals)
    # Define objective function which only takes arguments (AMTx1, AMTx2) and
    #   has access to parsIndiv via the environment
    obj <- function(AMTx1, AMTx2) {
      objectiveFunction(AMTx1,AMTx2, parsIndiv)
    }
    # Run isobole algorithm for this population
    fastIsoboles(objfun = obj, objvalue = objectiveValue,
                    gridmin = gridInfo$gridmin, gridmax = gridInfo$gridmax,
                    imin = itermin, imax = itermax, FLAGverbose = FLAGverbose,
                    k = k, .outputFolder = .outputFolder)
    # Zip each population to avoid too many files
    PI_zip(PI_files(.outputFolder, k))
  }
  return("done")
}




# Exit ----
