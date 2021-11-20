# -------------------------------------------------------------------------#
# Read Simulation results ----
# -------------------------------------------------------------------------#

#' This function is an extension base::list.files. In addition to searching
#' in .inputFolder, it will search in all xxxx.zip files found in this input
#' folder, where xxxx are four decimal digits.
#' Generate a character vector of filenames matching a pattern, where
#' these files may be stored in .inputFolder or inside xxxx.zip files
#' inside .inputFolder.
#'
#' Returns a character vector, like base::list.files. Matched files in
#' .zip files are denoted as
#' "path/to/file/xxxx.zip@#@matched.file".
#'
#' The above format is understood by the function readRDS.zipped
#'
#' @param path,pattern,all.files,full.names,recursive,ignore.case,include.dirs,no.. See [base::list.files()]
#'
#'
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::list.files.zipped(.inputFolder, "", full.names = TRUE)
list.files.zipped <- function(
  path,
  pattern      = NULL,
  all.files    = FALSE,
  full.names   = FALSE,
  recursive    = FALSE,
  ignore.case  = FALSE,
  include.dirs = FALSE,
  no..         = FALSE) {

  files.unzipped <- list.files(
    path,
    pattern = pattern,
    all.files = all.files,
    full.names = full.names,
    recursive = recursive,
    ignore.case = ignore.case,
    include.dirs = include.dirs,
    no.. = no..)

  files.zip <- list.files(
    path,
    pattern = "^[0-9]+.zip$",
    all.files = all.files,
    full.names = full.names,
    recursive = recursive,
    ignore.case = ignore.case,
    include.dirs =include.dirs,
    no.. = no..)

  files.zipped <- sapply(files.zip, function(zipfile) {
    found <- grep(pattern, zip::zip_list(zipfile)$filename, value = TRUE)
    if(length(found)) {
      paste0(zipfile, "@#@", found)
    } else {
      ""
    }
  }, USE.NAMES = FALSE)
  unlist(c(files.unzipped, files.zipped[!files.zipped == ""]))
}

#' Reads an RDS file which can be stored in a xxxx.zip file
#' The filename argument specifies the path to the RDS file.
#' If the files is inside a zip file, then "@#@" is used as
#' a delimiter.
#'
#' @param filename path to file
#'
#' @examples
#' # Read parts from zipped files
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' files <- populationIsoboles:::list.files.zipped(.inputFolder, "", full.names = TRUE)
#' file <- grep("@#@", files, value =TRUE)[1]
#' print(file)
#' populationIsoboles:::readRDS.zipped(file)
#'
#' # Read non-zipped files
#' tf <- tempfile(fileext = ".rds")
#' saveRDS(mtcars, tf)
#' populationIsoboles:::readRDS.zipped(tf)
readRDS.zipped <- function(filename) {
  if(length(grep(".zip@#@", filename, fixed = TRUE)) > 0) {
    names <- strsplit(filename, split = "@#@")[[1]]
    con <- unz(names[1], names[2])
    con2 <- gzcon(con)
    obj <- readRDS(con2)
    close(con2)
    obj
  } else {
    readRDS(filename)
  }
}


#' Get the identifier of a population from the filenam
#'
#' @param files vector of filenames
#'
#' @return Numeric vector of k's
#'
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#'
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' files <- list.files(.inputFolder)
#' files <- grep("^(\\d+)[._]", files, value = TRUE)
#' try(populationIsoboles:::extract_popkFromFilename(files))
extract_popkFromFilename <- function(files) {
  # Extract k from filenames
  k <- basename(files)
  if (any(!grepl("(\\d+)[._].*", k)))
    stop("All filenames processed should start with digits. Bad files: \n",
         paste0(grep("(\\d+)[._].*", k, value = TRUE, invert = TRUE), collapse = "\n"))
  k <- gsub("(\\d+)[._].*", "\\1", k)
  k <- as.numeric(k)
  k
}



#' Functions to read simulation results
#'
#' @param .inputFolder Folder where the results are stored.
#'   Is searched recursively
#' @param itermax_dt Output from [PI_readItermax()]
#' @param allPaths Output from [PI_readAllPaths()]
#'
#' @return data.table(itermax, k) containing the maximum iteration for
#'   each population
#'
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @importFrom data.table data.table rbindlist fwrite fread
#'
#' @examples
#' library(populationIsoboles)
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::PI_readItermax(.inputFolder)
PI_readItermax <- function(.inputFolder) {
  resultfile <- file.path(.inputFolder, "RESULT-01-itermax.csv")
  recalcNecessary <- determineRecalcNeccessary(.inputFolder, resultfile)
  if (!recalcNecessary) return(data.table::fread(resultfile))

  itermax_dt <- list.files.zipped(.inputFolder, pattern = "itermax", full.names = TRUE, recursive = TRUE)
  itermax_dt <- grep("(log|csv)$", itermax_dt, value = TRUE, invert = TRUE)
  k <- extract_popkFromFilename(itermax_dt)
  itermax_dt <- mapply(.x = itermax_dt, .y = k, FUN = function(.x,.y) data.table::data.table(itermax = readRDS.zipped(.x), k = .y), SIMPLIFY = FALSE)
  itermax_dt <- data.table::rbindlist(itermax_dt, use.names = TRUE, fill = TRUE)
  data.table::fwrite(itermax_dt, resultfile)
  itermax_dt
}


#' Functions to read simulation results
#'
#' @param .inputFolder Folder where the results are stored. Is searched recursively
#'
#' @return data.table(outsideGrid) containing if the isobole is outside the defined grid for each population
#' @importFrom data.table data.table rbindlist fwrite fread
#'
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::PI_readOutsideGrid(.inputFolder)
PI_readOutsideGrid <- function(.inputFolder) {
  resultfile <- file.path(.inputFolder, "RESULT-02-OutsideGrid.csv")
  recalcNecessary <- determineRecalcNeccessary(.inputFolder, resultfile)
  if (!recalcNecessary) return(data.table::fread(resultfile))
  outsideGrid_dt <- list.files.zipped(.inputFolder, pattern = "outsideGrid", full.names = TRUE, recursive = TRUE)
  outsideGrid_dt <- grep("(log|csv)$", outsideGrid_dt, value = TRUE, invert = TRUE)
  k <- extract_popkFromFilename(outsideGrid_dt)
  outsideGrid_dt <- mapply(.x = outsideGrid_dt, .y = k, FUN = function(.x,.y) data.table::data.table(outsideGrid = readRDS.zipped(.x), k = .y), SIMPLIFY = FALSE)
  outsideGrid_dt <- data.table::rbindlist(outsideGrid_dt, use.names = TRUE, fill = TRUE)
  data.table::fwrite(outsideGrid_dt, resultfile)
  outsideGrid_dt
}

#' @rdname PI_readItermax
#' @importFrom data.table rbindlist fwrite fread
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::PI_readAllPaths(.inputFolder)
PI_readAllPaths <- function(.inputFolder, itermax_dt = PI_readItermax(.inputFolder)) {
  resultfile <- file.path(.inputFolder, "RESULT-03-AllPaths.csv")
  recalcNecessary <- determineRecalcNeccessary(.inputFolder, resultfile)
  if (!recalcNecessary) return(data.table::fread(resultfile))

  allPaths <- list.files.zipped(.inputFolder, pattern = "pathlist", full.names = TRUE, recursive = TRUE)
  allPaths <- grep("(log|csv)$", allPaths, value = TRUE, invert = TRUE)
  k <- extract_popkFromFilename(allPaths)
  allPaths <- lapply(allPaths, readRDS.zipped)
  allPaths <- mapply(pathlist = allPaths, popk = k, FUN = function(pathlist,popk) {
    out <- mapply(.x = pathlist, .y = seq_along(pathlist), FUN = function(.x,.y) .x[,`:=`(iteration = .y, k = popk)], SIMPLIFY = FALSE)
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  }, SIMPLIFY = FALSE)
  allPaths <- data.table::rbindlist(allPaths, use.names = TRUE, fill = TRUE)
  allPaths <- itermax_dt[allPaths, on = "k"]
  data.table::fwrite(allPaths, resultfile)
  allPaths
}

#' @rdname PI_readItermax
#' @importFrom data.table copy
#' @importFrom stats predict loess
PI_getFinalPaths <- function(allPaths) {
  finalPaths <- data.table::copy(allPaths)
  finalPaths <- finalPaths[iteration == itermax]
  finalPaths
}


#' @rdname PI_readItermax
#' @importFrom data.table data.table rbindlist
PI_readGridLens <- function(.inputFolder) {
  gridInfo <- PI_readGridInfo(.inputFolder)
  gridlens <- mapply(.x = gridInfo$gridlens, .y = seq_along(gridInfo$gridlens), FUN = function(.x,.y) data.table::data.table(iteration = as.numeric(.y), gridlenx = .x[1], gridleny = .x[2]), SIMPLIFY = FALSE)
  gridlens <- data.table::rbindlist(gridlens, use.names = TRUE, fill = TRUE)
  gridlens
}

#' Save reading of gridInfo
#'
#' I moved the output of 001-gridInfo.rds from .outputFolder to .outputFolder/Simulations.
#' This function searches both folders to allow backwards compatibility
#' and that .inputFolder now means the "Simulations"-folder for all PI_read* functions
#'
#' @rdname PI_readItermax
#' @importFrom data.table data.table rbindlist
PI_readGridInfo <- function(.inputFolder) {

  # Search for gridInfo in Simulations or in Folder above
  if (basename(.inputFolder) == "Simulations") rootdir <- ".." else rootdir <- "."
  filepath <- list.files.zipped(file.path(.inputFolder, rootdir), pattern = "001-GridInfo.rds$", full.names = TRUE, recursive = TRUE)
  if (length(filepath) > 1) warning("More than one GridInfo was found. This one was used: ", filepath[1])
  filepath <- filepath[1]

  readRDS.zipped(filepath)
}


#' @rdname PI_readItermax
#' @param gridlens output of [PI_readGridLens]
#'
#' @author Daniel Lill (daniel.lill@intiquan.com), Venelin Metiv (venelin.mitov@intiquan.com)
#' @md
#'
#' @importFrom data.table rbindlist melt fwrite fread
#'
#' @examples
#' library(populationIsoboles)
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::PI_readGridValues(.inputFolder)
PI_readGridValues <- function(.inputFolder, itermax_dt = PI_readItermax(.inputFolder), gridlens = PI_readGridLens(.inputFolder)) {
  resultfile <- file.path(.inputFolder, "RESULT-04-Grid.csv")
  recalcNecessary <- determineRecalcNeccessary(.inputFolder, resultfile)
  if (!recalcNecessary) return(data.table::fread(resultfile))

  gridfiles <- list.files.zipped(.inputFolder, pattern = "grid", full.names = TRUE, recursive = TRUE)
  gridfiles <- grep("(log|csv)$"  , gridfiles, value = TRUE, invert = TRUE)
  gridfiles <- grep("GLOBAL", gridfiles, value = TRUE, invert = TRUE)
  k <- extract_popkFromFilename(gridfiles)
  grid <- mapply(.x = gridfiles, .y = k, FUN = function(.x,.y) cbind(readRDS.zipped(.x), k = .y), SIMPLIFY = FALSE)
  grid <- data.table::rbindlist(grid, use.names = TRUE, fill = TRUE)
  grid[,`:=`(objvalueCum = NULL, objvalue0 = NULL, evaluated0 = NULL, evaluatedCum = NULL)]
  if ("objvaluei" %in% names(grid)) grid[,`:=`(objvaluei = NULL)] # circumvent warning
  grid <- data.table::melt(grid, id.vars = c("x", "y", "k"),
                           measure.vars = lapply(c("^evaluated", "objvalue"), grep, names(grid)),
                           variable.name = "iteration",
                           value.name = c("evaluated","objvalue"),
                           variable.factor = FALSE)
  grid[,`:=`(iteration = as.numeric(as.character(iteration)))]
  grid <- itermax_dt[grid, on = c("k")]
  grid <- gridlens[grid, on = "iteration"]
  data.table::fwrite(grid, resultfile)
  grid
}

#' @rdname PI_readItermax
#' @importFrom data.table data.table rbindlist fwrite fread
PI_readAreas <- function(.inputFolder, itermax_dt = PI_readItermax(.inputFolder)) {
  resultfile <- file.path(.inputFolder, "RESULT-05-Areas.csv")
  recalcNecessary <- determineRecalcNeccessary(.inputFolder, resultfile)
  if (!recalcNecessary) return(data.table::fread(resultfile))
  areas <- list.files.zipped(.inputFolder, pattern = "areas", full.names = TRUE, recursive = TRUE)
  areas <- grep("(log|csv)$", areas, value = TRUE, invert = TRUE)
  k <- extract_popkFromFilename(areas)
  areas <- mapply(.x = areas, .y = k, FUN = function(.x,.y) data.table::data.table(areas = readRDS.zipped(.x), iteration = seq_along(readRDS.zipped(.x)), k = .y), SIMPLIFY = FALSE)
  areas <- data.table::rbindlist(areas, use.names = TRUE, fill = TRUE)
  areas <- areas[itermax_dt, on = "k"]
  data.table::fwrite(areas, resultfile)
  areas
}


# -------------------------------------------------------------------------#
# Post processing ----
# -------------------------------------------------------------------------#



#' Calculate Confidence Level Response Surface
#'
#' For each point on the full grid and each population,
#' check if it lies left or right of the population isobole and
#' assign 0 or 1 to it.
#' Then, calculate the mean for each point.
#' This is equivalent to the confidence level response surface
#'
#' This function will look in .inputFolder for a file called "GLOBAL-101-percentages-grid_fine.rds".
#' If it is found and FLAGoverwriteQuantiles = FALSE, it will load this result
#'
#' @inheritParams .plotDocumentation
#' @param .inputFolder path/to/Simulations
#' @param gridInfo Output of [PI_readGridInfo]
#'
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#'
#'
#' @return data.table(x,y, percentages)
#' @family Isobole simulation post processing
#' @export
#' @importFrom data.table data.table
#' @importFrom sp point.in.polygon
#' @importFrom stats setNames
aggregateIsoboles <- function(.inputFolder,
                          gridInfo = NULL,
                          finalPaths = NULL,
                          FLAGoverwriteQuantiles = determineRecalcNeccessary(.inputFolder = .inputFolder,
                                                                             fileToRecalc = "GLOBAL-101-percentages-grid_fine.rds")) {
  # Load old if existing
  if (!FLAGoverwriteQuantiles && file.exists(file.path(.inputFolder, "GLOBAL-101-percentages-grid_fine.rds"))){
    cat("Previously computed quantiles were loaded. To update, set FLAGoverwriteQuantiles.\n")
    return(readRDS.zipped(file.path(.inputFolder, "GLOBAL-101-percentages-grid_fine.rds")))
  }

  # Get results from inputFolder
  if (!is.null(.inputFolder) && is.null(gridInfo))
    gridInfo <- PI_readGridInfo(.inputFolder)
  if (!is.null(.inputFolder) && is.null(finalPaths))
    finalPaths <- PI_getFinalPaths(PI_readAllPaths(.inputFolder))

  # Process gridInfo:
  # * disallow non-integer value for gridlens
  # * Sample only up the maximum value of the isoboles, but not further
  gridmin  <- round(gridInfo$gridmin)  # For inPoly(), it's better to have actual 0s, not values very close to it.
  imax <- length(gridInfo$gridlens)
  if (any(gridInfo$gridlens[[imax]] < 1))
    cat("Minimal mesh length should not be smaller than one: 2d quantile calculation has not been implemented yet for this case")
  gridlen  <- round(gridInfo$gridlens[[imax]])
  gridmax <- gridInfo$gridmax
  gridmax <- vapply(stats::setNames(nm = names(gridmax)), function(.x) max(seq(gridmin[.x],gridmax[.x],round(gridlen[.x]))), FUN.VALUE = 1.0)
  kmax <- max(finalPaths$k)

  grid_fine <- data.table::data.table(expand.grid(x = seqminmax(gridmin["x"],gridmax["x"],gridlen["x"]),
                                                  y = seqminmax(gridmin["y"],gridmax["y"],gridlen["y"]),
                                                  popk = 1:kmax))
  # Get values on fine grid for each population
  grid_fine[popk %in% unique(finalPaths$k),`:=`(objvaluei = as.numeric(as.logical(
    sp::point.in.polygon(x, y,
                         c(finalPaths[k == popk,xp], max(x) + 1, max(x) + 1, 0     , finalPaths[k == popk,xp][1]),
                         c(finalPaths[k == popk,yp], 0,     max(y) + 1,  max(y) + 1, finalPaths[k == popk,yp][1]),
                         TRUE)))), by = "popk"]
  grid_fine[,`:=`(percentages = mean(objvaluei, na.rm = TRUE)), by = c("x", "y")]

  saveRDS(grid_fine, file.path(.inputFolder, "GLOBAL-101-percentages-grid_fine.rds"))

  grid_fine
}




# -------------------------------------------------------------------------#
# Plots ----
# -------------------------------------------------------------------------#

#' Documentation for plotting functions
#'
#' This function is only there to document parameters of plotting functions.
#' This is necessary because inheritParams can only inherit from a single function
#' and no plot uses all arguments of all plots
#'
#' @param finalPaths Result from [PI_getFinalPaths]
#' @param grid Result from [PI_readGridValues]
#' @param allPaths Result from [PI_readAllPaths]
#' @param itermax_dt Result from [PI_readItermax]
#' @param grid_fine Result from [aggregateIsoboles]
#' @param Compound1,Compound2 Character. Compound name, for cases when it's not available in the Malaria project
#' @param levels Vector of confidence levels, e.g. c(0.5,0.95). For plotRibbon, named vector with min and max, e.g. c(min = 0.025, max = 0.975)
#' @param FLAGoverwriteQuantiles Redo the calculation of the quantiles?
#'
#' @param .inputFolder Folder where simulations are stored. \cr
#' Example: file.path(.outputFolder, "Simulations")
#' @param filename Path to save plot
#' @param ... Arguments goint to \code{\link[ggplot2]{ggsave}}
#' @md
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' plotSpaghetti(.inputFolder)
#' plotQuantiles(.inputFolder)
#' plotCLRS(.inputFolder)
.plotDocumentation <- function(finalPaths, grid, allPaths,
                               itermax_dt, grid_fine,
                               levels,
                               Compound1 = getCompoundInfo("Compound1"),
                               Compound2 = getCompoundInfo("Compound2"),
                               filename = NULL, .inputFolder,
                               FLAGoverwriteQuantiles =
                                 determineRecalcNeccessary(.inputFolder), ...) {
  NULL
}

#' Plot isobole spaghetti plot
#'
#' @inheritParams .plotDocumentation
#'
#' @return ggplot
#'
#' @family Internal isobole plotting functions
#'
#' @importFrom data.table copy
#' @importFrom ggplot2 ggsave
PI_plotInternalSpaghetti <- function(finalPaths, Compound1 = getCompoundInfo("Compound1"),
                                     Compound2 = getCompoundInfo("Compound2"), filename = NULL, ...) {
  plot_df1 <- data.frame(data.table::copy(finalPaths[,list(xp,yp,k)]))
  pl <- ggplot_PI(plot_df1, aes(xp,yp, group = k)) + geom_path(alpha = 0.1) +
    guides(color = FALSE) + labs(x = .axlab(Compound1), y = .axlab(Compound2)) + geom_blank()
  if (!is.null(filename)) ggplot2::ggsave(plot = pl, filename = filename, ...)
  pl
}


#' Plot one population isobole only
#'
#' @inheritParams .plotDocumentation
#' @param Popk Integer denoting the population to plot
#'
#' @return ggplot
#'
#' @family Internal isobole plotting functions
#' @importFrom ggplot2 ggsave
PI_plotInternalSingleExample <- function(grid, allPaths,  Popk = 1, Compound1 = getCompoundInfo("Compound1"), Compound2 = getCompoundInfo("Compound2"), filename = NULL, ...) {
  # plot gridpoints and curve of single population
  pl <- ggplot_PI() +
    geom_point(data = grid[evaluated == 1 & k == Popk], mapping = aes(x,y, shape = factor(iteration))) +
    geom_path(data = allPaths[k == Popk & iteration == itermax], mapping = aes(xp,yp), color = "red", size = 2) +
    labs(x = .axlab(Compound1), y = .axlab(Compound2)) +
    ggtitle(paste0("Example: Population ", Popk))
  if (!is.null(filename)) ggplot2::ggsave(plot = pl, filename = filename, ...)
  pl
}

#' Plot multiple plots to demonstrate the algorithm
#'
#' @inheritParams .plotDocumentation
#' @param Popk Integer denoting the population to plot
#' @param .outputFolder Folder to save the plots
#'
#' @return ggplot
#'
#' @family Internal isobole plotting functions
#' @importFrom ggplot2 ggsave
PI_plotInternalAlgoExample <- function(grid, allPaths, itermax_dt, Popk = 1, Compound1 = getCompoundInfo("Compound1"), Compound2 = getCompoundInfo("Compound2"), .outputFolder = NULL, ...) {
  for (i in 1:itermax_dt[k == Popk]$itermax){
    pl42 <- ggplot_PI() +
      geom_path(data = allPaths[k==Popk & iteration %in% 1:i], mapping = aes(xp,yp, color = factor(iteration), linetype = factor(iteration))) +
      labs(x = .axlab(Compound1), y = .axlab(Compound2)) +
      ggtitle(paste0("Example: Population ", Popk))
    ggplot2::ggsave(plot = pl42, filename = file.path(.outputFolder, "Plots", paste0(sprintf("%03i", 10+i), "-Refine_Grid-Algo.png")))
    pl43 <- ggplot_PI() +
      geom_point(data = grid[evaluated == 1 & k == Popk & iteration %in% 1:i], mapping = aes(x,y, shape = factor(iteration))) +
      geom_path(data = allPaths[k==Popk & iteration %in% 1:i], mapping = aes(xp,yp, color = factor(iteration)), size = 1)
    ggplot2::ggsave(plot = pl43, filename = file.path(.outputFolder, "Plots", paste0(sprintf("%03i", 20+i), "-Refine_Grid-Algo.png")))
    pl44 <- ggplot_PI() +
      geom_path(data = allPaths[k==Popk & iteration == i], mapping = aes(xp,yp, color = factor(iteration), linetype = factor(iteration))) +
      geom_point(data = grid[evaluated == 1 & k==Popk & iteration == i], mapping = aes(x,y, shape = factor(iteration))) +
      labs(x = "Dose1 [a.u.]", y = "Dose2 [a.u.]") +
      ggtitle(paste0("Example: Population ", Popk))
    ggplot2::ggsave(plot = pl44, filename = file.path(.outputFolder, "Plots", paste0(sprintf("%03i", 30+i), "-Refine_Grid-Algo.png")))
  }
  NULL
}



#' Plot spaghetti with quantiles
#'
#' @inheritParams .plotDocumentation
#'
#' @return ggplot
#' @family Internal isobole plotting functions
#' @importFrom ggplot2 ggsave
PI_plotInternalQuantiles <- function(
  grid_fine,
  finalPaths,
  levels = c(.025,.5,.975),
  Compound1 = getCompoundInfo("Compound1"),
  Compound2 = getCompoundInfo("Compound2"),
  filename = NULL,
  ...) {

  pl <- ggplot_PI(unique(grid_fine[,list(x,y,percentages)]), aes(x,y,z = percentages)) +
    geom_path(data = finalPaths, aes(x = xp, y = yp, group = k, z = NULL), alpha = .02,color = "black") +
    geom_contour(aes(colour = factor(..level..)), breaks = levels, size = 1.2) +
    scale_color_manual(values = c(#"darkorange", "navyblue", "firebrick",
      PIcolors)) +
    labs(x = .axlab(Compound1), y = .axlab(Compound2), color = "Confidence\nLevel")
  if (!is.null(filename)) ggplot2::ggsave(plot = pl, filename = filename, ...)
  pl
}




#' Plot ribbon which is robust with respect to shape
#'
#' Deprecated: Plot ribbon of values which are compatible with 95% cure.
#' This is based on a confidence interval closed on both sides which is not
#' what is of interest...
#'
#' @inheritParams .plotDocumentation
#'
#' @return ggplot object
PI_plotInternalRibbon2 <- function(grid_fine, levels = c(min = 0.025, max = 0.975),
                                   Compound1 = getCompoundInfo("Compound1"),
                                   Compound2 = getCompoundInfo("Compound2"),
                                   filename = NULL, ...) {

  grid_fine <- grid_fine[,list(x,y,percentages)]
  grid_fine[,`:=`(x = round(x), y = round(y))]
  grid_fine <- unique(grid_fine)

  gf0 <- copy(grid_fine)

  grid_fine[percentages < levels["min"],`:=`(percentages = NA)]
  grid_fine[percentages > levels["max"],`:=`(percentages = NA)]
  grid_fine[,`:=`(percentages = abs(0.5 - percentages))]
  grid_fine[,`:=`(alpha = as.numeric(!is.na(percentages)))]

  pl <-   ggplot_PI(grid_fine,aes(x,y)) +
    geom_raster(aes(fill = percentages, alpha = alpha)) +
    scale_fill_gradient(low = "grey10", high = "grey80", na.value = "white") +
    geom_contour(data = gf0, mapping = aes(z = percentages), breaks = 0.5, color = "black", size = 2) +
    guides(color = FALSE, fill = FALSE, alpha = FALSE) +
    labs(x = .axlab(Compound1), y = .axlab(Compound2))
  if (!is.null(filename)) ggplot2::ggsave(plot = pl, filename = filename, ...)
  pl
}




#' Plot the distributions of the mono doses
#'
#' @inheritParams .plotDocumentation
#'
#' @return ggplot
#'
#' @family Internal isobole plotting functions
#' @importFrom data.table copy rbindlist melt
#' @importFrom ggplot2 ggsave
PI_plotInternalMonoDistributions <- function(grid_fine, filename = NULL, ...) {
  quantiles_Drug1 <- data.table::copy(grid_fine)
  quantiles_Drug1 <- unique(quantiles_Drug1[y == 0,list(x, percentages)])
  quantiles_Drug1 <- quantiles_Drug1[x != 1280]
  quantiles_Drug1[,`:=`(Drug = "Drug 1")]

  quantiles_Drug2 <- data.table::copy(grid_fine)
  quantiles_Drug2 <- unique(quantiles_Drug2[x == 0,list(x = y, percentages)])
  quantiles_Drug2 <- quantiles_Drug2[x != 1280]
  quantiles_Drug2[,`:=`(Drug = "Drug 2")]

  quantiles_mono <- data.table::rbindlist(list(quantiles_Drug1, quantiles_Drug2), use.names = TRUE, fill = TRUE)
  pl <- ggplot_PI(data.table::melt(quantiles_mono, id.vars = c("x", "Drug")), aes(x, value, linetype = variable)) + geom_line() +
    facet_wrap(~Drug, scales = "free_x") + labs(x = "Dose [mg]", y = "% of populations meeting objvalue")
  if (is.null(filename)) ggplot2::ggsave(plot = pl, filename = filename, ...)
  pl
}

#' See PI_plotInternalMonoDistributions
#'
#' @param grid_fine Output from aggregateIsoboles
#' @param cross_x,cross_y Doses at which the cross section should be evaluated
#' @param filename,... going to [ggplot2::ggsave()]
#'
#' @return ggplot
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
PI_plotInternalCrossSections <- function(grid_fine, cross_x = NULL, cross_y = NULL, filename = NULL, ...) {

  quantiles_x1 <- quantiles_x2 <- NULL
  if (is.null(cross_x) & is.null(cross_y)) stop("cross_x and cross_y can't both be NULL")
  if (!is.null(cross_x)){
    cross_x <- ensure_cross_in_vec(vec  = grid_fine$x, cross = cross_x)
    quantiles_x2 <- data.table::copy(grid_fine)
    quantiles_x2 <- unique(quantiles_x2[x %in% cross_x,list(x = y, cross = paste0(getCompoundInfo("Compound2")$Name, "@\n", x ," mg ", getCompoundInfo("Compound1")$Name), percentages)])
    quantiles_x2[,`:=`(Drug = getCompoundInfo("Compound2")$Name)]}
  if (!is.null(cross_y)){
    cross_y <- ensure_cross_in_vec(vec = grid_fine$y, cross = cross_y)
    quantiles_x1 <- data.table::copy(grid_fine)
    quantiles_x1 <- unique(quantiles_x1[y %in% cross_y,list(x = x, cross = paste0(getCompoundInfo("Compound1")$Name, "@\n", y, " mg ", getCompoundInfo("Compound2")$Name), percentages)])
    quantiles_x1[,`:=`(Drug = getCompoundInfo("Compound1")$Name)]}

  quantiles_mono <- data.table::rbindlist(list(quantiles_x1, quantiles_x2), use.names = TRUE, fill = TRUE)
  quantiles_mono[,`:=`(cross = factor(cross, levels = unique(cross)))]

  pl <- ggplot_PI(quantiles_mono, aes(x, percentages, color = Drug)) + geom_line() +
    facet_wrap(~cross, scales = "free_x") + labs(x = "Dose [mg]", y = "% Patients cured") +
    scale_color_PI(aesthetics = c("fill", "color"))
  if (!is.null(filename)) ggplot2::ggsave(plot = pl, filename = filename, ...)
  pl
}





#' Quick plotting of simulation results
#'
#' This plot gives a feeling of the spread of the population isoboles.
#' Its interpretation is difficult, therefore the plot is deprecated.
#'
#' @inheritParams .plotDocumentation
#'
#' @md
#'
#' @return ggplot
#' @family Isobole plotting
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::plotRibbon(.inputFolder)
plotRibbon <- function(.inputFolder, filename = NULL, levels = c(min = 0.025, max = 0.975),
                       FLAGoverwriteQuantiles = determineRecalcNeccessary(.inputFolder = .inputFolder,
                                                                          fileToRecalc = "GLOBAL-101-percentages-grid_fine.rds"),
                       ...) {

  grid_fine  <- aggregateIsoboles(.inputFolder, FLAGoverwriteQuantiles = FLAGoverwriteQuantiles)

  pl <- PI_plotInternalRibbon2(grid_fine, filename = filename, levels = levels, ...)

  pl
}




# -------------------------------------------------------------------------#
# Quick Plots ----
# -------------------------------------------------------------------------#



#' Quick plotting of simulation results
#'
#' @inheritParams .plotDocumentation
#'
#' @return A ggplot with the popIsoboles plotted as spaghetti
#' @export
#' @family Isobole plotting
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' plotSpaghetti(.inputFolder)
plotSpaghetti <- function(.inputFolder, filename = NULL, ...) {
  finalPaths <- PI_getFinalPaths(PI_readAllPaths(.inputFolder, itermax_dt = PI_readItermax(.inputFolder)))
  PI_plotInternalSpaghetti(finalPaths, filename = filename, ...)
}


#' Quick plotting of simulation results
#'
#' @inheritParams .plotDocumentation
#'
#' @return ggplot
#' @export
#' @family Isobole plotting
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' plotSpaghetti(.inputFolder)
plotQuantiles <- function(.inputFolder, filename = NULL,
                          FLAGoverwriteQuantiles = determineRecalcNeccessary(.inputFolder = .inputFolder,
                                                                             fileToRecalc = "GLOBAL-101-percentages-grid_fine.rds"),
                          levels = c(.025,.5,.975),
                          ...) {
  finalPaths <- PI_getFinalPaths(PI_readAllPaths(.inputFolder, itermax_dt = PI_readItermax(.inputFolder)))
  grid_fine  <- aggregateIsoboles(.inputFolder, finalPaths = finalPaths, FLAGoverwriteQuantiles = FLAGoverwriteQuantiles)

  PI_plotInternalQuantiles(grid_fine, finalPaths, filename = filename,levels = levels, ...)
}


#' Plot confidence level response surface
#'
#' @inheritParams .plotDocumentation
#'
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @export
#' @importFrom ggplot2 ggsave
#' @family Isobole plotting
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' plotCLRS(.inputFolder)
plotCLRS <- function(.inputFolder,
                     filename = NULL,
                     levels = c(0.5,0.95),
                     Compound1 = getCompoundInfo("Compound1"),
                     Compound2 = getCompoundInfo("Compound2"),
                     FLAGoverwriteQuantiles = determineRecalcNeccessary(.inputFolder = .inputFolder,
                                                                        fileToRecalc = "GLOBAL-101-percentages-grid_fine.rds"),
                     ...
) {

  finalPaths <- PI_getFinalPaths(PI_readAllPaths(.inputFolder, itermax_dt = PI_readItermax(.inputFolder)))
  grid_fine  <- aggregateIsoboles(.inputFolder, finalPaths = finalPaths, FLAGoverwriteQuantiles = FLAGoverwriteQuantiles)
  clrs <- unique(grid_fine[, list(x,y, percentages)])

  pl <- ggplot_PI(clrs, aes(x,y, fill= percentages)) +
    geom_tile() +
    scale_fill_viridis_c() +
    geom_contour(aes(z = percentages, color = factor(..level..)), breaks = levels, size = 2) +
    scale_color_PI() +
    labs(x = .axlab(Compound1), y = .axlab(Compound2), fill = "Confidence\nLevel", color = "Confidence\nLevel") +
    geom_blank()
  if (!is.null(filename)) ggplot2::ggsave(plot = pl, filename = filename, ...)
  pl
}



# -------------------------------------------------------------------------#
# Quick Tables ----
# -------------------------------------------------------------------------#

#' Summarize the number of iterations needed across populations
#'
#' @param .inputFolder path/to/Simulations
#'
#' @return data.table(Iterations = number of iterations, N = how many populations took that many iterations)
#'
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::PI_quickTableItermax(.inputFolder)
PI_quickTableItermax <- function(.inputFolder) {
  d <- PI_readItermax(.inputFolder)
  d <- d[,list(N = .N), by = "itermax"]
  d <- d[,list(Iterations = itermax, N)]
  d
}



# -------------------------------------------------------------------------#
# Helpers ----
# -------------------------------------------------------------------------#

#' Try getting the default Compound info
#'
#' If not existing, write dummy Compound Info
#'
#' @param Compound Example: "Compund1"
#'
#' @return IQRmalarie Compound info
getCompoundInfo <- function(Compound = "Compound1") {
  if (exists(Compound)) return(get(Compound))
  list(Name = Compound)
}


#' Make nice axis label
#'
#' @param Compound list with compound info.
#'
#' @return a nice label
#'
#' @examples
#' populationIsoboles:::.axlab(list(Name = "Drug 1"))
.axlab <- function(Compound) {
  paste0(Compound$Name, " [mg]")
}


#' Determine if calculated quantities need to be updated
#'
#' Check if some simulations are more recent than the file to recalculate.
#' This is useful for calculations which take some time, e.g. the 2d-quantiles ("percentages")
#'
#'
#' @param .inputFolder Path to "Simulations" with simulation results
#' @param fileToRecalc e.g. "GLOBAL-101-percentages-grid_fine.rds"
#'
#'
#' @examples
#' .inputFolder <- system.file("examples/Example-1/Outputs/Simulations",
#'   package = "populationIsoboles")
#' populationIsoboles:::determineRecalcNeccessary(.inputFolder)
determineRecalcNeccessary <- function(.inputFolder, fileToRecalc = file.path(.inputFolder, "GLOBAL-101-percentages-grid_fine.rds")) {

  if (!file.exists(fileToRecalc)) return(TRUE)
  fileToRecalc_mtime <- file.info(fileToRecalc)[,"mtime"]

  simresults   <- list.files(path = .inputFolder, pattern = "^\\d{3}",    full.names = TRUE)
  simresults_mtime   <- file.info(simresults)[,"mtime"]

  any(simresults_mtime > fileToRecalc_mtime)
}

#' Ensure a number is in a vector
#'
#' If cross not present in vec, choose the next bigger number of vec
#'
#' @param vec vector to test against
#' @param cross numbers which should be elements of vec
#'
#' @return numeric(length(cross))
#'
#' @examples
#' populationIsoboles:::ensure_cross_in_vec(vec = 1:10, cross = c(1,5,6.5))
ensure_cross_in_vec <- function(vec, cross) {
  vec <- sort(unique(vec))
  not_in_vec <- !cross %in% vec
  next_bigger <- vapply(cross[not_in_vec], function(x) {
    vec[which.min(x > vec)]
  }, 1)
  if (length(next_bigger)) {
    message(paste0(cross[not_in_vec], collapse = ","), "are not in grid and are replaced by ", paste0(next_bigger, collapse = ", "))
    cross[not_in_vec] <- next_bigger
  }
  cross
}

# -------------------------------------------------------------------------#
# Plotting with theme ----
# -------------------------------------------------------------------------#

#' Nice ggplot theme
#'
#' Taken from dMod
#'
#' @param base_size,base_family see ?theme_bw
#'
theme_PI <- function(base_size = 12, base_family = "") {
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]

  theme_bw(base_size = base_size, base_family = base_family) +
    theme(line = element_line(colour = "black"),
          rect = element_rect(fill = "white", colour = NA),
          text = element_text(colour = "black"),
          axis.text = element_text(size = rel(1.0), colour = "black", face = "bold"),
          axis.text.x = element_text(margin=unit(c(4, 4, 0, 4), "mm")),
          axis.text.y = element_text(margin=unit(c(4, 4, 4, 0), "mm")),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(-2, "mm"),
          legend.key = element_rect(colour = NA),
          panel.border = element_rect(colour = "black"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = rel(1.0)))
}


#' ggplot_PI
#'
#' @param data,mapping see ?ggplot
#'
#' @return ggplot
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
ggplot_PI <- function(data = NULL, mapping = aes()) {
  ggplot(data,mapping) +
    theme_PI()
}

#' scalecolorPI
#' Copied from dMod
#' @param ... see scale_color_manual
#'
#' @return added to plots
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
scale_color_PI <- function(...) {
  scale_color_manual(..., values = PIcolors)
}

#' scalefillPI
#' Copied from dMod
#' @param ... see scale_color_manual
#'
#' @return added to plots
#' @author Daniel Lill (daniel.lill@physik.uni-freiburg.de)
#' @md
scale_fill_PI <- function(...) {
  scale_fill_manual(..., values = PIcolors)
}


#' Nice colors
#'
#' Copied from dMod
#'
PIcolors <- c(rep(c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E",
                    "#4B1F6F", "#CC79A7","#006400", "#F0E442", "#8B4513"),2),
              rep("gray", 100))



