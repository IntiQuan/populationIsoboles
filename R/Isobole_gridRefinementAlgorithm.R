# -------------------------------------------------------------------------#
# Algorithm specific functions and helpers ----
# -------------------------------------------------------------------------#

#' Generate the default characteristics of the simulation grid
#'
#' @param AMTx1Max,AMTx2Max maximum dose for drug 1 and 2
#' @param imax Maximum number of iterations
#' @param offset small delta to avoid doses with value 0, which can cause problems if DOSELEVEL is a covariate
#'
#' @return list(gridmin, gridmax, gridlens)
#'
#' @examples
#' populationIsoboles:::gridInfo_default()
gridInfo_default <- function(AMTx1Max = 1000, AMTx2Max = 1000, imax = 7, offset = 1e-12) {
  offset <- if(length(offset)) {offset} else 1e-12

  gridmin = c(x = offset, y = offset)
  gridmax = c(x = AMTx1Max, y = AMTx2Max)
  gridlens = get_gridlens(gridmin = gridmin, gridmax = gridmax, imax = imax)
  list(gridmin = gridmin, gridmax = gridmax, gridlens = gridlens)
}


#' seq including min and max
#'
#' @param from,to,by as in \code{\link{seq}}
#' @author Daniel Lill, IntiQuan \email{daniel.lill@@intiquan.com}
#' @examples
#' seq(0,1,0.3)
#' populationIsoboles:::seqminmax(0,1,0.3)
seqminmax <- function(from, to, by) {
  unique(sort(c(from, to, seq(from, to, by))))
}


#' Construct a list of grid lenghts
#'
#' @param gridmin,gridmax named vectors
#' @param imax Maximum number of iterations
#'
#' @return list of vectors with decreasing gridlengths in powers of two.
#' @author Daniel Lill, IntiQuan \email{daniel.lill@@intiquan.com}
#' @family Grid refinement functions for isoboles
#'
#' @examples
#' gridmin <- c(x = 0, y = 0) + 0.0001
#' gridmax <- c(x = 10*2^7, y = 0.1*2^7) + 0.0001
#' populationIsoboles:::get_gridlens(gridmin, gridmax)
get_gridlens <- function(gridmin, gridmax, imax = 10) {
  lapply(seq_len(imax), function(.x) {(gridmax - gridmin) * 2 ^ (-.x)})
}




#' Interpolate grid values and run contourLines to obtain isobole path
#'
#' @param grid            data.table with colums (x,y,evaluatedCum, objvalueCum) and others
#' @param gridmin,gridmax General grid parameters: Vectors with names x and y c(x = 0, y = 0)
#' @param gridlen_i       Grid lengths of iteration i: c(x = 40, y = 40)
#' @param objvalue      Height of the countour line, for 95\% objvalue 95
#'
#' @return list with a data.table with columns (xp,yp), the nodes of the isobole path, and a logical variable outsideGrid to indicate if the isobole is out of the defined grid
#' @author Daniel Lill, IntiQuan \email{daniel.lill@@intiquan.com}
#'
#' @family Grid refinement functions for isoboles
#'
#' @importFrom akima interp
#' @importFrom grDevices contourLines
get_isobole_path <- function(grid, gridmin, gridmax, gridlen_i, objvalue = 0.95) {
  # 1 Interpolate the values on the grid
  grid_interp <- copy(grid)
  grid_interp <- grid_interp[evaluatedCum == 1]
  interpolated <- akima::interp(grid_interp$x, grid_interp$y, grid_interp$objvalueCum,
                                xo = seqminmax(gridmin["x"],gridmax["x"],min(gridlen_i["x"], (gridmax["x"] - gridmin["x"])/(2^6))),
                                yo = seqminmax(gridmin["y"],gridmax["y"],min(gridlen_i["y"], (gridmax["y"] - gridmin["y"])/(2^6))))
  # 2 Calculate contourlines
  path <- grDevices::contourLines(interpolated$z, levels = objvalue) # returns on scale [0,1]

  # Check if grid to small
  if(length(path)==0) {
    # Contourlines at max doses
    path <- data.table::data.table(xp = c(1,1,0),
                                   yp = c(0,1,1))

    # Logical variable for outer-grid
    outsideGrid <- TRUE

    # Little Warning:
    warning("Isobole for objvalue=", objvalue, " could not be calculated, and was set to the maximum doses: Please try increasing the grid size (AMTx1 and/or AMTx2).")
  }else{
    # Convert de data.table
    path <- data.table::data.table(xp = path[[1]]$x*1, yp = path[[1]]$y*1)

    # Logical variable for outer-grid
    outsideGrid <- FALSE
  }

  # When the first point of the isobole is at the maximum dose of drug 2 but not 0 for drug 1
  #   A segment on upper side is added
  if(path[1,xp]>0 && path[1,yp]==1){
    pathYFirst <- data.table::data.table(xp = c(0),
                                         yp = c(1))
    path <- data.table::rbindlist(list(pathYFirst,
                                       path
    ))
  }

  # When the last (down-est) point of the isobole is at the maximum dose of drug 1 but not 0 for drug 2
  #   A segment on the right side is added
  if(path[nrow(path),xp]==1 && path[nrow(path),yp]>0){
    pathXLast <- data.table::data.table(xp = c(1),
                                        yp = c(0))
    path <- data.table::rbindlist(list(path,
                                       pathXLast))
  }

  path <- arrange_isobole(path)
  path[,`:=`(xp = xp*gridmax["x"], yp = yp*gridmax["y"])] # rescale to [xmin,xmax]

  # Output:
  list(path=path, outsideGrid=outsideGrid)
}



#' Reorder rows of data.frame such that the shortest line is obtained
#'
#' Tests with TSP (traveling salesman problem) packages outputted nonsense,
#' therefore this function needed to be implemented anew.
#'
#' Order should be:
#' xp: 0         ... xp(yp=0)
#' yp: yp(xp=0)    ... 0
#'
#' Idea of the algorithm: Isoboles should be rather well behaved.
#' Start by arranging x and y coordinates and then reorder such that one always goes to the next closest point
#'
#' @param iso_unordered A data.frame with columns xp and yp possibly in unfavourable order
#' @param method handed to [stats::dist()] to calculate distance
#'
#' @return Arranged isobole
#'
#' @importFrom stats dist
#'
#' @author Daniel Lill, IntiQuan \email{daniel.lill@@intiquan.com}
#' @family Grid refinement functions for isoboles
#'
#' @examples
#' library(data.table)
#' library(ggplot2)
#'
#' iso_unordered <- data.frame(
#'   xp = c(00, 10, 10, 20, 30, 40, 50, 50, 50, 50, 60, 60, 60, 50, 50, 50),
#'   yp = c(80, 80, 70, 60, 60, 70, 70, 60, 50, 40, 30, 40, 20, 20, 10, 0)
#' )
#' iso_unordered <- iso_unordered[sample(1:nrow(iso_unordered), nrow(iso_unordered)),]
#' iso_unordered <- data.table(iso_unordered)
#' iso_ordered <- populationIsoboles:::arrange_isobole_internal(iso_unordered)
#' grid <- data.table(expand.grid(x = populationIsoboles:::seqminmax(0,100,10),
#'                                y = populationIsoboles:::seqminmax(0,100,10)))
#' grid[,`:=`(inpoly_unordered = as.numeric(as.logical(sp::point.in.polygon(
#'   x, y,
#'   c(iso_unordered[,xp], max(x), max(x), 0     , iso_unordered[,xp][1]),
#'   c(iso_unordered[,yp], 0,      max(y), max(y), iso_unordered[,yp][1]),
#'   TRUE))),
#'            inpoly_ordered = as.numeric(as.logical(sp::point.in.polygon(
#'   x, y,
#'   c(iso_ordered[,xp], max(x), max(x), 0     , iso_ordered[,xp][1]),
#'   c(iso_ordered[,yp], 0,      max(y), max(y), iso_ordered[,yp][1]),
#'   TRUE))))]
#' ggplot()  + geom_tile(aes(x,y, fill = inpoly_unordered), grid) +
#'   geom_path(aes(xp,yp),iso_unordered, color = "white")
#' ggplot()  + geom_tile(aes(x,y, fill = inpoly_ordered), grid) +
#'   geom_path(aes(xp,yp),iso_ordered, color = "white")
arrange_isobole_internal <- function(iso_unordered, method = "euclidean") {
  iso_arranged <- copy(iso_unordered)
  iso_arranged <- iso_arranged[order(xp,-yp)]
  mydist       <- stats::dist(iso_arranged[,list(xp,yp)], method = method)
  dmat                  <- matrix(0, nrow = nrow(iso_arranged), ncol = nrow(iso_arranged))
  dmat[lower.tri(dmat)] <- mydist
  dmat                  <- dmat + t(dmat)
  diag(dmat)            <- 10000
  heuristic_order <- c(1)
  i               <- 1
  for(i in 1:(nrow(iso_arranged)-1)) {
    j <- heuristic_order[i]
    dmat[heuristic_order[i],] <- 10000
    heuristic_order[i+1] <- which.min(dmat[,j])
  }
  out <- iso_arranged[heuristic_order]
  out <- data.table(out)
  attr(out, "order") <- heuristic_order
  out
}


#' Calculate the isobole path length
#'
#' @param iso data.table(xp,yp, ...)
#' @param FLAGreturnDT FALSE: return total path length (numeric)
#'   TRUE: Return data.table(xp,yp, distToNext, cumDist, totDist, ...)
#'
#' @return Depends on FLAGreturnDT: TRUE: data.table, FALSE: numeric(1L)
#' @examples
#' iso_unordered <- data.frame(
#' xp = c(00, 10, 10, 20, 30, 40, 50, 50, 50, 50, 60, 60, 60, 50, 50, 50),
#' yp = c(80, 80, 70, 60, 60, 70, 70, 60, 50, 40, 30, 40, 20, 20, 10, 0)
#' )
#' set.seed(1)
#' iso_unordered <- iso_unordered[sample(1:nrow(iso_unordered), nrow(iso_unordered)),]
#' iso_unordered <- data.table(iso_unordered)
#' iso_ordered <- populationIsoboles:::arrange_isobole_internal(iso_unordered)
#' populationIsoboles:::getPathlength(iso_unordered)
#' populationIsoboles:::getPathlength(iso_ordered)
getPathlength <- function(iso, FLAGreturnDT = FALSE) {
  iso <- copy(iso)
  iso[,`:=`(distToNext = c(sqrt(diff(xp)^2 + diff(yp)^2), 0))]
  if(!FLAGreturnDT) return(sum(iso$distToNext))
  iso[,`:=`(cumDist = cumsum(distToNext))]
  iso[,`:=`(totDist = sum(distToNext))]
  iso
}



#' Remove points which are causing trouble
#'
#' Sometimes arrange_isobole misses a point and then this point causes trouble
#' Identify them by finding points which have two large steps directly after each other
#'
#' @param iso data.table(xp,yp) arranged with [arrange_isobole]
#'
#' @return data.table(xp,yp) with potential outliers removed
#' @md
#'
#' @importFrom stats median
#'
#' @examples
#' # library(data.table)
#' # library(ggplot2)
#' #
#' # iso_outliers <- data.table(DATAPIOutlierExample1)
#' # iso_noOutliers <- removeOutliers(iso_outliers)
#' # PIggplot(mapping = aes(xp,yp)) +
#' # geom_path(data = iso_noOutliers, color = "black", linetype = 1) +
#' # geom_path(data = iso_outliers, color = "red", linetype = 2)
#' #
#' # iso_outliers <- data.table(DATAPIOutlierExample2)
#' # iso_noOutliers <- removeOutliers(iso_outliers)
#' # PIggplot(mapping = aes(xp,yp)) +
#' # geom_path(data = iso_noOutliers, color = "black", linetype = 1) +
#' # geom_path(data = iso_outliers, color = "red", linetype = 2)
removeOutliers <- function(iso) {

  iso <- getPathlength(iso, TRUE)
  # Identify all points with large steps
  largesteps <- (1:nrow(iso))[iso$distToNext > 4*stats::median(iso$distToNext)]
  if ((nrow(iso)-1) %in% largesteps) largesteps <- c(largesteps, nrow(iso))
  singlepoints <- largesteps[c(0, diff(largesteps)) == 1]

  if (!length(singlepoints)) return(iso)
  cat(length(singlepoints), " points are removed\n")
  iso_reduced <- iso[-singlepoints]
  iso_reduced <- getPathlength(iso_reduced, TRUE)
  if (unique(iso_reduced$totDist) < unique(iso$totDist)) iso <- iso_reduced
  iso[,list(xp,yp)]
}


#' Arrange the isobole more robustly
#'
#' Call arrange_isobole_internal in multiple ways to be a bit more robust.
#' Use either different metric or a different scale.
#' For all scale trafos, reverse them at the end to get the original isobole back
#'
#' 1. normal from x = 0 to y = 0
#' 2. log scale
#' 3. order from y = 0 to x = 0
#' 4. log scale from y = 0 to x = 0
#' 5. normal, L1 metric
#' 6. normal, canberra metrix (see \code{\link[stats]{dist}})
#'
#' Then, take the path with minimal path length, see \code{\link{getPathlength}}
#'
#' @param iso_unordered unordered data.table(xp,yp)
#' @param yps an offset for zero-values
#'
#' @return an ordered isobole
arrange_isobole <- function(iso_unordered, yps = 1e-9) {

  if (iso_unordered[1,xp]>iso_unordered[1,yp]&iso_unordered[.N,xp]<iso_unordered[.N,yp])
    iso_unordered <- iso_unordered[.N:1]

  # x to y
  iso1 <- copy(iso_unordered)
  iso1 <- arrange_isobole_internal(iso1)
  # y to x
  iso2 <- copy(iso_unordered)
  iso2[,`:=`(xp = yp, yp = xp)]
  iso2 <- arrange_isobole_internal(iso2)
  iso2[,`:=`(xp = yp, yp = xp)]
  iso2 <- iso2[.N:1]


  isolist <- list(iso_unordered, iso1,iso2)
  pathlengths <- vapply(isolist, function(.x) getPathlength(.x), FUN.VALUE = 1.)
  iso <- isolist[[which.min(pathlengths)]]

  iso <- removeOutliers(iso)
  iso[,list(xp,yp)]
}

#' Get area enclosed by isobole and right part of the grid
#'
#' This function calculates the fraction of the grid which is on the right hand side of the isobole, i.e.
#' outside of the polygon enclosed by the isobole and the coordinate axes.
#'
#' @param path data.table(xp,yp), e.g. output from [get_isobole_path()]
#' @param gridmin,gridmax,gridlens,imax See [fastIsoboles_iteration()] or [gridInfo_default()]
#' @param i Current iteration to subset gridlens
#'
#' @return A number between 0 and 1
#'
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
#' @importFrom sp point.in.polygon
#'
#' @examples
#' # Only approximate, gets better with higher imax
#' path <- data.table(xp = c(0,0.5,0.5), yp = c(0.5,.5,0))
#' imax <- 10
#' gridInfo <- populationIsoboles:::gridInfo_default(1,1,imax = imax)
#' i <- 9
#' do.call(populationIsoboles:::get_area, c(list(path =path, imax = imax, i = i), gridInfo))
get_area <- function(path, i, gridmin, gridmax, gridlens, imax){
  # Looks a bit complicated, was implemented like this to do a good interpolation for small iterations
  finegrid <- data.table::data.table(expand.grid(
    x = seqminmax(gridmin["x"], gridmax["x"], min(gridlens[[i]]["x"], (gridmax["x"] - gridmin["x"])/(2^max(7, imax)))),
    y = seqminmax(gridmin["y"], gridmax["y"], min(gridlens[[i]]["y"], (gridmax["y"] - gridmin["y"])/(2^max(7, imax))))))
  finegrid[,mean(as.numeric(as.logical(sp::point.in.polygon(
    point.x = x, point.y = y,
    pol.x = c(path[,xp], max(x), max(x), 0     , path[,xp][1]),
    pol.y = c(path[,yp], 0,      max(y), max(y), path[,yp][1]),
    mode.checked = TRUE))))]
}



#' Initialize the grid refinement algorithm
#'
#' In iteration 0, a dummy-grid and a dummy-path object are returned which will be used in iteration 1 to calculate objvalue values.
#' The algorithm evaluates only objvalue points which are close to the path of the isobole and haven't yet been evaluated.
#' The dummy path contains all grid points ensuring that all points of the grid in iteration 1 lie close to the path.
#' The dummy grid contains contains all grid points but evaluatedCum is set to zero, therefore the points will be evaluated next time.
#'
#' @param gridmin,gridmax as in fastIsoboles_iteration
#' @param gridlen_i,i as in fastIsoboles_iteration. gridlen_i and i are only there for backward compatibility
#'
#' @return list with grid and path as explained in the Arguments-section of fastIsoboles_iteration
#'
#' @author Daniel Lill, IntiQuan \email{daniel.lill@@intiquan.com}
#' @family Grid refinement functions for isoboles
#'
#' @importFrom data.table data.table copy ":="
#'
#' @examples
#' populationIsoboles:::fastIsoboles_iteration0(c(x = 0,y = 0), c(x = 16,y = 16), c(x = 8,y = 8))
fastIsoboles_iteration0 <- function(gridmin, gridmax, gridlen_i = (gridmax-gridmin)/2, i = 0) {
  grid <- data.table(expand.grid(x = seqminmax(gridmin["x"],gridmax["x"],gridlen_i["x"]), y = seqminmax(gridmin["y"],gridmax["y"],gridlen_i["y"])))
  grid[,`:=`(evaluatedCum = 0, objvalueCum = NA, objvalue0 = NA, evaluated0 = 0)]

  path <- data.table::copy(grid)
  path <- path[,list(xp = x, yp = y)]

  list(grid = grid, path = path, outsideGrid = FALSE)
}

#' Do one iteration of the grid refinement algorithm
#'
#' @param gridmin,gridmax         General grid parameters: Vectors with names x and y c(x = 0, y = 0)
#' @param gridlen_i               Grid lengths of iteration i: c(x = 40, y = 40)
#' @param path_j,grid_j Values from iteration i-1.
#'   path_j: data.table with columns (xp,yp) with nodes of iteration j's approximation of the isobole path
#'   grid_j: data.table with columns (x,y,evaluatedCum, objvalueCum, evaluated1, ..., evaluatedj, objvalue1,..., objvaluej)
#'   The *Cum columns contain all cumulated grid point evaluations of all iterations.
#'   evaluatedCum has a 1, if the objvalue has been calculated at that point, or a zero if this value has not been calculated
#'   The *j columns contain the information about iteration j.
#'   evaluated2 has a 1 for all grid points which were evaluated in iteration 2.
#' @param i                       Identifier for iteration i
#' @param objfun                Function(x,y) returning a vector of the same length with the corresponding objvalue values.
#'  Combined with objvalue, this could be called the objective function
#' @param objvalue              The value for which the contour lines are calculated. E.g. 0.95
#'
#' @return list(grid, path), as specified in the arguments `grid_j` and `path_j`
#'
#' @author Daniel Lill, IntiQuan \email{daniel.lill@@intiquan.com}
#' @family Grid refinement functions for isoboles
#' @importFrom data.table data.table copy ":="
#'
#'
#' @examples
#' # see fastIsoboles()
fastIsoboles_iteration <- function(gridmin, gridmax, gridlen_i,
                                      path_j, grid_j,
                                      i, objfun, objvalue = 0.95) {

  # .. 0 Deep copy so the external stuff isn't affected ----- #
  grid_j <- data.table::copy(grid_j)
  path_j <- data.table::copy(path_j)
  # .. 1 Calculate new values based on old path ----- #
  # .... 1 Create new grid ------ #
  grid <- data.table::data.table(expand.grid(x = seqminmax(gridmin["x"],gridmax["x"],gridlen_i["x"]), y = seqminmax(gridmin["y"],gridmax["y"],gridlen_i["y"])))
  grid[,`:=`(gridrowid = .I)]
  # .... 2 Determine points close to old isobole ------ #
  # Test every combination of (path-node, grid-point): Merge the two data.frames with x and y coordinates and calculate norm
  # This might appear very inefficient, but is not the bottleneck
  path_j[,`:=`(pathrowid = .I)]
  pathgrid_grid <- data.table::data.table(expand.grid(pathrowid = path_j$pathrowid, gridrowid = grid$gridrowid))
  pathgrid_grid <- pathgrid_grid[path_j, on = "pathrowid"][grid, on = "gridrowid"]
  pathgrid_grid[,`:=`(dist = sqrt(((xp-x)/(gridmax["x"] - gridmin["x"]))^2 + ((yp-y)/(gridmax["y"] - gridmin["y"]))^2))]
  pathgrid_grid[,`:=`(closetopath = dist < sqrt(sum((gridlen_i/(gridmax - gridmin))^2)))]

  close_points <- unique(pathgrid_grid[closetopath == TRUE, gridrowid])

  # .... 3 Join old and new grid ------ #
  grid <- grid_j[grid, on = c("x", "y")]
  grid[is.na(evaluatedCum),`:=`(evaluatedCum = 0)]
  # .... 4 Calculate new values if grid point was never evaluated ------ #
  grid[gridrowid %in% close_points & evaluatedCum == 0, c(paste0("objvalue",i), paste0("evaluated",i)) := list(objfun(x,y), 1)]

  # .... 5 Update objvalueCum and evaluatedCum ------ #
  # Need to sum rowwise over different columns => need apply
  evalCum <- apply(as.matrix(grid[,.SD,.SDcols = grep("evaluated(?!Cum)", names(grid), value = TRUE, perl = TRUE)]), 1, function(.x) sum(.x, na.rm = TRUE))
  grid[,`:=`(evaluatedCum = evalCum)]
  ACum <- apply(as.matrix(grid[,.SD,.SDcols = grep("objvalue(?!Cum)", names(grid), value = TRUE, perl = TRUE)]), 1, function(.x) {
    if (all(is.na(.x))) return(NA)
    .x[!is.na(.x)]})
  grid[,`:=`(objvalueCum = ACum)]
  # .... 6 Delete obsolete rows and columns ----- #
  grid <- grid[evaluatedCum == 1]
  grid[,c("gridrowid") := NULL]

  # .. 2 Calculate new path ----- #
  isoPath <- get_isobole_path(grid, gridmin, gridmax, gridlen_i, objvalue)

  # ggplot() +
  #   geom_tile(data = grid, aes(x,y,fill = objvalueCum)) +
  #   geom_label(data = grid, aes(x,y,label = round(objvalueCum,2))) +
  #   geom_path(data = path, aes(xp,yp)) +
  #   geom_blank()
  list(grid = grid, path = isoPath$path, outsideGrid = isoPath$outsideGrid)
}






#' Grid refinement algorithm to calculate isoboles
#'
#' @inheritParams fastIsoboles_iteration
#' @param imin,imax minimum/maximum number of iterations
#' @param k index of population for outputting the current isobole iteration to the disk
#' @param areaterm Terminate when the area above the curve (bounded by the grid) doesn't change more than this value (on relative scale)
#' @param FLAGverbose Print messages to the console
#' @param .outputFolder Path to store results in. See \code{\link{output_gridrefinement_results}}
#'
#'
#' @return list with
#'   * grid: Exact objvalue values where the objvalue has been evaluated
#'   * path: The last isobole path
#'   * pathlist: The paths of all iterations
#'   * itermax: The number of iterations used
#'   * areas: The areas above the isobole (To analyze the convergence)
#'   * times: proc.time - stamps. The first entry is before iteration 0, the second is after iteration 1 etc.
#' @export
#' @author Daniel Lill, IntiQuan \email{daniel.lill@@intiquan.com}
#' @family Grid refinement functions for isoboles
#' @importFrom sp point.in.polygon
#' @importFrom data.table data.table
#'
#' @examples
#' library(populationIsoboles)
#' library(data.table)
#' library(ggplot2)
#'
#' # Two drugs on different scales
#' objfun1 <- function(x,y) {sqrt((0.001)^2 * x^2 + (0.1)^2 * y^2)}
#' gridmin <- c(x = 0, y = 0) + 0.0001
#' gridmax <- c(x = 10*2^7, y = 0.1*2^7) + 0.0001
#'
#' # Quick plot of function
#' testobjvalue = data.table(expand.grid(
#'   x = populationIsoboles:::seqminmax(
#'       gridmin["x"], gridmax["x"], (gridmax["x"] - gridmin["x"])/(2^6)),
#'   y = populationIsoboles:::seqminmax(
#'       gridmin["y"], gridmax["y"], (gridmax["y"] - gridmin["y"])/(2^6))))
#' testobjvalue[,`:=`(objvalue1 = objfun1(x,y))]
#' ggplot(testobjvalue, aes(x,y, z = objvalue1)) +
#'   geom_contour(aes(color = ..level..))
#'
#' # Run algorithm
#' isobole <- fastIsoboles( objfun1,0.95,gridmin, gridmax, 5,7)
#' # Plot algo results
#' ggplot() + geom_tile(data = isobole$grid, aes(x,y,fill = objvalueCum)) +
#'   geom_path(data = isobole$pathlist[[5]], aes(xp,yp)) +
#'   geom_path(data = isobole$pathlist[[1]], aes(xp,yp), linetype = 2) +
#'   scale_fill_viridis_c()
#'
#' # Plot grid iterations
#' grd <- copy(isobole$grid)
#' grd[,`:=`(objvalue0 = NULL, evaluated0 = NULL, objvalueCum = NULL, evaluatedCum = NULL)]
#' grd <- melt(grd,c("x", "y"), measure = patterns("^evaluated", "objvalue"),
#'             variable.name = "iteration", value.name = c("evaluated","objvalue"),
#'             variable.factor = FALSE)
#' grd[,`:=`(iteration = as.numeric(as.character(iteration)))]
#' grd <- grd[!is.na(objvalue)]
#' grd[evaluated == 0 ,`:=`(evaluated = NA)]
#' ggplot() +
#'   geom_tile(data = grd, aes(x,y,fill = factor(evaluated * iteration))) +
#'   geom_point(data = grd, aes(x,y,shape = factor(evaluated * iteration))) +
#'   geom_path(data = isobole$pathlist[[5]], aes(xp,yp))
#'
#' # Run algorithm - different objvalue - Level
#' isobole <- fastIsoboles(objfun1, 0.5,gridmin, gridmax, 5,7)
#' # Plot algo results
#' ggplot() + geom_tile(data = isobole$grid, aes(x,y,fill = objvalueCum)) +
#'   geom_path(data = isobole$pathlist[[5]], aes(xp,yp)) + scale_fill_viridis_c()
#'
#' # Run algorithm - lower imin
#' isobole <- fastIsoboles(objfun1,0.95, gridmin, gridmax, 3,7)
#' isobole$itermax
#' # Plot algo results
#' ggplot() + geom_tile(data = isobole$grid, aes(x,y,fill = objvalueCum)) +
#'   geom_path(data = isobole$pathlist[[3]], aes(xp,yp)) + scale_fill_viridis_c()
#'
#' # Run algorithm - lower imin & imax
#' isobole <- fastIsoboles(objfun1,0.95, gridmin, gridmax, 3,3)
#' isobole$itermax
#' # Plot algo results
#' ggplot() + geom_tile(data = isobole$grid, aes(x,y,fill = objvalueCum)) +
#'   geom_path(data = isobole$pathlist[[3]], aes(xp,yp)) + scale_fill_viridis_c()
#'
#' # Run algorithm - Different function
#' objfun2 <- function(x,y) {(x+500)*(y+5)/(80^2)}
#' isobole <- fastIsoboles(objfun2,0.95, gridmin, gridmax, 3,7)
#' isobole$itermax
#' # Plot algo results
#' ggplot() + geom_tile(data = isobole$grid, aes(x,y,fill = objvalueCum)) +
#'   geom_path(data = isobole$pathlist[[3]], aes(xp,yp)) + scale_fill_viridis_c()
#'
fastIsoboles <- function(objfun, objvalue = 0.95,
                            gridmin = c(x=0,y=0),
                            gridmax = c(x=1,y=1),
                            imin = 5, imax = 7,
                            FLAGverbose = FALSE,
                            k = NULL, .outputFolder = NULL,
                            areaterm = 1.01
) {
  # 1 Initialize
  gridlens    <- get_gridlens(gridmin, gridmax, imax)
  dummy       <- fastIsoboles_iteration0(gridmin, gridmax, gridlen_i = gridlens[[1]], 0)
  grid        <- dummy$grid
  path        <- dummy$path
  outsideGrid <- dummy$outsideGrid
  areas       <- c()
  pathlist    <- list()

  # 2 Iterate
  times <- list(proc.time())
  for (i  in  1:imax) {
    if (FLAGverbose) cat(paste(Sys.time()), " ---- Iteration", i, "----------")
    # .. 1 Calculate new values
    if (FLAGverbose) cat("[-   ]")
    dummy <- fastIsoboles_iteration(gridmin = gridmin, gridmax = gridmax, gridlen_i = gridlens[[i]],
                                       path_j = path, grid_j = grid,
                                       i = i, objfun  = objfun, objvalue = objvalue)
    grid <- dummy$grid
    pathlist[[i]] <- path <- dummy$path
    outsideGrid   <- dummy$outsideGrid

    # .. 2 Check if isobole has converged: Calculate fraction of grid whole grid which is above isobole and check for stationarity
    if (FLAGverbose) cat("\b\b\b\b\b\b[--  ]")
    areas[i] <- get_area(path = path, i = i,
                         gridmin = gridmin, gridmax = gridmax,
                         gridlens = gridlens, imax = imax)
    times[[i+1]] <- proc.time()
    if (FLAGverbose) cat("\b\b\b\b\b\b[----]\n")
    if(imax == 1 || outsideGrid ||
       (i >= imin && abs(log10(areas[i]/areas[i-1])) < log10(areaterm)))
      break
  }

  # Collect output
  isobole <- list(grid = grid, pathlist = pathlist, itermax = i, areas = areas,
                  times = times, outsideGrid = outsideGrid)
  # 4 Write isobole to disk if k and .outputFolder is supplied
  if (length(k) && length(.outputFolder)){
    if (FLAGverbose) cat("writing results ... ")
    gridInfo <- list(gridmin = gridmin, gridmax = gridmax,
                     gridlens = get_gridlens(gridmin = gridmin,
                                             gridmax = gridmax, imax = imax))
    output_gridrefinement_results(isobole, k, .outputFolder)
    if (FLAGverbose) cat(" done \n")
  }

  # 5 Return
  isobole
}

# .. Termination helpers -----

# > Functions for a better termination criterion than "areaterm".
# Not yet implemented


# #' Helper function to construct a fine grid
# #'
# #' @param gridmin,gridmax from gridInfo
# #'
# #' @return data.table(x,y)
# #' @export
# #'
# #' @importFrom data.table data.table
# #'
# #' @examples
# #' getFineGrid(c(x = 0, y = 0), c(x = 1, y = 1), 2)
# getFineGrid <- function(gridmin, gridmax, imax = 8) {
#     data.table::data.table(expand.grid(
#       x = seqminmax(gridmin["x"], gridmax["x"], (gridmax["x"] - gridmin["x"])/(2^imax)),
#       y = seqminmax(gridmin["y"], gridmax["y"], (gridmax["y"] - gridmin["y"])/(2^imax))))
# }
#
#
# #' Join two isoboles together to obtain a polygon
# #'
# #' Note that the second isobole is appended in reverse order to go around in a loop
# #'
# #' @param pathlist list of isoboles
# #' @param i pathlist[i] and [i-1\] are going to be used
# #'
# #' @return data.table
# #' @export
# #'
# #' @importFrom data.table rbindlist
# #'
# #' @examples
# #' pathlist <- list(
# #'   data.table(xp = seq(0,1, length.out = 10), yp = seq(0,1, length.out = 10)^2),
# #'   data.table(xp = seq(0,1, length.out = 10), yp = seq(0,1, length.out = 10)^4))
# #' # right
# #' ggplot(getJoinedPath(pathlist, 2), aes(xp,yp)) + geom_path()
# #' # wrong
# #' ggplot(rbindlist(pathlist), aes(xp,yp)) + geom_path()
# getJoinedPath <- function(pathlist, i) {
#   data.table::rbindlist(list(pathlist[[i]], pathlist[[i-1]][.N:1]))
# }
#
#
# #' Get area between two isoboles
# #'
# #' Gets area between isoboles i and i-1 of pathlist
# #'
# #' @param finegrid Output of [getFinegrid()]
# #' @param pathlist List of isoboles
# #' @param i Iteration number
# #'
# #' @return value between 0 and 1, relative area of grid covered between both isoboles
# #' @export
# #'
# #' @importFrom sp point.in.polygon
# #'
# #' @examples
# #' finegrid <- getFineGrid(c(x = 0, y = 0), c(x = 1, y = 1), 8)
# #' pathlist <- list(
# #'   data.table(xp = seq(0,1, length.out = 10), yp = seq(0,1, length.out = 10)^2),
# #'   data.table(xp = seq(0,1, length.out = 10), yp = seq(0,1, length.out = 10)^4))
# #' get_areaDiff(finegrid, pathlist, 2)
# get_areaDiff <- function(finegrid, pathlist, i) {
#   path <- getJoinedPath(pathlist = pathlist, i = i)
#   areabetween <- finegrid[,mean(as.numeric(as.logical(sp::point.in.polygon(x, y,
#                                                                            c(path[,xp], path[,xp][1]),
#                                                                            c(path[,yp], path[,yp][1]),
#                                                                            TRUE))))]
#   areabetween
# }
#
#
# #' Get area left of isobole
# #'
# #' Given relative to grid size.
# #'
# #' It worked more reliably to determine the area right of the isobole, because of the little offset which always needed to be added.
# #' Therefore, the steps are:
# #' 1 Augment isobole path by the the three corners apart from the origin and connect to beginning of isobole
# #' 2 Determine area enclosed by this polygon
# #' 3 Subtract from 1 to obtain area left of isobole
# #'
# #' @param finegrid Output of [getFinegrid()]
# #' @param path An isobole path ordered from x = 0 to y = 0
# #'
# #' @return value between 0 and 1
# #' @export
# #'
# #' @examples
# #' finegrid <- getFineGrid(c(x = 0, y = 0), c(x = 1, y = 1), 8)
# #' pathlist <- list(
# #'   data.table(xp = seq(0,1, length.out = 10), yp = seq(1,0, length.out = 10)^2),
# #'   data.table(xp = seq(0,1, length.out = 10), yp = seq(1,0, length.out = 10)^4))
# #' get_area(finegrid, pathlist[[1]]) # 1/3 = integral of x^2 from 0 to 1
# #' get_area(finegrid, pathlist[[2]]) # 1/5
# get_area <- function(finegrid, path) {
#   arearight <- finegrid[,mean(as.numeric(as.logical(sp::point.in.polygon(x, y,
#                                                             c(path[,xp], max(x), max(x), 0     , path[,xp][1]),
#                                                             c(path[,yp], 0,      max(y), max(y), path[,yp][1]),
#                                                             TRUE))))]
#   arealeft <- 1-arearight
#   arealeft
# }
#
#
# #' Determine if the isobole has converged
# #'
# #' There are two criteria which are connected via OR
# #'
# #' 1. Areadiff / Pathlength < term_dArea_Pathlen
# #' 2. Areadiff / Area       < term_dArea_Area
# #'
# #' Areadiff, Pathlength and Area are calculated on the scale free grids, i.e. gridmax(x and y) = 1.
# #' Therefore, it is insensitive to the scale of the drugs.
# #'
# #' The default values have been found out by a small simulation study.
# #' They reduce the number of early terminations to ~ 1/8, while the number of late terminations is kept at ~ 5/16.
# #' Note that most early terminations also capture the final isobole shape reasonably well.
# #' Therefore, having some early terminations is not critical
# #'
# #' @param gridmin,gridmax From gridInfo
# #' @param pathlist List of isobole paths
# #' @param i Current iteration number
# #' @param term_dArea_Pathlen,term_dArea_Area Values to test the criteria agains
# #'
# #' @return TRUE = isobole has converged, FALSE = isobole has not converged
# #' @export
# #'
# #' @examples
# #' pathlist <- list(
# #'   data.table(xp = seq(0,1, length.out = 20), yp = seq(1,0, length.out = 20)^2),
# #'   data.table(xp = seq(0,1, length.out = 20), yp = seq(1,0, length.out = 20)^3),
# #'   data.table(xp = seq(0,1, length.out = 20), yp = seq(1,0, length.out = 20)^4),
# #'   data.table(xp = seq(0,1, length.out = 20), yp = seq(1,0, length.out = 20)^5))
# #' ggplot(rbindlist(pathlist, idcol = "i"), aes(xp,yp, color = factor(i))) + geom_line()
# #' get_termination(c(x = 0, y = 0), c(x = 1, y = 1), pathlist, 2)
# #' get_termination(c(x = 0, y = 0), c(x = 1, y = 1), pathlist, 3)
# #' get_termination(c(x = 0, y = 0), c(x = 1, y = 1), pathlist, 4)
# get_termination <- function(gridmin, gridmax, pathlist, i, term_dArea_Pathlen = 0.004, term_dArea_Area = 0.04) {
#   finegrid <- getFineGrid(gridmin = gridmin, gridmax = gridmax)
#   area <- get_area(finegrid = finegrid, path = pathlist[[i]])
#   areadiff <- get_areaDiff(finegrid = finegrid, pathlist = pathlist, i = i)
#   adiffByLeftArea <- areadiff/area
#   adiffByRightArea <- areadiff/(1-area)
#
#   iso_scalefree <- copy(pathlist[[i]])
#   iso_scalefree[,`:=`(xp = xp/(gridmax["x"] - gridmin["x"]), yp = yp/(gridmax["y"] - gridmin["y"]))]
#   totDist_scalefree <- getPathlength(iso_scalefree)
#   adiffByPathlength <- areadiff/totDist_scalefree
#
#   tAdiffPL <- adiffByPathlength < term_dArea_Pathlen
#   tAdiffALeft  <- adiffByLeftArea < term_dArea_Area
#   tAdiffARight  <- adiffByRightArea < term_dArea_Area
#
#   terminated <- tAdiffPL | tAdiffALeft | tAdiffARight
#   terminated
# }

# -------------------------------------------------------------------------#
# File interactions ----
# -------------------------------------------------------------------------#

#' Output the results of the grid-refinement algorithm in a standardized manner
#'
#' @md
#' @param isobole result of [fastIsoboles()]
#' @param k Identifier of the current population
#' @param .outputFolder output folder where the results are stored. Will save the results in a subdirectory called "Simulate
#'
#' @return Null. This function is called for its side-effect.
output_gridrefinement_results <- function(isobole, k, .outputFolder) {
  if (!dir.exists(file.path(.outputFolder, "Simulations")))
    dir.create(file.path(.outputFolder, "Simulations"), recursive = TRUE)
  files <- PI_files(.outputFolder, k)
  saveRDS(isobole$pathlist,    files$paths)
  saveRDS(isobole$grid,        files$grid)
  saveRDS(isobole$times,       files$times)
  saveRDS(isobole$itermax,     files$imax)
  saveRDS(isobole$areas,       files$areas)
  saveRDS(isobole$outsideGrid, files$outsideGrid)
  NULL
}


#' Generate default output file paths
#'
#' @param .outputFolder results will be saved in file.path(.outputFolder, "Simulations")
#' @param k Index of current population
#'
#' @return list of file paths
#' @examples
#' populationIsoboles:::PI_files("Output", 1)
PI_files <- function(.outputFolder, k) {
  list(
    paths           = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), "_01_pathlist.rds")) ,
    grid            = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), "_02_grid.rds")),
    times           = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), "_03_times.rds")),
    imax            = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), "_04_itermax.rds")),
    areas           = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), "_05_areas.rds")),
    outsideGrid     = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), "_10_outsideGrid.rds")),
    parameters      = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), "_11_parameters.rds")),
    zipFile         = file.path(.outputFolder, "Simulations", paste0(sprintf("%04i", k), ".zip"))
  )
}


#' Zip files belonging to an isobole simulation
#'
#' @param files List of filenames, result from [PI_files()]
#'
#' @return NULL (called for side-effect)
#'
#' @author Daniel Lill (daniel.lill@intiquan.com)
#' @md
PI_zip <- function(files) {
  filesToZip <- as.character(c(files, lapply(files, paste0, ".log")))
  filesToZip <- filesToZip[sapply(filesToZip, function(fName) {fName != files$zipFile && file.exists(fName)})]
  zip::zipr(files$zipFile, files = filesToZip, recurse = FALSE)
  unlink(filesToZip)
  NULL
}

#' Unzip files with pre-specified unzipping options
#'
#' @param files output from [PI_files()]
#'
#' @return NULL, called for side-effect
PI_unzip <- function(files) {
  if (!file.exists(files$zipFile)) return(NULL)
  zip::unzip(files$zipFile, exdir = dirname(files$zipFile))
  NULL
}

