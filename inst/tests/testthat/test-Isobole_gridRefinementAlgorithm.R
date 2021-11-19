test_that("grid_refinement works", {
  library(populationIsoboles)

  # Example 1 from ?grid_refinement
  objfun1 <- function(x,y) {sqrt((0.001)^2 * x^2 + (0.1)^2 * y^2)}
  gridmin <- c(x = 0, y = 0) + 0.0001
  gridmax <- c(x = 10*2^7, y = 0.1*2^7) + 0.0001
  isobole <- grid_refinement(objfun1, 0.95, gridmin, gridmax, 5,7)
  # Compare to previous run
  isobole_expectation <- readRDS("testdata/Expectations/test-grid_refinement.rds")
  expect_equivalent(isobole$grid, isobole_expectation$grid)

})
