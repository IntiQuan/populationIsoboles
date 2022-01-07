# MMVIsoboles


## Installation
```{r}
if (!require("devtools")) install.packages("devtools")
# Install without vignettes accessible via utils::vignette()
devtools::install_github("IntiQuan/populationIsoboles")
# Install with vignettes accessible via utils::vignette()
devtools::install_github("IntiQuan/populationIsoboles", build_vignettes = TRUE)
```

## Vignettes

* Long vignette `Population-Isoboles-DetailedSetup.Rmd`: Set up a simulation of confidence level response surfaces from scratch by manually defining the parameter sampling functions and model simulation function.
* Short vignette `Population-Isoboles.Rmd`: Use the simulation capabilites of the R-package *IQRtools* for a simplified setup.

## Main algorithms

* `fastIsoboles` See function documentation for fast examples.
* `aggregateIsoboles`


