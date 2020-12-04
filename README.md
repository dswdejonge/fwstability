# fwstability
R-package to calculate food web (fw) stability from an energy-flux model.

# Installation
```r
# To install the development version of the fwstability package:
devtools::install_github("dswdejonge/fwstability", build_vignettes = TRUE)

# You can also install the fwmodel package with 
# published energy-flux models to use as examples:
devtools::install_github("dswdejonge/fwmodels")

# Read vignette:
browseVignettes("fwstability") # click HTML
```

# Quick start
Please review the vignette for elaborate documentation and examples. 
The most important functions in this package are getJacobian and getStability:
```r
JM <- getJacobian(model)
s <- getStability(JM)
```
