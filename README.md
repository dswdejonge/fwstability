# fwstability
R-package to calculate food web (fw) stability from an energy-flux model.

NOTE: Package is in development, please be aware when using code in its current form.

# Installation
```r
# Install the fwstability package
devtools::install_github("dswdejonge/fwstability", build_vignettes = TRUE)

# Install the fwmodel package with published energy-flux models
# to use as examples.
devtools::install_github("dswdejonge/fwmodels")
```

# Quick start
Please review the vignette for elaborate documentation and examples. 
The most important functions in this package are getJacobian and getStability:
```r
JM <- getJacobian(model)
s <- getStability(JM)
```
