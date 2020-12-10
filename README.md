# fwstability
To infer food-web (fw) stability from energy-flux food-web models.  

Energy-flux food-web models comprise of compartments of biomass and flows of energy or material between those compartments. This type of food-web model is often employed for ecosystem wide approaches to an ecological problem. This package encapsulates published methods to infer stability from such food-web models. These methods translate the steady-state models to the dynamic domain, and then infer interaction strengths between all compartments. Subsequently, the mathematical concept of stability of differential systems is used to infer food-web stability. This package is useful if you are comparing different ecosystems, or want to study the effect of disturbances on the food-web.  

Please cite if you are using this package:  

Daniëlle S.W. de Jonge, Peter C. de Ruiter, Johan van de Koppel, Dick van Oevelen. Inferring stability of energy-flux food-web models. In preparation.  

# Important literature
The following list of references are central in this package. For a full overview of relevant literature, please review the vignette.

Arnoldi, J.F., Loreau, M., Haegeman, B., 2016. Resilience, reactivity and variability: A mathematical comparison of ecological stability measures. J. Theor. Biol. 389, 47–59. https://doi.org/10.1016/j.jtbi.2015.10.012

Boit, A.; Gaedke, U. (2014). Benchmarking successional progress in a quantitative food web. Plos One 9(2):e90404. doi:10.1371/journal.pone.0090404

de Ruiter, P. C., Neutel, A. M., & Moore, J. C. (1995). Energetics, Patterns of Interaction Strengths, and Stability in Real Ecosystems. Science, 269(5228), 1257–1260. https://doi.org/10.1126/science.269.5228.1257  

May, R. M. (1972). Will a large network be stable? Nature, 238, 37–38. https://doi.org/10.1038/238413a0  

Moore, J. C., Berlow, E. L., Coleman, D. C., De Ruiter, P. C., Dong, Q., Hastings, A., … Wall, D. H. (2004). Detritus, trophic dynamics and biodiversity. Ecology Letters, 7(7), 584–600. https://doi.org/10.1111/j.1461-0248.2004.00606.x  

Neutel, A. M., Heesterbeek, J. A. P., & Ruiter, P. C. De. (2002). Stability in Real Food Webs: Weak Links in Long Loops. Science, 296(5570), 1120–1123. https://doi.org/10.1126/science.1068326  

Neutel, A. M., Heesterbeek, J. A. P., Van De Koppel, J., Hoenderboom, G., Vos, A., Kaldeway, C., … De Ruiter, P. C. (2007). Reconciling complexity with stability in naturally assembling food webs. Nature, 449(7162), 599–602. https://doi.org/10.1038/nature06154  

Neutel, A. M., & Thorne, M. A. S. (2014). Interaction strengths in balanced carbon cycles and the absence of a relation between ecosystem complexity and stability. Ecology Letters, 17(6), 651–661. https://doi.org/10.1111/ele.12266  

Ulanowicz, R.E., 1997. Limitations on the connectivity of ecosystem flow networks. In: Rinaldo, A.; Marani, A. (eds) Biological models. Instituto Veneto de Scienze, Lettre ed Arti, Venica, pp 125-143.

van Altena, C., Hemerik, L., de Ruiter, P.C. (2016). Food web stability and weighted connectance: the complexity-stability debate revisited. Theor. Ecol. 9, 49–58. https://doi.org/10.1007/s12080-015-0291-7

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

The model input has to be a fully resolved energy-flux model, with biomasses, a flow matrix, energy conversion efficiencies, and information on any dead (e.g. detrital) and external compartments.  

The Jacobian matrix inferred from the food-web model will obtain values interpreted as interaction strengths between different compartments.  

The stability is the asymptotic rate back to equilibrium state after a small disturbance, but can also be used to infer return time and biological tipping points of the system.
