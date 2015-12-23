# flywatch
Code to analyse movies of flies. 

For the time being principally based on Yoshi Aso's (Janelia) behaviour rig as originally described in 

*Mushroom body output neurons encode valence and guide
memory-based action selection in Drosophila.*
Aso Y, Sitaraman D, Ichinose
T, Kaun KR, Vogt K, Belliart-Guerin G, Placais PY, Robie AA, Yamagata N,
Schnaitmann C, Rowell WJ, Johnston RM, Ngo TT, Chen N, Korff W, Nitabach
MN, Heberlein U, Preat T, Branson KM, Tanimoto H, Rubin GM.. **eLIFE**
2014 Dec 23;3:e04580. [doi:10.7554/eLife.04580](http://dx.doi.org/10.7554/eLife.04580)

## Quick Start

For the impatient ...

```r
# install
if (!require("devtools")) install.packages("devtools")
devtools::install_github(c("jefferis/nat", "jefferislab/flywatch"))

# use
library(flywatch)

# run examples for plotting data from Yoshi Behaviour Rig
example("plot_ybr")
example("plot_smoothed_displacement")

# get overview help for package
?flywatch

# help for some interesting functions
?read_ybr_summary
?read_ybr_xy
?plot_ybr
```

Note: Windows users need [Rtools](http://www.murdoch-sutherland.com/Rtools/) to 
install with devtools.
