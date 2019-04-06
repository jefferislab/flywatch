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
devtools::install_github("jefferislab/flywatch")

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

## Bugs, Feature Requests
If it doesn't work or is missing a feature you can check if anyone else has 
complained and file an issue at github

* https://github.com/jefferislab/flywatch/issues


## Dolan analysis scripts

This repository also wraps a number of analysis scripts written by 
Mike Dolan for the paper originally preprinted as 

Neurogenetic dissection of the Drosophila innate olfactory processing center
Michael-John Dolan, Shahar Frechter, Alexander Shakeel Bates, Chuntao Dan, Paavo Huoviala, Ruairi J.V. Roberts, Philipp Schlegel, Serene Dhawan, Remy Tabano, Heather Dionne, Christina Christoforou, Kari Close, Ben Sutcliffe, Bianca Giuliani, Feng Li, Marta Costa, Gudrun Ihrke, Geoffrey Meissner, Davi Bock, Yoshinori Aso, Gerald Rubin, Gregory Jefferis
bioRxiv 404277; doi: https://doi.org/10.1101/404277

most of these scripts don't actually depend on the package itself, but
can be found in the [demo](demo) folder.
