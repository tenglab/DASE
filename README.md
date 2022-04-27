# DASE: Differnetial Analysis of Super Enhancers
An R package to categorize and rank super-enhancer (SE) differences

## Introduction
Super enhancers (SEs) are broad enhancer domains usually containing multiple constituent enhancers that 
hold elevated activities in gene regulation. Disruption in one or more constituent enhancers causes 
aberrant SE activities that lead to gene dysregulation in diseases. To quantify SE aberrations, 
differential analysis is performed to compare SE activities between cell conditions. The state-of-art 
strategy in estimating differential SEs relies on overall activities and neglect the changes in length 
and structure of SEs. DASE uses a weighted spline model to identify differential SEs between two conditions 
by accounting for the combinatorial effects of constituent enhancers weighted with their activities
and locations (internal dynamics). In addition to overall changes, our method finds four novel types 
(*shortened*, *shifted*, *hollowed* and *other complex scenarios*) of differential SEs pointing to the 
structural differences within SEs.

## Installation

DASE is an R package with its source code documented in this [repository](https://github.com/tenglab/DASE).

Please follow the steps below to install DASE in R. Dependency packages noted in the 
[DESCRIPTION](https://github.com/tenglab/DASE/blob/master/DESCRIPTION) file are required to ensure successful
installation. These packages can be found from CRAN or Bioconductor repositories. In addition, 
R package *devtools* is required to start the installation.

```R
# Install DASE
library(devtools)
devtools::install_github("https://github.com/tenglab/DASE.git")
```

Or, Install with vignettes and dependencies.

```R
# Install DASE with vignettes and dependencies
devtools::install_github("https://github.com/tenglab/DASE.git",build_vignettes = TRUE)
```

## Using DASE
First load DASE,
```R
library(DASE)
```

Then, follow the [**User's Guide**](https://github.com/tenglab/DASE/blob/master/DASE_guide.pdf) 
to perform SE differential analysis. In the guide, we detail how to prepare input files, switch input 
formats, interpret outputs files as well as some spline fitting examples.

The users are also encouraged to refer to the help pages of R functions in this package. 

## Citation
If you use DASE, please cite our paper at (https://doi.org/10.1093/nar/gkac141).

"Liu X, Zhao B, Shaw TI, Fridley BL, Duckett DR, Tan AC, Teng M. Summarizing internal dynamics boosts differential analysis and functional interpretation of super enhancers. Nucleic Acids Res. 2022 Mar 2:gkac141. doi: 10.1093/nar/gkac141."

## Help
Feel free to leave any questions and bugs at [GitHub issues](https://github.com/tenglab/DASE/issues).

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Ftenglab%2FDASE&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
