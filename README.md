# DASE: Differnetial Analysis of Super Enhancers
R package used to categorize and rank super-enhancer (SE) internal differences in two conditions.

## Introduction
Super enhancers (SEs) were proposed as broad regulatory domains on genome, usually spanning a minimum
of thousands of base pairs and consisting of multiple constitute enhancers. The constitute enhancers work together as a unit, instead of separately, to facilitate high enhancer activity. Aberrant SE activities, 
which are critical to understand disease mechanisms, could be raised by the alterations of one or more of 
their constitute enhancers. However, the state-of-art binary strategy in calling differential SEs only relies on overall activity changes, neglecting the local dynamics of constitute enhancers within SEs. DASE uses a weighted spline model to identify differential SEs from two conditions by accounting for the combinatorial effects of constitute enhancers weighted with their activities and locations (internal dynamics). In addition to overall changes, our medthod finds four novel types (*Shortened/lengthened*, *shifted*, *hollowed* and *other complex scenarios*) of differential SEs pointing to the structural differences within SEs.

## Install

DASE is a R package, it can be installed with source code documented in [GitHub](https://github.com/tenglab/DASE).

To install DASE from GitHub, using the following command. It should also install all the dependencies DASE required. A complete list of dependencies can be found in the DESCRIPTION file.

```R
# Public package
devtools::install_github("https://github.com/tenglab/DASE.git")
```

## Using DASE
First load DASE,
```R
library(DASE)
```

Then follow the [DASE user's guild](http://127.0.0.1:13324/library/DASE/doc/DASE.html) to preform differneital analysis of SEs.

## Citation
If you use DASE, please cite our paper at (https://www.biorxiv.org/content/10.1101/2021.09.25.461810v1)

## Help
Feel free to leave any questions and bugs at [GitHub issues](https://github.com/tenglab/DASE/issues).
