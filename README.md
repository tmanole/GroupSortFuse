# GroupSortFuse
Implementation of the Group-Sort-Fuse (GSF) procedure for estimating the number of components in finite mixture models. Four families of mixture models are currently implemented:
 + `normalLocOrder`: Multidimensional Gaussian mixture models in location, with common but possibly unknown scale parameter;
 + `multinomialOrder`: Multinomial mixture models;
 + `poissonOrder`: Univariate Poisson mixture models;
 + `exponentialOrder`: Exponential distribution mixture models.

This package continues to be under development, and has only been tested on Ubuntu 16.04 with R  3.4.4. This release should not be considered stable. 

# Installation
This package may be installed as follows, using the `devtools` R package:
```{r}
library(devtools)
devtools::install_github("tmanole/GroupSortFuse")
```


