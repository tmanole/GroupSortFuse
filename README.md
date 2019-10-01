# GroupSortFuse
Implementation of the Group-Sort-Fuse (GSF) procedure for estimating the number of components in finite mixture models. Four families of mixture models are currently implemented:
 + `normalLocOrder`: Multidimensional Gaussian mixture models in location, with common but possibly unknown scale parameter;
 + `multinomialOrder`: Multinomial mixture models;
 + `poissonOrder`: Univariate Poisson mixture models;
 + `exponentialOrder`: Exponential distribution mixture models.

This package continues to be under development, and has only been tested on Ubuntu 16.04 with R  3.4.4. This release should not be considered stable. 

# Installation
This package may be installed as follows, using the `devtools` R package:
```r
library(devtools)
devtools::install_github("tmanole/GroupSortFuse")
```

# Examples
We provide a usage example for Gaussian mixtures in location, with unknown but common covariance matrix, based on the Old Faithful Geyser dataset [1]. 

```r
library(GroupSortFuse)

data(faithful) 
set.seed(1) 
out <- normalLocOrder(faithful, K=10, lambdas=c(0.1, 0.25, 0.5, 0.75, 1.0, 2), penalty="MCP-LLA", a=2, maxPgd=200, maxMem=500, verbose=FALSE) 
```


```r
plot(out, gg=FALSE)
```
Visuals are also available using `ggplot2`. To install this package, run `install.packages("ggplot2")`. 

```{r}
plot(out, gg=TRUE)
```


# References 
[1] Azzalini, A. and Bowman, A. W. (1990). A look at some data on the
Old Faithful geyser. Applied Statistics 39, 357-365.
