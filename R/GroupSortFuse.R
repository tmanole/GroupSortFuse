#' The GroupSortFuse Package.
#'
#' Implementation of the Group-Sort-Fuse (GSF) method (Manole and Khalili 2019), 
#' which estimates the number of components in finite mixture models via penalized 
#' maximum likelihood estimation. The method is implemented for the following classes of mixture models,
#' \itemize{
#'   \item{Location-Gaussian mixture models, with unknown or known covariance (\code{normalLocOrder})}
#'   \item{Multinomial mixture models (\code{multinomialOrder})}
#'   \item{Poisson mixture models (\code{poissonOrder})}
#'   \item{Exponential distribution mixture models (\code{exponentialOrder})}}
#' This package also provides coefficient path plotting functionality (\code{plot.gsf}).
#' Tuning parameter selection can be performed using the Bayesian Information
#' Criterion (\code{bicTuning}) or V-fold Cross Validation (\code{cvTuning}). 
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
#' 
#' @name GroupSortFuse-package
NULL
