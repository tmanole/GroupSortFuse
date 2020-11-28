#' Estimate the Number of Components in a Multivariate Normal Location Mixture Model.
#'
#' Estimate the order of a finite mixture of multivariate normal densities
#' with respect to the mean parameter, whose variance-covariance matrices
#' are common but potentially unknown. 
#'
#' @param y n by D matrix consisting of the data, where n is the sample size 
#'          and D is the dimension.
#' @param K Upper bound on the true number of components. 
#'          If \code{K} is \code{NULL}, at least one of \code{theta} and
#'          \code{pii} must be non-\code{NULL}, and K is inferred from their
#'          number of columns.
#' @param lambdas Vector of tuning parameter values. 
#' @param sigma   D by D matrix, which is the starting value for the common
#'                variance-covariance matrix. If \code{NULL}, \code{arbSigma}
#'                must be \code{TRUE}, and in this case the starting value is
#'                set to be the sample variance-covariance of the data.
#' @param arbSigma Equals \code{TRUE} if the common variance-covariance matrix should
#'                 be estimated, and FALSE if it should stay fixed, and equal to
#'                 \code{sigma}.
#' @param ... Additional control parameters. See the \strong{Details} section.
#'
#' @return An object with S3 classes \code{gsf} and \code{normalLocGsf}, 
#'         consisting of a list with the estimates produced for every tuning
#'         parameter in \code{lambdas}. 
#'
#' @details The following is a list of additional control parameters.
#'
#'  \describe{
#'   \item{\code{mu}}{D by K matrix of starting values where each column is the mean
#'                      vector for one component. If \code{theta=NULL}, the starting 
#'                      values are chosen using the procedure of Benaglia et al. (2009).}
#'   \item{\code{pii}}{Vector of size K whose elements must sum to 1, consisting of 
#'                      the starting values for the mixing proportions.  
#'                      If \code{NULL}, it will be set to a discrete 
#'                      uniform distribution with K support points.}
#'   \item{\code{penalty}}{Choice of penalty, which may be "\code{SCAD}", "\code{MCP}", 
#'                          "\code{SCAD-LLA}", "\code{MCP-LLA}" or "\code{ADAPTIVE-LASSO}". 
#'                          Default is "\code{SCAD}".}
#'   \item{\code{uBound}}{Upper bound on the tuning parameter of the proximal gradient descent algorithm.}
#'   \item{\code{C}}{Tuning parameter for penalizing the mixing proportions.}
#'   \item{\code{a}}{Tuning parameter for the SCAD or MCP penalty. Default is \code{3.7}.}
#'   \item{\code{convMem}}{Convergence criterion for the modified EM algorithm.}
#'   \item{\code{convPgd}}{Convergence criterion for the proximal gradient descent algorithm.}
#'   \item{\code{maxMem}}{Maximum number of iterations of the modified EM algorithm.}
#'   \item{\code{maxPgd}}{Maximum number of iterations of the proximal gradient descent algorithm.}
#'   \item{\code{verbose}}{If \code{TRUE}, print updates while the function is running.}}
#' 
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
#'
#'  Benaglia, T., Chauveau, D., Hunter, D., Young, D. 2009. "mixtools: An R package for analyzing
#'  finite mixture models". Journal of Statistical Software. 32(6): 1-29.
#'
#' @examples
#'  # Example 1: Seeds Data.
#'    data(seeds) 
#'    y <- cbind(seeds[,2], seeds[,6])
#'    set.seed(1)
#'    out <- normalLocOrder(y, K=12, lambdas=seq(0.1, 1.7, 0.3), maxPgd=200, maxMem=500)
#'    tuning <- bicTuning(y, out)
#'    plot(out, gg=FALSE, eta=TRUE, vlines=TRUE, opt=tuning$result$lambda)
#'   
#'  # Example 2: Old Faithful Data.
#'    data(faithful)
#'    set.seed(1)
#'    out <- normalLocOrder(faithful, K=10, 
#'              lambdas=c(0.1, 0.25, 0.5, 0.75, 1.0, 2), penalty="MCP-LLA", a=2, maxPgd=200, maxMem=500)
#' 
#'    # Requires ggplot2.
#'    plot(out, gg=TRUE, eta=FALSE)
normalLocOrder <- function (y, lambdas, K = NULL, sigma = NULL, arbSigma = TRUE, ...) {
  input <- list(...)

  mu  <- input$mu
  pii <- input$pii

  if (is.null(input[["penalty"]])) penalty <- "SCAD"
  else penalty <- input[["penalty"]]

  if (is.null(input[["uBound"]])) uBound <- 0.1
  else uBound <- input[["uBound"]]
  
  if (is.null(input[["C"]])) C <- 3
  else C <- input[["C"]]

  if (is.null(input[["a"]])) a <- 3.7
  else a <- input[["a"]]

  if (is.null(input[["convPgd"]])) delta <- 1e-5
  else delta <- input[["convPgd"]]

  if (is.null(input[["convMem"]])) epsilon <- 1e-8
  else epsilon <- input[["convMem"]]

  if (is.null(input[["maxMem"]])) maxMem <- 2500
  else maxMem <- input[["maxMem"]]

  if (is.null(input[["maxPgd"]])) maxPgd <- 1500
  else maxPgd <- input[["maxPgd"]]

  if (is.null(input[["verbose"]])) verbose <- T
  else verbose <- input[["verbose"]]

  if (is.data.frame(y)) {      
    y <- as.matrix(y)
  }

  .validateParams(y, mu, "mu", pii, K)
  .validateOther(maxMem, maxPgd, lambdas, penalty, a, C, epsilon, delta)
  .validateNormalLoc(y, mu, sigma, arbSigma)

  lambdas <- sort(lambdas)

  if (!is.null(mu) && !is.null(pii)) {
    if (arbSigma) sigma <- cov(y)

    out <- .myEm(y, mu, sigma, pii, arbSigma, -1, 1, C, a,
                .penCode(penalty), lambdas, epsilon, delta, maxMem, maxPgd, uBound, verbose, 0)

  } else if (!is.null(mu) && is.null(pii)) {
    K <- ncol(mu)    
    if (arbSigma) sigma <- cov(y)

    out <- .myEm(y, mu, sigma, rep(1.0/K, K), arbSigma, -1, 1, C, a,
                .penCode(penalty), lambdas, epsilon, delta, maxMem, maxPgd, uBound, verbose, 0)

  } else {
    if (is.null(pii)) {
      pii <- rep(1.0/K, K)

    } else {
      K <- length(pii)
    }

    if (arbSigma) sigma <- cov(y)

    # Compute starting values.
    means <- apply(y, 1, mean)
    yOrdered <- y[order(means), ]
    yBinned  <- list()

    n <- nrow(y)

    for (j in 1:K) {
      yBinned[[j]] <- yOrdered[max(1, floor((j - 1) * n/K)):ceiling(j * n/K), ]
    }
 
    hypTheta <- lapply(1:K, function(i) apply(yBinned[[i]], 2, mean))
    
    theta <- matrix(NA, ncol(y), K)
    for(k in 1:K) {
      theta[,k] <- as.vector(mixtools::rmvnorm(1, mu = as.vector(hypTheta[[k]]), sigma = sigma))
    }    

    out <- .myEm(y, theta, sigma, pii, arbSigma, -1, 1, C, a,
                .penCode(penalty), lambdas, epsilon, delta, maxMem, maxPgd, uBound, verbose, 0)
  }

  class(out) <- c("normalLocGsf", "gsf")     
  out
}

#' Estimate the Number of Components in a Multinomial Mixture Model.
#'
#' Estimate the order of a finite mixture of multinomial models with fixed and 
#' known number of trials.
#'
#' @param y n by D matrix consisting of the data, where n is the sample size 
#'          and D is the number of categories.
#'          The rows of \code{y} must sum to a constant value, 
#'          taken to be the number of trials. 
#' @param K Upper bound on the true number of components. 
#'          If \code{K} is \code{NULL}, at least one of the control parameters \code{theta} and
#'          \code{pii} must be non-\code{NULL}, and K is inferred from their
#'          number of columns.
#' @inheritParams normalLocOrder
#'
#' @return An object with S3 classes \code{gsf} and \code{multinomialGsf}, 
#'         consisting of a list with the estimates produced for every tuning
#'         parameter in \code{lambdas}. 
#'
#' @details The following is a list of additional control parameters.
#'
#' \itemize{
#'   \item{\code{theta}}    {D by K matrix of starting values where each column is the vector
#'                      of multinomial probabilities for one mixture component.              
#'                      The columns of \code{theta} should therefore 
#'                      sum to 1. If \code{theta=NULL}, the starting values are
#'                      chosen using the MCMC algorithm described by Grenier (2016).} 
#'   \item{\code{pii}}      {Vector of size K whose elements must sum to 1, consisting of 
#'                      the starting values for the mixing proportions.  
#'                      If \code{NULL}, it will be set to a discrete 
#'                      uniform distribution with K support points.}
#'   \item{\code{penalty}} {Choice of penalty, which may be \code{"SCAD"}, \code{"MCP"}, 
#'                          \code{"SCAD-LLA"}, \code{"MCP-LLA"} or \code{"ADAPTIVE-LASSO"}. 
#'                          Default is \code{"SCAD"}.}
#'   \item{\code{mcmcIter}} {Number of iterations for the starting value algorithm described 
#'                           by Grenier (2016).}
#'   \item{\code{uBound}}   {Upper bound on the tuning parameter of the proximal gradient algorithm.}
#'   \item{\code{C}}       {Tuning parameter for penalizing the mixing proportions.}
#'   \item{\code{a}}        {Tuning parameter for the SCAD or MCP penalty. Default is \code{3.7}.}
#'   \item{\code{convMem}}  {Convergence criterion for the modified EM algorithm.}
#'   \item{\code{convPgd}}    {Convergence criterion for the proximal gradient descent algorithm.}
#'   \item{\code{maxMem}}   {Maximum number of iterations of the Modified EM algorithm.}
#'   \item{\code{maxPgd}}   {Maximum number of iterations of the proximal gradient descent algorithm.}
#'   \item{\code{verbose}}  {If \code{TRUE}, print updates while the function is running.}
#' }
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
#'
#'  Grenier, I. (2016) Bayesian Model Selection for Deep Exponential Families. 
#'  M.Sc. dissertation, McGill University Libraries.
#'
#' @examples 
#'  require(MM)
#'  data(pollen)
#'  set.seed(1)
#'  out <- multinomialOrder(pollen, K=12, lambdas=seq(0.1, 1.6, 0.2))
#'  tuning <- bicTuning(pollen, out)
#'  plot(out, eta=TRUE, gg=FALSE, opt=tuning$result$lambda)
multinomialOrder <- function(y, lambdas, K = NULL, ...) {
  input <- list(...)

  theta <- input$theta
  pii <- input$pii

  if (is.null(input[["penalty"]])) penalty <- "SCAD"
  else penalty <- input[["penalty"]]

  if (is.null(input[["uBound"]])) uBound <- 0.1
  else uBound <- input[["uBound"]]
  
  if (is.null(input[["C"]])) C <- 3
  else C <- input[["C"]]

  if (is.null(input[["a"]])) a <- 3.7
  else a <- input[["a"]]

  if (is.null(input[["mcmcIter"]])) mcmcIter <- 50
  else mcmcIter <- input[["mcmcIter"]]

  if (is.null(input[["convPgd"]])) delta <- 1e-5
  else delta <- input[["convPgd"]]

  if (is.null(input[["convMem"]])) epsilon <- 1e-8
  else epsilon <- input[["convMem"]]

  if (is.null(input[["maxMem"]])) maxMem <- 2500
  else maxMem <- input[["maxMem"]]

  if (is.null(input[["maxPgd"]])) maxPgd <- 1500
  else maxPgd <- input[["maxPgd"]]

  if (is.null(input[["verbose"]])) verbose <- T
  else verbose <- input[["verbose"]]

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }

  .validateParams(y, theta, "theta", pii, K)
  .validateOther(maxMem, maxPgd, lambdas, penalty, a, C, epsilon, delta)
  .validateMultinomial(y, theta, mcmcIter)

  lambdas <- sort(lambdas)

  D <- ncol(y) - 1

  if(is.null(K)) {
    out <- .completeMultinomialCols(.myEm(y[, 1:D], theta[1:D, ], diag(D), pii, FALSE, sum(y[1,]), 3, C, a,
                .penCode(penalty), lambdas, epsilon, delta, maxMem, maxPgd, uBound, verbose, 0))
  } else {
    n <- nrow(y)

    newPii <- array(1 / K, dim = c(n, K))
    newZ    <- t(apply(newPii, 1, function(x) rmultinom(1,1,x)))
    piiSamp <- matrix(0, mcmcIter, K)

    for (d in 1:mcmcIter) {
      nkWords  <- apply(newZ, 2, function(x) colSums(x*y))
      newTheta <- apply(nkWords, 2, function(x) gtools::rdirichlet(1,1 + x))

      pxGivTheta <- 0
      for (k in 1:K) {
        pxGivTheta <- pxGivTheta + t(apply(y, 1, function(x) dmultinom(x, sum(x), newTheta[, k])))
      }

      for (k in 1:K) {
        newPii[, k] <- t(apply(y, 1, function(x) dmultinom(x, sum(x), newTheta[, k]))) / sum(pxGivTheta)
      }

      newZ <- t(apply(newPii, 1, function(x) rmultinom(1, 1, x)))
      piiSamp[d, ] <- colSums(newZ)/n
    }

    if (is.null(pii)) {
      newPii <- as.matrix(piiSamp[mcmcIter,])
  
    } else {
      newPii <- as.matrix(pii)
    }

    out <- .completeMultinomialCols(.myEm(y[,1:D], newTheta[1:D, ], diag(D), newPii, FALSE, sum(y[1,]), 3, C, a, 
                .penCode(penalty), lambdas, epsilon, delta, maxMem, maxPgd, uBound, verbose, 0))
  }
   
  class(out) <- c("multinomialGsf", "gsf")
  out
}

#' Estimate the Number of Components in a Poisson Finite Mixture Model.
#'
#' Estimate the order of a finite mixture of Poisson distributions.
#'
#' @param y Vector consisting of the data. 
#' @param K Upper bound on the true number of components. 
#'          If \code{K} is \code{NULL}, at least one of \code{theta} and
#'          \code{pii} must be non-\code{NULL}, and K is inferred from their
#'          number of columns.
#' @inheritParams normalLocOrder
#'
#' @return An object with S3 classes \code{gsf} and \code{poissonGsf}, 
#'         consisting of a list with the estimates produced for every tuning
#'         parameter in \code{lambdas}. 
#'
#' @details  The following is a list of additional control parameters.
#'
#' \itemize{
#'   \item{\code{theta}}  {Vector of starting values for the Poisson parameters, of length K.
#'                         If \code{NULL},  the starting values are chosen to be a discrete 
#'                         uniform distribution on the 100(k- 1/2)/K \% sample quantiles.}
#'   \item{\code{pii}}    {Vector of size K whose elements must sum to 1, consisting of 
#'                         the starting values for the mixing proportions.  
#'                         If \code{NULL}, it will be set to a discrete 
#'                         uniform distribution with K support points.}
#'   \item{\code{penalty}} {Choice of penalty, which may be \code{"SCAD"}, \code{"MCP"}, 
#'                          \code{"SCAD-LLA"}, \code{"MCP-LLA"} or \code{"ADAPTIVE-LASSO"}. 
#'                          Default is \code{"SCAD"}.}
#'   \item{\code{a}}        {Tuning parameter for the SCAD or MCP penalty. Default is \code{3.7}.}
#'   \item{\code{uBound}}   {Upper bound on the tuning parameter of the proximal gradient descent algorithm.}
#'   \item{\code{C}}       {Tuning parameter for penalizing the mixing proportions.}
#'   \item{\code{convMem}}  {Convergence criterion for the modified EM algorithm.}
#'   \item{\code{convPgd}}    {Convergence criterion for the proximal gradient descent algorithm.}
#'   \item{\code{maxMem}}   {Maximum number of iterations of the modified EM algorithm.}
#'   \item{\code{maxPgd}}   {Maximum number of iterations of the proximal gradient descent algorithm.}
#'   \item{\code{verbose}}  {If \code{TRUE}, print updates while the function is running.}
#' }
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
#'
#' @examples 
#'  data(notices)
#'
#'  # Run the GSF with the Adaptive Lasso penalty.
#'  set.seed(1)
#'  out <- poissonOrder(notices, lambdas=c(0, .001, .01, .5, 1, 2), 
#'                      K=12, penalty="ADAPTIVE-LASSO", maxMem=1000)
#'  # Requires ggplot2.
#'  plot(out, eta=FALSE)
#' 
#'  # Run the GSF with the SCAD penalty.
#'  set.seed(1)
#'  out <- poissonOrder(notices, lambdas=c(.00005, .001, .005, .01, .1, .5, 1, 2, 5), 
#'                      K=12, uBound=0.01, penalty="SCAD", maxMem=1000)
#'  plot(out, eta=FALSE)
#'
#'  # Select a tuning parameter using the BIC. 
#'  bicTuning(notices, out)
poissonOrder <- function (y, lambdas, K = NULL, ...) {
  input <- list(...)

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }

  theta <- input$theta
  pii   <- input$pii

  if (is.null(input[["penalty"]])) penalty <- "SCAD"
  else penalty <- input[["penalty"]]

  if (is.null(input[["uBound"]])) uBound <- 0.1
  else uBound <- input[["uBound"]]
  
  if (is.null(input[["C"]])) C <- 3
  else C <- input[["C"]]

  if (is.null(input[["a"]])) a <- 3.7
  else a <- input[["a"]]

  if (is.null(input[["convPgd"]])) delta <- 1e-5
  else delta <- input[["convPgd"]]

  if (is.null(input[["convMem"]])) epsilon <- 1e-8
  else epsilon <- input[["convMem"]]

  if (is.null(input[["maxMem"]])) maxMem <- 2500
  else maxMem <- input[["maxMem"]]

  if (is.null(input[["maxPgd"]])) maxPgd <- 1500
  else maxPgd <- input[["maxPgd"]]

  if (is.null(input[["verbose"]])) verbose <- T
  else verbose <- input[["verbose"]]

  if (is.data.frame(y)) {      
    y <- as.matrix(y)
  }
  
  y <- as.matrix(y, ncol=1)

  .validateParams(y, theta, "theta", pii, K)
  .validateOther(maxMem, maxPgd, lambdas, penalty, a, C, epsilon, delta)

  lambdas <- sort(lambdas)

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }

  if (!is.null(theta) && min(theta) < 0) 
    stop("Error: The elements of 'theta' must be positive.")

  if (is.null(theta)) {
    if (!is.null(K)) {
      theta <- .initialTheta(y, K)
    
    } else {
      theta <- .initialTheta(y, length(pii))
    }
  }

  if (is.null(pii) && !is.null(K))
    pii <- rep(1.0/K, K)

  if (is.null(pii) && is.null(K)) 
    pii <- rep(1.0/length(theta), length(theta))

  theta[which(theta==0)] <- 0.1

  out <- .myEm(matrix(y, ncol=1), matrix(theta, nrow=1), diag(1), pii, FALSE, -1, 5, C, a,
              .penCode(penalty), lambdas, epsilon, delta, maxMem, maxPgd, uBound, verbose, 0)
  
  class(out) <- c("poissonGsf", "gsf")
  out 
}

#' Estimate the Number of Components in a Finite Mixture of Exponential Distributions.
#'
#' Estimate the order of a finite mixture of Exponential distributions.
#'
#' @param y Vector consisting of the data. 
#' @param K Upper bound on the true number of components. 
#'          If \code{K} is \code{NULL}, at least one of \code{theta} and
#'          \code{pii} must be non-\code{NULL}, and K is inferred from their
#'          number of columns.
#' @param theta  Vector of starting values for the Exponential distribution parameters, of length K.
#' @inheritParams normalLocOrder
#'
#' @return An object with S3 classes \code{gsf} and \code{exponentialGsf}, 
#'         consisting of a list with the estimates produced for every tuning
#'         parameter in \code{lambdas}. 
#'
#' @details The following is a list of additional control parameters.
#'
#' \itemize{
#'   \item{\code{pii}}    {Vector of size K whose elements must sum to 1, consisting of 
#'                         the starting values for the mixing proportions.  
#'                         If \code{NULL}, it will be set to a discrete 
#'                         uniform distribution with K support points.}
#'   \item{\code{arbSigma}} {\code{TRUE} if the common variance-covariance matrix should
#'                      be estimated, and FALSE if it should stay fixed, and equal to
#'                      \code{sigma}.}
#'   \item{\code{penalty}} {Choice of penalty, which may be \code{"SCAD"}, \code{"MCP"} or 
#'                          \code{"ADAPTIVE-LASSO"}. Default is \code{"SCAD"}.}
#'   \item{\code{a}}        {Tuning parameter for the SCAD or MCP penalty. Default is \code{3.7}.}
#'   \item{\code{mcmcIter}} {Number of iterations for the starting value algorithm described 
#'                      in the details below). This value is ignored when \code{theta} 
#'                      is not \code{NULL}. }
#'   \item{\code{uBound}}   {Upper bound on the tuning parameter of the PGD algorithm.}
#'   \item{\code{C}}       {Tuning parameter for penalizing the mixing proportions.}
#'   \item{\code{convMem}}  {Convergence criterion for the modified EM algorithm.}
#'   \item{\code{convPgd}}    {Convergence criterion for the proximal gradient descent algorithm.}
#'   \item{\code{maxMem}}   {Maximum number of iterations of the modified EM algorithm.}
#'   \item{\code{maxPgd}}   {Maximum number of iterations of the proximal gradient descent algorithm.}
#'   \item{\code{verbose}}  {If \code{TRUE}, print updates while the function is running.}
#' }
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
exponentialOrder <- function (y, lambdas, K = NULL, theta, ...) {
  input <- list(...)

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }

  theta <- input$theta
  pii   <- input$pii

  if (is.null(input[["penalty"]])) penalty <- "SCAD"
  else penalty <- input[["penalty"]]

  if (is.null(input[["uBound"]])) uBound <- 0.1
  else uBound <- input[["uBound"]]
  
  if (is.null(input[["C"]])) C <- 3
  else C <- input[["C"]]

  if (is.null(input[["a"]])) a <- 3.7
  else a <- input[["a"]]

  if (is.null(input[["convPgd"]])) delta <- 1e-5
  else delta <- input[["convPgd"]]

  if (is.null(input[["convMem"]])) epsilon <- 1e-8
  else epsilon <- input[["convMem"]]

  if (is.null(input[["maxMem"]])) maxMem <- 2500
  else maxMem <- input[["maxMem"]]

  if (is.null(input[["maxPgd"]])) maxPgd <- 1500
  else maxPgd <- input[["maxPgd"]]

  if (is.null(input[["verbose"]])) verbose <- T
  else verbose <- input[["verbose"]]

  if (is.data.frame(y)) {      
    y <- as.matrix(y)
  }
  
  y <- as.matrix(y, ncol=1)

  .validateParams(y, theta, "theta", pii, K)
  .validateOther(maxMem, maxPgd, lambdas, penalty, a, C, epsilon, delta)

  lambdas <- sort(lambdas)

  if (!is.null(theta) && min(theta) < 0) 
    stop("Error: The elements of 'theta' must be positive.")

  if (is.null(pii) && !is.null(K))
    pii <- rep(1.0/K, K)

  if (is.null(pii) && is.null(K)) 
    pii <- rep(1.0/length(theta), length(theta))

  if (!is.null(theta)) {
    out <- .myEm(as.matrix(y, ncol=1), matrix(theta, nrow=1), diag(1), 
                pii, FALSE, -1, 6, C, a, .penCode(penalty), lambdas, 
                epsilon, delta, maxMem, maxPgd, uBound, verbose, 0)

  # TODO: Add starting values.
  } else {
     stop("Starting values for exponential mixtures are not yet supported.")
#    out <- .myEm(as.matrix(y, ncol=1), matrix(initialTheta(y, K), nrow=1), 
#                diag(1), as.matrix(pii), FALSE, -1, 6, ck, a, .penCode(penalty), 
#                lambdas, epsilon, delta, maxMem, maxPgd, uBound, verbose, 0)
  }
  
  class(out) <- c("exponentialGsf", "gsf")
  out 
}

#' Tuning Parameter Selection via the Bayesian Information Criterion 
#'
#' Bayesian Information Criterion (BIC) for tuning parameter selection. 
#' @param y n by D matrix consisting of the data.
#' @param result A \code{gsf} object.
#'
#' @details The BIC selects the best tuning parameter out of the ones used in a
#'  \code{gsf} object by minimizing the following function
#'
#' \deqn{\textrm{BIC}(\lambda) = -2 l_n(\hat{\mathbf{\Psi}}_{\lambda}) + 
#' \textrm{df}(\hat{\mathbf{\Psi}}_{\lambda}) \log n}
#'
#' where \eqn{l_n} is the log-likelihood function, and \eqn{
#' \hat{\mathbf{\Psi}}_{\lambda}} is the set of estimated 
#' parameters \code{theta} and \code{pii} corresponding to the tuning parameter \eqn{\lambda}.
#' \eqn{\textrm{df}(\hat{\mathbf{\Psi}}_{\lambda})} denotes the degrees
#' of freedom of the estimates.
#'
#' @return A list containing the selected tuning parameter and 
#'         corresponding estimates, as well as a summary of the
#'         computed RBIC values, log-likelihood values, and corresponding
#'         orders of the estimates.
#'
#' @examples 
#'  require(MM)
#'  data(pollen)
#'  set.seed(1)
#'  out <- multinomialOrder(pollen, K=12, lambdas=seq(0.1, 1.6, 0.2))
#'  tuning <- bicTuning(pollen, out)
#'  plot(out, eta=TRUE, gg=FALSE, opt=tuning$result$lambda)
bicTuning <- function(y, result) {
  if (is.vector(y)) {
    n <- length(y)

  } else {
    n <- nrow(y)
  }

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }
 
  if (is.vector(y)) {
    y <- matrix(y, ncol=1)
  }
 
  l <- length(result)
  thetas <- list()

  for (i in 1:l) {
    thetas[[i]] <- switch(class(result)[1], 
                  "normalLocGsf" = result[[i]]$mu,
                  "multinomialGsf" = result[[i]]$theta,
                  "poissonGsf" = result[[i]]$theta,
                  "exponentialGsf" = result[[i]]$theta,
                   stop("Error: Class not recognized."))
  }  

  if (is.vector(thetas[[1]])) {
    D <- length(thetas[[1]])

  } else {
    D <- nrow(thetas[[1]])
  }

  out      <- list()
  estimate <- list() 

  lambdaVals <- c()
  order      <- c()
  df         <- c()
  loglik     <- c()
  bic        <- c()  

  minRbic  <- 9999999999999999 #.Machine$double.max
  minIndex <- -1 
  
  arbSigma <- ("sigma" %in% names(result[[1]])) && (class(result)[1] == "normalLocGsf")

  index <- switch(class(result)[1], 
            "normalLocGsf" = 1,
            "multinomialGsf" = 3,
            "poissonGsf" = 5,
            "exponentialGsf" = 6,
            stop("Error: Class not recognized."))

  for (i in 1:l) {
    lambdaVals[i] <- result[[i]]$lambda

    if(lambdaVals[i] == 0) next
  
    order[i]      <- result[[i]]$order
    df[i] <- switch(class(result)[1], 
            "normalLocGsf" = {
              if(arbSigma) {
                order[i] * (D+1) + 0.5 * D * (D+1) - 1
              
              } else {
                order[i] * (D+1) - 1
              }
            }, 
            "multinomialGsf" = D * order[i] - 1,
            "poissonGsf" = 2 * order[i] - 1,
            "exponentialGsf" = 2 * order[i] - 1,
            stop("Error: Class not recognized."))
  
    if (arbSigma) {
      loglik[i] <- .bicLogLik(y, as.matrix(thetas[[i]], nrow=D), result[[i]]$pii, result[[i]]$sigma, index)

    } else {
      loglik[i] <- .bicLogLik(y, as.matrix(thetas[[i]], nrow=D), result[[i]]$pii, diag(D), index)
    }

    if (!is.finite(loglik[i])) {
      stop("Error: Log-likelihood is infinite.")
    }

    bic[i] <- -2 * loglik[i] + df[i] * log(n)

    if (bic[i] < minRbic) {
      minRbic <- bic[i]
      minIndex  <- i
    } 
  }

  if (index == 1) {
    estimate[["mu"]] <- result[[minIndex]]$mu
    estimate[["sigma"]] <- result[[minIndex]]$sigma

  } else {
    estimate[["theta"]] <- result[[minIndex]]$theta
  }

  estimate[["pii"]]   <- result[[minIndex]]$pii
  estimate[["order"]] <- result[[minIndex]]$order
  estimate[["lambda"]]<- result[[minIndex]]$lambda 

  out[["summary"]] = data.frame(lambdaVals, order, bic, loglik, df)
  out[["result"]] = estimate

  out 
}
 
#' Plot the coefficient paths of an \code{gsf} object. 
#'
#' @param x An \code{gsf} object.
#' @param gg A ggplot2 figure is produced if \code{TRUE}, and a Base R figure is produced otherwise.
#' @param eta If \code{TRUE}, a coefficient plot for the distance between consecutive coefficient is produced, and a 
#'            plot with all coefficients is produced otherwise.
#' @param vlines If \code{TRUE}, vertical dashed lines are plotted whenever the number of components decreases, 
#'               \code{FALSE} otherwise.
#' @param points If \code{TRUE}, points are plotted along the lines, and are not plotted otherwise.
#' @param opt Optimal tuning parameter value. If not \code{NULL}, a vertical red line is drawn at this point.
#' @param sqrtScale \code{TRUE} if lambda values are to be plotted on the square root scale, \code{FALSE} otherwise.
#' @param ... Additional Base R plotting parameters.
#'
#' @return Returns one or more plots depending on the value of \code{eta}.
#'
#' @examples 
#'  data(faithful)
#'  out <- normalLocOrder(faithful, K=12, 
#'                        lambdas=c(0.1, 0.25, 0.5, 0.75, 1.0), penalty="MCP", a=2)
#'  plot(out)
#'  bicTuning(faithful, out)
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
plot.gsf <- function (x, gg=FALSE, eta=TRUE, vlines=TRUE, points=FALSE, opt=NULL, sqrtScale=FALSE, ...) {
  result <- x
  l <- length(result)
  thetas <- list()

  for (i in 1:l) {
    thetas[[i]] <- switch(class(result)[1],
                  "normalLocGsf" = result[[i]]$mu,
                  "multinomialGsf" = result[[i]]$theta,
                  "poissonGsf" = result[[i]]$theta,
                  "exponentialGsf" = result[[i]]$theta,
                   stop("Error: Class not recognized."))
  }

  if (is.vector(thetas[[1]])) {
    K <- length(thetas[[1]])
    D <- 1

  } else {
    K <- ncol(thetas[[1]])
    D <- nrow(thetas[[1]])
  }

  lambdas <- c()

  for (j in 1:length(result)) {
    lambdas <- c(lambdas, result[[j]]$lambda)
  }

  if (!is.null(opt) && sqrtScale) {
    opt <- sqrt(opt)
  }

  # Base R version.
  if (!gg) {
    L <- length(lambdas)

    mat <- matrix(NA, L, K-1)
    for (i in 1:L) {
      th <- as.matrix(thetas[[i]][, order(thetas[[i]][1,])]) 

      if (ncol(th) == 1) th <- t(th)

      for (k in 1:(K-1)) {
        mat[i, k] <- sqrt(sum((th[, k+1] - th[, k])^2 ))
      }

      mat[i,] <- sort(mat[i,])
    }
  
    if (sqrtScale) { 
      lambdas <- sqrt(lambdas)
    }

    # Plot for eta variables.
    if (eta) {
      plot(lambdas, mat[,1], type="l", ylim=c(min(mat), max(mat)), col="navyblue", 
           xlab=if(sqrtScale) expression(sqrt(lambda)) else expression(lambda), ylab=expression(paste("||", hat(eta)[j], "||")), ...)

      if (points==TRUE) {
        points(lambdas, mat[,j], pch=19, col="navyblue")
      }
      
      for (j in 1:(K-1)) {
        points(lambdas, mat[,j], type="l", col="navyblue")
        
        if (points==TRUE) {
          points(lambdas, mat[,j], pch=19, col="navyblue")
        }
      }

      if (vlines) {
        lines <- c()
        
        for (i in 2:L) {
          if (nrow(unique(t(thetas[[i]]))) != nrow(unique(t(thetas[[i-1]])))) {
            abline(v=lambdas[i], lty="dashed")            
          }
        }
      }

      if (!is.null(opt)) {
        abline(v=opt, col="red", lwd=2)
      }
    
    # Plot for theta variables.
    } else {
      if (round(sqrt(D)) != floor(sqrt(D))) {
        parRows = ceiling(sqrt(D))
        parCols = round(sqrt(D))

      } else {
        parRows = round(sqrt(D))
        parCols = ceiling(sqrt(D))
      }

      par(mfrow=c(parRows, parCols), pty="s")

      for(j in 1:D) {
        out <- matrix(NA, l, K)

        for (i in 1:l) {
          for (k in 1:K) {
            out[i, k] <- thetas[[i]][j, k]
          }

           out[i,] <- sort(out[i,])
        }

        ch <- if (points) 'b' else 'l'

        matplot(out,
                 type = ch,
                 pch = 20,
                 xlab = if(sqrtScale) expression(sqrt(lambda)) else expression(lambda),
                 ylab = expression(theta[lambda]),
                 xaxt = "n",
                 yaxt = "n", ...
               )
        title(paste("Coordinate ", j))
        axis(side = 1, at = 1:length(result), labels = round(lambdas, 3), las=1)
        axis(2, las=1)

  
        if (vlines) {
          lines <- c()
          
          for (i in 2:L) {
            if (nrow(unique(t(thetas[[i]]))) != nrow(unique(t(thetas[[i-1]])))) {
              abline(v=i, lty="dashed")            
            }
          }
        }

        if (!is.null(opt)) {
          abline(v=which.min(abs(opt - lambdas))[1], col="red", lwd=2)
        } 
      }
    } 

  # GGplot version.
  } else {

    if (!requireNamespace ("ggplot2", quietly = TRUE)) {
      stop ("Error: Please load or install the ggplot package.", call. = FALSE)
    }

    if (sqrtScale) {
      lambdas <- sqrt(lambdas)
    }

    # Plot for eta variables.
    if (eta) {
      L <- length(lambdas)

      mat <- matrix(NA, L, K-1)
      for (i in 1:L) {
        th <- thetas[[i]][, order(thetas[[i]][1,])]

        for (k in 1:(K-1)) {
          mat[i, k] <- sqrt(sum((th[, k+1] - th[, k])^2 ))
        }

        mat[i,] <- sort(mat[i,])
      }

      out <- matrix(NA, L*(K-1), 3)

      for (k in 1:(K-1)) {
        for (i in 1:L) {
          out[(k-1)*l + i, 1] <- lambdas[i]
          out[(k-1)*l + i, 2] <- mat[i,k]
          out[(k-1)*l + i, 3] <- k
        }
      }

      df <- data.frame(out)

      df$X3 <- as.factor(df$X3)

      p <- ggplot2::ggplot(df, ggplot2::aes(X1, X2, group=X3))+#, linetype=X3)) +
        ggplot2::geom_line(color="navyblue") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
              panel.grid.major.y = ggplot2::element_blank(),
              panel.grid.minor.x = ggplot2::element_blank(),
              panel.grid.minor.y = ggplot2::element_blank(),
              aspect.ratio = 1,
              legend.position="none",
              strip.background=ggplot2::element_blank(),
              panel.spacing=ggplot2::unit(1, "centimetres")) +
        ggplot2::ylab(expression(paste("||", hat(eta), "||")))+
        ggplot2::xlab(if(sqrtScale) expression(sqrt(lambda)) else expression(lambda))

      # Add points.
      if (points) {
        p <- p + ggplot2::geom_point(color="navyblue")
      }

      # Add vertical line for optimal tuning parameter value.
      if (!is.null(opt)) {
        p <- p + ggplot2::geom_vline(xintercept=opt, color="red", size=1.5)
      }

      # Add vertical lines.
      if (vlines) {
        lines <- c()
        for (i in 2:L) {
          if (sum(mat[i, ] == 0) != sum(mat[i-1, ] == 0)) {
            lines <- c(lines, lambdas[i])
          }
        }

        p <- p + ggplot2::geom_vline(xintercept=lines, linetype=2)
      }

      return(p)

    # Plot for theta variables.
    } else {

      L <- length(lambdas)

      out <- matrix(NA, L*K*D, 4)

      labels <- c()

      for(j in 1:D) {

        labels <- c(labels, paste("Coordinate ", j))

        for (i in 1:L) {
          th <- thetas[[i]][, order(thetas[[i]][1,])]

          for (k in 1:K) {
            out[(j-1)*K*L + (k-1)*l + i, 1] <- lambdas[i]
            out[(j-1)*K*L + (k-1)*l + i, 2] <- th[j, k]
            out[(j-1)*K*L + (k-1)*l + i, 3] <- k
            out[(j-1)*K*L + (k-1)*l + i, 4] <- j #paste("Coordinate ", j)
          }
        }
      }

      df <- data.frame(out)

      df$X3 <- as.factor(df$X3)
      df$X4 <- as.factor(df$X4)

      levels(df$X4) <- labels

      # The following color palette is inspired by RColorBrewer. 
      # Neuwirth, Erich. Package 'RColorBrewer'. CRAN 2011-06-17 08: 34: 00. Apache License 2.0, 2011.
      colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
                  "#FF7F00", "#CAB2D6", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
                  "#E6AB02", "#A6761D", "#666666")

      p <- ggplot2::ggplot(df, ggplot2::aes(X1, X2, group=X3))+#, linetype=X3)) +
        ggplot2::geom_line(ggplot2::aes(color=X3)) +
        ggplot2::scale_color_manual(values = colors)+
        ggplot2::facet_wrap(~X4, ncol=2, scales="free_y") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                    panel.grid.major.y = ggplot2::element_blank(),
                    panel.grid.minor.x = ggplot2::element_blank(),
                    panel.grid.minor.y = ggplot2::element_blank(),
                    aspect.ratio = 1,
                    legend.position="none",
                    strip.background=ggplot2::element_blank(),
                    panel.spacing=ggplot2::unit(1, "centimetres")) +
        ggplot2::xlab(if(sqrtScale) expression(sqrt(lambda)) else expression(lambda)) +
        ggplot2::ylab(expression(hat(theta)))+
        ggplot2::annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
        ggplot2::annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)


      # Draw optimal vertical line.
      if (!is.null(opt)) {
        p <- p + ggplot2::geom_vline(xintercept=opt, color="red", size=1.5)
      }

      # Add vertical lines
      if (vlines) {
        lines <- c()
        for (i in 2:L) {
          if (nrow(unique(t(thetas[[i]]))) != nrow(unique(t(thetas[[i-1]])))) {
            lines <- c(lines, lambdas[i])
          }
        }

        p <- p + ggplot2::geom_vline(xintercept=lines, linetype=2)
      }

      # Draw points. 
      if (points) {
        p <- p + ggplot2::geom_point(ggplot2::aes(color=X3))
      }

      return(p)
    }
  }
}


