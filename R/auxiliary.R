.penCode <- function(pen) {
  if (pen == "SCAD") 1
  else if (pen == "MCP") 2
  else if (pen == "ADAPTIVE-LASSO") 3
  else if (pen == "SCAD-LLA") 4
  else if (pen == "MCP-LLA") 5
  else 0
}

.initialTheta <- function(y, K) {
  result <- c()

  for(k in 1:K){
    result[k] = quantile(y[k], (k-0.5)/K)
  }

  result
}

.frequency <- function(theta) {
  remList <- c()
 
  for (i in 1:ncol(theta)){
    if (i == ncol(theta)) next
      
    for (j in (i+1):ncol(theta)){
      if (dist(rbind(theta[,i], theta[,j])) < 1e-8){
        remList <- c(remList, j)
      }
    }
  }

  ncol(theta)-length(unique(remList))
}

.completeMultinomialCols <- function (res) {
  D <- nrow(res[[1]]$theta)
  K <- ncol(res[[1]]$theta)

  for (i in 1:length(res)) {
    temp <- c()
    for(j in 1:K) {
      temp[j] <- 1 - sum(res[[i]]$theta[1:D,j])
    }

    res[[i]]$theta <- rbind(res[[i]]$theta, temp)
    
  }

  res
}

.validateParams <- function (y, theta, str, pii, K) {
  if (!is.null(K) && (K%%1 != 0 || K < 1)) 
    stop("Error: 'K' must be a strictly positive integer.")

  if (!is.null(pii) && sum(pii) != 1)
    stop ("Error: The sum of the elements of pii must be equal to 1.")

  if (is.null(K) && is.null(theta) && is.null(pii)) 
    stop (paste("Error: At least one of 'K', '", str, "' and 'pii' must be non-NULL.", sep=""))

  if (!is.null(K) && !is.null(pii) && K != length(pii)) 
    stop("Error: 'pii' must be of length 'K'.")

  if (!is.null(K) && !is.null(theta)) {
    if ((is.vector(theta) && length(theta) != K) ||
        (!is.vector(theta) && ncol(theta) != K))
      stop(paste("Error: '", str, "' must be of length 'K'.", sep=""))
  }

  if (!is.null(theta) && !is.null(pii)) {
    if ((is.vector(theta) && length(theta) != length(pii)) ||
        (!is.vector(theta) && ncol(theta)   != length(pii)))
      stop(paste("Error: 'pii' must be of same length as the columns of '", str, "'.", sep=""))
  }

  if (nrow(y) <= ncol(y)) 
    stop("Error: Number of columns of 'y' must be less than its number of rows.")

  if (any(is.na(y)) || any(is.infinite(y))) 
    stop("Error: 'y' cannot contain missing or infinite values.")

  if (!is.null(pii) && (any(is.na(pii)) || any(is.infinite(pii))))
    stop("Error: 'pii' cannot contain missing or infinite values.")

  if (!is.null(pii) && (min(pii) < 0 || max(pii) > 1)) 
    stop("Error: 'pii' must only contain values between 0 and 1.")
}


.validateMultiParams <- function (y, mu, sigma, str1, str2, pii, K) {
  if (!is.null(pii) && sum(pii) != 1)
    stop ("Error: The sum of the elements of pii must be equal to 1.")

  if (is.null(K) && is.null(mu) && is.null(sigma) && is.null(pii)) 
    stop (paste("Error: At least one of 'K', '", str1, "', '", str2, "' and 'pii' must be non-NULL."))

  if (!is.null(K) && !is.null(pii) && K != length(pii)) 
    stop("Error: 'pii' must be of length K.")

  if (!is.null(K) && !is.null(mu) && K != length(mu)) 
    stop(paste("Error: '", str1, "' must be of length K."))
  
  if (!is.null(K) && !is.null(sigma) && K != length(sigma)) 
    stop(paste("Error: '", str2, "' must be of length K."))
  
  if (!is.null(mu) && !is.null(pii) && length(pii) != length(mu)) 
    stop(paste("Error: 'pii' must be of same length as '", str1, "'."))

  if (!is.null(mu) && !is.null(pii) && length(pii) != length(mu)) 
    stop(paste("Error: 'pii' must be of same length as '", str2, "'."))

  if (!is.vector(y) && nrow(y) <= ncol(y)) 
    stop("Error: Number of columns of 'y' must be less than its number of rows.")

  if (any(is.na(y)) || any(is.infinite(y))) 
    stop("Error: 'y' cannot contain missing or infinite values.")
}

.validateOther <- function (maxMem, maxItd, lambdaList, penalty, a, ck, epsilon, delta) {
  if (maxMem != floor(maxMem) || maxItd != floor(maxItd) || maxMem < 1 || maxItd < 1)
    stop ("Error: 'maxMem' and maxItd must be positive integers.")

  if (length(lambdaList) < 1 || min(lambdaList) < 0) 
    stop ("Error: 'lambdaList' must be non-empty and may only contain non-negative numbers.")

  if (!identical(penalty, "SCAD") && !identical(penalty, "MCP") && !identical(penalty, "ADAPTIVE-LASSO") &&
      !identical(penalty, "SCAD-LLA") && !identical(penalty, "MCP-LLA")) 
    stop ("Error: 'penalty' must be equal to 'SCAD', 'MCP', 'ADAPTIVE-LASSO', 'SCAD-LLA' or 'MCP-LLA'.")

  if ((identical(penalty, "SCAD") || identical(penalty, "SCAD-LLA")) && a <= 2) 
    stop ("Error: 'a' must be greater than 2 for the SCAD penalty.")

  if ((identical(penalty, "MCP") || identical(penalty, "MCP-LLA")) && a <= 1) 
    stop ("Error: 'a' must be greater than 1 for the MCP penalty.")

  if (ck < 0) 
    stop("Error: 'ck' must be nonnegative.")

  if (epsilon <= 0) 
    stop("Error: 'epsilon' must be positive.")

  if (delta <= 0) 
    stop("Error: 'delta' must be nonnegative.")
}

.validateMultinomial <- function (y, theta, mcmcIter) {
  if (length(unique(rowSums(y))) != 1) 
    stop ("Error: The sum of each row of 'y' must be equal to a constant.")

  if (!is.null(theta) && ((length(unique(colSums(theta))) != 1 || unique(colSums(theta)) != 1)))
    stop ("Error: The sum of each column of 'theta' must equal 1.")

  if (!(all(y == floor(y)) && all(y >= 0)))
    stop ("Error: Every element of 'y' must be a nonnegative integer.")

  if (!is.null(theta) && (mcmcIter != floor(mcmcIter) || mcmcIter < 1)) 
    stop ("Error: 'mcmcIter' must be an integer greater or equal to 1.")
}

.validateNormalLoc <- function (y, mu, sigma, arbSigma) {
  if (!is.null(sigma) && (!isSymmetric(sigma) || any(eigen(sigma, symmetric = T)$values <= 0))) 
    stop ("Error: 'sigma' must be a positive definite symmetric matrix.")

  if (is.null(sigma) && !arbSigma) 
    stop("Error: 'sigma' cannot be NULL when 'arbSigma' is FALSE.")

  if (!is.null(sigma) && (nrow(sigma) != ncol(y))) 
    stop("Error: 'sigma' must be a DxD dimensional matrix, where D is the number of columns of 'y'.")
  
  if (!is.null(mu) && (nrow(mu) != ncol(y)))
    stop("Error: 'mu' must have as many rows as the columns of 'y'.")

  if (!is.null(sigma) && any(is.na(sigma)) || any(is.infinite(sigma))) 
    stop("Error: 'sigma' cannot contain missing or infinite values.")

  if (!is.null(mu) && (any(is.na(mu)) || any(is.infinite(mu))))
    stop("Error: 'mu' cannot contain missing or infinite values.")
}

