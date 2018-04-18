# Utility Functions
logist <- function(eta) return(exp(eta) / (1 + exp(eta)))



newtonRaphson <- function(y, X, precision = precision, iterMax = iterMax) {
  beta <- rep(0, times = ncol(X))
  i <- 1
  convergence <- FALSE

  while (i <= iterMax & convergence == FALSE) {

    # calculate probabilities (mu)
    eta <- X %*% beta
    mu <- as.numeric(logist(X %*% beta))

    # init / update Matrix W
    W <- diag(mu * (1 - mu))

    # calculate and update beta (eq. 23 from Czepiel, 2002)
    deltaBeta <- solve(t(X) %*% W %*% X) %*% t(X) %*% (y - mu)
    beta <- beta + deltaBeta

    # check convergence condition
    if (sum(abs(deltaBeta)) < precision) {
      convergence <- TRUE
      break
    }

    # update i
    i <- i + 1
  }

  if (convergence == FALSE) {
    stop(paste(i, "iterations reached without convergence. Increase iterMax?"))
  }


  return(list(beta = beta,
              mu = mu,
              W = W,
              eta = eta,
              iterCount = i))
}
