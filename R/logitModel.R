#' Maximum Likelihood Estimation
#'
#' This function calculates the maximum likelihood for binary logistic regression
#' based on the Newton Raphson Method.
#'
#' This function is the calculation behind the \code{\link{logitMod}}
#' formula interface to the user.
#'
#' @param y Matrix/vector containing the response variable.
#' @param X Matrix/vector containing the design matrix (including intercept).
#' @param precision Degree of precision of algorithm.
#' @param iterMax Maximum number of iterations of optimisation algorithm.
#'
#' @return A list with maximum likelihood estimation results.
#'
#' @references Agresti, A (2013) Categorical Data Analysis,
#' John Wiley & Sons, 2013, Volume 792 of Wiley Series in Probability
#' and Statistics, ISSN 1940-6517
#'
#' @references Czepiel, S.A. (2002) Maximum Likelihood Estimation of Logistic
#' Regression Models: Theory and Implementation. Available at
#' \href{https://czep.net/stat/mlelr.pdf}{czep.net/stat/mlelr.pdf}. [visited:
#' 25.03.2018]
#'
#' @examples
#' set.seed(42)
#' X1 <- rnorm(100);    X2 <- rnorm(100);
#' points <- 2 * X1 - 3 * X2
#' y <- rbinom(100, 1, exp(points) / (1 + exp(points)))
#' X <- cbind(b1 = 1, b2 = X1, b3 = X2)
#' fit <- MLE(y = y, X = X)
#'
#' @export
MLE <- function(y, X, precision = 1e-10, iterMax = 25) {

  # Call Newton Raphson Method
  NRMRes <- newtonRaphson(y, X, precision, iterMax)

  # Degrees of Freedom
  dfRes <- nrow(X) - ncol(X)  # DF Residual Model
  dfNull <- nrow(X) - 1       # DF Null Model

  # Compute covariace matrix (eq 22 from Czepiel, 2002)
  W <- NRMRes$W
  vcov <- solve(t(X) %*% W %*% X)

  # Compute residual deviance TODO: learn this thing here
  eta <- NRMRes$eta
  s <- ifelse(y == 0, -1, y)
  Residuals <- as.numeric(s * sqrt(-2*((y * eta) - (log(1 + exp(eta))))))

  # calculate max value of logLikelihood
  beta <- NRMRes$beta
  maxLL <- (sum((y * X %*% beta) - (log(1 + exp(X %*% beta)))))

  # Return list of results
  list(coefficients = beta[,],
       fittedValues = NRMRes$mu,
       vcov = vcov,
       Residuals = Residuals,
       dfRes = dfRes,
       dfNull = dfNull,
       maxLL = maxLL,
       W = W,
       iterCount = NRMRes$iterCount)
}



#' Alternative Logistic Regression Implementation
#'
#' \code{logitMod} computes a logistic regression for a binary response variable
#'    via \code{\link{MLE}} estimation method.
#'
#' @param formula an object of class "formula".
#' @param data an optional data frame or coercible to one. If empty, the
#' variables are taken from the Global Environment.
#'
#' @param precision degree of precision of optimisation algorithm.
#' @param iterMax maximum number of iterations of optimisation algorithm.
#'
#' @return \code{logitMod} returns an object of \code{\link[base]{class}}
#'     \emph{logitMod}.
#'
#' @examples
#' set.seed(42)
#' X1 <- rnorm(100);    X2 <- rnorm(100);
#' points <- 2 * X1 - 3 * X2
#' y <- rbinom(100, 1, exp(points) / (1 + exp(points)))
#' X <- cbind(b1 = 1, b2 = X1, b3 = X2)
#' fit <- logitMod(y ~ X1 + X2)
#' @export
logitMod <- function(formula, data, precision = 1e-10, iterMax = 25) {

  ##############################################################################
  #  initialisation                                                            #
  ##############################################################################

  # generate model.frame. If data unused, search variables in Parent Env
  if(missing(data)) {
    modelFrame <- model.frame(formula, data = parent.frame())
  } else {
    modelFrame <- model.frame(formula, as.data.frame(data))
  }

  # Extract design matrix (X) and response var (y)
  X <- model.matrix(formula, modelFrame)
  y <- model.response(modelFrame)

  # make sure response variable is coded with 0s and 1s
  if (!(0 %in% y && 1 %in% y)) {
    y <- factor(y, labels = c(0, 1))
  }
  y <- as.numeric(as.character(y))

  # sanity check
  if (length(unique(y)) != 2) {
    stop("Response variable is expected to be binary")
  }

  ##############################################################################
  #  computation                                                               #
  ##############################################################################


  # conduct ML-Estimation for restricted Model
  restricModel <- MLE(y = y,
                      X = as.matrix(rep(1, length(y))),
                      precision = precision,
                      iterMax = iterMax)

  # conduct ML-Estimation for Full Model
  result <- MLE(y, X, precision = precision, iterMax = iterMax)

  idx <- match(c("coefficients",
                 "Residuals",
                 "dfRes",
                 "dfNull",
                 "maxLL",
                 "iterCount"),
               names(restricModel))


  # Results
  result$restricModel <- restricModel[idx]
  result$devianceNull <- -2 * restricModel$maxLL
  result$Residuals <- -2 * result$maxLL
  result$call <- match.call()
  result$formula <- formula
  result$X <- X
  result$y <- y

  # Assign Class
  class(result) <- "logitMod"

  result
}

#' S3 Print Method for logitMod Objects
#'
#' This internal function defines a \code{\link{print}} method for an object of
#' logitMod class.
#'
#' @param model a logitMod object
#' @param ... methods are required to include (at least) the same arguments
#' as their generic functions.
#'
#' @return a brief summary of the model results.
#'
#' @export
print.logitMod <- function(model, ...){

  # based on stats:::print.lm()
  cat("Call: ", paste0(deparse(model$call)), fill = TRUE)
  cat("\n\nCoefficients:\n")

  print.default(format(coef(model), digits = 4L),
                print.gap = 1L, quote = FALSE, right = TRUE)

  cat("\nDegrees of Freedom: ", model$dfNull, " Total (i.e. Null); ",
      model$dfRes, " Residual")

  # Akaike Information Criterium
  AIC <- (-2 * model$maxLL + 2 * ncol(model$X))

  cat("\nNull Deviance:\t", round(model$devianceNull, 2))

  cat("\nResidual Deviance:", round(model$devianceResid, 2),
      "\t","AIC: ", round(AIC,1), "\n")

  # invisibly return linMod object
  invisible(model)
}


#' S3 Summary Method for logitMod Objects
#'
#' This internal function defines a \code{\link{summary}} method for an object
#' of logitMod class, where also further statistics are computed for model
#' comparison and for the print.summary method.
#'
#' @param model a logitMod model
#' @param ... methods are required to include (at least) the same arguments
#' as their generic functions.
#'
#' @return a list with more detailed statistics and computations on the
#' estimated model
#'
#' @export
summary.logitMod <- function(model, ...) {

  # Standard Error
  betaSE <- sqrt(diag(model$vcov))
  model$betaSE <- betaSE

  # z-values
  zStat <- model$coefficients / betaSE
  model$zStat <- zStat

  # p-values
  pValue <- 2 * pnorm(-abs(zStat))
  model$pValue <- pValue

  # Berechnung von devianceNull, devianceResid & aic #################
  model$devianceNull <- sum(model$restricModel$devianceResid^2)
  model$devianceResid <- sum(model$devianceResid^2)
  model$AIC <- (-2 * model$maxLL + 2 * ncol(model$X))

  class(model) <- "summary.logitMod"

  return(model)
}


#' S3 Print Method for summary.logitMod Objects
#'
#' This internal function defines a \code{\link{print}} method for an object
#' of \code{\link{summary.logitMod}} class.
#'
#' @param model an object of class "summary.logitMod".
#' @param ... methods are required to include (at least) the same arguments
#' as their generic functions.
#'
#' @return a summary of the model results.
#'
#' @export
print.summary.logitMod <- function(model, ...) {

  cat("Call: ", deparse(model$call), fill = TRUE)

  cat("\nDeviance Residuals:\n")
  print.summaryDefault(summary(model$devianceResid), digits = 4L)

  cat("\nCoefficients:\n")

  model$coefTable <- cbind("Estimate" = model$coefficients,
                       "Std. error" = model$betaSE,
                       "z value" = model$zStat,
                       "Pr(>|z|)" = model$pValue)

  printCoefmat(model$coefTable, signif.legend = FALSE, digits = 4L)

  cat("\nNull deviance: ",
      round(model$devianceNull,2), " on ",
      model$dfNull, " degrees of freedom\n")

  cat("Residual deviance: ",
      round(model$devianceResid,2), " on ",
      model$dfRes, " degrees of freedom\n")

  cat("AIC: ", round(model$AIC,2))

  cat("\n\nIterations: ", model$iterCount, "\n")

  # invisibly return summary
  invisible(model)
}


#' S3 Plot Method for logitMod Objects
#'
#' This internal function defines a \code{\link{plot}} method for an object
#' of logitMod class.
#'
#' @param model an object of class logitMod
#' @param ... methods are required to include (at least) the same arguments
#' as their generic functions.
#'
#' @export
plot.logitMod <- function(model, which = c(1,2,3), ...) {

  if (1 %in% which) {

    plot(y = model$Residuals, model = (model$X %*% model$coefficients),
         main = "Residuals vs Fitted",
         ylab = "Residuals",
         xlab = paste("Predicted Values\n",
                      deparse(model$call)))
    abline(a = 0, b = 0, lty = 3)

  }
  if (2 %in% which) {

    qqnorm(model$Residuals,
           main = "Normal Q-Q",
           ylab = "Std. deviance resid.",
           xlab = paste("Theoretical Quantiles\n", deparse(model$call)))
    qqline(model$Residuals, lty = 3)
  }
  if (3 %in% which) {

    plot(y = sqrt(abs(model$Residuals)), model = (model$X %*% model$coefficients),
         main = "Scale Location",
         ylab = expression(sqrt("|Std. deviance resid.|")),
         xlab = paste("Predicted Values\n", deparse(model$call)))
  }
}
