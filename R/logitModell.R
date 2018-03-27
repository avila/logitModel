#' Maximum Likelihood Estimation
#'
#' This function calculates the maximum likelihood for binary logistic regression
#' based on the Newton Raphson Method.
#'
#' This function is the calculation behind the \code{\link{logitMod}}
#' formula interface to the user.
#'
#' @param y Matrix/vector containing the response variable
#' @param X Matrix/vector containing the design matrix (including intercept)
#' @param precision Degree of precision of algorithm
#' @param iterMax Maximum number of iterations of optimisation algorithm
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
#' y <- rbinom(n, 1, exp(points) / (1 + exp(points)))
#' X <- cbind(b1 = 1, b2 = X1, b3 = X2)
#' fit <- MLE(y = y, X = X)
MLE <- function(y, X, precision = 1e-14, iterMax = 100) {

  # Newton Raphson Method

  # follows Czepiel (2002) notation (page 8)
  beta <- rep(0, times = ncol(X))
  i <- 1
  convergence <- FALSE

  while (i <= iterMax & convergence == FALSE) { # epsilon > precision &

    # calculate probabilities (mu)
    eta <- X %*% beta
    mu <- as.numeric(logist(X %*% beta))

    # init / update ... Matrix W
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
    stop(paste(i, "iterations reached without convergence"))
  }

  # Degrees of Freedom #TODO: understand this
  dfRes <- nrow(X) - ncol(X)  # DF Residual Modell
  dfNull <- nrow(X) - 1       # DF Null Modell

  # Compute covariace matrix (eq 22 from Czepiel, 2002)
  vcov <- solve(t(X) %*% W %*% X)

  # Compute residual deviance TODO: learn this thing here
  s <- ifelse(y == 0, -1, y)
  devianceResidual <- as.numeric(s * sqrt(-2*((y * eta) - (log(1 + exp(eta))))))

  # calculate max value of logLikelihood
  maxLL <- (sum((y * X %*% beta) - (log(1 + exp(X %*% beta)))))

  # Liste der zurückgegebenen Werte
  result <- list(coefficients = beta[,],
                 ## residuals??
                 fittedValues = mu,
                 vcov = vcov,
                 devianceResidual = devianceResidual,
                 dfRes = dfRes,
                 dfNull = dfNull,
                 maxLL = maxLL,
                 W = W,

                 iterCount = i)
  return(result)
}

#' Alternative Interface for Logistic Regression Modelling
#'
#' \code{logitMod} computes a logistic regression for a binary response variable
#'    via \code{\link{MLE}} estimation method.
#'
#' @param formula an object of class "formula".
#' @param data a data frame or coercible to one
#' @param precision degree of precision of optimisation algorithm
#' @param iterMax maximum number of iterations of optimisation algorithm
#'
#' @return \code{logitMod} returns an object of \code{\link[base]{class}}
#'     \emph{logitMod}.
#'
#' @examples
#'
#' @export
logitMod <- function(formula, data, precision = 1e-14, iterMax = 100) {

  # generate model.frame, extract design matrix (X) and response var (y)
  modelFrame <- model.frame(formula, as.data.frame(data))
  X <- model.matrix(formula, modelFrame)
  y <- model.response(modelFrame)


  # TODO: sanity checks
  if (FALSE) {
    TRUE
  }

  # conduct ML-Estimation for Full Model
  result <- MLE(y, X, precision = precision, iterMax = iterMax)

  # conduct ML-Estimation for Null Model
  nullModell <- MLE(y = y, X = matrix(rep(1, times = nrow(X)),ncol = 1),
                    precision = precision, iterMax = iterMax)
  result$nullModell <- nullModell

  # Results
  result$formula <- formula
  result$call <- match.call()
  result$X <- X
  result$y <- y

  # Assign Class
  class(result) <- "logitMod"

  result
}


#' Printing method for logitMod estimations
#'
#' Printing method for class "logitMod"
#'
#' @param x an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return a standard print output equivalent to the built in binary logistic regression
#' @examples
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' print(logm)
#' @export
print.logitMod <- function(x, ...){

  cat("Call: ", paste0(deparse(x$call)), fill = TRUE)
  cat("\n\nCoefficients:\n")

  print.default(format(coef(x)[,1], digits = 4L),
                print.gap = 1L, quote = FALSE, right = TRUE)

  cat("\nDegrees of Freedom: ", x$dfNull, " Total (i.e. Null); ",
      x$dfRes, " Residual")

  # Berechnung von null deviance, residual deviance & aic
  # TODO: take it out of here
  nullDeviance <- sum(x$nullModell$devianceResidual^2)
  x$nullDeviance <- nullDeviance
  devianceResidual <- sum(x$devianceResidual^2)
  x$devianceResidual <- devianceResidual
  x_AIC <- (-2*x$maxLL + 2*ncol(x$X))
  x$AIC <- x_AIC

  cat("\nNull Deviance:\t", round(nullDeviance,1))
  cat("\nResidual Deviance:", round(devianceResidual,1), "\t",
      "AIC: ", round(x_AIC,1), "\n")

  # invisibly return linMod object
  invisible(x)

}


#' Summary method for "logitMod" estimations
#'
#' Summary method for class "logitMod"
#'
#' @param object an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return a list of all necessary values equivalent to the summary output of the
#' built in binary logistic regression
#' @examples
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' summary(logm)
summary.logitMod <- function(model, ...) {

  # Koeffizienten Standardfehler
  betaStandardError <- as.matrix(sqrt(diag(model$vcov)))
  model$betaStandardError <- betaStandardError

  # z-Statistik
  zStat <- model$coefficients / betaStandardError
  model$zStat <- zStat

  # p-Werte
  pValue <- 2 * pnorm(-abs(zStat))
  model$pValue <- pValue

  # Zusammenfassung der Werte für die Koeffizienten [,]
  model$coefficients <- cbind("Estimate" = model$coefficients,
                               "Std. error" = model$betaStandardError,
                               "z value" = model$zStat,
                               "Pr(>|z|)" = model$pValue)

  # Berechnung von nullDeviance, residualDeviance & aic
  nullDeviance <- sum(model$nullModell$devianceResidual^2)
  model$nullDeviance <- nullDeviance
  residualDeviance <- sum(model$devianceResidual^2)
  model$residualDeviance <- residualDeviance
  model$AIC <- (-2*model$maxLL + 2*ncol(model$X))

  class(model) <- "summary.logitMod"

  return(model)

}


#' Printing method for the summary of "logitMod" estimations
#'
#' Printing method for summary of class "logitMod"
#'
#' @param x an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return an equivalent output to the summary of the built in binary logistic regression
#' @examples
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' summary(logm)
#' @export
print.summary.logitMod <- function(x, ...) {

  cat("Call: ", deparse(x$call), fill = TRUE)

  cat("\nDeviance Residuals:\n")
  print.summaryDefault(summary(x$devianceResidual)[-4], digits = 4L)

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, signif.legend = TRUE, digits = 4L)

  cat("\n    Null deviance: ",
      round(x$nullDeviance,2), " on ", x$dfNull, " degrees of freedom\n")
  cat("Residual deviance: ",
      round(x$residualDeviance,2), " on ", x$dfRes, " degrees of freedom\n")
  cat("AIC: ", round(x$AIC,2))

  cat("\n\nNumber of Fisher Scoring iterations: ", "\n")

  # invisibly return summary
  invisible(x)
}


#' @title Plotting method for "logitMod" objects
#' @description Plotting method for objects of class "logitMod"
#' @param x an object of class "logitMod"
#' @param ... unused parameter, methods are required to have same
#' arguments as their generic functions
#' @return equivalent plots to those of the built in binary logistic regression
#' @examples
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logm <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' plot(logm)
#' @export
plot.logitMod <- function(x, which = 1, ...) {

  if (which == 1) {
    plot(y = x$devianceResidual, x = (x$X %*% x$coefficients),
         main = "Residuals vs Fitted",
         ylab = "Residuals",
         xlab = paste("Predicted Values\n",
                      deparse(x$call)))
    abline(a = 0, b = 0, lty = 3)
  } else if (which == 2) {
    qqnorm(x$devianceResidual,
           main = "Normal Q-Q",
           ylab = "Std. deviance resid.",
           xlab = paste("Theoretical Quantiles\n", deparse(x$call)))
    qqline(x$devianceResidual, lty = 3)
  } else {
    plot(y = sqrt(abs(x$devianceResidual)), x = (x$X %*% x$coefficients),
         main = "Scale Location",
         ylab = expression(sqrt("|Std. deviance resid.|")),
         xlab = paste("Predicted Values\n", deparse(x$call)))
  }

  # #4 shittt
  # pearsonResidual <-
  #     (x$y - x$fittedValues)/sqrt(x$fittedValues*(1 - x$fittedValues))
  # plot(y = pearsonResidual, x = diag(x$W))
}


## Helper Functions
logist <- function(eta) return(exp(eta) / (1 + exp(eta)))
