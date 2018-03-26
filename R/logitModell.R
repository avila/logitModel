#' @title Maximum Likelihood
#' @description This function calculates the maximum likelihood for binary logistic regression
#' @param y a matrix/vector containing the dependent variable
#' @param X a matrix containing the intercept and all independent variables
#' @return a list with maximum likelihood estimation results
#' @examples
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' MLE(X = model.matrix(testModell, testModelFrame), y = model.response(testModelFrame))
#' @export
MLE <- function(y, X, tolerance = 1e-14, epsilon = 1e-10, iterMax = 1000) {

  #  Newton Raphson Method
  beta <- rep(0, times = ncol(X))
  i <- 0
  convergence <- FALSE
  while (i <= iterMax & convergence == FALSE) { # epsilon > tolerance &

    # p := probabilities
    eta <- X %*% beta
    p <- as.numeric(logist(X %*% beta))

    # aktualisiere ...
    M <- diag(p * (1 - p))

    # calculate and update beta
    deltaBeta <- solve(t(X) %*% M %*% X) %*% t(X) %*% (y - p)
    beta <- beta + deltaBeta

    # update epsilon and i
    epsilon <- sum(abs(deltaBeta))
    i <- i + 1

    # check non-convergence condition
    if (epsilon < tolerance) {
      convergence <- TRUE
      break
    }
  }
  if (convergence == FALSE) {
    stop(paste(i, "iterations reached: Lack of convergence"))
  }
  # DF := Nr. Observations - Nr. Parameters
  dfRes <- nrow(X) - ncol(X)  # DF Residual Modell
  dfNull <- nrow(X) - 1       # DF Null Modell

  # Kovarianzmatrix
  vcov <- solve(t(X) %*% M %*% X)

  # Devianz Residual
  s <- y
  s[s == 0] = -1
  devianceResidual = as.numeric(s * sqrt(-2*((y * eta) - (log(1 + exp(eta))))))

  #https://onlinecourses.science.psu.edu/stat504/node/62

  # Maximumwert der Log Likelihood Funktion
  maxLogLikeValue <- (sum((y * X %*% beta) - (log(1 + exp(X %*% beta)))))

  # Fitted Values
  fittedWerte <- p
  M <- M # TODO: whats is this?

  # Liste der zurückgegebenen Werte
  result <- list(coefficients = beta,
                 vcov = vcov,
                 devianceResidual = devianceResidual,
                 dfRes = dfRes,
                 dfNull = dfNull,
                 maxLogLikeValue = maxLogLikeValue,
                 fittedWerte = fittedWerte,
                 M = M,

                 iterCount = i)
  return(result)

}

## Helper Functions
logist <- function(eta) return(exp(eta) / (1 + exp(eta)))


#' @title Interface for an alternative logistic regression implementation
#' @description This function computes coefficients of a binary logistic regression
#' and contructs an object of "logitMod" class
#' @param formula a formula object. Matrix or array are not accepted
#' @param data a data frame contains all variables
#' @return a list of class "logitMod" containing the coefficients, a null model, formula, call,
#' the dependent variable and all independent variables.
#' @examples
#' testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
#' testData$rank <- factor(testData$rank)
#' testModell <- as.formula("admit ~ gre + gpa + rank")
#' testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
#' logitMod(formula = admit ~ gre + gpa + rank, data = testData)
#' @export
logitMod <- function(formula, data) {

  # generate model.frame, extract design matrix (X) and response var (y)
  modelFrame <- model.frame(formula, data)
  X <- model.matrix(formula, modelFrame)
  y <- model.response(modelFrame)

  # TODO: sanity checks

  # falls y nicht als 0-1/Variable eingegeben wird
  if (!(0 %in% y && 1 %in% y)) {
    y <- factor(y, labels = c(0,1))
  }
  y <- as.numeric(as.character(y))

  # conduct ML-Estimation for Full Model
  result <- MLE(y, X)

  # conduct ML-Estimation for Null Model
  nullModell <- MLE(y = y, X = matrix(rep(1, times = nrow(X)), ncol = 1))
  result$nullModell <- nullModell

  # speichere notwendige Parameter in die Ergebnisliste
  result$formula <- formula
  result$call <- match.call()
  result$X <- X
  result$y <- y

  # ordne die Ergebnisliste der Klasse "logitMod" zu
  class(result) <- "logitMod"

  return(result)

}


#' @title Printing method for logitMod estimations
#' @description Printing method for class "logitMod"
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
  nullDeviance <- sum(x$nullModell$devianceResidual^2)
  x$nullDeviance <- nullDeviance
  devianceResidual <- sum(x$devianceResidual^2)
  x$devianceResidual <- devianceResidual
  x_AIC <- (-2*x$maxLogLikeValue + 2*ncol(x$X))
  x$AIC <- x_AIC

  cat("\nNull Deviance:\t", round(nullDeviance,1))
  cat("\nResidual Deviance:", round(devianceResidual,1), "\t",
      "AIC: ", round(x_AIC,1), "\n")

  # invisibly return linMod object
  invisible(x)

}


#' @title Summary method for "logitMod" estimations
#' @description Summary method for class "logitMod"
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
#' @export
summary.logitMod <- function(object, ...) {

  # Koeffizienten Standardfehler
  betaStandardError <- as.matrix(sqrt(diag(object$vcov)))
  object$betaStandardError <- betaStandardError

  # z-Statistik
  zStat <- object$coefficients / betaStandardError
  object$zStat <- zStat

  # p-Werte
  pValue <- 2 * pnorm(-abs(zStat))
  object$pValue <- pValue

  # Zusammenfassung der Werte für die Koeffizienten
  object$coefficients <- cbind("Estimate" = object$coefficients[,],
                               "Std. error" = object$betaStandardError[,],
                               "z value" = object$zStat[,],
                               "Pr(>|z|)" = object$pValue[,])

  # Berechnung von nullDeviance, residualDeviance & aic
  nullDeviance <- sum(object$nullModell$devianceResidual^2)
  object$nullDeviance <- nullDeviance
  residualDeviance <- sum(object$devianceResidual^2)
  object$residualDeviance <- residualDeviance
  x_AIC <- (-2*object$maxLogLikeValue + 2*ncol(object$X))
  object$AIC <- x_AIC

  class(object) <- "summary.logitMod"

  return(object)

}


#' @title Printing method for the summary of "logitMod" estimations
#' @description Printing method for summary of class "logitMod"
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
plot.logitMod <- function(x, ...) {

  #1
  plot(y = x$devianceResidual, x = (x$X %*% x$coefficients),
       main = "Residuals vs Fitted",
       ylab = "Residuals",
       xlab = paste("Predicted Values\n",
                    deparse(x$call)))
  abline(a = 0, b = 0, lty = 3)

  #2
  qqnorm(x$devianceResidual,
         main = "Normal Q-Q",
         ylab = "Std. deviance resid.",
         xlab = paste("Theoretical Quantiles\n", deparse(x$call)))
  qqline(x$devianceResidual, lty = 3)

  #3
  plot(y = sqrt(abs(x$devianceResidual)), x = (x$X %*% x$coefficients),
       main = "Scale Location",
       ylab = expression(sqrt("|Std. deviance resid.|")),
       xlab = paste("Predicted Values\n", deparse(x$call)))

  # #4 shittt
  # pearsonResidual <-
  #     (x$y - x$fittedWerte)/sqrt(x$fittedWerte*(1 - x$fittedWerte))
  # plot(y = pearsonResidual, x = diag(x$M))

}
