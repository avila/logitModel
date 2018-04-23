#' Alternative Logistic Regression Implementation
#'
#' \code{logitMod} computes a logistic regression for a binary response variable
#'    via \code{\link{newtonRaphson}} estimation method.
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

  ############################################################################
  # Initialisation                                                           #
  ############################################################################

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

  ############################################################################
  # Computation                                                              #
  ############################################################################

  # Restricted Model
  XNull <- as.matrix(rep(1, length(y)))
  restricModel <- newtonRaphson(y = y,
                                X = XNull,
                                precision = precision,
                                iterMax = iterMax)

  bNull <- restricModel$beta
  maxLLNull <- (sum((y * XNull %*% bNull) - (log(1 + exp(XNull %*% bNull)))))
  devianceNull <- -2 * maxLLNull
  dfNull <- nrow(X) - 1       # DF Null Model


  # Estimation Model
  MLEstimation <- newtonRaphson(y, X, precision = precision, iterMax = iterMax)

  # Degrees of Freedom
  dfRes <- nrow(X) - ncol(X)  # DF Residual Model

  # Compute covariace matrix (eq 22 from Czepiel, 2002)
  W <- MLEstimation$W
  vcov <- solve(t(X) %*% W %*% X)

  # Compute residual deviance TODO: learn this thing here
  eta <- MLEstimation$eta
  s <- ifelse(y == 0, -1, y)
  residuals <- as.numeric(s * sqrt(-2*((y * eta) - (log(1 + exp(eta))))))

  # calculate max value of logLikelihood
  beta <- MLEstimation$beta
  maxLL <- (sum((y * X %*% beta) - (log(1 + exp(X %*% beta)))))

  # Information Criteria
  AIC <- -2 * maxLL + 2 * ncol(X)

  # calculate Residual Deviance
  devianceResid <- -2 * maxLL

  # Return list of results
  result <- list(coefficients = beta[,],
                 fittedValues = MLEstimation$pi,
                 vcov = vcov,
                 residuals = residuals,
                 dfRes = dfRes,
                 maxLL = maxLL,
                 devianceResid = devianceResid,
                 AIC = AIC,

                 # Restricted Model
                 dfNull = dfNull,
                 devianceNull = devianceNull,

                 # Extras
                 formula = formula,
                 call = match.call(),
                 y = y,
                 X = X,
                 W = W,
                 iterCount = MLEstimation$iterCount)

  # Assign Class
  class(result) <- "logitMod"

  result
}


#' Newton-Raphson for Logistic Regression Estimation
#'
#' This function calculates the maximum likelihood for binary logistic
#' regression via the Newton Raphson algorithm.
#'
#' This algorithm is used on the \code{\link{logitMod}} regression estimation
#' function.
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
#' fit <- newtonRaphson(y = y, X = X)
#'
#' @export
newtonRaphson <- function(y, X, precision = 1e-10, iterMax = 25) {

  beta <- rep(0, times = ncol(X))
  i <- 1
  convergence <- FALSE

  while (i <= iterMax & convergence == FALSE) {

    # calculate probabilities (pi)
    eta <- X %*% beta
    pi <- as.numeric(logist(X %*% beta))

    # init / update Matrix W
    W <- diag(pi * (1 - pi))

    # calculate and update beta (eq. 23 from Czepiel, 2002)
    deltaBeta <- solve(t(X) %*% W %*% X) %*% t(X) %*% (y - pi)
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
              pi = pi,
              W = W,
              eta = eta,
              iterCount = i))
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

  # Extra Statistics
  model$devianceNull <- model$devianceNull
  model$devianceResid <- model$devianceResid
  model$AIC <- (-2 * model$maxLL + 2 * ncol(model$X))
  model$BIC <- -2 * model$maxLL + ncol(model$X) * log(length(model$y))

  class(model) <- "summary.logitMod"

  model
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

  cat("\nCoefficients:\n")

  print.default(format(coef(model),digits = 4L),
                print.gap = 2L, quote = FALSE, right = TRUE)

  # create named vector for printing two columns of summary stats. detail: cols
  # must be of even length for sprintf to properly iterate over the columns
  cols <- c("Obs."          = length(model$y),
            "DF"           = length(model$y) - ncol(model$X),
            "Deviance"     = round(model$devianceResid, 2),
            "AIC"          = model$AIC
  )

  # the following lines is just to make sure the output fits each
  # in its columns, no matter how many characters the variables have
  i2 <- max(nchar(format(cols[seq(1, length(cols), 2)], scientific=FALSE))) + 1
  i4 <- max(nchar(format(cols[seq(2, length(cols), 2)], scientific=FALSE))) + 1

  # concatenate over "four columns", as in: (names_l: 12.123 \t names_r: 45.678)
  cat(sprintf(paste0("\n%-9s:%",i2,".2f\t%-4s:%",i4,".2f"),
              names(cols[seq(1, length(cols), 2)]),    # first (names) col
              cols[seq(1, length(cols), 2)],           # second (number) col
              names(cols[seq(2, length(cols), 2)]),    # third (names) col
              cols[seq(2, length(cols), 2)])           # fourth (numbers) col
  )

  invisible(model)
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
  print.summaryDefault(summary(model$residuals), digits = 5L)

  cat("\nCoefficients:\n")

  model$coefTable <- cbind("Estimate" = model$coefficients,
                       "Std. error" = model$betaSE,
                       "z value" = model$zStat,
                       "Pr(>|z|)" = model$pValue)

  printCoefmat(model$coefTable, signif.legend = FALSE, digits = 5L,
               print.gap = 2)

  # create named vector for printing two columns of summary stats. detail: cols
  # must be of even length for sprintf to properly iterate over the columns
  cols <- c("Obs."            = length(model$y),
            "DF"              = length(model$y) - ncol(model$X),
            "Resid. Deviance" = round(model$devianceResid, 2),
            "AIC"             = model$AIC,
            "Null Devinance"  = model$devianceNull,
            "BIC"             = model$BIC,
            "Log Likelih."    = model$maxLL,
            "Iterations"      = model$iterCount
  )


  # the following lines is just to make sure the output fits each
  # in its columns, no matter how many characters each variable have
  i2 <- max(nchar(format(cols[seq(1, length(cols), 2)], scientific=FALSE))) + 1
  i4 <- max(nchar(format(cols[seq(2, length(cols), 2)], scientific=FALSE))) + 1

  # concatenate over "four columns", as in: (names_l: 12.123 \t names_r: 45.678)
  cat(sprintf(paste0("\n%-16s:%",i2,".2f\t%-11s:%",i4,".2f"),
              names(cols[seq(1, length(cols), 2)]),    # first (names) col
              cols[seq(1, length(cols), 2)],           # second (number) col
              names(cols[seq(2, length(cols), 2)]),    # third (names) col
              cols[seq(2, length(cols), 2)])           # fourth (numbers) col
  )


  invisible(model)
}


#' S3 Plot Method for logitMod Objects
#'
#' This internal function defines a \code{\link{plot}} method for an object
#' of logitMod class.
#'
#' @param model an object of class logitMod.
#' @param which numeric vector indicating which plot to be plotted.
#' 1: Residuals vs Fitted. 2: Normal Q-Q. 3: Scale Location. 4: Residual vs
#' Leverage
#' @param cookLevels vector indicating levels of Cook's Distance to be drawn in
#' plot 4 (Residual vs Leverage).
#' @param col optional string or numeric vector to indicate colours of data
#' points for response variable. If missing, assigned values for 0
#' and 1 will be c("#483D8BD0", "#B22222D0"), blue and red respectively.
#' @param ask optional boolean wether inteded to plot interactively or not.
#' Default value checks if graphics device are in interactive mode and
#' wheter the user previously defined a matrix plotting layout in order to
#' make a sensible guess about asking for each plot or not.
#' @param ... methods are required to include (at least) the same arguments
#' as their generic functions.
#'
#' @export
plot.logitMod <- function(model, which = 1:4, col, cookLevels = c(.5, 1), ask, ...) {

  # Common Parameters
  oldPar <- par(no.readonly = TRUE)
  on.exit(par(oldPar))
  par(mar = c(2, 0, 0, 0) + 2.2)
  if (missing(ask)) ask <- prod(par("mfcol")) < length(which) && dev.interactive()
  if (missing(col)) col <- c("#483D8BD0", "#B22222D0")
  else if (length(col) < 2) col <- rep(col, 2)

  # PLOT 01
  if (1 %in% which) {
    x <- model$X %*% model$coefficients
    y <- model$residuals
    plot(x = x, y = y,
         main = "Residuals vs Fitted",
         ylab = "Residuals",
         xlab = paste("Predicted Values\n", deparse(model$call)),
         col = col[1 + model$y], ...)
    abline(a = 0, b = 0, lty = 3, col = "gray")

    # add label of highest absolute residuals
    idx <- order(abs(y), decreasing = T)[1:5]
    text(x[idx], y[idx], model$residuals[idx], labels = idx, cex = .66, pos = 2,
         col = col[model$y[idx] + 1])

    # add smooth curve via LOWESS regression
    lines(lowess(x = x, y = y), col = "red")
  }

  # Plot 02
  if (2 %in% which) {
    if (length(which) > 1 && ask && nextPlot() == "no") {
      return(invisible())
    }
    # save plot as p for extracting (x,y) points for plotting labels
    p <- qqnorm(model$residuals,
                main = "Normal Q-Q",
                ylab = "Std. deviance resid.",
                xlab = paste("Theoretical Quantiles\n", deparse(model$call)),
                col = col[1 + model$y], ...)
    qqline(model$residuals, lty = 3, col = "gray")

    # add label of highest absolute residuals
    idx <- order(abs(model$residuals), decreasing = TRUE)[1:5]
    text(p$x[idx], p$y[idx],
         model$residuals[idx], labels = idx, cex = .66, pos = 2,
         col = col[model$y[idx] + 1])
  }

  # Plot 03
  if (3 %in% which) {
    if (length(which) > 1 && ask && nextPlot() == "no") {
      return(invisible())
    }

    x <- model$X %*% model$coefficients
    y <- sqrt(abs(model$residuals))
    ylab <- as.expression(substitute(sqrt(abs(x)),
                                     list(x = as.name("Std. Deviance Resid."))))
    plot(x = x, y = y,
         main = "Scale Location",
         ylab = ylab,
         xlab = paste("Predicted Values\n", deparse(model$call)),
         col = col[1 + model$y], ...)

    # Plot Smoothing via LOWESS regression
    lines(lowess(x = model$X %*% model$coefficients,
                  y = sqrt(abs(model$residuals))), col = "red")

    # add label of highest absolute residuals
    idx <- order(abs(y), decreasing = T)[1:5]
    text(x[idx], y[idx],
         model$residuals[idx], labels = idx, cex = .66, pos = 2,
         col = col[model$y[idx] + 1])

  }

  # Plot 04
  if (4 %in% which) {
    if (length(which) > 1 && ask && nextPlot() == "no") {
      return(invisible())
    }

    leverage <- diag(sqrt(model$W) %*% model$X %*%
                       (solve(t(model$X) %*%model$W %*% model$X)) %*%
                       t(model$X) %*% sqrt(model$W))
    pearsonResidual <- (model$y - model$fittedValues) /
                        sqrt(model$fittedValues*(1 - model$fittedValues))

    plot(x = leverage, y = pearsonResidual,
         main = "Residual vs Leverage",
         ylab = "Std. Pearson resid.",
         xlab = paste("Leverage\n", deparse(model$call)),
         col = col[1 + model$y],
         ylim = range(pearsonResidual) * 1.1, ...)

    abline(a = 0, b = 0, lty = 3, col = "gray")

    # Plot Smoothing via LOWESS regression
    lines(lowess(x = leverage,
                 y = pearsonResidual), col = "red")

    # PLOT COOKS DISTANCE
    rangeLev <- range(leverage)
    modelRank <- ncol(model$X)
    plotLim <- par("usr")
    seqPlot <- seq.int(min(rangeLev[1L], rangeLev[2L]/100),
                       plotLim[2L], length.out = 101)

    for (crit in cookLevels) {
      cl.h <- sqrt(crit * modelRank * (1 - seqPlot)/seqPlot)
      lines(seqPlot, cl.h, lty = 2, col = 2)
      lines(seqPlot, -cl.h, lty = 2, col = 2)
    }


    legend("bottomleft", legend = paste0("Cook's Distance [",
                                        paste(cookLevels, collapse = ", "),"]"),
           lty = 2, col = 2, bty = "n", cex = 0.75, text.col = "#000000AA")


    # add label of highest absolute residuals
    idx <- order(abs(model$residuals), decreasing = T)[1:5]
    text(leverage[idx], pearsonResidual[idx],
         model$residuals[idx], labels = idx, cex = .66, pos = 2,
         col = col[model$y[idx] + 1])
  }
}


#' S3 Pairs Method for logitMod Objects
#'
#' This internal function defines a \code{\link{pairs}} method for an object
#' of logitMod class. The main idea is to plot all explanatory variables
#' agaist the response variable of the model in order to quickly assess the
#' influence of each one.
#'
#' @param model an object of class logitMod.
#' @param nRow optional numeric vector indicating the number of rows the overall
#' plot should have.
#' @param single boolean to indicate whether each variable should be plotted
#' against the response variable individually.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dummies <- sample(LETTERS[1:5], n, replace = TRUE)
#' X1 <- rnorm(n); X2 <- rnorm(n)
#' points <- 2 * X1 - 3 * X2
#' y <- rbinom(n, 1, exp(points) / (1 + exp(points)))
#' df <- as.data.frame(cbind(b2 = X1, b3 = X2))
#' myFitD <- logitMod(y ~ X1 + X2 + dummies)
#' pairs(myFitD)
#' pairs(myFitD, single = TRUE)
#' @export
pairs.logitMod <- function(model, nRow, single = FALSE, col, ...) {

  # INIT
  if (missing(col)) col <- c("#483D8BD0", "#B22222D0")
  intercept <- attr(terms(model$formula), "intercept")
  if (intercept == 1) {
    X <- as.data.frame(model$X)[-1]
    betas <- model$coefficients[-1]
  } else {
    X <- as.data.frame(model$X)
    betas <- model$coefficients
  }

  # GRAPHICS PARAMETERS
  # keep previous parameters
  oldPar <- par(no.readonly = TRUE)
  on.exit(par(oldPar))

  if (missing(nRow)) nRow <- ceiling(sqrt(length(model$coefficients)))
  parCol <- ceiling(length(X) / nRow)
  parRow <- min(nRow, length(X))

  if (isTRUE(!single)) {
    par(mfrow = c(parRow, parCol), mar = c(0, 0, 0, -0.5) + 1, cex = 0.9,
        cex.axis = 0.9, tcl = -0.15, oma = c(0, 0, 0, 0), mgp = c(2, 0, 0))
  }

  # plot each
  for (i in seq_along(X)) {
    plot(X[, i], jitter(model$y, 1/10),
         xlab = "", ylab = "", main = names(X)[i],
         col = col[model$y + 1], ...)

    # col=rgb(0, 0, 0, 0.4))
    curve(logist(intercept + (betas[i] * x)),
          lwd=2, lty=3, col = "gray", add = TRUE)

    # points
    points(X[, i], logist(intercept + betas[i] * X[, i]),
           col = col[model$y + 1])

    if (isTRUE(single)) {
      prompt <- paste0("Plot ", i, "/", ncol(X),
                       ": Enter anything to continue or [q] to quit: ")
      UInput <- readline(prompt=prompt)
      if (UInput %in% c("q", "Q")) return(invisible())
    }
  }
}
