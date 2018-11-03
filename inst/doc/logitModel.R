## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#")

## ----library, include=FALSE----------------------------------------------
library("logitModel")

## ----shapeLogistic, echo=FALSE, fig.height=2.7, fig.cap=cap, fig.align='center'----
cap <- paste("Shape of logistic function for positive and negative values of",
              "$\\beta$ and $\\alpha = 0$.",
             "\\label{fig:logistic}")

op <- par(mfrow = c(1,2), mar = c(1.3, 1, -1, -0.5) + 2, cex = 0.9,
          cex.axis = 0.9, tcl = -0.15, oma = c(0, 0, 0, 0), mgp = c(2, 0.5, 0))

lim <- 5.5
curve(logitModel:::logist(x), from = -lim, to = lim, 
      ylab = "P(Y = 1)", main = "", xlab = "x")
legend(-2.5,.8, expression(beta[i] > 0), bty = "n")
curve(logitModel:::logist(-x), from = -lim, to = lim, 
      ylab = "P(Y = 1)", main = "", xlab = "x")
legend(0,0.4, expression(beta[i] < 0), bty = "n")

## ----NRExample, echo=FALSE, fig.cap=cap2, fig.height=2.5, fig.align='center', fig.width=5.5----
cap2 <- paste("Illustration of the Newton-Rahpson Method.",
             "\\label{fig:nr}")

cap = ""
par(mar = c(0,0,0,0))
curve(x^3, 4, 20, frame.plot = FALSE, xaxt='n', yaxt="n", 
      ylim = c(-300, 7000),xlab = "",ylab = "",main = "")
h<-200
abline(h=h)
text(6, y = 555, labels = expression(theta[real]), cex = 1)
dx2x <- deriv(~ x^3, "x")
x <- 16
pp <- eval(dx2x)
b <- attr(pp, "gradient")
a <- pp - x * attr(pp,"gradient") 
ad <- 1.25
lines(x =c(0, x*ad) , y = c(a, a + x*ad * b))
points(x, y = pp, pch=20)
text(x - 0.5, y = pp + 300, labels = expression("f'"[i]))
text(x+2, a + (x+2) * b + 1000, labels = "f(x)", cex = 0.8)
text(x = -a / b + .45, y = -270, labels = expression(hat(theta)[i]))
x <- (-a + h) / b
pp <- eval(dx2x)
lines(x =c(x, x) , y = c(h, pp), lty = 3)
points(x, h, pch=20)
points(5.8480, h, pch=20)
b <- attr(pp, "gradient")
a <- pp - x * attr(pp,"gradient") 
ad <- 1.22
lines(x =c(0, x*ad) , y = c(a, a + x*ad * b))
points(x, y = pp, pch=20)
text(x-0.5, y = pp + 300, labels = expression("f'"[1+i]))
text(x = -a / b + 0.95, y = -270, labels = expression(hat(theta)[i+1]))
points((-a + h)/b, h, pch=20)

## ------------------------------------------------------------------------
library(logitModel)
fit <- logitModel(survived ~ age + sex, data = DonnerData)

## ---- eval=FALSE---------------------------------------------------------
#  logitModel <- function(formula, data, precision = 1e-10, iterMax = 25) {
#  
#    # generate model.frame. If "data" missing, search variables in Parent Env
#    if(missing(data)) {
#      modelFrame <- model.frame(formula, data = parent.frame())
#    } else {
#      modelFrame <- model.frame(formula, as.data.frame(data))
#    }
#    ...

## ---- eval=FALSE---------------------------------------------------------
#    ...
#    # Extract design matrix (X) and response var (y)
#    X <- model.matrix(formula, modelFrame)
#    y <- model.response(modelFrame)
#  
#    # make sure response variable is coded with 0s and 1s
#    if (!(0 %in% y && 1 %in% y)) {
#      y <- factor(y, labels = c(0, 1))
#    }
#    y <- as.numeric(as.character(y))
#  
#    # sanity check
#    if (length(unique(y)) != 2) {
#      stop("Response variable is expected to be binary")
#    }
#    ...

## ---- eval=FALSE---------------------------------------------------------
#    ...
#    # Restricted Model
#    XNull <- as.matrix(rep(1, length(y)))
#    restricModel <- newtonRaphson(y = y,
#                                  X = XNull,
#                                  precision = precision,
#                                  iterMax = iterMax)
#  
#    bNull <- restricModel$beta
#    logLikelihoodNull <- (sum((y * XNull %*% bNull) -
#                                (log(1 + exp(XNull %*% bNull)))))
#  
#    devianceNull <- -2 * logLikelihoodNull
#    dfNull <- nrow(X) - 1       # DF Null Model
#  
#  
#    # Estimation Model
#    MLEstimation <- newtonRaphson(y, X, precision = precision, iterMax = iterMax)
#    ...

## ------------------------------------------------------------------------
fit

## ------------------------------------------------------------------------
summary(fit)

## ------------------------------------------------------------------------
op <- par(mfrow = c(2,2))
plot(fit)

## ---- eval=FALSE---------------------------------------------------------
#  plot.logitModel <- function(model, which = 1:4, col, ask,
#                              cookLevels = c(.5, 1), ...) {
#    ...
#    # Plot 02
#    if (2 %in% which) {
#      if (length(which) > 1 && ask && nextPlot() == "no") {
#        return(invisible())
#      }
#    ...
#  }
#  
#  nextPlot <- function() {
#    prompt <- "Enter anything to continue or [q] to quit: "
#    UInput <- readline(prompt=prompt)
#    if (UInput %in% c("q", "Q")) return(invisible("no"))
#    return("yes")
#  }

## ----include=FALSE-------------------------------------------------------
DonnerData$survived[3] <- 0

## ------------------------------------------------------------------------
fit2 <- logitModel(survived ~ age + sex + family, data = DonnerData)
pairs(fit2)

## ---- eval = FALSE-------------------------------------------------------
#  newtonRaphson <- function(y, X, precision = 1e-7, iterMax = 25) {
#    #01 initiate variables
#    beta <- rep(0, times = ncol(X))
#    i <- 0
#    convergence <- FALSE
#    ...

## ----eval=FALSE----------------------------------------------------------
#    ...
#    #02 compute coefficients until convergence or maximum iteration reached
#    while (i <= iterMax & convergence == FALSE) {
#      # update i
#      i <- i + 1
#  
#      # calculate probabilities (pi)
#      eta <- X %*% beta
#      pi <- as.numeric(logist(X %*% beta))
#  
#      # init / update Matrix W
#      W <- diag(pi * (1 - pi))
#  
#      # calculate and update beta (eq. 23 from Czepiel, 2002)
#      deltaBeta <- solve(t(X) %*% W %*% X) %*% t(X) %*% (y - pi)
#      beta <- beta + deltaBeta
#  
#      # check convergence condition
#      if (sum(abs(deltaBeta)) < precision) {
#        convergence <- TRUE
#        break
#      }
#    }
#  
#    if (convergence == FALSE) {
#      stop(paste(i, "iterations reached without convergence. Increase iterMax?"))
#    }
#  

## ---- eval=FALSE---------------------------------------------------------
#    ...
#    #03 return list of results
#    return(list(beta = beta,
#                pi = pi,
#                W = W,
#                eta = eta,
#                iterCount = i))
#  }

