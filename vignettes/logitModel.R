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

op <- par(mfrow = c(1,2), mar = c(-.2, 0, -1, -0.5) + 1.5, cex = 0.9,
          cex.axis = 0.9, tcl = -0.15, oma = c(0, 0, 0, 0), mgp = c(2, 0.5, 0))

lim <- 5.5
curve(logitModel:::logist(x), from = -lim, to = lim, ylab = "",
      main = "")
legend(-2.5,.8, expression(beta > 0), bty = "n")
curve(logitModel:::logist(-x), from = -lim, to = lim, 
      main = "")
legend(0,0.4, expression(beta < 0), bty = "n")

