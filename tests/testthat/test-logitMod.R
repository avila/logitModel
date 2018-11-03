## logitMod Function ----
context("logitModel Function: default")
set.seed(42); n <- 100
X1 <- rnorm(n); X2 <- rnorm(n)
points <- 2 * X1 - 3 * X2
y <- rbinom(n, 1, exp(points) / (1 + exp(points)))
df <- as.data.frame(cbind(b2 = X1, b3 = X2))
myFit <- logitModel(y ~ X1 + X2, data = df)
statFit <- glm(y ~ X1 + X2, family = binomial)
tolerance <- 1e-7
testthat::expect_equal(myFit$coefficients,
                       statFit$coefficients, tolerance = tolerance)



## logitModel Function: Sem 'data'" ----
context("logitModel Function: Sem 'data'")
myFit <- logitModel(y ~ X1 + X2)
testthat::expect_equal(myFit$coefficients,
                       statFit$coefficients, tolerance = tolerance)
all.equal(unname(myFit$coefficients),
          unname(statFit$coefficients), tolerance = tolerance)



## logitModel Function: With Dummies ----
context("logitModel Function: With Dummies")
set.seed(42)
n <- 100
dummies <- sample(LETTERS[1:5], n, replace = TRUE)
X1 <- rnorm(n); X2 <- rnorm(n)
points <- 2 * X1 - 3 * X2
y <- rbinom(n, 1, exp(points) / (1 + exp(points)))
df <- as.data.frame(cbind(b2 = X1, b3 = X2))
myFitD <- logitModel(y ~ X1 + X2 + dummies)
statFitD <- glm(y ~ X1 + X2 + dummies, family = binomial)
testthat::expect_equal(
  myFitD$coefficients, statFitD$coefficients, tolerance = tolerance
)

testthat::expect_equal(
  myFitD$AIC, statFitD$aic
)




## newtonRaphson function ----
context("newtonRaphson Function")
set.seed(42)
n <- 100
X1 <- rnorm(n);    X2 <- rnorm(n);
points <- 2 * X1 - 3 * X2
y <- rbinom(n, 1, exp(points) / (1 + exp(points)))
X <- cbind(b1 = 1, b2 = X1, b3 = X2)
myMLE <- newtonRaphson(y = y, X = X)
statFit <- glm(formula = y ~ X1 + X2, family = binomial)
tolerance <- 1e-7

# check if beta coef are the same as stats:::glm()
testthat::expect_equal(c(myMLE$beta), unname(c(statFit$coefficients)),
                       tolerance = tolerance)

all.equal(c(myMLE$beta), unname(c(statFit$coefficients)), tolerance = tolerance)





context("Donner Data")
DonnerData <- DonnerData

myfitDonner <- logitModel(survived ~ age + sex, data = DonnerData)
statFitDonner <- glm(survived ~ age + sex, data = DonnerData)
testthat::expect_equal(myFit$coefficients,
                       statFit$coefficients, tolerance = tolerance)




