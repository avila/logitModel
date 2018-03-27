####------------
testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
testData$rank <- factor(testData$rank)
####------------

# debugonce(MLE)
myFit <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
statFit <- glm(formula = admit ~ gre + gpa + rank, data = testData, family = binomial)

myFit
statFit
## -----------------------------------
rm(myFit, statFit)
##------------------------------------
set.seed(42)
X1 <- rnorm(100); X2 <- rnorm(100)
points <- 2 * X1 - 3 * X2
y <- rbinom(n, 1, exp(points) / (1 + exp(points)))
X <- cbind(b1 = 1, b2 = X1, b3 = X2)
fit <- MLE(y = y, X = X)
fit$coefficients

