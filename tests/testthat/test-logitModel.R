# testData <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")[1:100,]
# testData$rank <- factor(testData$rank)
# testModell <- as.formula("admit ~ gre + gpa + rank")
# testModelFrame <- model.frame(admit ~ gre + gpa + rank, testData)
myFit <- logitMod(formula = admit ~ gre + gpa + rank, data = testData)
statFit <- glm(formula = admit ~ gre + gpa + rank, data = testData, family = binomial)

myFit
statFit
