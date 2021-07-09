# logitModel

As part of the Statistical Programming with R course ministered in the Winter Term of 2017 / 2018, the 
final project was to write a R package from scratch. 
This pacakge implements an logistic regression estimation via Maximum Likelihood Estimator based on a
Newton-Rahpson algorithm. It also includes a formula interface for easier usage, a `print`, `summary` and `plot` 
S3-Methods that mimics the `glm(..., family = binomial)` implementation of the stats package as well
as a `pairs` method for an overview of the interaction between the explanatory and explained variables.

# Installation 

``` r
#install.packages("devtools")
devtools::install_github("avila/logitModel")
```

# Vignette

To read the full vignette explaining the package, run from R: 

``` r
vignette("logitModel")
```
