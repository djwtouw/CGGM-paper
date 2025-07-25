# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## data are obtained from https://openpsychometrics.org/_rawdata/HSQ.zip

## load packages
library("CGGMR")               # clusterpath Gaussian graphical modeling
source("simulations/utils.R")  # cross-validation folds


## preprocess data

# load data
HSQ <- read.csv("applications/HSQ/data.csv")

# keep only respondents who reported a high accuracy of their responses
keep <- HSQ$accuracy == 100
HSQ <- HSQ[keep, ]

# keep only rating-scale items
HSQ <- HSQ[, 1:32]

# remove observations with missing values (coded as -1)
HSQ[HSQ == -1] <- NA
keep <- !apply(is.na(HSQ), 1, any)
HSQ <- HSQ[keep, ]

# which variables belong to which scale
affiliative   <- c(1, 5,  9, 13, 17, 21, 25, 29)
selfenhancing <- c(2, 6, 10, 14, 18, 22, 26, 30)
aggressive <-    c(3, 7, 11, 15, 19, 23, 27, 31)
selfdefeating <- c(4, 8, 12, 16, 20, 24, 28, 32)

# reorder variables
HSQ <- HSQ[, c(affiliative, selfenhancing, aggressive, selfdefeating)]

# reverse scale of negatively keyed items
keys <- c(
  # affiliative
  -1,  1, -1,  1, -1,  1, -1, -1,
  # selfenhancing
   1,  1,  1,  1,  1, -1,  1,  1,
  # aggressive
   1, -1,  1, -1,  1, -1,  1, -1,
  # selfdefeating
   1,  1,  1, -1,  1,  1,  1,  1
)
offset <- ifelse(keys == 1, 0, 6)
HSQ <- sweep(HSQ, 2, keys, FUN = "*")
HSQ <- sweep(HSQ, 2, offset, FUN = "+")

# check that all items of a scale are in the same key
cov(HSQ[, 1:8])
cov(HSQ[, 9:16])
cov(HSQ[, 17:24])
cov(HSQ[, 25:32])


## CGGM

# tuning parameters (regularization parameter is determined automatically)
phi <- seq(from = 1, to = 3, by = 0.5) # tuning parameter for weights
k <- 1:5                               # number of nearest neighbors for weights
tune_grid <- expand.grid(phi = phi, k = k)

# generate folds for cross-validation (to use the same folds for each trait)
set.seed(20240526)
folds <- cv_folds(n = nrow(HSQ), K = 5)  # 5-fold CV

# fit CGGM
fit <- cggm_cv(HSQ, tune_grid = tune_grid, folds = folds, refit = TRUE,
               estimate_Sigma = TRUE)
clusters <- get_clusters(fit)
print(table(clusters))
print(clusters)

## save results to file
save(fit, file = "applications/HSQ/results.RData")


## check covariances of misclassified items

min(cov(HSQ[, "Q6", drop = FALSE], HSQ[, 1:8]))
min(cov(HSQ[, "Q6", drop = FALSE], HSQ[, setdiff(names(HSQ)[9:16], "Q6")]))

mean(cov(HSQ[, "Q30", drop = FALSE], HSQ[, 1:8]))
mean(cov(HSQ[, "Q30", drop = FALSE], HSQ[, setdiff(names(HSQ)[9:16], "Q30")]))

mean(cov(HSQ[, "Q28", drop = FALSE], HSQ[, 1:8]))
mean(cov(HSQ[, "Q28", drop = FALSE], HSQ[, 9:16]))
mean(cov(HSQ[, "Q28", drop = FALSE], HSQ[, 17:24]))
mean(cov(HSQ[, "Q28", drop = FALSE], HSQ[, setdiff(names(HSQ)[25:32], "Q28")]))


## correlations of scores
cor(rowMeans(HSQ[, 1:8]), rowMeans(HSQ[, clusters == 1]))
cor(rowMeans(HSQ[, 9:16]), rowMeans(HSQ[, clusters == 2]))
cor(rowMeans(HSQ[, 17:24]), rowMeans(HSQ[, clusters == 3]))
cor(rowMeans(HSQ[, 25:32]), rowMeans(HSQ[, clusters == 4]))
