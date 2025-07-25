# Libraries
# install.packages(c("Rcpp", "RcppEigen", "igraph", "rare"))
# install.packages("ComGGL/modified/ComGGL_1.0.tar.gz", repos = NULL)
# install.packages("clusterglasso/clusterglasso_0.1.0.tar.gz", repos = NULL)
library(Rcpp); library(RcppEigen); library(igraph); library(rare); library(ComGGL); library(clusterglasso)
# Source utility functions, this file contains the CV folds function
source("simulations/utils.R")

# Load cleaned data and sector info
load("applications/finance/data/rd_clean.RData")
sector_data = read.csv("applications/finance/data/sp100_sectors.csv")
sectors = unique(sector_data[,2])

# Find minimum standard deviation and scale entire data set.
X = scale(rd_clean)
# min_sd = min(apply(rd_clean, 2, sd))
# X = rd_clean / min_sd

# Folds for cross validation
seed = 20231116  # seed for the random number generator
set.seed(seed)
# folds = cv_folds(nrow(X), 4)
folds = cv_folds(nrow(X), 5)

# Create A matrix for taglasso based on sector info
A = diag(1, ncol(X))
for(s in 1:length(sectors)){
  get_sectors = which(sector_data[, 2]==sectors[s])
  new_col = rep(0, ncol(X))
  new_col[get_sectors] = 1
  A = cbind(A, new_col)
}
A = cbind(A, rep(1, ncol(X)))
colnames(A) = c(sector_data[, 1], sectors, "all")

# Initialize tuning parameters for taglasso
p = ncol(X)
ptm = proc.time()
lambda1_max_TAGL = binary_search(target = "lambda1", fun = taglasso,
                                  X = X, A = A, lambda2 = 0,
                                  refitting = TRUE, adaptive = FALSE)
timing1 = proc.time() - ptm # 100 seconds

S = cov(X)
S_offdiagonal = S - diag(p) * diag(S)
lambda_max_glasso = max(max(S_offdiagonal), -min(S_offdiagonal))
lambda_frac = lambda_grid(c(0, 1), n_points = 11, factor = 2)

# perform grid search to find regularization parameters
ptm = proc.time()
opt_lambda = cv_taglasso(X, A,
                          lambda1 = lambda_frac * lambda1_max_TAGL,
                          lambda2 = lambda_frac * lambda_max_glasso,
                          folds = folds, type = "grid",
                          refitting = TRUE, adaptive = FALSE)
timing2 = proc.time() - ptm # 5512.51 seconds

# Apply taglasso with fixed tuning parameter values
ptm = proc.time()
fit_opt <- taglasso(X, A, lambda1 = opt_lambda$lambda1_opt,
                    lambda2 = opt_lambda$lambda2_opt,
                    refitting = TRUE, adaptive = FALSE)
timing3 = proc.time() - ptm # 9 seconds
fit_opt$cluster
table(factor(fit_opt$cluster))

#
# save.image(file = "applications/finance/output/taglasso-finance-4fold.RData")
# save.image(file = "applications/finance/output/taglasso-finance-5fold.RData")
save.image(file = "applications/finance/output/taglasso-finance-5fold-std.RData")
