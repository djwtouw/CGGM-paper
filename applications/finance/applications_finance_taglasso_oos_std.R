# Source utility functions, this file contains the CV folds function
source("simulations/utils.R")

# Load libraries
# install.packages(c("Rcpp", "RcppEigen", "igraph"))
# install.packages("ComGGL/modified/ComGGL_1.0.tar.gz", repos = NULL)
# install.packages("clusterglasso/clusterglasso_0.1.0.tar.gz", repos = NULL)

# Load cleaned data and sector info
load("applications/finance/data/rd_clean.RData")
sector_data = read.csv("applications/finance/data/sp100_sectors.csv")
sectors = unique(sector_data[,2])

# Find minimum standard deviation and scale entire data set.
X = scale(rd_clean)
# min_sd = min(apply(rd_clean, 2, sd))
# X = rd_clean / min_sd

# Create A matrix for taglasso based on sector info
p = ncol(X)
A = diag(1, ncol(X))
for(s in 1:length(sectors)){
  get_sectors = which(sector_data[, 2]==sectors[s])
  new_col = rep(0, ncol(X))
  new_col[get_sectors] = 1
  A = cbind(A, new_col)
}
A = cbind(A, rep(1, ncol(X)))
colnames(A) = c(sector_data[, 1], sectors, "all")

# Folds for cross validation
seed = 20231116  # seed for the random number generator
set.seed(seed)
folds_outer = cv_folds(nrow(X), 10)
lambda_frac = lambda_grid(c(0, 1), n_points = 11, factor = 2)
cvscore_outer  = clusters_outer = rep(NA, 10)


for(i in 1:length(folds_outer)){
  cat("start outer iteration", i,'\n')

  Xtrain = X[-folds_outer[[i]], ]
  Xtest = X[folds_outer[[i]], ]
  Stest = cov(Xtest)

  set.seed(seed)
  # folds_train = cv_folds(nrow(Xtrain), 4)
  folds_train = cv_folds(nrow(Xtrain), 5)

  # Initialize tuning parameters for taglasso
  ptm = proc.time()
  lambda1_max_TAGL = binary_search(target = "lambda1", fun = taglasso,
                                   X = Xtrain, A = A, lambda2 = 0,
                                   refitting = TRUE, adaptive = FALSE)
  timing1 = proc.time() - ptm

  S = cov(Xtrain)
  S_offdiagonal = S - diag(p) * diag(S)
  lambda_max_glasso = max(max(S_offdiagonal), -min(S_offdiagonal))

  # perform grid search to find regularization parameters
  ptm = proc.time()
  opt_lambda = cv_taglasso(Xtrain, A,
                           lambda1 = lambda_frac * lambda1_max_TAGL,
                           lambda2 = lambda_frac * lambda_max_glasso,
                           folds = folds_train, type = "grid",
                           refitting = TRUE, adaptive = FALSE)
  timing2 = proc.time() - ptm

  # Apply taglasso with fixed tuning parameter values
  ptm = proc.time()
  fit_opt <- taglasso(Xtrain, A, lambda1 = opt_lambda$lambda1_opt,
                      lambda2 = opt_lambda$lambda2_opt,
                      refitting = TRUE, adaptive = FALSE)
  timing3 = proc.time() - ptm

  # Out of sample cv score
  cvscore_outer[i] = - log(det(fit_opt$omega_full)) + sum(diag(Stest%*%fit_opt$omega_full))
  clusters_outer[i] = length(unique(fit_opt$cluster))

  cat("end outer iteration", i, " timings:", timing1[3], timing2[3], timing3[3], "\n")
}

# save.image(file = "applications/finance/output/taglasso-finance-4fold-oos.RData")
# save.image(file = "applications/finance/output/taglasso-finance-5fold-oos.RData")
save.image(file = "applications/finance/output/taglasso-finance-5fold-oos-std.RData")
