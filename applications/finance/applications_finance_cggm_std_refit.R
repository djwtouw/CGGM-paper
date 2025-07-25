rm(list=ls())
library("CGGMR")

# Source utility functions, this file contains the CV folds function
source("simulations/utils.R")

# Load cleaned data
load("applications/finance/data/rd_clean.RData")

# Load sectors
sectors = read.csv("applications/finance/data/sp100_sectors.csv")
sectors = data.frame(sectors$Sector, row.names = sectors$Symbol)

# Standardize data
X = scale(rd_clean)

# Settings
kfold = 5

# Folds for cross validation
seed = 20231116  # seed for the random number generator
set.seed(seed)
folds = cv_folds(nrow(X), kfold)

# Inspect distribution of nonzero weights
og = par(mfrow = c(2, 3))
phis = 2 + 4 * c(0:5)
S = cov(X)
for (phi in phis) {
  W = cggm_weights(S, phi = phi, k = 10, connected = FALSE)
  hist(W[W>0], breaks = seq(0, 1, 0.025))
}
par(mfrow = og)


## Parallel cross validation for CGGM
# Make cluster
cl = parallel::makePSOCKcluster(
  min(length(phis), parallel::detectCores(logical = FALSE))
)

# Export variables and functions to the workers
parallel::clusterExport(cl, c("X", "folds"),
                        envir=environment())
parallel::clusterEvalQ(cl, library(CGGMR))

# Perform cross validation
fit_opts = parallel::parLapply(cl, phis, function(phi) {
  # Cross validation for CGGM
  fit_opt = cggm_cv(
    X = X,
    folds = folds,
    tune_grid = expand.grid(phi = phi,
                            k = seq(1, 15, 2)),
    refit = TRUE
  )

  return(fit_opt)
})

# Stop cluster
parallel::stopCluster(cl)


## Combine the results from the parallel cross validation and select the best
# Store best scores from objects
best_cv_scores = do.call(rbind, lapply(1:length(fit_opts), function(i) {
  for (j in 1:nrow(fit_opts[[i]]$scores)) {
    if (sum(fit_opts[[i]]$scores[j, c("phi", "k", "lambda")] ==
            fit_opts[[i]]$opt_tune) == 3)
      index = j
  }

  return(fit_opts[[i]]$scores[index, ])
}))

# Sort scores, take into account the interval length for lambda
best_cv_scores_sorted =
  dplyr::arrange(cbind(1:nrow(best_cv_scores), best_cv_scores),
                 score, dplyr::desc(lambda_intv_length))

# Select index with lowest score
fit_opt_index = best_cv_scores_sorted[1, 1]

# Select fit_opt
fit_opt = fit_opts[[fit_opt_index]]

# Store scores from other objects in the same one
fit_opt$scores = do.call(rbind, lapply(1:length(fit_opts), function(i) {
  return(fit_opts[[i]]$scores)
}))
fit_opt$scores = dplyr::arrange(fit_opt$scores, k, phi)

# Do some cleanup
rm(best_cv_scores, best_cv_scores_sorted, fit_opt_index, fit_opts)

# Save (and load) results
save.image("applications/finance/output/cggm-finance-5fold-std-refit.RData")
load("applications/finance/output/cggm-finance-5fold-std-refit.RData")

fit_opt$opt_tune
get_clusters(fit_opt)
table(get_clusters(fit_opt))


## Comparison with taglasso, need taglasso results for matching data set
# Load cggm results
load("applications/finance/output/cggm-finance-5fold-std-refit.RData")
# Load taglasso results
fit_opt_cggm = fit_opt
load("applications/finance/output/taglasso-finance-5fold-std.RData")
fit_opt_tagl = fit_opt

# Visualize a heatmap of the results
library(pheatmap)
T_cggm = get_Theta(fit_opt_cggm)
T_tagl = fit_opt_tagl$omega_full
colnames(T_tagl) = colnames(T_cggm)
rownames(T_tagl) = rownames(T_cggm)
breaks = seq(-0.1, 0.1, length.out = 101)
pheatmap(T_cggm, cluster_rows = FALSE, cluster_cols = FALSE, breaks = breaks)
pheatmap(T_tagl, cluster_rows = FALSE, cluster_cols = FALSE, breaks = breaks)
