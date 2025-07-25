rm(list=ls())
library("CGGMR")

tr <- function(M)
{
  return(sum(diag(M)))
}

# Load cleaned data
load("applications/finance/data/rd_clean.RData")

# Load sectors
sectors = read.csv("applications/finance/data/sp100_sectors.csv")
sectors = data.frame(sectors$Sector, row.names = sectors$Symbol)

# Standardize data
X = scale(rd_clean)

# Inspect distribution of nonzero weights
og = par(mfrow = c(2, 3))
phis = 2 + 4 * c(0:5)
S = cov(X)
for (phi in phis) {
  W = cggm_weights(S, phi = phi, k = round(sqrt(ncol(X))), connected = FALSE)
  hist(W[W>0], breaks = seq(0, 1, 0.025))
}
par(mfrow = og)
rm(S, W, og)

# Set sequence for k
ks = seq(1, 15, 2)

# Folds for cross validation
seed = 20231116  # seed for the random number generator
set.seed(seed)
folds_outer = cv_folds(nrow(X), 10)
outer_weights = sapply(folds_outer, length) / sum(sapply(folds_outer, length))

# Track statistics of results
clusters_outer = rep(NA, length(folds_outer))
cvscore_outer = matrix(NA, nrow = length(folds_outer), ncol = 3)
colnames(cvscore_outer) = c("NLL", "KLD", "norm")

# Store optimal models for each of the folds
outer_fit_opts = list()

# Filename for output
file_name = paste0("applications/finance/output/cggm-finance-5fold-oos-std",
                   "-refit.RData")

# Start time
t_start = Sys.time()

# Perform nested cv
for(outer_i in 1:length(folds_outer)) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S\n"))
  cat("Start outer iteration", outer_i, '\n')

  # Train and test sets
  X_train = X[-folds_outer[[outer_i]], ]
  X_test = X[folds_outer[[outer_i]], ]
  S_test = cov(X_test)

  # Folds for inner cv
  set.seed(seed)
  folds_inner = cv_folds(nrow(X_train), 5)

  ## Parallel cross validation for CGGM
  # Make cluster
  cl = parallel::makePSOCKcluster(
    min(length(phis), parallel::detectCores(logical = FALSE))
  )

  # Export variables and functions to the workers
  parallel::clusterExport(cl, c("X_train", "folds_inner", "ks"),
                          envir=environment())
  parallel::clusterEvalQ(cl, library(CGGMR))

  # Perform cross validation
  fit_opts = parallel::parLapply(cl, phis, function(phi) {
    # Cross validation for CGGM
    fit_opt = cggm_cv(
      X = X_train,
      folds = folds_inner,
      tune_grid = expand.grid(phi = phi, k = ks),
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

  # Compute cv scores
  Theta = get_Theta(fit_opt)
  cvscore_outer[outer_i, "NLL"] = -log(det(Theta)) + tr(S_test %*% Theta)
  cvscore_outer[outer_i, "KLD"] = -log(det(Theta %*% S_test)) +
    tr(Theta %*% S_test)
  cvscore_outer[outer_i, "norm"] = norm(solve(Theta) - S_test, "F")

  # Number of clusters
  clusters_outer[outer_i] = length(unique(get_clusters(fit_opt)))

  # Store fit_opt
  outer_fit_opts[[outer_i]] = fit_opt

  # Do some time computations to estimate finish time
  t_intermediary = Sys.time()
  t_avg_per_fold = (t_intermediary - t_start) / outer_i
  t_est_total = t_avg_per_fold * length(folds_outer)

  cat(format(Sys.time(), "\n%Y-%m-%d %H:%M:%S\n"))
  cat("Scores:\n")
  print(cvscore_outer)
  cat("\nMean scores (weighted with fold size):\n")
  print(apply(cvscore_outer * outer_weights, 2, sum, na.rm = TRUE) /
          sum(outer_weights[1:outer_i]))
  cat("\nNumber of clusters =", clusters_outer[outer_i], "\n")
  if (outer_i < length(folds_outer)) {
    cat("\nETA:", format(t_start + t_est_total, "\n%Y-%m-%d %H:%M:%S\n"))
  }
  cat("\n-------------------------------\n\n")

  # Save after each fold
  save.image(file = file_name)
}

# Cleanup
rm(cl)

# Save
save.image(file = file_name)

# Load
load(file = file_name)
