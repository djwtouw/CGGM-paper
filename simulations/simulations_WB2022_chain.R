# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## install Eugen's package ComGGL, which is available from
## https://perso.uclouvain.be/eugen.pircalabelu/software/ComGGL_1.0.tar.gz
## (even though 'igraph' is not listed in DESCRIPTION as a required package,
## it is imported in NAMESPACE)
# install.packages(c("Rcpp", "RcppEigen", "igraph"))
# install.packages("ComGGL/modified/ComGGL_1.0.tar.gz", repos = NULL)

## Here a version of ComGGL with the following modifications is used:
##  - No information is printed on the status of the algorithm or
##    cross-validation.
##  - The seed of the random number generator is not set internally
## That is, the lines corresponding to the above actions are commented out,
## but no other modifications were made to the code.

# install package clusterglasso for TAGL, which is available in the private
# GitHub repository https://github.com/ineswilms/clusterglasso
# install.packages("clusterglasso/clusterglasso_0.1.0.tar.gz", repos = NULL)

# install our package for CGGM, which is available in the private GitHub
# repository https://github.com/djwtouw/CGGMR
# install.packages("CGGMR/CGGMR_0.1.0.tar.gz", repos = NULL)

## load packages
library("CGGMR")          # clusterpath Gaussian graphical modeling
library("ComGGL")         # community-based group graphical lasso
library("clusterglasso")  # our old package (for tag-lasso)
library("mvtnorm")
source("simulations/evaluation_criteria.R")
source("simulations/utils.R")


## control parameters for simulation
seed <- 20231116  # seed for the random number generator
n <- 120          # number of observations
p <- 15           # number of variables
B <- 3            # number of blocks
theta_w <- 0.5    # (off-diagonal) elements within blocks
theta_b <- 0.25   # non-zero elements between blocks
R <- 100          # number of simulation runs

## additional parameters
blocks <- sort(rep(1:B, length.out = p))  # works if p is not a multiple of B
p_b <- tabulate(blocks)                   # number of variables per block


## control parameters for CGGM
# tuning parameters (regularization parameter is determined automatically)
phi <- seq(from = 1, to = 3, by = 0.5) # tuning parameter for weights
k <- 1:5                               # number of nearest neighbors for weights
tune_grid <- expand.grid(phi = phi, k = k)
# other parameters to loop over
args_CGGM <- list(
  "CGGM-raw" = list(refit = FALSE),
  "CGGM-refit" = list(refit = TRUE)
)

## control parameters for other methods
# sparsity parameter will be set based on smallest lambda that sets all
# elements to 0 in graphical lasso
lambda_frac <- lambda_grid(c(0, 1), n_points = 10, factor = 2)
# ComGGL
lambda2_max_ComGGL <- 1  # largest grouping parameter
lambda3_ComGGL <- 1      # balance parameter
# default values except for tolerance for convergence
control_ComGGL <- list(penalize.diagonal = FALSE, plotConvergence = FALSE,
                       rho = 1, maxiterADMM = 500, maxiterKMEANS = 500,
                       tol = 1e-06, alpha = 0.7)

## other control parameters for cross-validation
K <- 3  # number of folds


## utility functions for data generation

# function to generate diagonal blocks of precision matrix
generate_diagonal_block <- function(d, theta) {
  mat <- matrix(theta, nrow = d, ncol = d)
  diag(mat) <- 1
  mat
}

# function to generate offdiagonal blocks of precision matrix
generate_offdiagonal_block <- function(nrow, ncol, theta) {
  matrix(theta, nrow = nrow, ncol = ncol)
}


## store R session info with results
session_info <- sessionInfo()


## run simulation
set.seed(seed)
cat(paste(Sys.time(), ": starting ...\n"))
results_list <- lapply(1:R, function(r) {

  ## print simulation run
  cat(paste(Sys.time(), sprintf(":     run = %d\n", r)))

  # initialize precision matrix as a list of lists:
  # (outer list corresponds to columns, inner list to rows of the given column)
  Theta <- replicate(B, replicate(B, NULL, simplify = FALSE),
                     simplify = FALSE)
  # loop over indices of blocks in the final precision matrix, and generate
  # those blocks following the R convention of building a matrix by column
  # (i is the index of the row, j is the index of the column)
  for (j in seq_len(B)) {
    for (i in seq_len(B)) {
      if (i == j) Theta[[j]][[i]] <- generate_diagonal_block(p_b[i], theta_w)
      else if (i == j+1) {
        Theta[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j], theta_b)
      } else if (i > j+1) {
        Theta[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j], 0)
      } else Theta[[j]][[i]] <- t(Theta[[i]][[j]])  # upper diagonal blocks
    }
  }
  # put precision matrix together
  Theta <- do.call(cbind, lapply(Theta, function(column) do.call(rbind, column)))

  ## generate covariance matrix by inverting
  Sigma <- solve(Theta)

  ## generate training data
  X <- rmvnorm(n, sigma = Sigma)

  ## compute sample covariance matrix
  S <- cov(X)

  ## estimate the smallest lambda1 that sets everything to 0 in glasso
  S_offdiagonal <- S - diag(p) * diag(S)
  lambda_max_glasso <- max(max(S_offdiagonal), -min(S_offdiagonal))

  ## generate folds for cross-validation
  folds <- cv_folds(n, K = K)

  ## tree aggregation matrices for TAG-lasso
  A <- get_tree_matrices(blocks)

  ## apply CGGM with and without refitting
  df_CGGM <- mapply(function(args, label) {
    tryCatch({
      # perform grid search to find tuning parameters and re-fit with optimal
      # tuning parameters
      fit_cv <- cggm_cv(X, tune_grid = tune_grid, folds = folds,
                        refit = args$refit)
      Theta_hat <- get_Theta(fit_cv)
      clusters <- get_clusters(fit_cv)
      # compute evaluation criteria and combine into data frame
      data.frame(Run = r, n = n, p = p, B = B, Method = label,
                 Lambda = fit_cv$opt_tune$lambda,
                 Tuning1 = fit_cv$opt_tune$phi,
                 Tuning2 = fit_cv$opt_tune$k,
                 Error = Frobenius(Theta, Theta_hat),
                 NLL = neg_log_lik(Sigma, Theta_hat),
                 B_hat = length(unique(clusters)),
                 ARI = adjusted_Rand_index(blocks, clusters),
                 stringsAsFactors = FALSE)
    }, error = function(condition) NULL)
  }, args = args_CGGM, label = names(args_CGGM),
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  df_CGGM <- do.call(rbind, df_CGGM)

  ## apply tag-lasso with ideal, realistic, and misspecified tree hierarchies
  df_TAGL <- mapply(function(current_A, label) {
    tryCatch({
      # find aggregation parameter that yields maximal aggregation
      lambda1_max_TAGL <- binary_search(target = "lambda1", fun = taglasso,
                                        X = X, A = current_A, lambda2 = 0,
                                        refitting = TRUE, adaptive = FALSE)
      # perform grid search to find regularization parameters
      opt_lambda <- cv_taglasso(X, current_A,
                                lambda1 = lambda_frac * lambda1_max_TAGL,
                                lambda2 = lambda_frac * lambda_max_glasso,
                                folds = folds, type = "grid",
                                refitting = TRUE, adaptive = FALSE)
      # re-fit with optimal regularization parameters
      fit_opt <- taglasso(X, current_A, lambda1 = opt_lambda$lambda1_opt,
                          lambda2 = opt_lambda$lambda2_opt,
                          refitting = TRUE, adaptive = FALSE)
      # compute evaluation criteria and combine into data frame
      data.frame(Run = r, n = n, p = p, B = B, Method = label,
                 Lambda = opt_lambda$lambda1_opt,
                 Tuning1 = opt_lambda$lambda2_opt,
                 Tuning2 = NA_real_,
                 Error = Frobenius(Theta, fit_opt$omega_full),
                 NLL = neg_log_lik(Sigma, fit_opt$omega_full),
                 B_hat = length(unique(fit_opt$cluster)),
                 ARI = adjusted_Rand_index(blocks, fit_opt$cluster),
                 stringsAsFactors = FALSE)
    }, error = function(condition) NULL)
  }, current_A = A, label = paste("TAGL", names(A), sep = "-"),
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  df_TAGL <- do.call(rbind, df_TAGL)

  ## apply ComGGL
  df_ComGGL <- tryCatch({
    # perform grid search to find regularization parameters
    opt_lambda <- cv_ComGGL(X, lambda1 = lambda_frac * lambda_max_glasso,
                            lambda2 = lambda_frac * lambda2_max_ComGGL,
                            lambda3 = lambda3_ComGGL, folds = folds,
                            type = "grid", list.controls = control_ComGGL)
    # re-fit with optimal regularization parameters
    fit_opt <- ComGGL(X, lambda1 = opt_lambda$lambda1_opt,
                      lambda2 = opt_lambda$lambda2_opt,
                      lambda3 = lambda3_ComGGL,
                      list.controls = control_ComGGL)
    # compute evaluation criteria and combine into data frame
    data.frame(Run = r, n = n, p = p, B = B, Method = "ComGGL",
               Lambda = opt_lambda$lambda2_opt,
               Tuning1 = opt_lambda$lambda1_opt,
               Tuning2 = lambda3_ComGGL,
               Error = Frobenius(Theta, fit_opt$Theta),
               NLL = neg_log_lik(Sigma, fit_opt$Theta),
               B_hat = length(unique(fit_opt$Components)),
               ARI = adjusted_Rand_index(blocks, fit_opt$Components),
               stringsAsFactors = FALSE)
  }, error = function(condition) NULL)

  # inverse of sample covariance matrix
  df_Sinv <- tryCatch({
    # compute inverse of sample covariance matrix
    Sinv <- solve(S)
    # compute evaluation criteria and combine into data frame
    data.frame(Run = r, n = n, p = p, B = B, Method = "Sinv",
               Lambda = NA_real_, Tuning1 = NA_real_,
               Tuning2 = NA_real_, Error = Frobenius(Theta, Sinv),
               NLL = neg_log_lik(Sigma, Sinv), B_hat = p,
               ARI = adjusted_Rand_index(blocks, 1:p), stringsAsFactors = FALSE)
  })

  ## return results
  rbind(df_CGGM, df_TAGL, df_ComGGL, df_Sinv)

})

## combine results into data frame
results <- do.call(rbind, results_list)

## store results
file <- "simulations/results/results_WB2022_chain_n=%d_p=%d_B=%d_R=%d.RData"
save(results, n, p, B, theta_w, theta_b, session_info,
     file = sprintf(file, n, p, B, R))

## done with simulation
cat(paste(Sys.time(), ": finished.\n"))
