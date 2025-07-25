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
library("mvtnorm")
source("simulations/evaluation_criteria.R")
source("simulations/utils.R")


## control parameters for simulation
seed <- 20231116  # seed for the random number generator
n <- 120          # number of observations
p <- 15           # number of variables
B <- 3            # number of blocks
sigma_w <- 0.5    # (off-diagonal) elements within blocks
sigma_b <- 0.25   # non-zero elements between blocks
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
  "CGGM_Theta_inv" = list(estimate_Sigma = FALSE),
  "CGGM_Sigma" = list(estimate_Sigma = TRUE)
)

## other control parameters for cross-validation
K <- 3  # number of folds


## utility functions for data generation

# function to generate diagonal blocks of covariance matrix
generate_diagonal_block <- function(d, sigma) {
  mat <- matrix(sigma, nrow = d, ncol = d)
  diag(mat) <- 1
  mat
}

# function to generate offdiagonal blocks of covariance matrix
generate_offdiagonal_block <- function(nrow, ncol, sigma) {
  matrix(sigma, nrow = nrow, ncol = ncol)
}


## store R session info with results
session_info <- sessionInfo()


## run simulation
set.seed(seed)
cat(paste(Sys.time(), ": starting ...\n"))
results_list <- lapply(1:R, function(r) {

  ## print simulation run
  cat(paste(Sys.time(), sprintf(":     run = %d\n", r)))

  # initialize covariance matrix as a list of lists:
  # (outer list corresponds to columns, inner list to rows of the given column)
  Sigma <- replicate(B, replicate(B, NULL, simplify = FALSE),
                     simplify = FALSE)
  # loop over indices of blocks in the final covariance matrix, and generate
  # those blocks following the R convention of building a matrix by column
  # (i is the index of the row, j is the index of the column)
  for (j in seq_len(B)) {
    for (i in seq_len(B)) {
      if (i == j) Sigma[[j]][[i]] <- generate_diagonal_block(p_b[i], sigma_w)
      else if (i == j+1) {
        Sigma[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j], sigma_b)
      } else if (i > j+1) {
        Sigma[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j], 0)
      } else Sigma[[j]][[i]] <- t(Sigma[[i]][[j]])  # upper diagonal blocks
    }
  }
  # put covariance matrix together
  Sigma <- do.call(cbind, lapply(Sigma, function(column) do.call(rbind, column)))

  ## generate training data
  X <- rmvnorm(n, sigma = Sigma)

  ## generate folds for cross-validation
  folds <- cv_folds(n, K = K)

  ## apply CGGM for clustering precision matrix and take the inverse
  df_CGGM_Theta_inv <- tryCatch({
    # perform grid search to find tuning parameters and re-fit with optimal
    # tuning parameters
    fit_cv <- cggm_cv(X, tune_grid = tune_grid, folds = folds, refit = TRUE,
                      estimate_Sigma = FALSE)
    Sigma_hat <- solve(get_Theta(fit_cv))
    clusters <- get_clusters(fit_cv)
    # compute evaluation criteria and combine into data frame
    data.frame(Run = r, n = n, p = p, B = B, Method = "CGGM_Theta_inv",
               Lambda = fit_cv$opt_tune$lambda,
               Tuning1 = fit_cv$opt_tune$phi,
               Tuning2 = fit_cv$opt_tune$k,
               Error = Frobenius(Sigma, Sigma_hat),
               B_hat = length(unique(clusters)),
               ARI = adjusted_Rand_index(blocks, clusters),
               stringsAsFactors = FALSE)
  }, error = function(condition) NULL)

  ## apply modified CGGM for clustering covariance matrix
  df_CGGM_Sigma <- tryCatch({
    # perform grid search to find tuning parameters and re-fit with optimal
    # tuning parameters
    fit_cv <- cggm_cv(X, tune_grid = tune_grid, folds = folds, refit = TRUE,
                      estimate_Sigma = TRUE)
    Sigma_hat <- get_Theta(fit_cv)
    clusters <- get_clusters(fit_cv)
    # compute evaluation criteria and combine into data frame
    data.frame(Run = r, n = n, p = p, B = B, Method = "CGGM_Sigma",
               Lambda = fit_cv$opt_tune$lambda,
               Tuning1 = fit_cv$opt_tune$phi,
               Tuning2 = fit_cv$opt_tune$k,
               Error = Frobenius(Sigma, Sigma_hat),
               B_hat = length(unique(clusters)),
               ARI = adjusted_Rand_index(blocks, clusters),
               stringsAsFactors = FALSE)
  }, error = function(condition) NULL)

  # sample covariance matrix
  df_S <- tryCatch({
    # computesample covariance matrix
    S <- cov(X)
    # compute evaluation criteria and combine into data frame
    data.frame(Run = r, n = n, p = p, B = B, Method = "S",
               Lambda = NA_real_, Tuning1 = NA_real_,
               Tuning2 = NA_real_, Error = Frobenius(Sigma, S),
               B_hat = p, ARI = adjusted_Rand_index(blocks, 1:p),
               stringsAsFactors = FALSE)
  })

  ## return results
  rbind(df_CGGM_Theta_inv, df_CGGM_Sigma, df_S)

})

## combine results into data frame
results <- do.call(rbind, results_list)

## store results
file <- "simulations/results/results_Sigma_block_n=%d_p=%d_B=%d_R=%d.RData"
save(results, n, p, B, sigma_w, sigma_b, session_info,
     file = sprintf(file, n, p, B, R))

## done with simulation
cat(paste(Sys.time(), ": finished.\n"))
