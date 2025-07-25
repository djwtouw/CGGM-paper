rm(list=ls())
library("CGGMR")

source("applications/utils/to_hclust.R")


## SETTINGS for the analysis
estimate_Theta = TRUE
cluster = 1
std = TRUE
## SETTINGS for the analysis


# Load cleaned data
load("applications/oecd/data/OECDData.RData")

# Standardize data
if (cluster == 1) {
  X = df1
} else {
  X = df2
}
if (std) {
  X = scale(X)
} else {
  X = as.matrix(X)
}

# Compute sample covariance matrix
if (estimate_Theta) {
  S = cov(X)
} else {
  # Compute sample precision matrix
  S = solve(cov(X))
}


# Compute weights
W = cggm_weights(S, phi = 0.5, k = 3)
hist(W[W > 0] / max(W))

# Compute clusterpath
res = cggm(S, W, lambda = seq(0, 8.1, 0.005), expand = TRUE)
res = cggm_refit(res)

# Plot dendrogram
res_hclust = to_hclust(res, height_type = "index")
if (estimate_Theta) {
  plot_title = "Estimated Theta"
} else {
  plot_title = "Estimated Sigma"
}
plot(res_hclust, main = plot_title, xlab = "", sub = "", hang = -1, ylab = "")

# Save
save.image(paste0("applications/oecd/output/results_cluster_", cluster,
                  ".RData"))

# Load
load(paste0("applications/oecd/output/results_cluster_", cluster, ".RData"))

# Clusterings obtained by Cavicchia, Vichi, and Zaccaria
cvz_1 = c(2, 1, 2, 3, 3, 1, 1, 1, 2, 1, 1)
cvz_2 = c(2, 2, 1, 1, 1, 1, 1, 2, 1, 2, 3)

if (cluster == 1) {
  idx = which.max(sapply(1:nrow(res$clusters), function (i) {
    mclust::adjustedRandIndex(res$clusters[i, ], cvz_1)
  }))
} else {
  idx = which.max(sapply(1:nrow(res$clusters), function (i) {
    mclust::adjustedRandIndex(res$clusters[i, ], cvz_2)
  }))
}

# Most similar clustering according to ARI
clusters = get_clusters(res, idx)


# Plot dendrogram with most similar clustering
if (estimate_Theta) {
  plot_title = "Estimated Theta"
} else {
  plot_title = "Estimated Sigma"
}
plot(res_hclust, main = plot_title, xlab = "", sub = "", hang = -1, ylab = "")
rect.hclust(res_hclust, k = length(unique(clusters)))
