to_hclust <- function(res, height_type = "index")
{
  if (!(height_type %in% c("lambda", "index"))) {
    stop(paste0("height_type should be one of lambda or index. This determines",
                "the height at which clusters are fused in the dendrogram"))
  }

  # Get relevant solution indices
  indices = which(c(-1, diff(res$cluster_counts)) < 0)

  # Get relevant parts of the output
  clusters = res$clusters[indices, ]
  lambdas = res$lambdas[indices]

  # Refit the result, just in case, can be done much more efficient
  # Initialize result
  result = list()

  # Initialize merge matrix
  n = ncol(clusters)
  merge = matrix(0, nrow = n - 1, ncol = 2)

  # Initialize height vector
  height = rep(0, n - 1)

  # Introduce some helper variables
  cluster_id = -c(1:n)
  cluster_check = rep(1, n)
  n_fuse = 1

  # Re-used code from an older project, no idea how and why it works, but it
  # does
  for (i in 1:nrow(clusters)) {
    for (j in 1:(n - 1)) {
      if (cluster_check[j] > 0) {
        for (k in (j + 1):n) {
          if (cluster_check[k] > 0) {
            if (clusters[i, j] == clusters[i, k]) {
              cluster_check[k] = -1
              merge[n_fuse, 1] = cluster_id[j]
              merge[n_fuse, 2] = cluster_id[k]
              cluster_id[j] = n_fuse
              cluster_id[k] = n_fuse
              if (height_type == "lambda") {
                height[n_fuse] = lambdas[i]
              } else {
                height[n_fuse] = i
              }
              n_fuse = n_fuse + 1
            }
          }
        }
      }
    }
  }

  # Order vector
  order = list()

  order[[1]] = abs(merge[1, ])
  for (i in 2:nrow(merge)) {
    if (sum(merge[i, ] < 0) == 2) {
      order[[i]] = abs(merge[i, ])
    } else if (sum(merge[i, ] < 0) == 1) {
      idx_pos = c(1, 2)[(merge[i, ] > 0)]
      idx_neg = c(1, 2)[(merge[i, ] < 0)]

      order[[i]] = c(order[[merge[i, idx_pos]]], abs(merge[i, idx_neg]))

      # Reduce the size of irrelevant list indices
      order[[merge[i, idx_pos]]] = c(0)
    } else if (sum(merge[i, ] == c(0, 0)) == 2) {
      final_order = c()
      for (i2 in 1:length(order)) {
        if (length(order[[i2]]) > 1) {
          final_order = c(final_order, order[[i2]])
        }
      }

      if (length(final_order) < n) {
        ordered = final_order[order(final_order)]
        correction = 0
        missings = c()

        for (i3 in 1:length(ordered)) {
          if (ordered[i3] != (i3 + correction)) {
            missings = c(missings, i3)
            correction = correction + 1
          }
        }

        final_order = c(final_order, missings)
      }

      order[[i]] = final_order

      break
    } else {
      order[[i]] = c(order[[merge[i, 1]]], order[[merge[i, 2]]])

      # Reduce the size of irrelevant list indices
      order[[merge[i, 1]]] = c(0)
      order[[merge[i, 2]]] = c(0)
    }
  }

  # Fill result
  result$order = order[[length(order)]]
  result$merge = merge
  result$height = height

  # Labels
  result$labels = colnames(res$clusters)

  # Method
  result$method = "CGGM"

  # Call
  result$call = list()
  result$call[[1]] = "CGGM"
  result$call$d[[1]] = "dist"
  result$call$d[[2]] = "X"

  # Dist.method
  result$dist.method = "CGGM"

  # Set class
  class(result) = "hclust"

  return(result)
}
