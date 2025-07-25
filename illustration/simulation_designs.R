# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


# load packages
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggh4x")
library("ggnewscale")
library("scales")


# utility functions for generating block matrices -----

# function to generate diagonal blocks of covariance matrix
generate_diagonal_block <- function(d, theta) {
  if (length(theta) == 1) {
    mat <- matrix(theta, nrow = d, ncol = d)
    diag(mat) <- 1
  } else {
    mat <- diag(1, d)
    lower <- lower.tri(mat)
    upper <- upper.tri(mat)
    mat[lower] <- runif(d * (d-1) / 2, min = theta[1], max = theta[2])
    mat[upper] <- t(mat)[upper]
  }
  mat
}

# function to generate offdiagonal blocks of covariance matrix
generate_offdiagonal_block <- function(nrow, ncol, theta) {
  if (length(theta) == 1) {
    matrix(theta, nrow = nrow, ncol = ncol)
  } else {
    matrix(runif(nrow * ncol, min = theta[1], max = theta[2]),
           nrow = nrow, ncol = ncol)
  }
}

# function to generate blocks of zeros for covariance matrix
generate_zero_block <- function(nrow, ncol) {
  matrix(0, nrow = nrow, ncol = ncol)
}


# generate list of precision matrices -----

# simulation designs
designs_main <- c("baseline", "diagonal", "blockdiagonal")
designs_baseline <- c("random", "chain", "unbalanced", "unstructured")
designs_diagonal <- c("balanced", "unbalanced")
designs_blockdiagonal <- c("balanced", "unbalanced")
# nice labels
labels_main <- c("baseline" = "Baseline simulation designs",
                 "diagonal" = "Clustering structure on diagonal",
                 "blockdiagonal" = "Blockdiagonal structure")
labels_variation <- c("random" = "Random",
                      "chain" = "Chain",
                      "balanced" = "Balanced",
                      "unbalanced" = "Unbalanced",
                      "unstructured" = "Unstructured")
# construct data frame
designs <- data.frame(
  Main = rep.int(designs_main, times = c(length(designs_baseline),
                                         length(designs_diagonal),
                                         length(designs_blockdiagonal))),
  Variation = c(designs_baseline, designs_diagonal, designs_blockdiagonal)
) %>%
  mutate(Main = factor(labels_main[Main], levels = labels_main),
         Variation = factor(labels_variation[Variation],
                            levels = labels_variation))

# initializations
n_designs <- nrow(designs)
df_Theta <- df_lines <- vector("list", length = n_designs)

# baseline random design
# parameters
p <- 15           # number of variables
B <- 3            # number of blocks
theta_w <- 0.5    # (off-diagonal) elements within blocks
theta_b <- 0.25   # non-zero elements between blocks
blocks <- sort(rep(1:B, length.out = p))  # works if p is not a multiple of B
p_b <- tabulate(blocks)                   # number of variables per block
# select one non-zero edge in the aggregated graph (in simulation at random)
connected <- c(1, 3)
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
    else if (i > j) {
      # lower diagonal blocks
      if (i == connected[2] && j == connected[1]) {
        Theta[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j], theta_b)
      } else Theta[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j], 0)
    } else Theta[[j]][[i]] <- t(Theta[[i]][[j]])  # upper diagonal blocks
  }
}
# put precision matrix together
Theta <- do.call(cbind, lapply(Theta, function(column) do.call(rbind, column)))
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[1]] <- data.frame(designs[1, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[1]] <- data.frame(designs[1, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL) %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))

## baseline chain design
# parameters
p <- 15           # number of variables
B <- 3            # number of blocks
theta_w <- 0.5    # (off-diagonal) elements within blocks
theta_b <- 0.25   # non-zero elements between blocks
blocks <- sort(rep(1:B, length.out = p))  # works if p is not a multiple of B
p_b <- tabulate(blocks)                   # number of variables per block
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
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[2]] <- data.frame(designs[2, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[2]] <- data.frame(designs[2, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL) %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))

## baseline unbalanced design
# parameters
p_b <- c(3, 5, 7)  # number of variables per block
theta_w <- 0.5     # (off-diagonal) elements within blocks
theta_b <- 0.25    # non-zero elements between blocks
p <- sum(p_b)      # overall number of variables
B <- length(p_b)   # number of blocks
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
      Theta[[j]][[i]] <- generate_zero_block(p_b[i], p_b[j])
    } else Theta[[j]][[i]] <- t(Theta[[i]][[j]])  # upper diagonal blocks
  }
}
# put precision matrix together
Theta <- do.call(cbind, lapply(Theta, function(column) do.call(rbind, column)))
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[3]] <- data.frame(designs[3, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[3]] <- data.frame(designs[3, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL) %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))

## baseline unstructured design
# parameters
seed <- 20240604  # seed for the random number generator
p <- 15           # number of variables
pi <- 0.1         # connection probability of variables
theta <- 0.25     # non-zero elements of precision matrix
# initialize as diagonal matrix
set.seed(seed)
Theta <- diag(p)
# index matrices for lower and upper triangle
lower <- lower.tri(Theta)
upper <- upper.tri(Theta)
# in the simulations, set nonzero offdiagonal elements with probability pi, but
# here we ensure that pi*100% of possible edges are drawn for nice illustration
candidates <- which(lower)
connected <- sample(candidates, ceiling(pi * length(candidates)))
Theta[connected] <- theta
Theta[upper] <- t(Theta)[upper]
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[4]] <- data.frame(designs[4, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[4]] <- data.frame(designs[4, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL)

## diagonal balanced design
# parameters
p <- 15           # number of variables
B <- 3            # number of blocks
theta <- 0.5      # (off-diagonal) elements
blocks <- sort(rep(1:B, length.out = p))  # works if p is not a multiple of B
p_b <- tabulate(blocks)                   # number of variables per block
# initialize precision matrix as a constant matrix
Theta <- matrix(theta, nrow = p, ncol = p)
# put diagonal elements equal to block number
diag(Theta) <- blocks
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[5]] <- data.frame(designs[5, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[5]] <- data.frame(designs[5, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL) %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))

## diagonal unbalanced design
# parameters
p_b <- c(3, 5, 7)  # number of variables per block
theta <- 0.5       # (off-diagonal) elements
p <- sum(p_b)                    # overall number of variables
B <- length(p_b)                 # number of blocks
blocks <- rep(1:B, times = p_b)  # true block membership
# initialize precision matrix as a constant matrix
Theta <- matrix(theta, nrow = p, ncol = p)
# put diagonal elements equal to block number
diag(Theta) <- blocks
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[6]] <- data.frame(designs[6, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[6]] <- data.frame(designs[6, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL) %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))

## blockdiagonal balanced design
# parameters
seed <- 20240527  # seed for the random number generator
p <- 15           # number of variables
B <- 3            # number of blocks
theta_w <- 0.5    # (off-diagonal) elements within blocks
theta_b <- 0.25   # non-zero elements in off-diagonal blocks
pi <- 0.1         # connection probability of variables in off-diagonal blocks
blocks <- sort(rep(1:B, length.out = p))  # works if p is not a multiple of B
p_b <- tabulate(blocks)                   # number of variables per block
# use rejection sampling to generate precision matrix Theta: if the
# resulting Theta is not positive semi-definite, generate a new Theta
set.seed(seed)
# initialize negative eigenvalue
eigenvalues <- -Inf
# loop until we obtain a Theta with nonnegative eigenvalues
while(any(eigenvalues < 0)) {
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
      else Theta[[j]][[i]] <- generate_zero_block(p_b[i], p_b[j])
    }
  }
  # put precision matrix together
  Theta <- do.call(cbind, lapply(Theta, function(column) do.call(rbind, column)))
  # in the simulations, we set nonzero offdiagonal elements with probability pi
  # in the loops above, but here we ensure that pi*100% of possible edges in
  # the offdiagonal blocks are drawn for nice illustration
  lower <- lower.tri(Theta)
  upper <- upper.tri(Theta)
  p_sum <- cumsum(p_b)
  indices <- data.frame(expand.grid(i = seq_len(p), j = seq_len(p)),
                        index = seq_len(p^2)) %>%
    # ATTENTION: this is hardcoded for three blocks on the diagonal
    filter((i > p_sum[1] & j <= p_sum[1]) | (i > p_sum[2] & j <= p_sum[2]))
  candidates <- indices$index
  connected <- sample(candidates, ceiling(pi * length(candidates)))
  Theta[connected] <- theta_b
  Theta[upper] <- t(Theta)[upper]
  # compute eigenvalues
  eigenvalues <- eigen(Theta)$values
}
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[7]] <- data.frame(designs[7, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[7]] <- data.frame(designs[7, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL) %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))

## blockdiagonal unbalanced design
# parameters
seed <- 20240527   # seed for the random number generator
p_b <- c(3, 5, 7)  # number of variables per block
theta_w <- 0.5     # (off-diagonal) elements within blocks
theta_b <- 0.25    # non-zero elements in off-diagonal blocks
pi <- 0.1          # connection probability of variables in off-diagonal blocks
p <- sum(p_b)                    # overall number of variables
B <- length(p_b)                 # number of blocks
blocks <- rep(1:B, times = p_b)  # true block membership
# use rejection sampling to generate precision matrix Theta: if the
# resulting Theta is not positive semi-definite, generate a new Theta
set.seed(seed)
# initialize negative eigenvalue
eigenvalues <- -Inf
# loop until we obtain a Theta with nonnegative eigenvalues
while(any(eigenvalues < 0)) {
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
      else Theta[[j]][[i]] <- generate_zero_block(p_b[i], p_b[j])
    }
  }
  # put precision matrix together
  Theta <- do.call(cbind, lapply(Theta, function(column) do.call(rbind, column)))
  # in the simulations, we set nonzero offdiagonal elements with probability pi
  # in the loops above, but here we ensure that pi*100% of possible edges in
  # the offdiagonal blocks are drawn for nice illustration
  lower <- lower.tri(Theta)
  upper <- upper.tri(Theta)
  p_sum <- cumsum(p_b)
  indices <- data.frame(expand.grid(i = seq_len(p), j = seq_len(p)),
                        index = seq_len(p^2)) %>%
    # ATTENTION: this is hardcoded for three blocks on the diagonal
    filter((i > p_sum[1] & j <= p_sum[1]) | (i > p_sum[2] & j <= p_sum[2]))
  candidates <- indices$index
  connected <- sample(candidates, floor(pi * length(candidates)))
  Theta[connected] <- theta_b
  Theta[upper] <- t(Theta)[upper]
  # compute eigenvalues
  eigenvalues <- eigen(Theta)$values
}
# add column names
names(Theta) <- seq_len(p)
# convert to data frame and add to list
df_Theta[[8]] <- data.frame(designs[8, ], i = seq_len(p), Theta,
                            check.names = FALSE, row.names = NULL)
# data frame for horizontal and vertical lines
df_lines[[8]] <- data.frame(designs[8, ],
                            Intercept = 1:(p-1) + 0.5,
                            Start = 0.5,
                            End = p + 0.5,
                            Width = "thin",
                            row.names = NULL) %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))


## prepare data frames for plotting
df_Theta <- bind_rows(df_Theta) %>%
  pivot_longer(cols = -(Main:i), names_to = "j", values_to = "Value") %>%
  mutate(#Main = factor(Main, levels = labels_main),
    #Variation = factor(Variation, levels = labels_variation),
    j = as.numeric(j))
df_lines <- bind_rows(df_lines)

## extract separate data frames for diagonal and offdiagonal
## (to have different color scales)
df_diagonal <- df_Theta %>% filter(i == j)
df_offdiagonal <- df_Theta %>% filter(i != j)

## data frame for rectangle around matrix
df_rect <- data.frame(Min = 0.5, Max = p + 0.5,
                      Width = "thin")

## create plot
p_heatmap <-
  ggplot() +
  # diagonal elements with their own color scale
  geom_tile(mapping = aes(x = j, y = i, fill = Value),
            data = df_diagonal, show.legend = FALSE) +
  scale_fill_gradient(name = "diagonal", low = "#F5BF30",
                      high = muted("#F5BF30")) +
  # off-diagonal elements with their own color scale
  new_scale("fill") +
  geom_tile(mapping = aes(x = j, y = i, fill = Value),
            data = df_offdiagonal, show.legend = FALSE) +
  scale_fill_gradient(name = "off-diagonal", limits = c(NA, 0.6),
                      low = "white", high = "#3596E1") +
  # rest of the plot
  geom_segment(mapping = aes(x = Intercept, y = Start, yend = End,
                             linewidth = Width),
               data = df_lines, color = "black", show.legend = FALSE) +
  geom_hline(mapping = aes(yintercept = Intercept, linewidth = Width),
             data = df_lines, color = "black", show.legend = FALSE) +
  geom_rect(mapping = aes(xmin = Min, xmax = Max, ymin = Min, ymax = Max,
                          linewidth = "thin"),
            data = df_rect, color = "black", fill = "transparent",
            show.legend = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)),
                     transform = "reverse") +
  scale_linewidth_manual(values = c(thin = 0.25, thick = 0.5)) +
  facet_wrap2(~ Main + Variation, #labeller = label_parsed,
              nrow = 2, strip = strip_nested(bleed = FALSE)) +
  coord_fixed(clip = "off") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 10))

## save plot to file
pdf(file = "illustration/simulation_designs.pdf",
    width = 6, height = 4)
print(p_heatmap)
dev.off()
