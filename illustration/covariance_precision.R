# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


## load packages
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggh4x")
library("ggnewscale")
library("scales")

## control parameters
# seed <- 20240531                  # seed for the random number generator
seed <- 20240604                  # seed for the random number generator
p <- 15                           # number of variables
B <- 3                            # number of blocks
sigma_w_block <- 0.5              # (off-diagonal) elements within blocks
sigma_b_block <- 0.25             # non-zero elements between blocks
sigma_w_imperfect <- c(0.4, 0.6)  # (off-diagonal) elements within blocks
sigma_b_imperfect <- c(0.2, 0.3)  # non-zero elements between blocks

## additional parameters
blocks <- sort(rep(1:B, length.out = p))  # works if p is not a multiple of B
p_b <- tabulate(blocks)                   # number of variables per block


## utility functions for generating block matrices

# function to generate diagonal blocks of covariance matrix
generate_diagonal_block <- function(d, sigma) {
  if (length(sigma) == 1) {
    mat <- matrix(sigma, nrow = d, ncol = d)
    diag(mat) <- 1
  } else {
    mat <- diag(1, d)
    lower <- lower.tri(mat)
    upper <- upper.tri(mat)
    mat[lower] <- runif(d * (d-1) / 2, min = sigma[1], max = sigma[2])
    mat[upper] <- t(mat)[upper]
  }
  mat
}

# function to generate offdiagonal blocks of covariance matrix
generate_offdiagonal_block <- function(nrow, ncol, sigma) {
  if (length(sigma) == 1) {
    matrix(sigma, nrow = nrow, ncol = ncol)
  } else {
    matrix(runif(nrow * ncol, min = sigma[1], max = sigma[2]),
           nrow = nrow, ncol = ncol)
  }
}

# function to generate blocks of zeros for covariance matrix
generate_zero_block <- function(nrow, ncol) {
  matrix(0, nrow = nrow, ncol = ncol)
}


## construct covariance matrix
set.seed(seed)
# initialize covariance matrix as a list of lists:
# (outer list corresponds to columns, inner list to rows of the given column)
Sigma_block <- Sigma_imperfect <-
  replicate(B, replicate(B, NULL, simplify = FALSE), simplify = FALSE)
# loop over indices of blocks in the final covariance matrix, and generate
# those blocks following the R convention of building a matrix by column
# (i is the index of the row, j is the index of the column)
for (j in seq_len(B)) {
  for (i in seq_len(B)) {
    if (i == j) {
      Sigma_block[[j]][[i]] <- generate_diagonal_block(p_b[i], sigma_w_block)
      Sigma_imperfect[[j]][[i]] <- generate_diagonal_block(p_b[i],
                                                           sigma_w_imperfect)
    } else if (i == j+1) {
      Sigma_block[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j],
                                                          sigma_b_block)
      Sigma_imperfect[[j]][[i]] <- generate_offdiagonal_block(p_b[i], p_b[j],
                                                              sigma_b_imperfect)
    } else if (i > j+1) {
      Sigma_block[[j]][[i]] <-     generate_zero_block(p_b[i], p_b[j])
      Sigma_imperfect[[j]][[i]] <- generate_zero_block(p_b[i], p_b[j])
    } else {
      # upper diagonal blocks
      Sigma_block[[j]][[i]] <- t(Sigma_block[[i]][[j]])
      Sigma_imperfect[[j]][[i]] <- t(Sigma_imperfect[[i]][[j]])
    }
  }
}
# put covariance matrix together
Sigma_block <- do.call(cbind,
                       lapply(Sigma_block,
                              function(column) do.call(rbind, column)))
Sigma_imperfect <- do.call(cbind,
                           lapply(Sigma_imperfect,
                                  function(column) do.call(rbind, column)))

## compute precision matrix
Theta_block <- solve(Sigma_block)
Theta_imperfect <- solve(Sigma_imperfect)

## add column names to matrices
colnames(Sigma_block) <- colnames(Theta_block) <- seq_len(p)
colnames(Sigma_imperfect) <- colnames(Theta_imperfect) <- seq_len(p)


## nice labels for plot
design_labels <- c("Exact~block~structure",
                   "Approximate~block~structure")
matrix_labels <- c("plain(Covariance~matrix)~bold(Sigma)",
                   "plain(Precision~matrix)~bold(Theta)")

## convert matrices to data frames for plotting
df_all <- bind_rows(data.frame(Design = design_labels[1],
                               Matrix = matrix_labels[1],
                               i = seq_len(p), Sigma_block,
                               check.names = FALSE),
                    data.frame(Design = design_labels[1],
                               Matrix = matrix_labels[2],
                               i = seq_len(p), Theta_block,
                               check.names = FALSE),
                    data.frame(Design = design_labels[2],
                               Matrix = matrix_labels[1],
                               i = seq_len(p), Sigma_imperfect,
                               check.names = FALSE),
                    data.frame(Design = design_labels[2],
                               Matrix = matrix_labels[2],
                               i = seq_len(p), Theta_imperfect,
                               check.names = FALSE)) %>%
  pivot_longer(cols = -(Design:i), names_to = "j", values_to = "Value") %>%
  mutate(Design = factor(Design, levels = design_labels),
         Matrix = factor(Matrix, levels = matrix_labels),
         j = as.numeric(j))

## extract separate data frames for diagonal and offdiagonal
## (to have different color scales)
df_diagonal <- df_all %>% filter(i == j)
df_offdiagonal <- df_all %>% filter(i != j)


## data frame for horizontal and vertical lines
df_lines <- data.frame(Intercept = 1:(p-1) + 0.5,
                       Start = 0.5,
                       End = p + 0.5,
                       Width = "thin") %>%
  mutate(Width = replace(Width, Intercept %in% (cumsum(p_b[-B])+0.5), "thick"))

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
  scale_fill_gradient2(name = "off-diagonal", limits = c(NA, 0.6),
                       low = "#A83828", mid = "white", high = "#3596E1") +
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
  facet_nested(~ Design + Matrix, labeller = label_parsed) +
  coord_fixed(clip = "off") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

## save plot to file
pdf(file = "illustration/covariance_precision.pdf",
    width = 6, height = 2)
print(p_heatmap)
dev.off()
