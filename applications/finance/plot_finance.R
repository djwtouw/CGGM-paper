# --------------------------------------
# Author: Daniel Touw
#         Erasmus Universiteit Rotterdam
# --------------------------------------

rm(list = ls())

## load packages
library("dplyr")
library("ggplot2")
library("ggh4x")
library("scales")
library("stringr")
library("tidyr")
library("CGGMR")
library("gridExtra")


# Function to load out of sample results
load_oos_results <- function() {
  load("applications/finance/output/cggm-finance-5fold-oos-std-refit.RData")
  cggm = cvscore_outer[, "NLL"]
  load("applications/finance/output/taglasso-finance-5fold-oos-std.RData")
  tagl = cvscore_outer

  return(data.frame(cbind(cggm, tagl)))
}

# Function to load cluster results
load_cluster_results <- function()
{
  load("applications/finance/output/cggm-finance-5fold-std-refit.RData")
  cggm = get_clusters(fit_opt)
  load("applications/finance/output/taglasso-finance-5fold-std.RData")
  tagl = fit_opt$cluster
  sect = sector_data$Sector

  res = data.frame(cbind(cggm, tagl))
  res$sector = sect

  return(res)
}

## Load out of sample results
oos_res = load_oos_results()
oos_res$overall = "Out-of-sample errors"

# Calculate axis limits
lim = c(52.5, 187.5)

## create plot
oos_plot = ggplot(oos_res, aes(x=cggm, y=tagl)) +
  # Scatter plot
  geom_point(shape = 20, size = 2, col = "black", show.legend = FALSE) +
  # 45-degree line
  geom_abline(slope = 1) +
  # White background theme
  theme_bw() +
  # Labels for the axes
  labs(x = "Negative log-lik. CGGM", y = "Negative log-lik. TAGL") +
  # Fixed aspect ratio
  coord_fixed(xlim = lim, ylim = lim) +
  # X and y tick labels
  scale_x_continuous(breaks = c(60, 90, 120, 150, 180)) +
  scale_y_continuous(breaks = c(60, 90, 120, 150, 180)) +
  # Facet plot title
  facet_nested(~ overall) +
  theme(axis.title = element_text(size = 11),
        strip.text.x = element_text(size = 10))

oos_plot

## Load cluster results
clust_res = load_cluster_results()

# Matrix for heatmap, 11 rows for the number of sectors. Columns is the number
# of clusters
cggm_mat = matrix(0, nrow = 11, ncol = length(unique(clust_res$cggm)))
tagl_mat = matrix(0, nrow = 11, ncol = length(unique(clust_res$tagl)))

# Sectors
sectors = unique(clust_res$sector)

# Fill heatmap matrices
for (i in 1:nrow(clust_res)) {
  # Sector of the ith company
  sector_id = which(clust_res[i, ]$sector == sectors)

  cggm_mat[sector_id, clust_res[i, ]$cggm] =
    cggm_mat[sector_id, clust_res[i, ]$cggm] + 1

  tagl_mat[sector_id, clust_res[i, ]$tagl] =
    tagl_mat[sector_id, clust_res[i, ]$tagl] + 1
}

# Reorder clusters of CGGM result
cggm_mat = cggm_mat[, c(1, 3, 2)]

# Reverse order of the sectors
cggm_mat = cggm_mat[11:1, ]
tagl_mat = tagl_mat[11:1, ]

# Not pretty. Results in a dataframe that is in a format that works with facets
df_all = rbind(data.frame(Composition = "Cluster composition", Method = "CGGM",
                          i = seq_len(11), cggm_mat, check.names = FALSE) %>%
                 pivot_longer(cols = -(Composition:i), names_to = "j",
                              values_to = "Value") %>%
                 mutate(Composition = factor(Composition,
                                             levels = c("Cluster composition")),
                        Method = factor(Method, levels = c("CGGM", "TAGL")),
                        j = as.numeric(j)),
               data.frame(Composition = "Cluster composition", Method = "TAGL",
                          i = seq_len(11), tagl_mat, check.names = FALSE) %>%
                 pivot_longer(cols = -(Composition:i), names_to = "j",
                              values_to = "Value") %>%
                 mutate(Composition = factor(Composition,
                                             levels = c("Cluster composition")),
                        Method = factor(Method, levels = c("CGGM", "TAGL")),
                        j = as.numeric(j)))

## data frame for rectangle around matrix
df_rect <- data.frame(xmin = c(0.5, 0.5), xmax = c(3.5, 9.5),
                      ymin = c(0.5, 0.5), ymax = c(11.5, 11.5),
                      Width = "thin")
df_rect = transform(df_rect, Method = c("CGGM", "TAGL"))
df_rect$Composition = "Cluster composition"

# Plot
clust_plot = ggplot() +
  # Colored squares, gradient between white and blue
  geom_tile(mapping = aes(x = j, y = i, fill = Value),
            data = df_all, show.legend = FALSE) +
  scale_fill_gradient(name = "diagonal", low = "#FFFFFF",
                      high = "#8388fc") +
  # Numbers for within the squares, size needs to be really small to prevent
  # overlapping
  geom_text(data = df_all, mapping = aes(x = j, y = i, label = Value),
            size = 2.8) +
  # Rectangle around the plots
  geom_rect(mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                          linewidth = "thin"),
            data = df_rect,
            color = "black", fill = "transparent", show.legend = FALSE,
            inherit.aes = FALSE) +
  # X tick labels are the cluster IDs
  scale_x_continuous(breaks = 1:9, expand = c(0, 0)) +
  xlab("Cluster index") +
  # Y tick labels are the sector names
  scale_y_continuous(breaks = 11:1, labels = sectors,
                     expand = expansion(mult = c(0, 0.04))) +
  # Not sure
  scale_linewidth_manual(values = c(thin = 0.25, thick = 0.5)) +
  # Facets
  facet_nested(~ Composition + Method, space = "free_x", scales = "free_x") +
  # No fixed aspect ratio because facets have the same width, looks not good
  # with fixed aspect ratio
  #coord_fixed(clip = "off") +
  # White background
  theme_bw() +
  #force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "in")) +
  # Remove ticks and titles. Make the padding between x tick labels smaller
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 10))

clust_plot

## save plot to file
pdf(file = "applications/finance/figures/applications_finance.pdf", width = 7.5,
    height = 2.35)
grid.arrange(clust_plot, oos_plot, nrow = 1, ncol = 2)
dev.off()
