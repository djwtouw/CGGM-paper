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
library("ggdendro")


load_hclusts <- function()
{
  load("applications/oecd/output/results_cluster_1.RData")
  clust1 = res_hclust
  load("applications/oecd/output/results_cluster_2.RData")
  clust2 = res_hclust

  res = list(clust1 = clust1, clust2 = clust2)

  return(res)
}

hclusts = load_hclusts()
hclusts[[1]]$height = hclusts[[1]]$height + 1
hclusts[[2]]$height = hclusts[[2]]$height + 1

# Convert to dendrogram data
dendro_data1 = dendro_data(hclusts[[1]], type = "rectangle")
dendro_data2 = dendro_data(hclusts[[2]], type = "rectangle")

# Add a label to differentiate the dendrograms
dendro_data1$segments$label = "Group 1"
dendro_data2$segments$label = "Group 2"

# Combine the data
combined_data = rbind(dendro_data1$segments, dendro_data2$segments)

# Extract and combine the label data
combined_labels = rbind(dendro_data1$labels, dendro_data2$labels)
colnames(combined_labels) = c("x", "y", "text")
combined_labels$label = combined_labels$label =
  rep(c("Group 1", "Group 2"), c(nrow(dendro_data1$labels),
                                 nrow(dendro_data2$labels)))

# Make labels pretty
combined_labels$text =
  sapply(combined_labels$text, function(x) gsub("\\.", " ", x))


# Create the ggplot object with facet_wrap
p = ggplot(combined_data) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = combined_labels, aes(x = x, y = y, label = text, hjust = 0),
            size = 3, nudge_x = 0.45) +
  coord_flip() +
  facet_wrap(~label, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#d9d9d9"),
        strip.text.x = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())


# Display the plot
print(p)

## save plot to file
pdf(file = "applications/oecd/figures/applications_oecd.pdf", width = 7.5,
    height = 2.35)
p
dev.off()
