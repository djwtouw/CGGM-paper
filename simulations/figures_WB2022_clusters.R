# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## load packages
library("dplyr")
library("ggplot2")
library("ggh4x")
library("scales")
library("stringr")
library("tidyr")

## functions to parse only method names that contain a superscript (for x-axis)
# function to be passed to scale_x_discrete()
label_custom <- function() {
  function(text) parse_custom(as.character(text))
}
# underlying function that does the work (parsing)
parse_custom <- function(text) {
  # find which labels contain a superscript
  to_parse <- grep("^", text, fixed = TRUE)
  # return labels without superscripts unparsed, otherwise loop over labels to
  # be parsed and return a vector of expressions
  if (length(to_parse) == 0) text
  else {
    out <- vector("expression", length(text))
    out[-to_parse] <- text[-to_parse]
    for (i in to_parse) {
      expr <- parse(text = text[[i]])
      out[[i]] <- if (length(expr) == 0) NA else expr[[1]]
    }
    out
  }
}


## different designs to be analyzed
design <- "chain"

## parameter settings
n <- 240             # number of observations
p <- 30              # number of variables
B <- c(3, 5, 6, 10)  # number of blocks
R <- 100             # number of simulation runs

## path and file names containing results
path <- "simulations/results"
file <- paste(path, "results_WB2022_%s_n=%d_p=%d_B=%d_R=%d.Rdata", sep = "/")

## nice labels for methods and evaluation criteria
method_labels <- c("CGGM-refit" = "CGGM-refit",
                   "TAGL-ideal" = "TAGL-ideal",
                   "TAGL-realistic" = "TAGL-realistic",
                   "TAGL-misspecified" = "TAGL-misspecified",
                   ComGGL = "ComGGL",
                   Sinv = "S^-1")
criterion_labels <- c(Error = "Frobenius norm",
                      B_hat = "No. of clusters",
                      ARI = "ARI")
B_label <- "K = %d"
B_levels <- sprintf(B_label, B)


## define colors and line types for methods
colors <- hue_pal()(6)[c(4, 6, 6, 6, 1, 5)]


##  load results
df_results <-
  lapply(B, function(B) {
    load(sprintf(file, design, n, p, B, R))
    data.frame(Design = design, results)
  }) %>%
  bind_rows() %>%
  filter(Method %in% names(method_labels)) %>%
  select(Design:Method, Error, B_hat, ARI) %>%
  pivot_longer(cols = Error:ARI, names_to = "Criterion",
               values_to = "Value") %>%
  mutate(Design = factor(Design, levels = design),
         B = factor(sprintf(B_label, B), levels = B_levels),
         Method = factor(method_labels[Method],
                         levels = method_labels),
         Criterion = factor(criterion_labels[Criterion],
                            levels = criterion_labels))


## prepare data frames for plot

# prepare data frame of estimation performance
df_estimation <- df_results %>%
  filter(Criterion == criterion_labels["Error"])

# prepare data frame of means over simulation runs
df_estimation_means <- df_estimation %>%
  group_by(Design, n, p, B, Method, Criterion) %>%
  summarize(Value = mean(Value),
            .groups = "drop")

# prepare data frame of aggregation performance
df_aggregation <- df_results %>%
  filter(Method != method_labels["Sinv"],
         Criterion %in% criterion_labels[c("B_hat", "ARI")])

# prepare data frame of means over simulation runs
df_aggregation_means <- df_aggregation %>%
  group_by(Design, n, p, B, Method, Criterion) %>%
  summarize(Mean = mean(Value),
            Nr_Runs = n(),
            .groups = "drop")

# prepare data frame of frequencies of individuals values over simulation runs
df_aggregation_frequencies <- df_aggregation %>%
  group_by(Design, n, p, B, Method, Criterion, Value) %>%
  summarize(Frequency = n(),
            .groups = "drop") %>%
  left_join(df_aggregation_means,
            by = c("Design", "n", "p", "B", "Method", "Criterion")) %>%
  mutate(Frequency = Frequency / Nr_Runs)

# construct data frame of true values for reference lines
df_true <- rbind(data.frame(Design = design, B, Criterion = "B_hat", Value = B),
                 data.frame(Design = design, B, Criterion = "ARI", Value = 1)) %>%
  mutate(Design = factor(Design, levels = design),
         B = factor(sprintf("K = %d", B), levels = B_levels),
         Criterion = factor(criterion_labels[Criterion],
                            levels = criterion_labels))

## create plot
plt <- ggplot() +
  # box plots in top row
  geom_boxplot(aes(x = Method, y = Value, fill = Method),
               data = df_estimation, outlier.color = gray(0.5),
               show.legend = FALSE) +
  geom_point(aes(x = Method, y = Value), data = df_estimation_means,
             shape = 23, size = 1.5, fill = "black") +
  # dot plots in other rows
  geom_hline(aes(yintercept = Value), data = df_true) +
  geom_point(aes(x = Method, y = Value, size = Frequency),
             data = df_aggregation_frequencies, shape = 21,
             color = "transparent", fill = "darkgray",
             show.legend = FALSE) +
  geom_point(aes(x = Method, y = Mean, fill = Method),
             data = df_aggregation_means, shape = 23,
             size = 2.5, show.legend = FALSE) +
  # facets and fine tuning
  facet_grid(Criterion ~ B, scales = "free_y") +
  scale_x_discrete(labels = label_custom()) +
  facetted_pos_scales(
    y = list(NULL,
             scale_y_continuous(breaks = seq(from = 1, to = p, by = 6),
                                expand = expansion(mult = c(0.12, 0.05))),
             scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                                expand = expansion(mult = 0.12)))
  ) +
  scale_fill_manual(values = colors) +
  scale_size(limits = c(0, 1), range = c(0, 5)) +
  labs(x = NULL, y = "Evaluation criterion") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        strip.text.x = element_text(size = 10))

## save plot to file
pdf(file = "simulations/figures/WB2022_clusters.pdf",
    width = 6.45, height = 4.3)
print(plt)
dev.off()
