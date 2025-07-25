# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

## load packages
library("dplyr")
library("ggplot2")
library("ggh4x")
library("gridExtra")
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
  to_parse <- grep("[", text, fixed = TRUE)
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
designs <- c("block", "imperfect")

## parameter settings
n <- 120  # number of observations
p <- 15   # number of variables
B <- 3    # number of blocks
R <- 100  # number of simulation runs

## path and file names containing results
path <- "simulations/results"
file <- paste(path, "results_Sigma_%s_n=%d_p=%d_B=%d_R=%d.Rdata", sep = "/")

## nice labels for methods and evaluation criteria
design_labels <- c(block = "Exact",
                   imperfect = "Approximate")
method_labels <- c("CGGM_Theta_inv" = "paste(\"CGGM-\", hat(Theta)^-1)",
                   "CGGM_Sigma" = "paste(\"CGGM-\", hat(Sigma))",
                   S = "S")
criterion_labels <- c(Error = "Frobenius norm",
                      B_hat = "Number of clusters",
                      ARI = "ARI")

## define colors and line types for methods
colors <- hue_pal()(6)[c(4, 4, 5)]


##  load results
df_results <-
  mapply(function(design, B) {
    load(sprintf(file, design, n, p, B, R))
    data.frame(Main = "Block structure", Design = design, results)
  }, design = designs, B = B, SIMPLIFY = FALSE, USE.NAMES = FALSE) %>%
  bind_rows() %>%
  select(Main:Method, Error, B_hat, ARI) %>%
  pivot_longer(cols = Error:ARI, names_to = "Criterion",
               values_to = "Value") %>%
  mutate(Design = factor(design_labels[Design],
                         levels = design_labels),
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
  group_by(Main, Design, n, p, B, Method, Criterion) %>%
  summarize(Value = mean(Value),
            .groups = "drop")

# prepare data frame of aggregation performance
df_aggregation <- df_results %>%
  filter(Method != method_labels["S"],
         Criterion %in% criterion_labels[c("B_hat", "ARI")])

# prepare data frame of means over simulation runs
df_aggregation_means <- df_aggregation %>%
  group_by(Main, Design, n, p, B, Method, Criterion) %>%
  summarize(Mean = mean(Value),
            Nr_Runs = n(),
            .groups = "drop")

# prepare data frame of frequencies of individuals values over simulation runs
df_aggregation_frequencies <- df_aggregation %>%
  group_by(Main, Design, n, p, B, Method, Criterion, Value) %>%
  summarize(Frequency = n(),
            .groups = "drop") %>%
  left_join(df_aggregation_means,
            by = c("Main", "Design", "n", "p", "B", "Method", "Criterion")) %>%
  mutate(Frequency = Frequency / Nr_Runs)

# construct data frame of true values for reference lines
df_true <- rbind(data.frame(Main = "Block structure", Design = designs,
                            Criterion = "B_hat", Value = B),
                 data.frame(Main = "Block structure", Design = designs,
                            Criterion = "ARI", Value = 1)) %>%
  mutate(Design = factor(design_labels[Design],
                         levels = design_labels),
         Criterion = factor(criterion_labels[Criterion],
                            levels = criterion_labels))


## create plots

# Frobenius norm
plt_error <- ggplot() +
  # box plots
  geom_boxplot(aes(x = Method, y = Value, fill = Method),
               data = df_estimation, outlier.color = gray(0.5),
               show.legend = FALSE) +
  geom_point(aes(x = Method, y = Value), data = df_estimation_means,
             shape = 23, size = 1.5, fill = "black") +
  # facets and fine tuning
  facet_nested(. ~ Main + Design) +
  scale_x_discrete(labels = label_parse()) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.07))) +
  scale_fill_manual(values = colors) +
  scale_size(limits = c(0, 1), range = c(0, 5)) +
  labs(x = NULL, y = criterion_labels["Error"]) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        strip.text.x = element_text(size = 10))

# number of clusters
plt_B_hat <- ggplot() +
  # dot plots in other rows
  geom_hline(aes(yintercept = Value),
             data = filter(df_true,
                           Criterion == criterion_labels["B_hat"])) +
  geom_point(aes(x = Method, y = Value, size = Frequency),
             data = filter(df_aggregation_frequencies,
                           Criterion == criterion_labels["B_hat"]),
             shape = 21, color = "transparent", fill = "darkgray",
             show.legend = FALSE) +
  geom_point(aes(x = Method, y = Mean, fill = Method),
             data = filter(df_aggregation_means,
                           Criterion == criterion_labels["B_hat"]),
             shape = 23, size = 2.5, show.legend = FALSE) +
  # facets and fine tuning
  facet_nested(. ~ Main + Design) +
  scale_x_discrete(labels = label_parse()) +
  scale_y_continuous(breaks = seq(from = 1, to = p, by = 4)) +
  scale_fill_manual(values = colors) +
  scale_size(limits = c(0, 1), range = c(0, 5)) +
  labs(x = NULL, y = criterion_labels["B_hat"]) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        strip.text.x = element_text(size = 10))

# ARI
plt_ARI <- ggplot() +
  # dot plots in other rows
  geom_hline(aes(yintercept = Value),
             data = filter(df_true,
                           Criterion == criterion_labels["ARI"])) +
  geom_point(aes(x = Method, y = Value, size = Frequency),
             data = filter(df_aggregation_frequencies,
                           Criterion == criterion_labels["ARI"]),
             shape = 21, color = "transparent", fill = "darkgray",
             show.legend = FALSE) +
  geom_point(aes(x = Method, y = Mean, fill = Method),
             data = filter(df_aggregation_means,
                           Criterion == criterion_labels["ARI"]),
             shape = 23, size = 2.5, show.legend = FALSE) +
  # facets and fine tuning
  facet_nested(. ~ Main + Design) +
  scale_x_discrete(labels = label_parse()) +
  scale_y_continuous(breaks = seq(from = 0.25, to = 1, by = 0.25),
                     expand = expansion(mult = c(0.05, 0.13))) +
  scale_fill_manual(values = colors) +
  scale_size(limits = c(0, 1), range = c(0, 5)) +
  labs(x = NULL, y = criterion_labels["ARI"]) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        strip.text.x = element_text(size = 10))


## save plot to file
pdf(file = "simulations/figures/Sigma.pdf", width = 7.5, height = 2.15)
grid.arrange(plt_error, plt_B_hat, plt_ARI, nrow = 1, ncol = 3)
dev.off()
