##########################################################################
# This R script generates plots for the distributions of rotational      #
# scores.                                                                #
##########################################################################

##### SETUP #####
require(ggplot2)
require(jsonlite)
require(stringr)
require(rstudioapi)

# Directories
setwd(dirname(getActiveDocumentContext()$path))

# Variables
colour_scale <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#004949", "#009999", "#22cf22", "#490092", "#006ddb",
  "#b66dff", "#ff6db6", "#920000", "#8f4e00", "#db6d00", "#ffdf4d",
  "#000000", "#252525", "#676767", "#ffffff", "#171723"
)

##### LOAD DATA #####
dirs <- c(
  "NucleosomePattern/cl_rt_log2",
  "NucleosomePattern/cl_ice_log2",
  "NucleosomePattern/cl_ice_abs_diff",
  "NucleosomePattern/composite_models/composite_ice"
)
dataset_names <- c(
  "Cellular vs. Room Temperature log2 Ratio",
  "Cellular vs. Ice log2 Ratio",
  "Cellular vs. Ice Absolute Difference",
  "Cellular vs. Ice Combined"
)

D <- data.frame(
  dataset = character(0),
  score = numeric(0)
)
for (i in 1:length(dirs)) {
  data <- read.delim(
    paste0(dirs[i], "/rotational_scores.tsv"),
    header = FALSE,
    sep = "\t"
  )
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep(dataset_names[i], nrow(data)),
      score = data$V3
    )
  )
}

##### PLOTS #####
ggplot() +
  geom_freqpoly(
    data = D,
    mapping = aes(
      color = dataset,
      x = score
    ),
    bins = 500
  ) +
  scale_x_continuous(name = "Rotational Score", limits = c(-100, 600)) +
  scale_y_continuous(name = "Count") +
  scale_color_discrete(name = "Dataset", type = colour_scale)

ggsave(
  "RotationalScoreDistribution.png",
  width = 4000,
  height = 1500,
  units = "px"
)
