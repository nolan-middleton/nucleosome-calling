##########################################################################
# This R script generates plots for the distributions of linker lengths  #
# in the models.                                                         #
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
  "#E69F00", "#56B4E9", "#009E73", "#000000", "#0072B2", "#D55E00",
  "#CC79A7", "#004949", "#009999", "#22cf22", "#490092", "#006ddb",
  "#b66dff", "#ff6db6", "#920000", "#8f4e00", "#db6d00", "#ffdf4d",
  "#000000", "#252525", "#676767", "#ffffff", "#171723"
)

##### LOAD DATA #####
dirs <- c(
  "NucleosomePattern/cl_rt_log2",
  "NucleosomePattern/cl_ice_log2",
  "NucleosomePattern/cl_ice_abs_diff",
  "NucleosomePattern/composite_models/composite_ice",
  "NucleosomePattern/composite_models/composite_ice/GreedyLinkerLens",
  "NucleosomePattern/ReferenceLinkerLengths/all",
  "NucleosomePattern/ReferenceLinkerLengths/weiner"
)
dataset_names <- c(
  "Cellular vs. Room Temperature log2 Ratio",
  "Cellular vs. Ice log2 Ratio",
  "Cellular vs. Ice Absolute Difference",
  "Cellular vs. Ice Combined",
  "Cellular vs. Ice Combined (Greedy)",
  "Brogaard et al. (All Nucleosomes)",
  "Weiner et al."
)

D <- data.frame(
  dataset = character(0),
  length = numeric(0)
)
for (i in 1:length(dirs)) {
  data <- read.delim(
    paste0(dirs[i], "/linker_lens.tsv"),
    header = FALSE,
    sep = "\t"
  )
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep(dataset_names[i], length(data$V1)),
      length = data$V1
    )
  )
}

##### PLOTS #####

points = seq(-146, 146, 10.1)
annotationData <- data.frame(
  x = rep(points, 2),
  y = rep(c(-Inf, Inf), each = length(points)),
  group = rep(1:length(points), 2)
)

ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linetype = "dotted"
  ) +
  geom_freqpoly(
    data = D,
    mapping = aes(
      color = dataset,
      x = length
    ),
    bins = 147 + 73
  ) +
  scale_x_continuous(
    name = "Linker Length",
    limits = c(-73, 146),
    breaks = points[points >= -73]
  ) +
  scale_y_continuous(name = "Count") +
  scale_color_discrete(name = "Dataset", type = colour_scale)

ggsave(
  "LinkerLengthsDistribution.png",
  width = 4000,
  height = 1500,
  units = "px"
)
