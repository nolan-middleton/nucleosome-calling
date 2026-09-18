##########################################################################
# This R script generates plots for the rotational categories to see if  #
# the Viterbi nucleosomes maintain rotational positining better than the #
# random nucleosome placements.                                          #
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
  "NucleosomePattern/composite_models/composite_ice",
  "NucleosomePattern/ReferenceRotationalCategories/all",
  "NucleosomePattern/ReferenceRotationalCategories/weiner"
)
dataset_names <- c(
  "Cellular vs. Room Temperature log2 Ratio",
  "Cellular vs. Ice log2 Ratio",
  "Cellular vs. Ice Absolute Difference",
  "Cellular vs. Ice Combined",
  "Brogaard et al. (All Nucleosomes)",
  "Weiner et al."
)
random_dirs <- c(
  paste0("RandomPattern/RotationalCategories/",unlist(lapply(0:9999,toString)))
)

D <- data.frame(
  dataset = character(0),
  L = numeric(0),
  freq = numeric(0)
)
for (i in 1:length(dirs)) {
  data <- read_json(paste0(dirs[i], "/rotational_categories.json"))
  
  Ls <- c()
  for (chrom in names(data)) {
    Ls <- c(Ls, unlist(lapply(data[[chrom]], length)))
  }
  tab <- as.data.frame(table(Ls))
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep(dataset_names[i], nrow(tab)),
      L = tab$Ls,
      freq = tab$Freq
    )
  )
}
D$L <- as.numeric(D$L)

for (i in 1:length(random_dirs)) {
  print(paste0("> ", toString(i), "..."))
  data <- read_json(paste0(random_dirs[i], "/rotational_categories.json"))
  
  Ls <- c()
  for (chrom in names(data)) {
    Ls <- c(Ls, unlist(lapply(data[[chrom]], length)))
  }
  tab <- as.data.frame(table(Ls))
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep(i, nrow(tab)),
      L = tab$Ls,
      freq = tab$Freq
    )
  )
}
D$L <- as.numeric(D$L)

means <- c()
U <- unique(D$dataset)
for (i in 1:length(U)) {
  print(paste0("> ", toString(i), "..."))
  subset <- D[D$dataset == U[i],]
  means <- c(means, sum(subset$L*subset$freq)/sum(subset$freq))
}

plotData <- data.frame(
  means = means,
  dataset = U
)

w = (max(plotData$means) - min(plotData$means))/251
ggplot(data = plotData) +
  geom_histogram(
    mapping = aes(x = means),
    bins = 250
  ) +
  annotate(
    geom = "line",
    x = rep(
      plotData$means[plotData$dataset %in% dataset_names] - w/2,
      each = 2
    ),
    y = rep(
      c(100, 10), sum(plotData$dataset %in% dataset_names)
    ),
    group = rep(plotData$dataset[plotData$dataset %in% dataset_names], each=2),
    arrow = arrow(length = unit(5, "pt")),
    linewidth = 1
  ) +
  annotate(
    geom = "text",
    x = plotData$means[plotData$dataset %in% dataset_names] - w/2,
    y = rep(100, sum(plotData$dataset %in% dataset_names)),
    label = plotData$dataset[plotData$dataset %in% dataset_names],
    angle = 90,
    hjust = 0,
    size = 3
  ) +
  scale_y_continuous(name = "Frequency") +
  scale_x_continuous(name = "Mean Length of Rotationally-Consistent Block")

ggsave(
  paste0("rotational_categories_histogram.png"),
  width = 2500,
  height = 1250,
  units = "px"
)
