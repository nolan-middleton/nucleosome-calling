##########################################################################
# This R script generates plots for the rotational classes to see if the #
# Viterbi nucleosomes maintain rotational positining better in genes or  #
# outside of genes as compared to the random nucleosome positions.       #
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
  length = numeric(0),
  in_gene = logical(0),
  threshold = numeric(0)
)
for (i in 1:length(dirs)) {
  data <- read_json(paste0(dirs[i], "/rotational_classes.json"))
  
  Ls <- c()
  Gs <- c()
  for (chrom in names(data)) {
    unlisted <- unlist(data[[chrom]])
    Ls <- c(Ls, unlisted[seq(1,length(unlisted),2)])
    Gs <- c(Gs, unlisted[seq(2,length(unlisted),2)])
  }
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep(dataset_names[i], length(Ls)),
      length = Ls,
      in_gene = Gs,
      threshold = rep(1.5, length(Ls))
    )
  )
  
  data <- read_json(paste0(dirs[i], "/rotational_classes_lenient.json"))
  
  Ls <- c()
  Gs <- c()
  for (chrom in names(data)) {
    unlisted <- unlist(data[[chrom]])
    Ls <- c(Ls, unlisted[seq(1,length(unlisted),2)])
    Gs <- c(Gs, unlisted[seq(2,length(unlisted),2)])
  }
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep(dataset_names[i], length(Ls)),
      length = Ls,
      in_gene = Gs,
      threshold = rep(2.5, length(Ls))
    )
  )
}

for (i in 1:length(random_dirs)) {
  print(paste0("> ", toString(i), "..."))
  data<-read_json(paste0(random_dirs[i],"/rotational_classes.json"))
  
  Ls <- c()
  Gs <- c()
  for (chrom in names(data)) {
    unlisted <- unlist(data[[chrom]])
    Ls <- c(Ls, unlisted[seq(1,length(unlisted),2)])
    Gs <- c(Gs, unlisted[seq(2,length(unlisted),2)])
  }
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep("Random", length(Ls)),
      length = Ls,
      in_gene = Gs,
      threshold = rep(1.5, length(Ls))
    )
  )
  
  data<-read_json(paste0(random_dirs[i],"/rotational_classes_lenient.json"))
  
  Ls <- c()
  Gs <- c()
  for (chrom in names(data)) {
    unlisted <- unlist(data[[chrom]])
    Ls <- c(Ls, unlisted[seq(1,length(unlisted),2)])
    Gs <- c(Gs, unlisted[seq(2,length(unlisted),2)])
  }
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep("Random", length(Ls)),
      length = Ls,
      in_gene = Gs,
      threshold = rep(2.5, length(Ls))
    )
  )
}
rm(Ls, Gs, i, data, chrom)

##### PLOTS #####

plotData <- data.frame(
  dataset = character(0),
  in_gene = logical(0),
  mean = numeric(0),
  sd = numeric(0),
  n = numeric(0),
  threshold = numeric(0)
)

for (dataset in unique(D$dataset)) {
  for (threshold in unique(D$threshold)) {
    subset = D[(D$dataset == dataset) & (D$threshold == threshold),]
    plotData <- rbind(
      plotData,
      data.frame(
        dataset = rep(dataset, 2),
        in_gene = c(TRUE, FALSE),
        mean = c(
          mean(subset$length[subset$in_gene == 1]),
          mean(subset$length[subset$in_gene == 0])
        ),
        sd = c(
          sd(subset$length[subset$in_gene == 1]),
          sd(subset$length[subset$in_gene == 0])
        ),
        n = c(
          nrow(subset[subset$in_gene == 1,]),
          nrow(subset[subset$in_gene == 0,])
        ),
        threshold = rep(threshold, 2)
      )
    )
  }
}

ggplot(data = plotData[plotData$threshold == 1.5,]) +
  geom_col(
    mapping = aes(
      x = dataset,
      y = mean,
      fill = in_gene
    ),
    position = position_dodge()
  ) +
  geom_errorbar(
    mapping = aes(
      x = dataset,
      ymin = mean - sd,
      ymax = mean + sd,
      group = in_gene
    ),
    position = position_dodge(),
    width = 0.9
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Dataset") +
  scale_y_continuous(name = "Mean Length of\nRotationally-Consistent Block") +
  scale_fill_discrete(name = "In Gene", type = colour_scale)

ggsave(
  "rotational_classes_bargraph.png",
  width = 2000,
  height = 1500,
  units = "px"
)

ggplot(data = plotData[plotData$threshold == 2.5,]) +
  geom_col(
    mapping = aes(
      x = dataset,
      y = mean,
      fill = in_gene
    ),
    position = position_dodge()
  ) +
  geom_errorbar(
    mapping = aes(
      x = dataset,
      ymin = mean - sd,
      ymax = mean + sd,
      group = in_gene
    ),
    position = position_dodge(),
    width = 0.9
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "Dataset") +
  scale_y_continuous(name = "Mean Length of\nRotationally-Consistent Block") +
  scale_fill_discrete(name = "In Gene", type = colour_scale)

ggsave(
  "rotational_classes_lenient_bargraph.png",
  width = 2000,
  height = 1500,
  units = "px"
)
