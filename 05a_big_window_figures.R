##########################################################################
# This R script generates plots for the big window scores for the        #
# models.                                                                #
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

models <- c(
  "unfloored_cl",
  "cl_rt_log2",
  "cl_ice_log2",
  "cl_rt_diff",
  "cl_ice_diff",
  "cl_rt_abs_diff",
  "cl_ice_abs_diff"
)

datanames <- list(
  bquote("CPD counts"),
  bquote("log"[2]("CPD Ratio")),
  bquote("log"[2]("CPD Ratio")),
  bquote("difference in CPDs"),
  bquote("difference in CPDs"),
  bquote("|difference in CPDs|"),
  bquote("|difference in CPDs|")
)

dataset_names <- c(
  "Cellular",
  "Cellular / Room Temperature",
  "Cellular / Ice",
  "Cellular - Room Temperature",
  "Cellular - Ice",
  "|Cellular - Room Temperature|",
  "|Cellular - Ice|"
)

W <- 80

omega = 10.3

##### LOAD DATA #####
test_points <- read.delim(
  "evaluationPositions.tsv",
  sep = "\t",
  header = TRUE
)
test_points$answer <- test_points$answer == "True"
M <- nrow(test_points)

##### MAIN LOOP #####
for (i in 1:length(models)) {
  model <- models[i]
  print(paste0("> ", model, "..."))
  
  scores <- read.delim(
    paste0("NucleosomePattern/",model,"/models/bigWindowEvalBayes.tsv"),
    header = TRUE,
    sep = "\t"
  )
  colnames(scores) <- paste0("X", unlist(lapply(-80:80, toString)))
  
  results <- read.delim(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/bigWindowEvalResults.tsv"
    ),
    header = TRUE,
    sep = "\t"
  )
  
  H <- hist(results[test_points$answer,]$offset, breaks = 2*W+1)
  
  plot <- ggplot(data = results[test_points$answer,]) +
    geom_histogram(mapping = aes(x = offset), bins = 2*W+1) +
    scale_y_continuous(name = "Counts") +
    scale_x_continuous(
      name = "Nucleosome Placement",
      breaks = -floor(W/10):ceiling(W/10) * 10
    ) +
    ggtitle(paste0("Nucleosome Position Calling (",dataset_names[i],")"))
  
  for (x in omega*(-round(W/omega):round(W/omega))) {
    plot = plot + annotate(
      "line",
      x = c(x, x),
      y = c(0, max(H$counts)),
      linetype = "dashed"
    )
  }
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/bigWindowPlacementHist.png"
    ),
    plot,
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  plotData <- data.frame(score = numeric(0), offset = numeric(0))
  for (offset in -W:W) {
    plotData <- rbind(
      plotData,
      data.frame(
        score = mean(
          scores[[paste0("X", toString(offset))]][test_points$answer]
        ),
        offset = offset
      )
    )
  }
  ggplot(data = plotData) +
    geom_line(mapping = aes(x = offset, y = score)) +
    scale_x_continuous(
      "Position (Relative to Dyad)",
      breaks = -round(W/10):round(W/10) * 10
    ) +
    scale_y_continuous("Mean ln(Score)") +
    ggtitle(paste0("Average Score (", dataset_names[i], ")"))
  
  ggsave(
    paste0("NucleosomePattern/",model,"/models/bigWindowMeanScore.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
}