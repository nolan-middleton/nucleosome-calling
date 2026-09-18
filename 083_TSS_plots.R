##########################################################################
# This R script generates plots for the TSS evaluations of the various   #
# nucleosome placements.                                                 #
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
dataset_names <- c(
  "Cellular vs. Room Temperature log2 Ratio",
  "Cellular vs. Ice log2 Ratio",
  "Cellular vs. Ice Absolute Difference",
  "Cellular vs. Ice Combined"
)
dirs <- c(
  "NucleosomePattern/cl_rt_log2/TSS_evals",
  "NucleosomePattern/cl_ice_log2/TSS_evals",
  "NucleosomePattern/cl_ice_abs_diff/TSS_evals",
  "NucleosomePattern/composite_models/composite_ice/TSS_evals"
)

greedy <- list()
viterbi <- list()
scores <- list()
median_scores <- list()
probs <- list()
for (i in 1:length(dirs)) {
  greedy[[i]] <- read.delim(
    paste0(dirs[i], "/greedy_frequencies.tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(greedy[[i]]) <- c("position", "freq")
  viterbi[[i]] <- read.delim(
    paste0(dirs[i], "/viterbi_frequencies.tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(viterbi[[i]]) <- c("position", "freq")
  scores[[i]] <- read.delim(
    paste0(dirs[i], "/mean_score.tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[i]]) <- c("position", "mean_score")
  median_scores[[i]] <- read.delim(
    paste0(dirs[i], "/median_score.tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(median_scores[[i]]) <- c("position", "median_score")
  probs[[i]] <- read.delim(
    paste0(dirs[i], "/mean_prob.tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(probs[[i]]) <- c("position", "mean_prob")
}

reference_datasets <- c("strong", "decent", "all", "Weiner")
reference_names <- c(
  "Brogaard et al. (Score > 5)",
  "Brogaard et al. (Score > 1)",
  "Brogaard et al. (Complete)",
  "Weiner et al."
)
reference <- list()
for (i in 1:length(reference_datasets)) {
  reference[[i]] <- read.delim(
    paste0(
      "NucleosomePattern/TSS_reference_datasets/",
      reference_datasets[i],
      "_frequencies.tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(reference[[i]]) <- c("position", "freq")
}

##### MAIN LOOP #####
for (i in 1:length(dirs)) {
  print(paste0("> ", dataset_names[i], "..."))
  for (j in 1:length(reference)) {
    print(paste0(">> ", reference_names[j], "..."))
    R = max(reference[[j]]$freq) / sum(reference[[j]]$freq)
    
    M = max(greedy[[i]]$freq) / sum(greedy[[i]]$freq)
    m = min(greedy[[i]]$freq) / sum(greedy[[i]]$freq)
    plotData <- greedy[[i]]
    plotData$freq = (plotData$freq/sum(plotData$freq) - m) * R/(M-m)
    ggplot() +
      geom_rect(
        data = reference[[j]],
        mapping = aes(
          xmin = position - 0.5,
          xmax = position + 0.5,
          ymin = 0,
          ymax = freq / sum(freq)
        )
      ) +
      geom_line(
        data = plotData,
        mapping = aes(x = position, y = freq),
        color = colour_scale[1],
        linewidth = 0.25
      ) +
      annotate(
        geom = "line",
        x = c(0,0),
        y = c(-Inf, Inf),
        linetype = "dashed"
      ) +
      scale_y_continuous(
        name = "Proportion (Reference)",
        limits = c(0,max(reference[[j]]$freq / sum(reference[[j]]$freq))),
        sec.axis = sec_axis(
          transform = ~.*(M-m)/R + m,
          name = "Proportion"
        )
      ) +
      scale_x_continuous(
        name = "Position (Relative to TSS)",
        breaks = (-10:13)*50
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(
        paste0(
          dataset_names[i],
          " Greedy Dyad Placements\nCompared to ",
          reference_names[j]
        )
      )
    
    ggsave(
      paste0(dirs[i], "/", reference_datasets[j], "_greedy_plot.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    M = max(viterbi[[i]]$freq) / sum(viterbi[[i]]$freq)
    m = min(viterbi[[i]]$freq) / sum(viterbi[[i]]$freq)
    plotData <- viterbi[[i]]
    plotData$freq = (plotData$freq/sum(plotData$freq) - m) * R/(M-m)
    ggplot() +
      geom_rect(
        data = reference[[j]],
        mapping = aes(
          xmin = position - 0.5,
          xmax = position + 0.5,
          ymin = 0,
          ymax = freq / sum(freq)
        )
      ) +
      geom_line(
        data = plotData,
        mapping = aes(x = position, y = freq),
        color = colour_scale[1],
        linewidth = 0.25
      ) +
      annotate(
        geom = "line",
        x = c(0,0),
        y = c(-Inf, Inf),
        linetype = "dashed"
      ) +
      scale_y_continuous(
        name = "Proportion (Reference)",
        limits = c(0,max(reference[[j]]$freq / sum(reference[[j]]$freq))),
        sec.axis = sec_axis(
          transform = ~.*(M-m)/R + m,
          name = "Proportion"
        )
      ) +
      scale_x_continuous(
        name = "Position (Relative to TSS)",
        breaks = (-10:13)*50
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(
        paste0(
          dataset_names[i],
          " Dynamic Dyad Placements\nCompared to ",
          reference_names[j]
        )
      )
    
    ggsave(
      paste0(dirs[i], "/", reference_datasets[j], "_viterbi_plot.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    M = max(scores[[i]]$mean_score)
    m = min(scores[[i]]$mean_score)
    plotData <- scores[[i]]
    plotData$mean_score = (plotData$mean_score - m) * R/(M-m)
    ggplot() +
      geom_rect(
        data = reference[[j]],
        mapping = aes(
          xmin = position - 0.5,
          xmax = position + 0.5,
          ymin = 0,
          ymax = freq / sum(freq)
        )
      ) +
      geom_line(
        data = plotData,
        mapping = aes(
          x = position,
          y = mean_score
        ),
        color = colour_scale[1],
        linewidth = 0.25
      ) +
      annotate(
        geom = "line",
        x = c(0,0),
        y = c(-Inf, Inf),
        linetype = "dashed"
      ) +
      scale_x_continuous(
        name = "Position (Relative to TSS)",
        breaks = (-10:13)*50
      ) +
      scale_y_continuous(
        name = "Proportion (Reference)",
        limits = c(0,max(reference[[j]]$freq / sum(reference[[j]]$freq))),
        sec.axis = sec_axis(
          transform = ~.*(M-m)/R + m,
          name = "Mean Score"
        )
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(
        paste0(
          dataset_names[i],
          " Score\nCompared to ",
          reference_names[j]
        )
      )
    
    ggsave(
      paste0(dirs[i], "/", reference_datasets[j], "_score_plot.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    M = max(median_scores[[i]]$median_score)
    m = min(median_scores[[i]]$median_score)
    plotData <- median_scores[[i]]
    plotData$median_score = (plotData$median_score - m) * R/(M-m)
    ggplot() +
      geom_rect(
        data = reference[[j]],
        mapping = aes(
          xmin = position - 0.5,
          xmax = position + 0.5,
          ymin = 0,
          ymax = freq / sum(freq)
        )
      ) +
      geom_line(
        data = plotData,
        mapping = aes(
          x = position,
          y = median_score
        ),
        color = colour_scale[1],
        linewidth = 0.25
      ) +
      annotate(
        geom = "line",
        x = c(0,0),
        y = c(-Inf, Inf),
        linetype = "dashed"
      ) +
      scale_x_continuous(
        name = "Position (Relative to TSS)",
        breaks = (-10:13)*50
      ) +
      scale_y_continuous(
        name = "Proportion (Reference)",
        limits = c(0,max(reference[[j]]$freq / sum(reference[[j]]$freq))),
        sec.axis = sec_axis(
          transform = ~.*(M-m)/R + m,
          name = "Mean Score"
        )
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(
        paste0(
          dataset_names[i],
          " Median Score\nCompared to ",
          reference_names[j]
        )
      )
    
    ggsave(
      paste0(dirs[i], "/", reference_datasets[j], "_median_score_plot.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    M = max(probs[[i]]$mean_prob)
    m = min(probs[[i]]$mean_prob)
    plotData <- probs[[i]]
    plotData$mean_prob = (plotData$mean_prob - m) * R/(M-m)
    ggplot() +
      geom_rect(
        data = reference[[j]],
        mapping = aes(
          xmin = position - 0.5,
          xmax = position + 0.5,
          ymin = 0,
          ymax = freq / sum(freq)
        )
      ) +
      geom_line(
        data = plotData,
        mapping = aes(
          x = position,
          y = mean_prob
        ),
        color = colour_scale[1],
        linewidth = 0.25
      ) +
      annotate(
        geom = "line",
        x = c(0,0),
        y = c(-Inf, Inf),
        linetype = "dashed"
      ) +
      scale_x_continuous(
        name = "Position (Relative to TSS)",
        breaks = (-10:13)*50
      ) +
      scale_y_continuous(
        name = "Proportion (Reference)",
        limits = c(0,max(reference[[j]]$freq / sum(reference[[j]]$freq))),
        sec.axis = sec_axis(
          transform = ~.*(M-m)/R + m,
          name = "Mean Probability"
        )
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(
        paste0(
          dataset_names[i],
          " Probability\nCompared to ",
          reference_names[j]
        )
      )
    
    ggsave(
      paste0(dirs[i], "/", reference_datasets[j], "_prob_plot.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
  }
}
