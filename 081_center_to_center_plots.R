##########################################################################
# This R script generates plots for the distance evaluations of the      #
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
  "Cellular vs. Ice Combined",
  "Weiner et al."
)
dirs <- c(
  "NucleosomePattern/cl_rt_log2/distance_evals",
  "NucleosomePattern/cl_ice_log2/distance_evals",
  "NucleosomePattern/cl_ice_abs_diff/distance_evals",
  "NucleosomePattern/composite_models/composite_ice/distance_evals",
  "NucleosomePattern/weiner_brogaard_distance_evals"
)

reference_datasets <- c("strong", "decent", "all", "Weiner")
reference_names <- c(
  "Brogaard et al. (Score > 5)",
  "Brogaard et al. (Score > 1)",
  "Brogaard et al. (Complete)",
  "Weiner et al."
)

algorithms <- list(
  c("greedy", "viterbi"),
  c("greedy", "viterbi"),
  c("greedy", "viterbi"),
  c("greedy", "viterbi"),
  c("Weiner")
)
algorithm_names <- list(
  c("Greedy", "Dynamic"),
  c("Greedy", "Dynamic"),
  c("Greedy", "Dynamic"),
  c("Greedy", "Dynamic"),
  c("")
)

distance_frequencies <- data.frame(
  distance = numeric(0),
  frequency = numeric(0),
  total = numeric(0),
  reference = character(0),
  putative = character(0),
  algorithm = character(0)
)
periodic <- data.frame(
  proportion = numeric(0),
  reference = character(0),
  putative = character(0),
  algorithm = character(0)
)
periodic_E <- c()
anti_periodic <- data.frame(
  proportion = numeric(0),
  reference = character(0),
  putative = character(0),
  algorithm = character(0)
)
anti_periodic_E <- c()
for (i in 1:length(dataset_names)) {
  print(paste0("> ", dataset_names[i], "..."))
  thisPeriodicTable <- read.delim(
    paste0(dirs[i], "/periodic_proportions.txt"),
    header = FALSE,
    sep = "\t"
  )
  periodic_E <- c(
    periodic_E,
    as.numeric(unlist(strsplit(thisPeriodicTable$V1[1],"="))[2])
  )
  thisPeriodicTable <- thisPeriodicTable[2:nrow(thisPeriodicTable),]
  colnames(thisPeriodicTable) <- c(
    "reference",
    "algorithm",
    "successes",
    "total"
  )
  
  periodic <- rbind(
    periodic,
    data.frame(
      proportion = thisPeriodicTable$successes/thisPeriodicTable$total,
      reference = reference_names[
        unlist(
          lapply(thisPeriodicTable$reference, grep, reference_datasets)
        )
      ],
      putative = rep(dataset_names[i], nrow(thisPeriodicTable)),
      algorithm = algorithm_names[[i]][
        unlist(lapply(thisPeriodicTable$algorithm,grep,algorithms[[i]]))
      ]
    )
  )
  
  thisAntiPeriodicTable <- read.delim(
    paste0(dirs[i], "/anti_periodic_proportions.txt"),
    header = FALSE,
    sep = "\t"
  )
  anti_periodic_E <- c(
    anti_periodic_E,
    as.numeric(unlist(strsplit(thisAntiPeriodicTable$V1[1],"="))[2])
  )
  thisAntiPeriodicTable <- thisAntiPeriodicTable[
    2:nrow(thisAntiPeriodicTable),
  ]
  colnames(thisAntiPeriodicTable) <- c(
    "reference",
    "algorithm",
    "successes",
    "total"
  )
  
  anti_periodic <- rbind(
    anti_periodic,
    data.frame(
      proportion = thisAntiPeriodicTable$successes/
        thisAntiPeriodicTable$total,
      reference = reference_names[
        unlist(
          lapply(thisAntiPeriodicTable$reference,grep,reference_datasets)
        )
      ],
      putative = rep(dataset_names[i], nrow(thisAntiPeriodicTable)),
      algorithm = algorithm_names[[i]][
        unlist(
          lapply(thisAntiPeriodicTable$algorithm,grep,algorithms[[i]])
        )
      ]
    )
  )
  
  for (j in 1:length(reference_names)) {
    print(paste0(">> ", reference_names[j], "..."))
    for (k in 1:length(algorithm_names[[i]])) {
      if (reference_datasets[j] != algorithms[[i]][k]) {
        print(paste0(">>> ", algorithm_names[[i]][k], "..."))
        thisFreqTable <- read.delim(
          paste0(
            dirs[i],
            "/",
            reference_datasets[j],
            "_",
            algorithms[[i]][k],
            "_freq_table.tsv"
          ),
          header = FALSE,
          sep = "\t"
        )
        names(thisFreqTable) <- c("distance", "frequency")
        thisFreqTable <- thisFreqTable[thisFreqTable$frequency > 0,]
        N <- nrow(thisFreqTable)
        
        distance_frequencies <- rbind(
          distance_frequencies,
          data.frame(
            distance = thisFreqTable$distance,
            frequency = thisFreqTable$frequency,
            total = rep(sum(thisFreqTable$frequency), N),
            reference = rep(reference_names[j], N),
            putative = rep(dataset_names[i], N),
            algorithm = rep(algorithm_names[[i]][k], N)
          )
        )
      }
    }
  }
}
rm(thisFreqTable, thisPeriodicTable, N, i, j, k)

##### EXTRA SETUP #####
W <- 100 # The width of distances to consider
weinerPoints <- data.frame(
  mean = numeric(0),
  std = numeric(0),
  n = numeric(0),
  reference = character(0)
)

for (j in 1:length(reference_names)) {
  if (reference_names[j] != "Weiner et al.") {
    thisDataSubset <- distance_frequencies[
      (distance_frequencies$putative == "Weiner et al.") &
        (distance_frequencies$reference == reference_names[j]) &
        (abs(distance_frequencies$distance) <= W),
    ]
    
    N <- sum(thisDataSubset$frequency)
    M <- sum(
      abs(thisDataSubset$distance)*thisDataSubset$frequency
    ) / N
    S <- sum(
      thisDataSubset$frequency*(abs(thisDataSubset$distance)-M)^2
    ) / N
    weinerPoints <- rbind(
      weinerPoints,
      data.frame(mean = M, std=S, n=N, reference=reference_names[j])
    )
    
    thisDataSubset <- periodic[
      (periodic$putative == "Weiner et al.") &
        (periodic$reference == reference_names[j]),
    ]
  }
}

##### MAIN LOOPS #####
B <- 15 # The number of bins for the periodic freqpoly graph.
for (i in 1:length(dataset_names)) {
  print(paste0("> ", dataset_names[i], "..."))
  if (dataset_names[i] != "Weiner et al.")
  {
    ##### REGULAR DATASETS #####
    pointPlotData <- data.frame(
      mean = numeric(0),
      std = numeric(0),
      n = numeric(0),
      reference = character(0),
      algorithm = character(0)
    )
    for (j in 1:length(reference_names)) {
      print(paste0(">> ", reference_names[j], "..."))
      for (k in 1:length(algorithm_names[[i]])) {
        thisDataSubset <- distance_frequencies[
          (distance_frequencies$putative == dataset_names[i]) &
            (distance_frequencies$reference == reference_names[j]) &
            (distance_frequencies$algorithm==algorithm_names[[i]][k]) &
            (abs(distance_frequencies$distance) <= W),
        ]
        
        N <- sum(thisDataSubset$frequency)
        M <- sum(
          abs(thisDataSubset$distance)*thisDataSubset$frequency
        ) / N
        S <- sum(
          thisDataSubset$frequency*(abs(thisDataSubset$distance)-M)^2
        ) / N
        pointPlotData <- rbind(
          pointPlotData,
          data.frame(
            mean = M,
            std = S,
            n = N,
            reference = reference_names[j],
            algorithm = algorithm_names[[i]][k]
          )
        )
      }
      
      if (reference_names[j] != "Weiner et al.")
      {
        ggplot() +
          geom_line(
            data = distance_frequencies[
              (distance_frequencies$putative == dataset_names[i]) &
                (distance_frequencies$reference == reference_names[j]),
            ],
            mapping = aes(
              x = distance,
              y = frequency / total,
              color = algorithm
            )
          ) +
          scale_color_discrete(
            type = colour_scale,
            name = "Placement\nAlgorithm"
          ) +
          scale_x_continuous(
            name = "Distance to Reference Dyad",
            limits = c(-W, W)
          ) +
          scale_y_continuous(name = "Proportion") +
          ggtitle(
            paste0(
              dataset_names[i],
              " Distances to\n",
              reference_names[j]
            )
          ) +
          geom_line(
            data = distance_frequencies[
              (distance_frequencies$putative == "Weiner et al.") &
                (distance_frequencies$reference == reference_names[j]),
            ],
            mapping = aes(
              x = distance,
              y = frequency / total,
              linetype = "Weiner et al."
            ),
          ) +
          guides(
            linetype = guide_legend(title = NULL)
          )
        
        ggsave(
          paste0(dirs[i], "/", reference_datasets[j], "_distances.png"),
          width = 2000,
          height = 1000,
          units = "px"
        )
      } else {
        ggplot() +
          geom_line(
            data = distance_frequencies[
              (distance_frequencies$putative == dataset_names[i]) &
                (distance_frequencies$reference == reference_names[j]),
            ],
            mapping = aes(
              x = distance,
              y = frequency / total,
              color = algorithm
            )
          ) +
          scale_color_discrete(
            type = colour_scale,
            name = "Placement\nAlgorithm"
          ) +
          scale_x_continuous(
            name = "Distance to Reference Dyad",
            limits = c(-W, W)
          ) +
          scale_y_continuous(name = "Proportion") +
          ggtitle(
            paste0(
              dataset_names[i],
              " Distances to\n",
              reference_names[j]
            )
          )
        
        ggsave(
          paste0(dirs[i], "/", reference_datasets[j], "_distances.png"),
          width = 2000,
          height = 1000,
          units = "px"
        )
      }
    }
    
    # Point plot
    ggplot() +
      geom_point(
        data = pointPlotData,
        mapping = aes(
          x = algorithm,
          y = mean,
          color = reference
        ),
        position = position_dodge(width = 0.1)
      ) +
      geom_errorbar(
        data = pointPlotData,
        mapping = aes(
          x = algorithm,
          ymin = mean - std/sqrt(n),
          ymax = mean + std/sqrt(n),
          color = reference
        ),
        position = position_dodge(width = 0.1),
        width = 0.1
      ) +
      geom_point(
        data = weinerPoints,
        mapping = aes(
          x = "Weiner et al.",
          y = mean,
          color = reference
        ),
        position = position_dodge(width = 0.1)
      ) +
      geom_errorbar(
        data = weinerPoints,
        mapping = aes(
          x = "Weiner et al.",
          ymin = mean - std/sqrt(n),
          ymax = mean + std/sqrt(n),
          color = reference
        ),
        position = position_dodge(width = 0.1),
        width = 0.1
      ) +
      scale_color_discrete(type=colour_scale, name="Reference Dyads") +
      scale_x_discrete(name = "Putative Dyad Placement") +
      scale_y_continuous(name = "Mean Distance to True Dyad") +
      ggtitle(
        paste0(
          "Average Distances for ",
          dataset_names[i],
          "\nDyad Placements within ",
          toString(W),
          "bp of True Dyad"
        )
      )
    
    ggsave(
      paste0(dirs[i], "/mean_distances.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
  }
  
  # Periodic and Anti-Periodic Plot
  if (dataset_names[i] != "Weiner et al.") {
    ggplot(
      data = periodic[
        periodic$putative == dataset_names[i] |
          periodic$putative == "Weiner et al.",
      ]
    ) +
      geom_col(
        mapping = aes(x = reference, y = proportion, fill = algorithm),
        position = position_dodge()
      ) +
      scale_fill_discrete(
        name = "Placement\nAlgorithm",
        type = colour_scale,
        limits = c("Dynamic", "Greedy"),
        na.value = "#777777"
      ) +
      scale_y_continuous(name = "Proportion") +
      scale_x_discrete(name = "Reference Dataset") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle(
        paste0(
          dataset_names[i],
          " Proportion of Nucleosomes\n",
          "Placed at Minor Out Positions"
        )
      ) +
      guides(
        custom = guide_custom(
          grob = grid::polygonGrob(
            x = c(0,1,1,0),
            y = c(0,0,1,1),
            gp = grid::gpar(fill = "#777777", col = "grey92")
          ),
          width = unit(0.6, "cm"),
          height = unit(0.6, "cm"),
          title = "Weiner et al."
        )
      ) +
      annotate(
        geom = "line",
        x = c(-Inf, Inf),
        y = c(periodic_E[i], periodic_E[i]),
        linetype = "dashed"
      )
    
    ggsave(
      paste0(dirs[i], "/periodic_proportion_barplot.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    ggplot(
      data = anti_periodic[
        anti_periodic$putative == dataset_names[i] |
          anti_periodic$putative == "Weiner et al.",
      ]
    ) +
      geom_col(
        mapping = aes(x = reference, y = proportion, fill = algorithm),
        position = position_dodge()
      ) +
      scale_fill_discrete(
        name = "Placement\nAlgorithm",
        type = colour_scale,
        limits = c("Dynamic", "Greedy"),
        na.value = "#777777"
      ) +
      scale_y_continuous(name = "Proportion") +
      scale_x_discrete(name = "Reference Dataset") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle(
        paste0(
          dataset_names[i],
          " Proportion of Nucleosomes\n",
          "Placed at Minor In Positions"
        )
      ) +
      guides(
        custom = guide_custom(
          grob = grid::polygonGrob(
            x = c(0,1,1,0),
            y = c(0,0,1,1),
            gp = grid::gpar(fill = "#777777", col = "grey92")
          ),
          width = unit(0.6, "cm"),
          height = unit(0.6, "cm"),
          title = "Weiner et al."
        )
      ) +
      annotate(
        geom = "line",
        x = c(-Inf, Inf),
        y = c(anti_periodic_E[i], anti_periodic_E[i]),
        linetype = "dashed"
      )
    
    ggsave(
      paste0(dirs[i], "/anti_periodic_proportion_barplot.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
  }
}
