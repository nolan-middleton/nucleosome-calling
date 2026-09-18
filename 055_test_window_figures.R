#########################################################################
# This R script generates plots for the window scores for the models.   #
#########################################################################

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

N <- as.numeric(readLines("RandomPattern/datasets.txt"))

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

W <- 10

##### LOAD DATA #####
test_points <- read.delim(
  "evaluationPositions.tsv",
  sep = "\t",
  header = TRUE
)
test_points$answer <- test_points$answer == "True"
M <- nrow(test_points)

R <- 10
random_vals <- ceiling(runif(R) * M/2)

for (i in 1:length(models)) {
  model <- models[i]
  print(paste0("> ", model, "..."))
  
  mean_scores <- data.frame(
    mean_fun = character(0),
    background = numeric(0),
    offset = numeric(0),
    mean_score = numeric(0),
    sd_score = numeric(0),
    answer = logical(0)
  )
  random_scores <- data.frame(
    mean_fun = character(0),
    background = numeric(0),
    offset = numeric(0),
    score = numeric(0),
    answer = logical(0),
    index = numeric(0)
  )
  for (mean_fun in c("flat", "linear", "quadratic")) {
    print(paste0(">> ", mean_fun, "..."))
    for (w in -W:W) {
      thisData <- read.delim(
        paste0(
          "NucleosomePattern/",
          model,
          "/models/",
          mean_fun,
          "/windowEvalScores_",
          toString(w),
          ".tsv"
        )
      )
      
      for (n in 0:(N-1)) {
        theseScores <- thisData[[paste0("X", toString(n))]]
        mean_scores <- rbind(
          mean_scores,
          data.frame(
            mean_fun = mean_fun,
            background = n,
            offset = w,
            mean_score = mean(theseScores[test_points$answer]),
            sd_score = sd(theseScores[test_points$answer]),
            answer = TRUE
          ),
          data.frame(
            mean_fun = mean_fun,
            background = n,
            offset = w,
            mean_score = mean(theseScores[!test_points$answer]),
            sd_score = sd(theseScores[!test_points$answer]),
            answer = FALSE
          )
        )
        
        random_scores <- rbind(
          random_scores,
          data.frame(
            mean_fun = rep(mean_fun, R),
            background = rep(n, R),
            offset = rep(w, R),
            score = theseScores[test_points$answer][random_vals],
            answer = rep(TRUE, R),
            index = random_vals
          ),
          data.frame(
            mean_fun = rep(mean_fun, R),
            background = rep(n, R),
            offset = rep(w, R),
            score = theseScores[!test_points$answer][random_vals],
            answer = rep(FALSE, R),
            index = random_vals
          )
        )
      }
    }
  }
  
  # Messy plots with all backgrounds
  ggplot(data = mean_scores[mean_scores$answer,]) +
    geom_line(
      mapping = aes(
        x = offset,
        y = mean_score,
        color = mean_fun,
        group = interaction(mean_fun, background)
      )
    ) +
    scale_y_continuous(name = "Mean ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    ggtitle(paste0("Nucleosome Positions (", dataset_names[i], ")"))
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_meanNucleosomeScores.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ggplot(data = mean_scores[!mean_scores$answer,]) +
    geom_line(
      mapping = aes(
        x = offset,
        y = mean_score,
        color = mean_fun,
        group = interaction(mean_fun, background)
      )
    ) +
    scale_y_continuous(name = "Mean ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    ggtitle(paste0("Random Positions (", dataset_names[i], ")"))
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_meanRandomScores.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Clean plots with only one background
  ggplot(
    data = mean_scores[mean_scores$answer & mean_scores$background == 0,]
  )+
    geom_line(
      mapping = aes(
        x = offset,
        y = mean_score,
        color = mean_fun
      )
    ) +
    scale_y_continuous(name = "Mean ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    ggtitle(paste0("Nucleosome Positions (", dataset_names[i], ")"))
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_oneMeanNucleosomeScores.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ggplot(
    data = mean_scores[!mean_scores$answer & mean_scores$background == 0,]
  ) +
    geom_line(
      mapping = aes(
        x = offset,
        y = mean_score,
        color = mean_fun,
        group = interaction(mean_fun, background)
      )
    ) +
    scale_y_continuous(name = "Mean ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    ggtitle(paste0("Random Positions (", dataset_names[i], ")"))
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_oneMeanRandomScores.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Point plots with error bars
  ggplot(
    data = mean_scores[mean_scores$answer & mean_scores$background == 0,]
  )+
    geom_line(mapping = aes(x=offset, y=mean_score, color=mean_fun)) +
    geom_point(mapping = aes(x=offset, y=mean_score, color=mean_fun)) +
    geom_errorbar(
      mapping = aes(
        x = offset,
        ymin = mean_score - sd_score/sqrt(M/2),
        ymax = mean_score + sd_score/sqrt(M/2),
        color = mean_fun
      )
    ) +
    scale_y_continuous(name = "Mean ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    ggtitle(paste0("Nucleosome Positions (", dataset_names[i], ")"))
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_nucleosomeScores_stderr.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ggplot(
    data = mean_scores[!mean_scores$answer & mean_scores$background == 0,]
  )+
    geom_line(mapping = aes(x=offset, y=mean_score, color=mean_fun)) +
    geom_point(mapping = aes(x=offset, y=mean_score, color=mean_fun)) +
    geom_errorbar(
      mapping = aes(
        x = offset,
        ymin = mean_score - sd_score/sqrt(M/2),
        ymax = mean_score + sd_score/sqrt(M/2),
        color = mean_fun
      )
    ) +
    scale_y_continuous(name = "Mean ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    ggtitle(paste0("Random Positions (", dataset_names[i], ")"))
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_randomScores_stderr.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Combined plot
  ggplot(data = mean_scores[mean_scores$background == 0,]) +
    geom_line(
      mapping = aes(
        x = offset,
        y = mean_score,
        color = mean_fun,
        linetype = answer,
        group = interaction(mean_fun, background, answer)
      )
    ) +
    scale_y_continuous(name = "Mean ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    guides(linetype = guide_legend(title = "Nucleosome")) +
    ggtitle(paste0("Average Window Scores (", dataset_names[i], ")"))
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_meanScores.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Individual window scores plots
  ggplot(
    data = random_scores[
      random_scores$answer &
      random_scores$background == 0 &
      random_scores$mean_fun == "flat",
    ]
  ) +
    geom_line(mapping = aes(x = offset, y = score, group = index)) +
    scale_y_continuous(name = "ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    ggtitle(
      paste0(
        "Individual Nucleosome Window Scores (",
        dataset_names[i],
        ", Flat Mean)"
      )
    )
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_individualNucleosomeScores.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ggplot(
    data = random_scores[
      !random_scores$answer &
      random_scores$background == 0 &
      random_scores$mean_fun == "flat",
    ]
  ) +
    geom_line(mapping = aes(x = offset, y = score, group = index)) +
    scale_y_continuous(name = "ln(Score)") +
    scale_x_continuous(name = "Offset", breaks = -W:W) +
    ggtitle(
      paste0(
        "Individual Random Window Scores (",
        dataset_names[i],
        ", Flat Mean)"
      )
    )
  
  ggsave(
    paste0(
      "NucleosomePattern/",
      model,
      "/models/window_individualRandomScores.png"
    ),
    width = 2000,
    height = 1000,
    units = "px"
  )
}
