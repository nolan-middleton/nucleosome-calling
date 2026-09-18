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

class_windows = list(
  near = c(0,1),
  mid = c(2,2),
  far = c(3,3),
  distant = c(4,10)
)

for (i in 1:length(models)) {
  model <- models[i]
  print(paste0("> ", model, "..."))
  
  for (mean_fun in c("flat", "linear", "quadratic")) {
    data = data.frame(
      score = numeric(0),
      offset = numeric(0),
      class = numeric(0),
      index = numeric(0),
      category = character(0),
      answer = logical(0)
    )
    
    thisClass <- read.delim(
      paste0(
        "NucleosomePattern/",
        model,
        "/models/windowClasses.tsv"
      ),
      header = TRUE,
      sep = "\t"
    )[[mean_fun]]
    
    thisIndex <- rep(0, nrow(data))
    thisCategory <- rep("None", nrow(data))
    
    for (j in 1:length(class_windows)) {
      window <- class_windows[[j]]
      selected <- (thisClass >= window[1]) & (thisClass <= window[2])
      
      thisIndex[selected] <- 1:sum(selected)
      thisCategory[selected]<-rep(names(class_windows)[j],sum(selected))
    }
    
    for (w in -W:W) {
      data = rbind(
        data,
        data.frame(
          score = read.delim(
            paste0(
              "NucleosomePattern/",
              model,
              "/models/",
              mean_fun,
              "/windowEvalScores_",
              toString(w),
              ".tsv"
            ),
            header = TRUE,
            sep = "\t"
          )$X0,
          offset = rep(w, M),
          class = thisClass,
          index = thisIndex,
          category = thisCategory,
          answer = test_points$answer
        )
      )
    }
    
    # Heatmap
    ggplot() +
      geom_rect(
        data = data[data$category == "near" & data$answer,],
        mapping = aes(
          xmin = offset - 0.5,
          xmax = offset + 0.5,
          ymin = index - 0.5,
          ymax = index + 0.5,
          fill = score
        )
      ) +
      geom_rect(
        data = data[(data$category == "mid") & (data$answer),],
        mapping = aes(
          xmin = offset - 0.5,
          xmax = offset + 0.5,
          ymin = max(
            data$index[(data$category == "near") & (data$answer)]
          ) + 10 + index - 0.5,
          ymax = max(
            data$index[(data$category == "near") & (data$answer)]
          ) + 10 + index + 0.5,
          fill = score
        )
      ) +
      geom_rect(
        data = data[(data$category == "far") & (data$answer),],
        mapping = aes(
          xmin = offset - 0.5,
          xmax = offset + 0.5,
          ymin = max(
            data$index[(data$category == "near") & (data$answer)]
          ) + max(
            data$index[(data$category == "mid") & (data$answer)]
          ) + 20 + index - 0.5,
          ymax = max(
            data$index[(data$category == "near") & (data$answer)]
          ) + max(
            data$index[(data$category == "mid") & (data$answer)]
          ) + 20 + index + 0.5,
          fill = score
        )
      ) +
      geom_rect(
        data = data[(data$category == "distant") & (data$answer),],
        mapping = aes(
          xmin = offset - 0.5,
          xmax = offset + 0.5,
          ymin = max(
            data$index[(data$category == "near") & (data$answer)]
          ) + max(
            data$index[(data$category == "mid") & (data$answer)]
          ) + max(
            data$index[(data$category == "far") & (data$answer)]
          ) + 30 + index - 0.5,
          ymax = max(
            data$index[(data$category == "near") & (data$answer)]
          ) + max(
            data$index[(data$category == "mid") & (data$answer)]
          ) + max(
            data$index[(data$category == "far") & (data$answer)]
          ) + 30 + index + 0.5,
          fill = score
        )
      ) +
      scale_fill_continuous(name = "ln(Score)", type = "viridis") +
      scale_x_continuous(name = "Offset", breaks = -W:W) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      ) +
      ggtitle(paste0("Peak Locations (", dataset_names[i], ")"))
    
    ggsave(
      paste0(
        "NucleosomePattern/",
        model,
        "/models/",
        mean_fun,
        "/windowClassesHeatmap.png"
      ),
      width = 2000,
      height = 1000,
      units = "px"
    )
  }
}
