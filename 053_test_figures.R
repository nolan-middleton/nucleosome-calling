#########################################################################
# This R script generates plots for the evaluations of the models.      #
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

##### LOAD DATA #####
test_points <- read.delim(
  file = "evaluationPositions.tsv",
  header = TRUE,
  sep = "\t"
)
test_points$answer <- test_points$answer == "True"
M <- nrow(test_points)

modelComparisons <- cbind(
  test_points[0,],
  data.frame(
    numerator = character(0),
    denominator = character(0),
    score = numeric(0),
    model = character(0)
  )
)

backgroundComparisons <- cbind(
  test_points[0,],
  data.frame(
    numerator = numeric(0),
    denominator = numeric(0),
    score = numeric(0),
    model = character(0)
  )
)

for (model in models) {
  print(paste0(">> ", model, "..."))
  thisModelComparison <- read.delim(
    file = paste0(
      "NucleosomePattern/",
      model,
      "/models/comparisonScores.tsv"
    ),
    header = TRUE,
    sep = "\t"
  )
  
  modelComparisons <- rbind(
    modelComparisons,
    cbind(
      test_points,
      data.frame(
        numerator = rep("flat", M),
        denominator = rep("linear", M),
        score = thisModelComparison$flat_linear,
        model = rep(model, M)
      )
    ),
    cbind(
      test_points,
      data.frame(
        numerator = rep("flat", M),
        denominator = rep("quadratic", M),
        score = thisModelComparison$flat_quadratic,
        model = rep(model, M)
      )
    ),
    cbind(
      test_points,
      data.frame(
        numerator = rep("linear", M),
        denominator = rep("quadratic", M),
        score = thisModelComparison$linear_quadratic,
        model = rep(model, M)
      )
    )
  )
  
  thisBackgroundTable <- read.delim(
    file = paste0("RandomPattern/", model, "/comparisonScores.tsv"),
    header = TRUE,
    sep = "\t"
  )
  for (n in 0:(N-1)) {
    for (d in 0:(N-1)) {
      backgroundComparisons <- rbind(
        backgroundComparisons,
        cbind(
          test_points,
          data.frame(
            numerator = rep(n, M),
            denominator = rep(d, M),
            score = thisBackgroundTable[[
              paste0("X", toString(n), "_", toString(d))
            ]],
            model = rep(model, M)
          )
        )
      )
    }
  }
}
rm(thisBackgroundTable, thisModelComparison, n, d)

evaluationTable <- read.delim(
  "evaluationTable.tsv",
  header = TRUE,
  sep = "\t"
)
curveTable <- read.delim(
  "ROC_PR_curveTable.tsv",
  header = TRUE,
  sep = "\t"
)

##### BACKGROUND COMPARISONS #####
print("> Background plots...")
for (i in 1:length(models)) {
  model <- models[i]
  print(paste0(">> ", model, "..."))
  
  plotData <- backgroundComparisons[
    (
      (backgroundComparisons$model == model) &
      (!backgroundComparisons$answer)
    ),
  ]
  
  # Boxplots
  ggplot() +
    geom_boxplot(
      data = plotData,
      mapping=aes(y=score,x=factor(numerator),fill=factor(denominator)),
      outlier.size = 0.1
    ) +
    scale_fill_discrete(type = colour_scale, name = "Denominator") +
    scale_x_discrete(name = "Numerator") +
    scale_y_continuous("ln(Score)") +
    ggtitle(paste0("Background Comparisons (", dataset_names[i], ")"))
  
  ggsave(
    paste0("RandomPattern/", model, "/boxplot.png"),
    width = 5000,
    height = 2000,
    units = "px"
  )
  
  meanData <- data.frame(
    numerator = numeric(0),
    denominator = numeric(0),
    mean_score = numeric(0),
    std_dev = numeric(0)
  )
  for (n in 0:(N-1)) {
    for (d in 0:(N-1)) {
      meanData <- rbind(
        meanData,
        data.frame(
          numerator = n,
          denominator = d,
          mean_score = mean(
            plotData$score[
              (plotData$numerator == n) & (plotData$denominator == d)
            ]
          ),
          std_dev = sd(
            plotData$score[
              (plotData$numerator == n) & (plotData$denominator == d)
            ]
          )
        )
      )
    }
  }
  
  ggplot(data = meanData) +
    geom_rect(
      mapping = aes(
        xmin = numerator - 0.5,
        xmax = numerator + 0.5,
        ymin = denominator - 0.5,
        ymax = denominator + 0.5,
        fill = mean_score
      )
    ) +
    scale_y_continuous(name = "Denominator", trans = "reverse") +
    geom_text(
      mapping = aes(
        x = numerator,
        y = denominator,
        label = paste0(round(mean_score, 3), "±", round(std_dev, 3)),
        color = mean_score < -0.002
      ),
      size = 1.4
    ) +
    scale_color_discrete(type = c("#000000", "#FFFFFF")) +
    guides(color = guide_none()) +
    scale_x_continuous("Numerator") +
    scale_fill_continuous(name = "Mean ln(Score)", type = "viridis") +
    ggtitle(paste0("Background Comparisons (", dataset_names[i], ")"))
    
  ggsave(
    paste0("RandomPattern/", model, "/heatmap.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
}

##### MODEL COMPARISONS #####
print("> Model comparison plots...")
for (i in 1:length(models)) {
  model = models[i]
  print(paste0(">> ", model, "..."))
  
  plotData <- modelComparisons[
    (modelComparisons$model == model) & (modelComparisons$answer),
  ][,c(6,7,8)]
  
  plotData <- rbind(
    plotData,
    data.frame(
      numerator = rep("linear", M),
      denominator = rep("flat", M),
      score = -plotData$score[
        (
          (plotData$numerator == "flat") &
          (plotData$denominator == "linear")
        )
      ]
    ),
    data.frame(
      numerator = rep("quadratic", M),
      denominator = rep("flat", M),
      score = -plotData$score[
        (
          (plotData$numerator == "flat") &
          (plotData$denominator == "quadratic")
        )
      ]
    ),
    data.frame(
      numerator = rep("quadratic", M),
      denominator = rep("linear", M),
      score = -plotData$score[
        (
          (plotData$numerator == "linear") &
          (plotData$denominator == "quadratic")
        )
      ]
    ),
    data.frame(
      numerator = rep("flat", M),
      denominator = rep("flat", M),
      score = rep(0,M)
    ),
    data.frame(
      numerator = rep("linear", M),
      denominator = rep("linear", M),
      score = rep(0,M)
    ),
    data.frame(
      numerator = rep("quadratic", M),
      denominator = rep("quadratic", M),
      score = rep(0,M)
    )
  )
  
  # Boxplots
  ggplot() +
    geom_boxplot(
      data = plotData,
      mapping=aes(y = score, x = numerator,fill = denominator),
      outlier.size = 0.5
    ) +
    scale_fill_discrete(type = colour_scale, name = "Denominator") +
    scale_x_discrete(name = "Numerator") +
    scale_y_continuous("ln(Score)") +
    ggtitle(paste0("Model Comparisons (", dataset_names[i], ")"))
  
  ggsave(
    paste0("NucleosomePattern/", model, "/models/boxplot.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  meanData <- data.frame(
    numerator = character(0),
    denominator = character(0),
    mean_score = numeric(0),
    std_dev = numeric(0)
  )
  for (n in unique(plotData$numerator)) {
    for (d in unique(plotData[plotData$numerator == n,]$denominator)) {
      meanData <- rbind(
        meanData,
        data.frame(
          numerator = n,
          denominator = d,
          mean_score = mean(
            plotData$score[
              (plotData$numerator == n) & (plotData$denominator == d)
            ]
          ),
          std_dev = sd(
            plotData$score[
              (plotData$numerator == n) & (plotData$denominator == d)
            ]
          )
        )
      )
    }
  }
  
  ggplot(data = meanData) +
    geom_rect(
      mapping = aes(
        xmin = stage(numerator, after_scale = xmin - 0.5),
        xmax = stage(numerator, after_scale = xmax + 0.5),
        ymin = stage(denominator, after_scale = ymin - 0.5),
        ymax = stage(denominator, after_scale = ymax + 0.5),
        fill = mean_score
      )
    ) +
    scale_y_discrete(
      name = "Denominator",
      limits = c("quadratic", "linear", "flat")
    ) +
    geom_text(
      mapping = aes(
        x = numerator,
        y = denominator,
        label = paste0(round(mean_score, 3), "±", round(std_dev, 3)),
        color = mean_score < -0.002
      ),
      size = 4.5
    ) +
    scale_color_discrete(type = c("#000000", "#FFFFFF")) +
    guides(color = guide_none()) +
    scale_x_discrete("Numerator") +
    scale_fill_continuous(name = "Mean ln(Score)", type = "viridis") +
    ggtitle(paste0("Model Comparisons (", dataset_names[i], ")"))
  
  ggsave(
    paste0("NucleosomePattern/", model, "/models/heatmap.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
}

##### MODEL STATISTICS #####
for (threshold in unique(evaluationTable$threshold)) {
  for (i in 1:length(models)) {
    model = models[i]
    ggplot(
      data = evaluationTable[
        (
          (evaluationTable$threshold == threshold) &
          (evaluationTable$dataset == model) &
          !(evaluationTable$statistic %in% c("TP", "TN", "FP", "FN"))
        ),
      ]
    ) +
      geom_col(
        mapping = aes(
          x = statistic,
          y = value,
          fill = model
        ),
        position = position_dodge2()
      ) +
      scale_x_discrete(name = "Statistic") +
      scale_y_continuous(name = "Value") +
      scale_fill_discrete(name = "Model", type = colour_scale) +
      ggtitle(
        paste0(
          "Model Evaluation, Threshold = ",
          toString(threshold),
          " (",
          dataset_names[i],
          ")"
        )
      )
    
    ggsave(
      paste0(
        "NucleosomePattern/",
        model,
        "/models/evaluation",
        str_replace(toString(threshold), "[.]", "_"),
        ".png"
      ),
      width = 2500,
      height = 1000,
      units = "px"
    )
  }
  
  stats <- c(
    "accuracy",
    "precision",
    "sensitivity",
    "NPV",
    "phi",
    "specificity"
  )
  for (stat in stats) {
    ggplot(
      data = evaluationTable[
        (
          (evaluationTable$statistic == stat) &
          (evaluationTable$threshold == threshold)
        ),
      ]
    ) +
      geom_col(
        mapping = aes(x = dataset, y = value, fill = model),
        position = position_dodge2()
      ) +
      scale_x_discrete(
        name = "Dataset",
        breaks = models,
        labels = dataset_names
      ) +
      scale_y_continuous(name = "Value") +
      scale_fill_discrete(name = "Model", type = colour_scale) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      ggtitle(
        paste0(
          str_to_sentence(stat),
          ", Threshold = ",
          toString(threshold)
        )
      )
    
    ggsave(
      paste0(
        "NucleosomePattern/",
        stat,
        str_replace(toString(threshold), "[.]", "_"),
        ".png"
      ),
      width = 2500,
      height = 1000,
      units = "px"
    )
  }
}

##### ROC AND PR CURVES #####
ggplot(data = curveTable[curveTable$background == 0,]) +
  geom_step(
    mapping = aes(
      x = FPR,
      y = sensitivity,
      group = interaction(dataset, model),
      linetype = model,
      color = dataset
    ),
    linewidth = 0.5,
    direction = "hv"
  ) +
  scale_color_discrete(
    type = colour_scale,
    name = "Dataset",
    breaks = models,
    labels = dataset_names
  ) +
  guides(
    linetype = guide_legend(
      title = "Model"
    )
  ) +
  scale_x_continuous(
    name = "1 - Specificity",
    limits = c(0, 1),
    breaks = 0:10 / 10
  ) +
  scale_y_continuous(name="Sensitivity", limits=c(0,1), breaks=1:10/10) +
  ggtitle("Receiver-Operator Characteristic Curve") +
  annotate(
    "segment",
    x = 0,
    xend = 1,
    y = 0,
    yend = 1
  )

ggsave(
  "NucleosomePattern/ROC_curve.png",
  width = 4000,
  height = 3000,
  units = "px"
)

ggplot(data = curveTable[curveTable$background == 0,]) +
  geom_step(
    mapping = aes(
      x = sensitivity,
      y = precision,
      group = interaction(dataset, model),
      linetype = model,
      color = dataset
    ),
    linewidth = 0.5,
    direction = "hv"
  ) +
  scale_color_discrete(
    type = colour_scale,
    name = "Dataset",
    breaks = models,
    labels = dataset_names
  ) +
  guides(
    linetype = guide_legend(
      title = "Model"
    )
  ) +
  scale_x_continuous(name="Sensitivity", limits=c(0, 1), breaks=0:10/10) +
  scale_y_continuous(name="Precision", limits=c(0.5,1), breaks=1:10/10) +
  ggtitle("Precision-Recall Curve") +
  annotate(
    "segment",
    x = 0,
    xend = 1,
    y = 0.5,
    yend = 0.5
  )

ggsave(
  "NucleosomePattern/PR_curve.png",
  width = 4000,
  height = 3000,
  units = "px"
)
