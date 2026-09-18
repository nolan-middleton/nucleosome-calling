##########################################################################
# This R script generates plots for the distributions of rotational      #
# scores stratified by nucleosome class.                                 #
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
  chrom = character(0),
  dyad = numeric(0),
  score = numeric(0)
)
C <- data.frame()
for (i in 1:length(dirs)) {
  data <- read.delim(
    paste0(dirs[i], "/rotational_scores.tsv"),
    header = FALSE,
    sep = "\t"
  )
  
  classes <- read.delim(
    paste0(dirs[i], "/nucleosome_classification.tsv"),
    header = TRUE,
    sep = "\t"
  )
  
  if (nrow(C) == 0) {
    C <- cbind(
      dataset = rep(dataset_names[i], nrow(classes)),
      classes
    )
  } else {
    C <- rbind(
      C,
      cbind(
        dataset = rep(dataset_names[i], nrow(classes)),
        classes
      )
    )
  }
  
  D <- rbind(
    D,
    data.frame(
      dataset = rep(dataset_names[i], nrow(data)),
      chrom = data$V1,
      dyad = data$V2,
      score = data$V3
    )
  )
}

N <- c(colnames(D), colnames(C)[4:ncol(C)])
for (i in 4:15) {
  D <- cbind(D, rep(NA, nrow(D)))
  C[[i]] <- C[[i]] == "True"
}
for (i in 16:ncol(C)) {
  C[[i]][is.nan(C[[i]])] <- NA
  D <- cbind(D, rep(NA, nrow(D)))
}
colnames(D) <- N

for (dnm in unique(D$dataset)) {
  print(paste0("> ", dnm, "..."))
  for (chrom in unique(D$chrom[D$dataset == dnm])) {
    print(paste0(">> ", chrom, "..."))
    subC <- C[(C$dataset == dnm) & (C$chrom == chrom),]
    for (j in 1:nrow(subC)) {
      selection <- (D$dataset==dnm) & (D$chrom==chrom) & (D$dyad==subC$loc[j])
      D[selection,5:ncol(D)] <- subC[j,4:ncol(C)]
    }
  }
}; save(D, file = "nucleosomeClassDataTable.rda")
load("nucleosomeClassDataTable.rda")

##### PLOTS #####
for (i in 1:length(dataset_names)) {
  print(paste0("> ", dataset_names[i], "..."))
  plotData <- data.frame(
    score = numeric(0),
    class = character(0)
  )
  for (cnm in colnames(D)[5:16]) {
    theseScores <- D$score[(D$dataset == dataset_names[i]) & D[[cnm]]]
    plotData <- rbind(
      plotData,
      data.frame(
        score = theseScores,
        class = rep(cnm, length(theseScores))
      )
    )
    theseScores <- D$score[(D$dataset == dataset_names[i]) & (!D[[cnm]])]
    plotData <- rbind(
      plotData,
      data.frame(
        score = theseScores,
        class = rep(paste0("!", cnm), length(theseScores))
      )
    )
  }
  theseScores <- D$score[(D$dataset == dataset_names[i])&(D$X.1&(!D$in_gene))]
  plotData <- rbind(
    plotData,
    data.frame(
      score = theseScores,
      class = rep("distal", length(theseScores))
    )
  )
  
  plotData$class <- factor(
    plotData$class,
    levels = c(
      "in_gene", "!in_gene",
      "X.1", "!X.1", "distal",
      "ARS", "!ARS",
      "LTR", "!LTR",
      "misc_mRNA", "!misc_mRNA",
      "tRNA", "!tRNA",
      "snoRNA", "!snoRNA",
      "ncRNA", "!ncRNA",
      "telomeric", "!telomeric",
      "centromeric", "!centromeric",
      "in_low_expr_gene", "!in_low_expr_gene",
      "in_high_expr_gene", "!in_high_expr_gene"
    )
  )
  
  ggplot(data = plotData) +
    geom_boxplot(mapping = aes(x = class, y = score), outlier.size = 0.2) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_continuous(name = "Score") +
    scale_x_discrete(
      name = "Nucleosome Classification",
      labels = c(
        "Genic", "Intergenic", "+1", "Not +1", "Distal",
        "ARS", "Not ARS",
        "LTR", "Not LTR",
        "Misc. mRNA", "Not Misc. mRNA",
        "tRNA", "Not tRNA",
        "snoRNA", "Not snoRNA",
        "ncRNA", "Not ncRNA",
        "Telomeric", "Not Telomeric",
        "Centromeric", "Not Centromeric",
        "In Low Expr. Gene", "Not in Low Expr. Gene",
        "In High Expr. Gene", "Not in High Expr. Gene"
      )
    )
  
  ggsave(
    paste0(dirs[i], "/RotationalScoreDistribution_byclass.png"),
    width = 4000,
    height = 1500,
    units = "px"
  )
  
  plotData2 <- D[D$dataset == dataset_names[i], c(4,17:ncol(D))]
  modData <- data.frame(score = numeric(0), class = character(0))
  for (j in 2:ncol(plotData2)) {
    thisData <- plotData2[!is.na(plotData2[,j]),c(1,j)]
    colnames(thisData) <- c("y", "x")
    model = lm(y ~ x, data = thisData)
    
    C <- cor.test(thisData$x, thisData$y)
    r <- C$estimate
    p <- C$p.value
    
    ggplot(data = thisData) +
      geom_point(mapping = aes(x = x, y = y), size = 0.2) +
      scale_x_continuous(name = "Modification Value") +
      scale_y_continuous(name = "Rotational Score") +
      ggtitle(
        paste0(
          colnames(plotData2)[j],
          ", r=",
          toString(signif(r, 5)),
          ", p=",
          toString(signif(p, 5))
        )
      )
    
    ggsave(
      paste0(dirs[i], "/", colnames(plotData2)[j], ".png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    bottomScores <- thisData$y[thisData$x < quantile(thisData$x, 0.25)]
    topScores <- thisData$y[thisData$x > quantile(thisData$x, 0.75)]
    modData <- rbind(
      modData,
      data.frame(
        score = c(bottomScores, topScores),
        class = c(
          rep(paste0("!", colnames(plotData2)[j]), length(bottomScores)),
          rep(colnames(plotData2)[j], length(topScores))
        )
      )
    )
  }
  
  modData$class <- factor(
    modData$class,
    levels = paste0(
      c("!", ""),
      rep(colnames(plotData2)[2:ncol(plotData2)], each = 2)
    )
  )
  
  ggplot(data = modData) +
    geom_boxplot(mapping = aes(x = class, y = score)) +
    scale_y_continuous(name = "Score") +
    scale_x_discrete(name = "Class") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    paste0(dirs[i], "/RotationalScoreDistribution_bymod.png"),
    width = 4000,
    height = 1500,
    units = "px"
  )
}
