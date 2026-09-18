#########################################################################
# This R script generates plots for the nucleosome CPD damage pattern   #
# for easy visualization.                                               #
#########################################################################

##### SETUP #####

# Imports
require(ggplot2)
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
  "Cellular - Room Temperature",
  "Cellular - Ice"
)

##### MAIN LOOP #####
for (i in 1:length(models)) {
  model <- models[i]
  dataname <- datanames[[i]]
  dataset_name <- dataset_names[i]
  print(paste0("> ", toString(model), "..."))
  
  data <- read.delim(
    paste0("NucleosomePattern_all/", model, "/damagePattern.tsv"),
    header = FALSE,
    sep = "\t"
  )
  
  minus_values <- data.frame(chrom = character(0))
  minus_selected <- data.frame(chrom = character(0))
  plus_values <- data.frame(chrom = character(0))
  plus_selected <- data.frame(chrom = character(0))
  for (j in 1:ncol(data)) {
    thisCol <- data.frame(x = numeric(0))
    colnames(thisCol) <- paste0("V", toString(j))
    minus_values <- cbind(minus_values, thisCol)
    
    thisCol <- data.frame(x = numeric(0))
    colnames(thisCol) <- paste0("V", toString(j))
    plus_values <- cbind(plus_values, thisCol)
    
    thisCol <- data.frame(x = logical(0))
    colnames(thisCol) <- paste0("V", toString(j))
    minus_selected <- cbind(minus_selected, thisCol)
    
    thisCol <- data.frame(x = logical(0))
    colnames(thisCol) <- paste0("V", toString(j))
    plus_selected <- cbind(plus_selected, thisCol)
  }
  
  thisDir <- paste0("NucleosomePattern_all/", model, "/minus_vals")
  for (file in list.files(thisDir)) {
    chrom <- unlist(strsplit(file, "[.]"))[1]
    thisData <- read.delim(
      paste0(thisDir, "/", file),
      sep = "\t",
      header = FALSE
    )
    minus_values <- rbind(
      minus_values,
      cbind(
        data.frame(
          chrom = rep(chrom, nrow(thisData))
        ),
        thisData
      )
    )
  }
  
  thisDir <- paste0("NucleosomePattern_all/", model, "/plus_vals")
  for (file in list.files(thisDir)) {
    chrom <- unlist(strsplit(file, "[.]"))[1]
    thisData <- read.delim(
      paste0(thisDir, "/", file),
      sep = "\t",
      header = FALSE
    )
    plus_values <- rbind(
      plus_values,
      cbind(
        data.frame(
          chrom = rep(chrom, nrow(thisData))
        ),
        thisData
      )
    )
  }
  
  thisDir <- paste0("NucleosomePattern_all/", model, "/minus_selected")
  for (file in list.files(thisDir)) {
    chrom <- unlist(strsplit(file, "[.]"))[1]
    thisData <- read.delim(
      paste0(thisDir, "/", file),
      sep = "\t",
      header = FALSE
    ) == 1
    minus_selected <- rbind(
      minus_selected,
      cbind(
        data.frame(
          chrom = rep(chrom, nrow(thisData))
        ),
        thisData
      )
    )
  }
  
  thisDir <- paste0("NucleosomePattern_all/", model, "/plus_selected")
  for (file in list.files(thisDir)) {
    chrom <- unlist(strsplit(file, "[.]"))[1]
    thisData <- read.delim(
      paste0(thisDir, "/", file),
      sep = "\t",
      header = FALSE
    ) == 1
    plus_selected <- rbind(
      plus_selected,
      cbind(
        data.frame(
          chrom = rep(chrom, nrow(thisData))
        ),
        thisData
      )
    )
  }
  
  plotData <- data.frame(
    x = rep(-73:72 + 0.5, 2),
    mean = c(unlist(data[3,]), unlist(data[4,])),
    var = c(unlist(data[5,]), unlist(data[6,])),
    n = c(unlist(data[7,]), unlist(data[8,])),
    strand = c(rep("-", ncol(data)), rep("+", ncol(data)))
  )
  
  # Mean Pattern
  ggplot(data = plotData, mapping = aes(x=x, y=mean, color=strand)) +
    geom_line() +
    scale_color_discrete(type = colour_scale, name = "Strand") +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(paste0("Nucleosome Damage Pattern (", dataset_name, ")"))
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/mean.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Variance
  ggplot(data = plotData, mapping = aes(x = x, y = var, color = strand)) +
    geom_line() +
    scale_color_discrete(type = colour_scale, name = "Strand") +
    scale_y_continuous(name=bquote("Variance of"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(paste0("Nucleosome Damage Pattern (", dataset_name, ")"))
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/var.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Align Pattern
  plotData$x[plotData$strand=="-"] <- -plotData$x[plotData$strand=="-"]
  
  # Aligned Mean Pattern
  ggplot(data = plotData, mapping = aes(x=x, y=mean, color=strand)) +
    geom_line() +
    scale_color_discrete(type = colour_scale, name = "Strand") +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(paste0("Nucleosome Damage Pattern (", dataset_name, ")"))
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/mean_aligned.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Get Trended data
  plotData <- data.frame(
    x = -73:72 + 0.5,
    detrended = unlist(data[1,]),
    trended = unlist(data[2,])
  )
  
  # Detrended Plot
  ggplot(data=plotData, mapping=aes(x=x,y=detrended)) +
    geom_line() +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(
      paste0("Detrended Nucleosome Damage Pattern\n(", dataset_name, ")")
    )
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/detrended_mean.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # Trended Plot
  ggplot(data=plotData, mapping=aes(x=x,y=trended)) +
    geom_line() +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(
      paste0("Overall Nucleosome Damage Pattern\n(", dataset_name, ")")
    )
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/overall_mean.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  # The trend itself
  ggplot(data=plotData, mapping=aes(x=x,y=trended - detrended)) +
    geom_line() +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(
      paste0("Overall Nucleosome Damage Pattern\n(", dataset_name, ")")
    )
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/trend.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
}
