#########################################################################
# This R script uses the Lomb-Scargle periodogram to analyze the        #
# periodicity of the datasets.                                          #
#########################################################################

##### SETUP #####

# Imports
require(ggplot2)
require(lomb)
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
  
  ##### INDIVIDUAL #####
  L_minus <- lsp(
    unlist(data[4,]),
    times = -73:72 + 0.5,
    type = "period",
    to = 15,
    ofac = 100
  )
  L_plus <- lsp(
    unlist(data[5,]),
    times = -73:72 + 0.5,
    type = "period",
    to = 15,
    ofac = 100
  )
  
  plotData <- data.frame(
    period = c(L_minus$scanned, L_plus$scanned),
    normalized_power = c(L_minus$power, L_plus$power),
    strand = c(
      rep("-", length(L_minus$scanned)),
      rep("+", length(L_plus$scanned))
    )
  )
  
  write.table(
    plotData,
    file = paste0("NucleosomePattern_all/", model, "/periodogram.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  peaks <- data.frame(
    period = c(L_minus$peak.at, L_plus$peak.at),
    strand = c(
      rep("-", length(L_minus$peak.at)),
      rep("+", length(L_plus$peak.at))
    )
  )
  
  write.table(
    peaks,
    file=paste0("NucleosomePattern_all/",model,"/periodogram_peaks.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # Periodicity
  ggplot(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = strand)
  ) +
    geom_line() +
    scale_color_discrete(type = colour_scale, name = "Strand") +
    scale_y_continuous(name = bquote("Normalized Power"~.(dataname))) +
    scale_x_continuous(name = "Period", breaks = 2:15) +
    ggtitle(paste0("Nucleosome Periodogram (", dataset_name, ")")) +
    annotate(
      "line",
      x = c(L_minus$peak.at[1], L_minus$peak.at[1]),
      y = c(-Inf, Inf),
      linetype = "dashed",
      color = colour_scale[1]
    ) +
    annotate(
      "line",
      x = c(L_plus$peak.at[1], L_plus$peak.at[1]),
      y = c(-Inf, Inf),
      linetype = 4,
      color = colour_scale[2]
    )
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/periodogram.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### COMBINED #####
  L <- lsp(
    unlist(data[1,]),
    times = -73:72 + 0.5,
    type = "period",
    to = 15,
    ofac = 100
  )
  
  plotData <- data.frame(
    period = L$scanned,
    normalized_power = L$power
  )
  
  write.table(
    plotData,
    file = paste0(
      "NucleosomePattern_all/",
      model,
      "/detrended_periodogram.tsv"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  peaks <- data.frame(period = L$peak.at)
  
  write.table(
    peaks,
    file=paste0(
      "NucleosomePattern_all/",
      model,
      "/detrended_periodogram_peaks.tsv"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  # Combined Periodicity
  ggplot(
    data = plotData,
    mapping = aes(x = period, y = normalized_power)
  ) +
    geom_line() +
    scale_y_continuous(name = bquote("Normalized Power"~.(dataname))) +
    scale_x_continuous(name = "Period", breaks = 2:15) +
    ggtitle(paste0("Nucleosome Periodogram (", dataset_name, ")")) +
    annotate(
      "line",
      x = c(L$peak.at[1], L$peak.at[1]),
      y = c(-Inf, Inf),
      linetype = "dashed"
    )
  
  ggsave(
    paste0("NucleosomePattern_all/", model, "/detrended_periodogram.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
}
