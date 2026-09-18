#########################################################################
# This R script generates plots for the maximum likelihood estimates    #
# for visualization and evaluation.                                     #
#########################################################################

##### SETUP #####

# Imports
require(ggplot2)
require(jsonlite)
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
  "Cellular - Room Temperature",
  "Cellular - Ice"
)

smooth_x <- -730:720 / 10 + 0.5
X <- -73:72 + 0.5

##### MAIN LOOP #####
for (i in 1:length(models)) {
  model <- models[i]
  dataname <- datanames[[i]]
  dataset_name <- dataset_names[i]
  print(paste0("> ", toString(model), "..."))
  for (n in 0:(N-1)) {
    print(paste0(">> ", toString(n), "..."))
    ##### LOAD DATA #####
    MLEs <- read_json(
      paste0("Randompattern/",model,"/",toString(n),"/MLE_values.json")
    )
    
    data <- read.delim(
      paste0(
        "RandomPattern/",
        model,
        "/",
        toString(n),
        "/damagePattern.tsv"
      ),
      header = FALSE,
      sep = "\t"
    )
    
    trendData <- data.frame(
      x = X,
      y = unlist(data[2,])
    )
    
    ##### FLAT GRAPH #####
    ggplot() + 
      geom_line(
        data = trendData,
        mapping = aes(x = x, y = y)
      ) +
      geom_line(
        data = data.frame(
          x = smooth_x,
          y = MLEs$a
        ),
        mapping = aes(x = x, y = y),
        color = colour_scale[1]
      ) +
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
        bquote(
          "Mean ~ "*.(round(MLEs$a, 3))~(.(dataset_name))
        )
      )
    
    ggsave(
      paste0(
        "RandomPattern/",
        model,
        "/",
        toString(n),
        "/underlying_mean.png"
      ),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    ##### FLAT HYPERPARAMETERS #####
    a_hats <- c()
    for (j in MLEs$a_hat) {
      a_hats <- c(a_hats, unlist(j))
    }
    hyperData <- data.frame(a_hat = a_hats)
    
    HyperX<-(
      floor(min(hyperData$a_hat))*100
    ):(ceiling(max(hyperData))*100)/100
    
    ggplot() +
      geom_histogram(
        data = hyperData,
        mapping = aes(x = a_hats, y = after_stat(count) / max(count)),
        bins = 200
      ) +
      geom_line(
        data = data.frame(
          x = HyperX,
          y = dnorm(HyperX, MLEs$mu, MLEs$tau)
            / max(dnorm(HyperX, MLEs$mu, MLEs$tau))
        ),
        mapping = aes(x = x, y = y)
      ) +
      scale_x_continuous(name = bquote(hat(a))) +
      scale_y_continuous(name = "Normalized Counts") +
      ggtitle(
        bquote(
          hat(a)~"~ N"(.(round(MLEs$mu,3)),.(round(MLEs$tau,3)))
        )
      )
    
    ggsave(
      paste0(
        "RandomPattern/",
        model,
        "/",
        toString(n),
        "/hyperparameters.png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
  }
}
