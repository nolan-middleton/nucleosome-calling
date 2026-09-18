#########################################################################
# This R script generates plots for the maximum likelihood estimates    #
# for visualization and evaluation.                                     #
#########################################################################

##### SETUP #####

# Imports
require(ggplot2)
require(jsonlite)

# Directories
setwd(
  "C:\\Users\\nolan\\Desktop\\MiddletonPhDProject\\NucleosomeCalling"
)

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

smooth_relu <- function(x) {
  return(x*exp(x)/(1+exp(x)))
}
std_relations <- c(
  function(MLEs, means) {
    return(MLEs$k*means + MLEs$l)
  },
  function(MLEs, means) {
    return(MLEs$k*means + MLEs$l)
  },
  function(MLEs, means) {
    return(MLEs$k*(means - MLEs$m)^2 + MLEs$l)
  },
  function(MLEs, means) {
    return(MLEs$k*smooth_relu(MLEs$n*(means - MLEs$m)) + MLEs$l)
  },
  function(MLEs, means) {
    return(MLEs$k*smooth_relu(MLEs$n*(means - MLEs$m)) + MLEs$l)
  },
  function(MLEs, means) {
    return(MLEs$k*means + MLEs$l)
  },
  function(MLEs, means) {
    return(MLEs$k*means + MLEs$l)
  }
)

smooth_x <- -730:720 / 10 + 0.5
X <- -73:72 + 0.5

##### MAIN LOOP #####
for (i in 1:length(models)) {
  ##### LOAD DATA #####
  model <- models[i]
  dataname <- datanames[[i]]
  dataset_name <- dataset_names[i]
  print(paste0("> ", toString(model), "..."))
  
  MLE_dir <- paste0("Nucleosomepattern/", model, "/models")
  
  data <- read.delim(
    paste0("NucleosomePattern/", model, "/damagePattern.tsv"),
    header = FALSE,
    sep = "\t"
  )
  
  minus_values <- data.frame(
    chrom = character(0)
  )
  minus_selected <- data.frame(
    chrom = character(0)
  )
  plus_values <- data.frame(
    chrom = character(0)
  )
  plus_selected <- data.frame(
    chrom = character(0)
  )
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
  
  thisDir <- paste0("NucleosomePattern/", model, "/minus_vals")
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
  
  thisDir <- paste0("NucleosomePattern/", model, "/plus_vals")
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
  
  thisDir <- paste0("NucleosomePattern/", model, "/minus_selected")
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
  
  thisDir <- paste0("NucleosomePattern/", model, "/plus_selected")
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
  
  flat_MLEs <- read_json(paste0(MLE_dir, "/flat/MLE_values.json"))
  linear_MLEs <- read_json(paste0(MLE_dir, "/linear/MLE_values.json"))
  quadratic_MLEs<-read_json(paste0(MLE_dir,"/quadratic/MLE_values.json"))
  
  flat <- rep(flat_MLEs$a, length(X))
  linear <- linear_MLEs$a + linear_MLEs$b*X
  quadratic <- quadratic_MLEs$a + quadratic_MLEs$p*(X-quadratic_MLEs$h)^2
  
  trendData <- data.frame(
    x = rep(X, 4),
    y = c(
      unlist(data[2,]),
      c(flat, linear, quadratic) + rep(unlist(data[1,]), 3)
    ),
    model = rep(c("Base", "Flat", "Linear", "Quadratic"), each=length(X))
  )
  
  ##### FLAT GRAPH #####
  ggplot() + 
    geom_line(
      data = trendData[trendData$model %in% c("Base", "Flat"),],
      mapping = aes(x = x, y = y, color = model)
    ) +
    geom_line(
      data = data.frame(
        x = smooth_x,
        y = flat_MLEs$a
      ),
      mapping = aes(x = x, y = y),
      color = colour_scale[2]
    ) +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(
      bquote(
        "Mean ~ "*.(round(flat_MLEs$a, 3))~(.(dataset_name))
      )
    )
  
  ggsave(
    paste0(MLE_dir, "/flat/underlying_mean.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### LINEAR GRAPH #####
  ggplot() + 
    geom_line(
      data = trendData[trendData$model %in% c("Base", "Linear"),],
      mapping = aes(x = x, y = y, color = model)
    ) +
    geom_line(
      data = data.frame(
        x = smooth_x,
        y = linear_MLEs$a + linear_MLEs$b * smooth_x
      ),
      mapping = aes(x = x, y = y),
      color = colour_scale[2]
    ) +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(
      bquote(
        "Mean ~ "*.(
          round(linear_MLEs$a,3)
        )*"+"*.(
          round(linear_MLEs$b,3)
        )*"*i"~(.(dataset_name))
      )
    )
  
  ggsave(
    paste0(MLE_dir, "/linear/underlying_mean.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### QUADRATIC GRAPH #####
  ggplot() + 
    geom_line(
      data = trendData[trendData$model %in% c("Base", "Quadratic"),],
      mapping = aes(x = x, y = y, color = model)
    ) +
    geom_line(
      data = data.frame(
        x = smooth_x,
        y=quadratic_MLEs$a+quadratic_MLEs$p*(smooth_x-quadratic_MLEs$h)^2
      ),
      mapping = aes(x = x, y = y),
      color = colour_scale[2]
    ) +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(
      bquote(
        "Mean ~ "*.(
          round(quadratic_MLEs$a,3)
        )*"+"*.(round(quadratic_MLEs$p,3))(
          "i-"*.(round(quadratic_MLEs$h,3))
        )^2~(.(dataset_name))
      )
    )
  
  ggsave(
    paste0(MLE_dir, "/quadratic/underlying_mean.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### COMBINED GRAPH #####
  ggplot() + 
    geom_line(
      data = trendData,
      mapping = aes(x = x, y = y, color = model)
    ) +
    scale_y_continuous(name = bquote("Mean"~.(dataname))) +
    scale_x_continuous(
      name = "Inbetween Position",
      breaks = c(seq(-73,0,5),0,seq(3,73,5)),
      minor_breaks = NULL
    ) +
    scale_color_discrete(name = "Model", type = colour_scale) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    ggtitle(paste0("Model Comparison (", dataset_name, ")"))
  
  ggsave(
    paste0(MLE_dir, "/overall_mean.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### FLAT HYPERPARAMETERS #####
  a_hats <- c()
  for (j in flat_MLEs$a_hat) {
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
        y = dnorm(HyperX, flat_MLEs$mu, flat_MLEs$tau)
          / max(dnorm(HyperX, flat_MLEs$mu, flat_MLEs$tau))
      ),
      mapping = aes(x = x, y = y)
    ) +
    scale_x_continuous(name = bquote(hat(a))) +
    scale_y_continuous(name = "Normalized Counts") +
    ggtitle(
      bquote(
        hat(a)~"~ N"(.(round(flat_MLEs$mu,3)),.(round(flat_MLEs$tau,3)))
      )
    )
  
  ggsave(
    paste0(MLE_dir, "/flat/hyperparameters.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### LINEAR HYPERPARAMETERS #####
  a_hats <- c()
  for (j in linear_MLEs$a_hat) {
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
        y = dnorm(HyperX, flat_MLEs$mu, flat_MLEs$tau)
        / max(dnorm(HyperX, flat_MLEs$mu, flat_MLEs$tau))
      ),
      mapping = aes(x = x, y = y)
    ) +
    scale_x_continuous(name = bquote(hat(a))) +
    scale_y_continuous(name = "Normalized Counts") +
    ggtitle(
      bquote(
        hat(a)~"~ N"(.(round(flat_MLEs$mu,3)),.(round(flat_MLEs$tau,3)))
      )
    )
  
  ggsave(
    paste0(MLE_dir, "/linear/hyperparameters.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### QUADRATIC HYPERPARAMETERS #####
  a_hats <- c()
  for (j in quadratic_MLEs$a_hat) {
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
        y = dnorm(HyperX, flat_MLEs$mu, flat_MLEs$tau)
        / max(dnorm(HyperX, flat_MLEs$mu, flat_MLEs$tau))
      ),
      mapping = aes(x = x, y = y)
    ) +
    scale_x_continuous(name = bquote(hat(a))) +
    scale_y_continuous(name = "Normalized Counts") +
    ggtitle(
      bquote(
        hat(a)~"~ N"(.(round(flat_MLEs$mu,3)),.(round(flat_MLEs$tau,3)))
      )
    )
  
  ggsave(
    paste0(MLE_dir, "/quadratic/hyperparameters.png"),
    width = 2000,
    height = 1000,
    units = "px"
  )
  
  ##### POSITION PLOTS #####
  if (!dir.exists(paste0(MLE_dir, "/flat/positions"))) {
    dir.create(paste0(MLE_dir, "/flat/positions"))
  }
  if (!dir.exists(paste0(MLE_dir, "/linear/positions"))) {
    dir.create(paste0(MLE_dir, "/linear/positions"))
  }
  if (!dir.exists(paste0(MLE_dir, "/quadratic/positions"))) {
    dir.create(paste0(MLE_dir, "/quadratic/positions"))
  }
  for (j in 1:ncol(data)) {
    print(paste0("> Position ", toString(j)))
    
    minus_vals <- minus_values[[paste0("V", toString(146 - j + 1))]]
    minus_dipys <- minus_selected[[paste0("V", toString(146 - j + 1))]]
    plus_vals <- plus_values[[paste0("V", toString(j))]]
    plus_dipys <- plus_selected[[paste0("V", toString(j))]]
    
    plotData <- data.frame(
      x = c(minus_vals[minus_dipys], plus_vals[plus_dipys])
    )
    posX <- (floor(min(plotData$x))*10):(ceiling(max(plotData$x))*10)/10
    
    means <-  flat_MLEs$a  + unlist(data[1,])[j]
    ggplot() +
      geom_histogram(
        data = plotData,
        mapping = aes(x = x, y = after_stat(count) / max(count)),
        bins = 100
      ) +
      geom_line(
        data = data.frame(
          x = posX,
          y = dnorm(
            posX,
            means,
            std_relations[[i]](flat_MLEs, means)
          )
          / max(
            dnorm(posX, means, std_relations[[i]](flat_MLEs, means))
          )
        ),
        mapping = aes(x = x, y = y)
      ) +
      scale_x_continuous(name = bquote(hat(a))) +
      scale_y_continuous(name = "Normalized Counts") +
      ggtitle(paste0("Position ", toString(j - 73.5)))
    
    ggsave(
      paste0(MLE_dir, "/flat/positions/", toString(j), ".png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    means <- linear_MLEs$a + linear_MLEs$b*(j-73.5) + unlist(data[1,])[j]
    ggplot() +
      geom_histogram(
        data = plotData,
        mapping = aes(x = x, y = after_stat(count) / max(count)),
        bins = 100
      ) +
      geom_line(
        data = data.frame(
          x = posX,
          y = dnorm(
            posX,
            means,
            std_relations[[i]](linear_MLEs, means)
          )
          / max(
            dnorm(posX, means, std_relations[[i]](linear_MLEs, means))
          )
        ),
        mapping = aes(x = x, y = y)
      ) +
      scale_x_continuous(name = bquote(hat(a))) +
      scale_y_continuous(name = "Normalized Counts") +
      ggtitle(paste0("Position ", toString(j - 73.5)))
    
    ggsave(
      paste0(MLE_dir, "/linear/positions/", toString(j), ".png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
    
    means <- quadratic_MLEs$a + quadratic_MLEs$p*(
      (j-73.5) - quadratic_MLEs$h
    )^2 + unlist(data[1,])[j]
    ggplot() +
      geom_histogram(
        data = plotData,
        mapping = aes(x = x, y = after_stat(count) / max(count)),
        bins = 100
      ) +
      geom_line(
        data = data.frame(
          x = posX,
          y = dnorm(
            posX,
            means,
            std_relations[[i]](quadratic_MLEs, means)
          )
          / max(
            dnorm(posX, means, std_relations[[i]](quadratic_MLEs, means))
          )
        ),
        mapping = aes(x = x, y = y)
      ) +
      scale_x_continuous(name = bquote(hat(a))) +
      scale_y_continuous(name = "Normalized Counts") +
      ggtitle(paste0("Position ", toString(j - 73.5)))
    
    ggsave(
      paste0(MLE_dir, "/quadratic/positions/", toString(j), ".png"),
      width = 2000,
      height = 1000,
      units = "px"
    )
  }
}
