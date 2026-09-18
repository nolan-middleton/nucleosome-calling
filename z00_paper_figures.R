########################################################################
# This R script compiles all the data for the plots used in the paper. #
########################################################################

##### SETUP #####
require(ggplot2)
require(ggforce)
require(jsonlite)
require(stringr)
require(rstudioapi)
require(grid)
require(lomb)
require(rlang)
require(eulerr)
require(ggseqlogo)
require(ggtext)

# Directories
setwd(dirname(getActiveDocumentContext()$path))

if (!dir.exists("PaperFigures")) {
  dir.create("PaperFigures")
}

if (!dir.exists("PaperFigures/FigureData")) {
  dir.create("PaperFigures/FigureData")
}

tinyTxtSize <- 12
smallTxtSize <- 14
largeTxtSize <- 16

# Plotting
plotTheme <- theme(
  plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
  plot.title = element_text(size = largeTxtSize, face = "bold", hjust = 0.5),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line(linewidth = 1.5, lineend = "square"),
  axis.ticks = element_line(linewidth = 1.5, color = "black"),
  axis.ticks.length = unit(7.5, "pt"),
  axis.title = element_text(face = "bold", size = smallTxtSize),
  axis.text = element_text(face = "bold", size = tinyTxtSize, color = "black"),
  legend.text = element_text(face="bold", size=tinyTxtSize, color="black"),
  legend.position = "inside",
  legend.position.inside = c(0.025,1),
  legend.justification = c("left", "top"),
  legend.direction = "vertical",
  legend.margin = margin()
)

yeast_bg <- c(0.31, 0.31, 0.19, 0.19)
names(yeast_bg) <- c("A", "T", "C", "G")

plot_sequence_logo <- function(
    data,
    bg_freqs = NULL,
    pseudocount = NULL,
    alpha = 0.05,
    doAnnotate = TRUE,
    method = "bits",
    seq_type = "auto",
    namespace = NULL,
    ...
) {
  ##### Defaults #####
  if (method != "custom") {
    #### Determine Sequence Type and Default Pseudocount #####
    # If we don't supply a namespace to ggseqlogo, assume DNA
    if (seq_type == "dna") {
      letters <- c("A", "C", "T", "G")
      eps <- 4
    } else if (seq_type == "rna") {
      letters <- c("A", "C", "U", "G")
      eps <- 4
    } else if (seq_type == "aa") {
      letters <- c("A", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                   "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W")
      eps <- 10
    } else if (seq_type == "auto" && is.null(namespace)) {
      if (is.matrix(data)) {
        letters <- rownames(data)
      } else {
        letters <- unique(unlist(strsplit(data, "")))
      }
      
      if ((length(letters) <= 4) &&
          !(FALSE %in% (letters %in% c("A","C","G","T")))) {
        letters <- c("A", "C", "T", "G")
        eps <- 4
      } else if ((length(letters) <= 4) &&
                 !(FALSE %in% (letters %in% c("A","C","G","U")))) {
        letters <- c("A", "C", "U", "G")
        eps <- 4
      } else if ((length(letters) <= 20) &&
                 !(FALSE %in% (letters %in%
                               c("A","C","D","E","F","G","H","I","J","K",
                                 "L","M","N","P","Q","R","S","T","V","W")))) {
        letters <- c("A", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                     "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W")
        eps <- 10
      } else {
        eps <- 1
      }
    } else if (!is.null(namespace)) {
      letters <- namespace
      eps <- 1
    }
    
    ##### Get Data Matrix #####
    # We're doing this manually, so I'm declaring the heights myself. I need the
    # data to be in a matrix format.
    if (!is.matrix(data)) {
      L <- strsplit(data, "")
      M <- matrix(unlist(L), nrow = length(L), byrow = TRUE)
      
      letter_freqs = matrix(nrow = length(letters), ncol = ncol(M))
      rownames(letter_freqs) <- letters
      colnames(letter_freqs) <- unlist(lapply(1:ncol(letter_freqs), toString))
      for (letter in letters) {
        letter_freqs[letter,] <- colSums(M == letter)
      }
    } else {
      letter_freqs <- data
      letters <- rownames(data)
    }
    
    ##### Get Pseudocount and Background Frequencies #####
    if (!is.null(pseudocount)) {
      eps <- rep(pseudocount, length(letters))
    } else {
      eps <- rep(eps, length(letters))
    }
    names(eps) <- letters
    
    # If we're not supplying a background model, just use a uniform distribution
    if (is.null(bg_freqs)) {
      BG <- rep(1/nrow(letter_freqs), nrow(letter_freqs))
      names(BG) <- letters
    } else {
      BG <- bg_freqs
      # If a background model is applied, scale the pseudocount
      for (letter in letters) {
        eps[letter] <- eps[letter] * BG[letter]
      }
    }
    
    ##### Get Logo Height #####
    # Apply pseudocount
    for (letter in letters) {
      letter_freqs[letter,] <- letter_freqs[letter,] + eps[letter]
    }
    
    # Get letter probabilities
    N <- colSums(letter_freqs)
    P <- letter_freqs / N
    
    # Set up letter height matrix
    letter_heights <- matrix(nrow=nrow(letter_freqs), ncol=ncol(letter_freqs))
    rownames(letter_heights) <- letters
    colnames(letter_heights) <- colnames(letter_freqs)
    
    # Find heights
    if (method == "bits") {
      # KL Divergence
      height <- rep(0, ncol(letter_freqs))
      for (letter in letters) {
        height <- height + P[letter,]*log2(P[letter,]/BG[letter])
      }
      
      if (is.null(bg_freqs)) {
        yname <- "Information Content (bits)"
        ylims <- c(0, log2(length(letters)))
      } else {
        yname <- "Relative Entropy (bits)"
        ylims <- c(0, -log2(min(BG)))
      }
      
      for (letter in letters) {
        letter_heights[letter,] <- height * P[letter,]
      }
    } else if (method == "probability") {
      # Probability logo
      yname <- "Probability"
      ylims <- c(0,1)
      
      for (letter in letters) {
        letter_heights[letter,] <- P[letter,]
      }
    } else if (method == "pLogo") {
      # pLogo (see doi.org/10.1038/nmeth.2646)
      for (letter in letters) {
        letter_heights[letter,] <- stats::pbinom(
          letter_freqs[letter,],
          N,
          BG[letter],
          log.p = TRUE
        ) - stats::pbinom(
          letter_freqs[letter,] - 1,
          N,
          BG[letter],
          lower.tail = FALSE,
          log.p = TRUE
        )
      }
      letter_heights <- letter_heights / log(10) # Change base
      
      alpha_prime <- alpha / (ncol(letter_heights)*length(letters))
      lineHeight <- log10(alpha_prime / (1 - alpha_prime))
      
      annotationData <- data.frame(
        x = c(-Inf, Inf, -Inf, Inf),
        y = c(lineHeight, lineHeight, -lineHeight, -lineHeight),
        group = c(0, 0, 1, 1)
      )
      
      ylims <- NULL
      yname <- "Log Odds"
    }
    
    ##### Get Plot #####
    if (method != "pLogo") {
      suppressMessages(
        result <- ggplot2::ggplot() +
          ggseqlogo::geom_logo(
            data = letter_heights,
            method = "custom",
            seq_type = seq_type,
            namespace = namespace,
            ...
          ) +
          ggplot2::scale_x_continuous(
            name = NULL,
            breaks = 1:ncol(letter_heights),
            labels = colnames(letter_heights),
            expand = c(0,0)
          ) +
          ggplot2::scale_y_continuous(
            name = yname,
            limits = ylims,
            expand = c(0,0)
          ) +
          plotTheme
      )
    } else {
      if (doAnnotate) {
        suppressMessages(
          result <- result <- ggplot2::ggplot() +
            ggplot2::geom_line(
              data = data.frame(x = c(-Inf, Inf), y = 0),
              mapping = ggplot2::aes(x = x, y = y),
              color = "black",
              linewidth = 1.5
            ) +
            ggplot2::geom_line(
              data = annotationData,
              mapping = ggplot2::aes(x = x, y = y, group = group),
              color = "red",
              linewidth = 0.75
            ) +
            ggseqlogo::geom_logo(
              data = letter_heights,
              method = "custom",
              seq_type = seq_type,
              namespace = namespace,
              ...
            ) +
            ggplot2::scale_x_continuous(
              name = NULL,
              breaks = 1:ncol(letter_heights),
              labels = colnames(letter_heights),
              expand = c(0,0)
            ) +
            ggplot2::scale_y_continuous(
              name = yname,
              limits = ylims,
              expand = c(0,0)
            ) +
            plotTheme +
            ggplot2::theme(
              axis.line.x = ggplot2::element_blank(),
              axis.ticks.x = ggplot2::element_blank()
            )
        )
      } else {
        suppressMessages(
          result <- result <- ggplot2::ggplot() +
            ggplot2::geom_line(
              data = data.frame(x = c(-Inf, Inf), y = 0),
              mapping = ggplot2::aes(x = x, y = y),
              color = "black",
              linewidth = 1.5
            ) +
            ggseqlogo::geom_logo(
              data = letter_heights,
              method = "custom",
              seq_type = seq_type,
              namespace = namespace,
              ...
            ) +
            ggplot2::scale_x_continuous(
              name = NULL,
              breaks = 1:ncol(letter_heights),
              labels = colnames(letter_heights),
              expand = c(0,0)
            ) +
            ggplot2::scale_y_continuous(
              name = yname,
              limits = ylims,
              expand = c(0,0)
            ) +
            plotTheme +
            ggplot2::theme(
              axis.line.x = ggplot2::element_blank(),
              axis.ticks.x = ggplot2::element_blank()
            )
        )
      }
    }
  } else {
    # With predetermined custom heights, we just apply the geom_logo. We won't
    # supply a y-axis label in this case because we have no idea what the
    # heights actually represent if the heights are predetermined custom.
    suppressMessages(
      result <- ggplot2::ggplot() +
        ggseqlogo::geom_logo(
          data = data,
          seq_type = seq_type,
          method = method,
          namespace = namespace,
          ...
        ) +
        plotTheme +
        ggplot2::scale_x_continuous(
          name = NULL,
          breaks = 1:ncol(data),
          labels = colnames(data),
          expand = c(0,0)
        ) +
        ggplot2::scale_y_continuous(expand = c(0,0))
    )
  }
  
  return(result)
}

# Reference TSS Datasets
R <- c(-500, 650) # The range of values around the TSS to plot
w <- 5 # The width of the window to average to Brogaard bars over
S = seq(R[1], R[2]-1, w)

brogaard_TSS <- read.table(
  "NucleosomePattern/TSS_reference_datasets/all_frequencies.tsv",
  header = FALSE,
  sep = "\t"
)
colnames(brogaard_TSS) <- c("position", "frequency")

brogaard_plotData <- data.frame(
  x = S,
  y = rep(0, length(S))
)
for (i in 1:length(S)) {
  brogaard_plotData$y[i] <- mean(
    brogaard_TSS$frequency[
      (brogaard_TSS$position >= S[i]) & (brogaard_TSS$position < S[i] + w)
    ]
  )
}
write.table(
  brogaard_plotData,
  "PaperFigures/FigureData/brogaard_TSS.tsv",
  quote = FALSE,
  row.names = FALSE
)

weiner_TSS <- read.table(
  "NucleosomePattern/TSS_reference_datasets/weiner_frequencies.tsv",
  header = FALSE,
  sep = "\t"
)
colnames(weiner_TSS) <- c("position", "frequency")

weiner_plotData <- data.frame(
  x = S,
  y = rep(0, length(S))
)
for (i in 1:length(S)) {
  weiner_plotData$y[i] <- mean(
    weiner_TSS$frequency[
      (weiner_TSS$position >= S[i]) & (weiner_TSS$position < S[i] + w)
    ]
  )
}
write.table(
  weiner_plotData,
  "PaperFigures/FigureData/weiner_TSS.tsv",
  quote = FALSE,
  row.names = FALSE
)
rm(i)

weiner_minorOutData <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/periodic_proportions.txt"
  ),
  sep = "\t",
  header = FALSE,
  comment.char = "#"
)
colnames(weiner_minorOutData) <- c("reference","algorithm","successes","total")
write.table(
  weiner_minorOutData,
  "PaperFigures/FigureData/weiner_minor_out.tsv",
  quote = FALSE,
  row.names = FALSE
)

weiner_minorInData <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/anti_periodic_proportions.txt"
  ),
  sep = "\t",
  header = FALSE,
  comment.char = "#"
)
colnames(weiner_minorInData) <- c("reference","algorithm","successes","total")

write.table(
  weiner_minorInData,
  "PaperFigures/FigureData/weiner_minor_in.tsv",
  quote = FALSE,
  row.names = FALSE
)

weiner_dist <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/all_Weiner_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(weiner_dist) <- c("position", "frequency")

write.table(
  weiner_dist,
  "PaperFigures/FigureData/weiner_dyad_distances.tsv",
  quote = FALSE,
  row.names = FALSE
)

weiner_offset <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/all_Weiner_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(weiner_offset) <- c("position", "frequency")

write.table(
  weiner_offset,
  "PaperFigures/FigureData/weiner_dyad_offsets.tsv",
  quote = FALSE,
  row.names = FALSE
)

occu <- rep(0, length(brogaard_TSS$position))
for (i in 1:length(brogaard_TSS$position)) {
  occu[i] <- sum(
    brogaard_TSS$frequency[
      abs(brogaard_TSS$position - brogaard_TSS$position[i]) <= 73
    ]
  )
}

brogaard_occupancyData <- data.frame(
  position = brogaard_TSS$position,
  occupancy = occu
)
write.table(
  brogaard_occupancyData,
  "PaperFigures/FigureData/brogaard_occupancy.tsv",
  quote = FALSE,
  row.names = FALSE
)

occu <- rep(0, length(weiner_TSS$position))
for (i in 1:length(weiner_TSS$position)) {
  occu[i] <- sum(
    weiner_TSS$frequency[
      abs(weiner_TSS$position - weiner_TSS$position[i]) <= 73
    ]
  )
}

weiner_occupancyData <- data.frame(
  position = weiner_TSS$position,
  occupancy = occu
)
write.table(
  weiner_occupancyData,
  "PaperFigures/FigureData/weiner_occupancy.tsv",
  quote = FALSE,
  row.names = FALSE
)
rm(occu, i)

weiner_readData <- read.delim(
  "NucleosomePattern/TSS_reference_datasets/read_coverage.tsv",
  sep = "\t",
  header = FALSE
)
colnames(weiner_readData) <- c("position", "coverage")

##### CPD-seq Schematic #####

# Nothing to generate here...

##### All Nucleosome log2 Ratio Pattern #####
pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,])[ncol(pattern):1], unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5, 0),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_pattern.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#-: ",
  toString(periodicity$period[periodicity$strand == "-"][1]),
  "\n#+: ",
  toString(periodicity$period[periodicity$strand == "+"][1])
)

write(string, "PaperFigures/FigureData/all_log2_pattern.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_log2_pattern.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(pattern, periodicity, plotData, annotationData, string, P, G)

##### All Nucleosome Periodogram of log2 Ratios #####
periodogram <- read.table(
  "NucleosomePattern_all/cl_ice_log2/periodogram.tsv",
  header = TRUE,
  sep = "\t"
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

annotationData <- data.frame(
  x = c(
    rep(periodicity$period[periodicity$strand == "-"][1], 2),
    rep(periodicity$period[periodicity$strand == "+"][1], 2)
  ),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("-","+"), each = 2)
)

labs <- rep("", 14)
labs[c(4,9,14)] = c(5,10,15)

periodicity_annotations <- c(
  as.character(round(annotationData$x[annotationData$group == "-"][1], 2)),
  as.character(round(annotationData$x[annotationData$group == "+"][1], 2))
)
for (i in 1:length(periodicity_annotations)) {
  item <- periodicity_annotations[i]
  if (str_length(unlist(strsplit(item, "[.]"))[2]) < 2) {
    periodicity_annotations[i] <- paste0(periodicity_annotations[i], "0")
  }
}

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group, linetype = group),
    linewidth = 0.75
  ) +
  geom_line(
    data = periodogram,
    mapping = aes(x = period, y = normalized_power, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (nt)",
    limits = c(2,15),
    breaks = 2:15,
    labels = labs,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,0.6),
    breaks = (0:6)/10,
    expand = c(0,0)
  ) +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand", "Plus Strand")
  ) +
  annotate(
    "text",
    x = 10.3,
    y = 0.5,
    label = paste0(periodicity_annotations[1], "nt"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "bottom",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.3,
    y = 0.5,
    label = paste0(periodicity_annotations[2], "nt"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "blue"
  ) +
  guides(linetype = guide_none())

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_periodogram.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#-: ",
  toString(periodicity$period[periodicity$strand == "-"][1]),
  "\n#+: ",
  toString(periodicity$period[periodicity$strand == "+"][1])
)

write(string, "PaperFigures/FigureData/all_log2_periodogram.tsv")

write.table(
  periodogram,
  file = "PaperFigures/FigureData/all_log2_periodogram.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(annotationData,G,P,periodicity,periodogram,labs,string,pattern,plotData)

##### TSS and Scatter Plot of Absolute Differences #####

# Load Data
chroms <- unlist(read.table("PreprocessedData/TSS/chromosomes.txt"))

TSS_data <- list()
N <- 0
abs_diff_minus <- list()
abs_diff_plus <- list()
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  TSS_data[[chrom]] <- read.table(
    paste0("PreprocessedData/TSS/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(TSS_data[[chrom]]) <- c(
    "position",
    "gene",
    "internal",
    "flag",
    "strand"
  )
  N <- N + nrow(TSS_data[[chrom]])
  
  abs_diff_minus[[chrom]] <- read.table(
    paste0("PreprocessedData/cl_ice_abs_diff_minus/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(abs_diff_minus[[chrom]]) <- c("position", "value")
  
  abs_diff_plus[[chrom]] <- read.table(
    paste0("PreprocessedData/cl_ice_abs_diff_plus/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(abs_diff_plus[[chrom]]) <- c("position", "value")
}

# Get TSS Spots
values <- matrix(nrow = N, ncol = R[2]-R[1])
inbetween_positions <- R[1]:(R[2]-1) + 0.5
I <- 1
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  for (i in 1:nrow(TSS_data[[chrom]])) {
    position <- TSS_data[[chrom]]$position[i]
    if (TSS_data[[chrom]]$strand[i] == "+") {
      theseAbsDiffs <- rbind(
        abs_diff_minus[[chrom]][
          (abs_diff_minus[[chrom]]$position >= position + R[1]) &
            (abs_diff_minus[[chrom]]$position <= position + R[2]),
        ],
        plus_vals <- abs_diff_plus[[chrom]][
          (abs_diff_plus[[chrom]]$position >= position + R[1]) &
            (abs_diff_plus[[chrom]]$position <= position + R[2]),
        ]
      )
      theseAbsDiffs <- rbind(
        theseAbsDiffs,
        data.frame(
          position = (position + inbetween_positions)[
            !((position+inbetween_positions) %in% theseAbsDiffs$position)
          ],
          value = NA
        )
      )
      
      theseAbsDiffs <- theseAbsDiffs[order(theseAbsDiffs$position),]
      
      values[I,] <- theseAbsDiffs$value
    } else {
      theseAbsDiffs <- rbind(
        abs_diff_minus[[chrom]][
          (abs_diff_minus[[chrom]]$position >= position - R[2]) &
            (abs_diff_minus[[chrom]]$position <= position - R[1]),
        ],
        plus_vals <- abs_diff_plus[[chrom]][
          (abs_diff_plus[[chrom]]$position >= position - R[2]) &
            (abs_diff_plus[[chrom]]$position <= position - R[1]),
        ]
      )
      theseAbsDiffs <- rbind(
        theseAbsDiffs,
        data.frame(
          position = (position - inbetween_positions)[
            !((position-inbetween_positions) %in% theseAbsDiffs$position)
          ],
          value = NA
        )
      )
      
      theseAbsDiffs <- theseAbsDiffs[
        order(theseAbsDiffs$position, decreasing = TRUE),
      ]
      
      values[I,] <- theseAbsDiffs$value
    }
    
    I <- I + 1
  }
}

plotData <- data.frame(
  position = inbetween_positions,
  abs_diff = colMeans(values, na.rm = TRUE)
)

regressionData <- data.frame(
  x = (
    brogaard_TSS$frequency[1:(nrow(brogaard_TSS)-1)] +
      brogaard_TSS$frequency[2:nrow(brogaard_TSS)]
  ) / 2 ,
  y = plotData$abs_diff
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(25, 41),
    breaks = c(25, 30, 35, 40),
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    x = 100,
    y = 27,
    size = smallTxtSize,
    size.unit = "pt",
    fontface = "bold"
  ) +
  plotTheme

P

png(
  "PaperFigures/abs_diff_scatter.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/abs_diff_scatter.tsv",
  quote = FALSE,
  row.names = FALSE
)

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = 25,
      ymax = y * 16/100 + 25
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(
      x = position,
      y = abs_diff
    ),
    linewidth = 0.75,
    color = "blue"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(25, 41),
    breaks = seq(25, 40, 5),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/16 - 25*100/16,
      name = "Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "blue")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) + annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    x = -100,
    y = 40,
    size = smallTxtSize,
    size.unit = "pt",
    fontface = "bold"
  )

P

png(
  "PaperFigures/abs_diff_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

data <- data.frame(
  abs_diff_position = c(inbetween_positions, ""),
  abs_diff_mean = c(plotData$abs_diff, "")
)

write.table(
  data,
  "PaperFigures/FigureData/abs_diff_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(abs_diff_minus, abs_diff_plus, data, P, plotData, plus_vals, theseAbsDiffs,
   TSS_data, values, chrom, chroms, i, I, inbetween_positions, N, position, C,
   model, regressionData, r, p, X)

##### All Nucleosome Absolute Difference Pattern #####
pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_abs_diff/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,])[ncol(pattern):1], unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

model <- lm(
  values ~ position+I(position^2),
  data = plotData[plotData$strand == "-",]
)

annotationData <- data.frame(
  x = seq(min(plotData$position), max(plotData$position), 0.1)
)
annotationData$y <- model$coefficients[1] +
  model$coefficients[2]*annotationData$x +
  model$coefficients[3]*annotationData$x^2

maxes = plotData[plotData$strand == "+",][
  order(plotData$values[plotData$strand == "+"], decreasing = TRUE)[1:2],
]

P <- ggplot() +
  annotate(
    "line",
    x = rep(maxes$position, 2),
    y = rep(c(-Inf, Inf), each = 2),
    group = rep(c(1,2), each = 2),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y),
    linewidth = 1.5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = c(-60 + 20*(0:6), maxes$position),
    labels = c(-60 + 20*(0:6), maxes$position),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(30, 46),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_pattern.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_abs_diff_pattern.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(pattern, G, P, plotData, model, annotationData, maxes)

##### Schematic of Algorithm #####

# Nothing here...

##### Average Score vs. Relative Dyad Position #####
periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

# Load Scores
chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "genome_scores/chromosomes.txt"
  ),
  header = FALSE
)$V1
scores <- list()
shuffled <- list()
for (chrom in chroms) {
  scores[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/genome_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[chrom]]) <- c("position", "score")
  
  shuffled[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/shuffled_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(shuffled[[chrom]]) <- c("position", "score")
}

chroms <- read.table(
  "PreprocessedData/all_nucleosomes/chromosomes.txt",
  header = FALSE
)$V1
nucleosomes <- list()
for (chrom in chroms) {
  nucleosomes[[chrom]] <- read.table(
    paste0("PreprocessedData/all_nucleosomes/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(nucleosomes[[chrom]]) <- c("start", "stop", "dyad")
}

W <- c(-100,100) # The window width to graph around the true dyad

total <- rep(0, W[2]-W[1] + 1)
T_shuff <- rep(0, W[2]-W[1] + 1)
n <- 0
m <- 0
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  for (i in 1:nrow(nucleosomes[[chrom]])) {
    dyad <- nucleosomes[[chrom]]$dyad[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      total <- total + theseScores
      n <- n + 1
    }
    
    theseShuffled <- shuffled[[chrom]]$score[
      (shuffled[[chrom]]$position >= dyad + W[1]) &
        (shuffled[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseShuffled) == W[2] - W[1] + 1) {
      T_shuff <- T_shuff + theseShuffled
      m <- m + 1
    }
  }
}

plotData <- data.frame(
  relative_position = rep(W[1]:W[2], 2),
  average_score = c(total / n, T_shuff / m),
  dataset = rep(c("score", "shuffled"), each = W[2]-W[1]+1)
)

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linetype = "dotted",
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = relative_position, y = average_score, color = dataset),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-1.75,0),
    breaks = seq(-1.5,0,0.5),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("black", "#777777"),
    labels = c("Score", "Shuffled")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/score_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))

write(string, "PaperFigures/FigureData/score_dyad_position.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/score_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(chrom, scores, chroms, periodicity, nucleosomes, W, total, n, plotData,
   annotationData, P, G, string, theseScores, dyad, i, shuffled, T_shuff, m)

##### Score TSS, Occupancy, and Scatter Plot for Brogaard and Weiner #####
score_TSS <- read.table(
  "NucleosomePattern/composite_models/composite_ice/TSS_evals/mean_score.tsv",
  header = FALSE,
  sep = "\t"
)
colnames(score_TSS) <- c("position", "mean_score")

plotData <- data.frame(
  position = score_TSS$position,
  mean_score = score_TSS$mean_score
)

# > Brogaard et al. dyad frequency vs. Scores ----
regressionData <- data.frame(
  x = brogaard_TSS$frequency[brogaard_TSS$position > 0],
  y = score_TSS$mean_score[score_TSS$position > 0]
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-1.3, -0.7),
    breaks = (-13):(-7)*0.1,
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    x = 100,
    y = -1.2,
    size = smallTxtSize,
    size.unit = "pt",
    fontface = "bold"
  ) +
  ggtitle(
    bquote(bold(underline("Post-TSS")))
  ) +
  plotTheme

P

png(
  "PaperFigures/score_scatter.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/score_scatter.tsv",
  quote = FALSE,
  row.names = FALSE
)

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = -2.6,
      ymax = y * 1.9/100 - 2.6
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(
      x = position,
      y = mean_score
    ),
    linewidth = 0.75,
    color = "red"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.5, -0.7, 0.5),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/1.9 + 2.6*100/1.9,
      name = "Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "red")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 500,
    y = -0.8,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/score_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Brogaard et al. occupancy vs. scores ----
in_plot <- (
  (brogaard_occupancyData$position >= R[1] + 73) &
    (brogaard_occupancyData$position <= R[2] - 73)
)
regressionData <- data.frame(
  x = brogaard_occupancyData$occupancy[in_plot],
  y = score_TSS$mean_score[in_plot]
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_occupancyData$occupancy), 1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Nucleosome Occupancy",
    limits = c(1800,5000),
    breaks = seq(1800, 5000, 400),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.6, -0.7, 0.3),
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 4400,
    y = -2.4,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme +
  theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt"))

P

png(
  "PaperFigures/score_scatter_occupancy.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/score_scatter_occupancy.tsv",
  quote = FALSE,
  row.names = FALSE
)

P <- ggplot() +
  geom_rect(
    data = brogaard_occupancyData,
    mapping = aes(
      xmin = position - 0.5,
      xmax = position + 0.5,
      ymin = -2.6,
      ymax = occupancy * 1.9/5000 - 2.6
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(
      x = position,
      y = mean_score
    ),
    linewidth = 0.75,
    color = "red"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = c(R[1] + 73, R[2] - 73),
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.5, -0.7, 0.5),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 5000/1.9 + 2.6*5000/1.9,
      name = "Nucleosome Occupancy"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "red")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = -100,
    y = -0.8,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/score_TSS_occupancy.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  "PaperFigures/FigureData/score_TSS_occupancy.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# > Weiner et al. dyad frequency vs. scores ----
regressionData <- data.frame(
  x = weiner_TSS$frequency[weiner_TSS$position > 0],
  y = score_TSS$mean_score[weiner_TSS$position > 0]
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(weiner_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-1.3, -0.7),
    breaks = (-13):(-7)*0.1,
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 100,
    y = -1.2,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  ggtitle(
    bquote(bold(underline("Post-TSS")))
  ) +
  plotTheme

P

png(
  "PaperFigures/score_scatter_weiner.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/score_scatter_weiner.tsv",
  quote = FALSE,
  row.names = FALSE
)

P <- ggplot() +
  geom_rect(
    data = weiner_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = -2.6,
      ymax = y * 1.9/120 - 2.6
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(
      x = position,
      y = mean_score
    ),
    linewidth = 0.75,
    color = "red"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.5, -0.7, 0.5),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 120/1.9 + 2.6*120/1.9,
      name = "Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "red")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 500,
    y = -0.8,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/score_TSS_weiner.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Weiner et al. occupancy vs. scores ----
in_plot <- (
  (weiner_occupancyData$position >= R[1] + 73) &
    (weiner_occupancyData$position <= R[2] - 73)
)
regressionData <- data.frame(
  x = weiner_occupancyData$occupancy[in_plot],
  y = score_TSS$mean_score[in_plot]
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(weiner_occupancyData$occupancy), 1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Nucleosome Occupancy",
    limits = c(750, 5100),
    breaks = seq(750, 5100, 500),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.6, -0.7, 0.3),
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 4400,
    y = -2.4,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme +
  theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt"))

P

png(
  "PaperFigures/score_scatter_weiner_occupancy.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/score_scatter_weiner_occupancy.tsv",
  quote = FALSE,
  row.names = FALSE
)

P <- ggplot() +
  geom_rect(
    data = weiner_occupancyData,
    mapping = aes(
      xmin = position - 0.5,
      xmax = position + 0.5,
      ymin = -2.6,
      ymax = occupancy * 1.9/5100 - 2.6
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(
      x = position,
      y = mean_score
    ),
    linewidth = 0.75,
    color = "red"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = c(R[1] + 73, R[2] - 73),
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.5, -0.7, 0.5),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 5100/1.9 + 2.6*5100/1.9,
      name = "Nucleosome Occupancy"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "red")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = -100,
    y = -0.8,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/score_TSS_weiner_occupancy.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  "PaperFigures/FigureData/score_TSS_weiner_occupancy.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# > Weiner et al. reads ----
regressionData <- data.frame(
  x = weiner_readData$coverage,
  y = score_TSS$mean_score
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(weiner_readData$coverage), 100)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Read Coverage",
    limits = c(700000,5700000),
    breaks = seq(1000000, 5700000, 500000),
    labels = paste0(unlist(lapply((2:11)/2, format, nsmall=1)), "E6"),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.6, -0.7, 0.3),
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 5000000,
    y = -2.3,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme

P

png(
  "PaperFigures/score_reads.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/score_reads.tsv",
  quote = FALSE,
  row.names = FALSE
)

P <- ggplot() +
  geom_rect(
    data = weiner_readData,
    mapping = aes(
      xmin = position,
      xmax = position + w,
      ymin = -2.6,
      ymax = coverage * 1.9/5700000 - 2.6
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(
      x = position,
      y = mean_score
    ),
    linewidth = 0.75,
    color = "red"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6, -0.7),
    breaks = seq(-2.5, -0.7, 0.5),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 5700000/1.9 + 2.6*5700000/1.9,
      name = "Read Coverage",
      breaks = seq(1000000, 5700000, 1000000),
      labels = paste0(unlist(lapply(1:5, format, nsmall=1)), "E6")
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "red")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = -125,
    y = -0.8,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/score_TSS_reads.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

rm(P, plotData, score_TSS, regressionData, model, C, p, r, X, in_plot)

##### Greedy Placement Explanation #####
X <- -500:500 / 100
N <- length(X)
W <- c(10, pi, sqrt(2), 1, 3, 8, 1.0258)
A <- c(1, 0.9, 0.95, 1.25, 1.1, 0.6, 1.23)
B <- c(0.1, 2, 4, sqrt(3), exp(1), 27, 0)
M <- 7 # length(W)
f <- function(x, W, A, B) {
  return(sum(A*sin(W*x+B)))
}
plotData <- data.frame(
  x = X,
  y = unlist(lapply(X, f, W, A, B))
)
plotData$valid <- as.logical(
  replace(
    replace(plotData$y, plotData$y < -2, FALSE),
    plotData$y >= -2,
    TRUE
  )
)
plotData$group <- cumsum(
  c(TRUE, !((plotData$valid)[2:N] == (plotData$valid)[1:(N-1)]))
)

highPt <- plotData$x[which.max(plotData$y)]
ggplot(data = plotData) +
  geom_line(
    mapping = aes(x = x, y = y, linetype = valid, group = group),
    linewidth = 1
  ) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(-2, -2),
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    geom = "point",
    x = highPt,
    y = max(plotData$y),
    size = 2
  ) +
  annotate(
    geom = "text",
    x = highPt,
    y = max(plotData$y) + 0.15,
    hjust = "center",
    vjust = "bottom",
    label = "230bp",
    size = 7,
    fontface = "bold",
    color = "#FFFFFF"
  ) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.border = element_blank()
  ) +
  guides(linetype = guide_none()) +
  scale_y_continuous(name = "Score", expand = expansion(0.05, c(0,0.5))) +
  scale_x_continuous(
    name = "Genome Position (bp)",
    limits = c(-5,5),
    expand = c(0,0)
  )

ggsave(
  "PaperFigures/GreedyExplanation_1.png",
  width = 2000,
  height = 1000,
  units = "px"
)

width <- 1
plotData2 <- plotData
plotData2$valid[
  (plotData2$x >= highPt - width) & (plotData2$x <= highPt + width)
] <- FALSE
plotData2$group <- cumsum(
  c(TRUE, !((plotData2$valid)[2:N] == (plotData2$valid)[1:(N-1)]))
)

ggplot(data = plotData2) +
  geom_line(
    mapping = aes(x = x, y = y, linetype = valid, group = group),
    linewidth = 1
  ) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(-2, -2),
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    geom = "point",
    x = highPt,
    y = max(plotData2$y),
    size = 2
  ) +
  annotate(
    geom = "line",
    x = c(highPt - width, highPt + width),
    y = c(max(plotData2$y), max(plotData2$y)),
    linewidth = 1,
    arrow = arrow(
      angle = 20,
      length = unit(5, "pt"),
      ends = "both",
      type = "closed"
    )
  ) +
  annotate(
    geom = "line",
    x = c(highPt - width, highPt - width),
    y = c(-Inf, Inf),
    linewidth = 1,
    linetype = "dashed"
  ) +
  annotate(
    geom = "line",
    x = c(highPt + width, highPt + width),
    y = c(-Inf, Inf),
    linewidth = 1,
    linetype = "dashed"
  ) +
  annotate(
    geom = "text",
    x = highPt,
    y = max(plotData2$y) + 0.15,
    hjust = "center",
    vjust = "bottom",
    label = "235bp",
    size = 7,
    fontface = "bold"
  ) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.border = element_blank()
  ) +
  guides(linetype = guide_none()) +
  scale_y_continuous(name = "Score", expand = expansion(0.05, c(0,0.5))) +
  scale_x_continuous(
    name = "Genome Position (bp)",
    limits = c(-5,5),
    expand = c(0,0)
  )

ggsave(
  "PaperFigures/GreedyExplanation_2.png",
  width = 2000,
  height = 1000,
  units = "px"
)

rm(plotData, plotData2, A, B, highPt, M, N, W, width, X, f)

##### Greedy Nucleosome Position Relative to Reference Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_greedy_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = rep(-80:80, 2),
  frequency = c(
    data$frequency[data$position %in% -80:80],
    weiner_dist$frequency[weiner_dist$position %in% -80:80]
  ),
  dataset = rep(c("greedy", "weiner"), each = 161)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "black"),
    labels = c("Greedy", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1800),
    breaks = seq(0, 1800, 300),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/greedy_dyad_position.tsv")
write.table(
  data,
  file = "PaperFigures/FigureData/greedy_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string)

##### Greedy Nucleosome Offset Relative to Reference Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_greedy_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = rep(-100:100, 2),
  frequency = c(
    data$frequency[data$position %in% -100:100],
    weiner_offset$frequency[weiner_offset$position %in% -100:100]
  ),
  dataset = rep(c("greedy", "weiner"), each = 201)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "black"),
    labels = c("Greedy", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1800),
    breaks = seq(0, 1800, 300),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/greedy_dyad_offset.tsv")
write.table(
  data,
  file = "PaperFigures/FigureData/greedy_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string)

##### Greedy Nucleosome Placement TSS and Scatter Plot #####
greedy_TSS <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/TSS_evals",
    "/greedy_frequencies.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(greedy_TSS) <- c("position", "frequency")

regressionData <- data.frame(
  x = brogaard_TSS$frequency,
  y = greedy_TSS$frequency
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Reference Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Greedy Dyad Frequency",
    limits = c(0, 60),
    breaks = 0:6*10,
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 100,
    y = 10,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme

P

png(
  "PaperFigures/greedy_scatter.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/greedy_scatter.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  position = greedy_TSS$position,
  frequency = greedy_TSS$frequency
)

s <- 5 # Smoothing window radius
smooth_plotData = data.frame(
  position = (R[1]+s):(R[2]-s),
  frequency = rep(0, R[2]-R[1]-(2*s-1))
)

for (i in 1:nrow(smooth_plotData)) {
  smooth_plotData$frequency[i] <- mean(
    plotData$frequency[
      (plotData$position >= smooth_plotData$position[i]-s) &
        (plotData$position <= smooth_plotData$position[i]+s)
    ]
  )
}

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = 0,
      ymax = y * 60/100
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 0.75,
    color = "#FFAAAA"
  ) +
  geom_line(
    data = smooth_plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "red"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Greedy Dyad Frequency",
    limits = c(0, 60),
    breaks = 0:6*10,
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/60,
      name = "Reference Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "red")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = -100,
    y = 55,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/greedy_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

data <- plotData
data$smooth_frequency <- rep(NA, nrow(data))
data$smooth_frequency[
  data$position %in% smooth_plotData$position
] <- smooth_plotData$frequency

write.table(
  data,
  "PaperFigures/FigureData/greedy_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

rm(C, data, greedy_TSS, model, P, plotData, regressionData, smooth_plotData, i,
   p, r, s, X)
##### Minor-Out vs. Minor-In Positions Bar Graphs #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/periodic_proportions.txt"
  ),
  sep = "\t",
  header = FALSE
)

E_out <- as.numeric(unlist(strsplit(data$V1[1], "="))[2])
data <- data[2:nrow(data),]
colnames(data) <- c("reference", "algorithm", "successes", "total")
data <- data[data$algorithm %in% c("greedy", "viterbi"),]

outData <- data.frame(
  reference = data$reference,
  algorithm = data$algorithm,
  diff = data$successes / data$total - E_out
)
outData = rbind(
  outData,
  data.frame(
    reference = weiner_minorOutData$reference,
    algorithm = weiner_minorOutData$algorithm,
    diff = weiner_minorOutData$successes / weiner_minorOutData$total - E_out
  )
)
outData <- cbind(outData, position = rep("minor_out", nrow(outData)))

data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/anti_periodic_proportions.txt"
  ),
  sep = "\t",
  header = FALSE
)

E_in <- as.numeric(unlist(strsplit(data$V1[1], "="))[2])
data <- data[2:nrow(data),]
colnames(data) <- c("reference", "algorithm", "successes", "total")
data <- data[data$algorithm %in% c("greedy", "viterbi"),]

inData <- data.frame(
  reference = data$reference,
  algorithm = data$algorithm,
  diff = data$successes / data$total - E_in
)
inData = rbind(
  inData,
  data.frame(
    reference = weiner_minorInData$reference,
    algorithm = weiner_minorInData$algorithm,
    diff = weiner_minorInData$successes / weiner_minorInData$total - E_in
  )
)
inData <- cbind(inData, position = rep("minor_in", nrow(inData)))

plotData <- rbind(outData, inData)
plotData <- plotData[plotData$reference %in% c("all", "strong"),]
plotData <- cbind(
  plotData,
  signif = c("***","***","***","***","*","*","***","***","***","***","","")
)

# > Greedy All ----

P <- ggplot(
  data=plotData[(plotData$reference=="all")&(plotData$algorithm!="viterbi"),]
) +
  geom_col(
    mapping = aes(x = algorithm, y = diff, group = position, fill = position),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Greedy", "MNase-seq")) +
  scale_y_continuous(
    name = "Deviation from Expectation",
    limits = c(-0.2, 0.2),
    breaks = seq(-0.2, 0.2, 0.05),
    expand = c(0,0),
    labels = paste0(unlist(lapply(-4:4*5, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("Minor-In", "Minor-Out")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.015 - 0.0275*(position == "minor_in"),
      group = position,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width = 0.9),
    vjust = 0
  ) +
  ggtitle(bquote(underline(bold("All Nucleosomes")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 30)
  )

P

png(
  "PaperFigures/greedy_rotational_barplot_all.png",
  width = 1000,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()


# > Greedy Strong ----

P <- ggplot(
  data=plotData[(plotData$reference=="strong")&(plotData$algorithm!="viterbi"),]
) +
  geom_col(
    mapping = aes(x = algorithm, y = diff, group = position, fill = position),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Greedy", "MNase-seq")) +
  scale_y_continuous(
    name = "Deviation from Expectation",
    limits = c(-0.3, 0.3),
    breaks = seq(-0.3, 0.3, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(-3:3*10, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("Minor-In", "Minor-Out")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.025 - 0.035*(position == "minor_in"),
      group = position,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width = 0.9),
    vjust = 0
  ) +
  ggtitle(bquote(underline(bold("Strongly-Positioned")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 30)
  )

P

png(
  "PaperFigures/greedy_rotational_barplot_strong.png",
  width = 1000,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Viterbi All ----

P <- ggplot(
  data=plotData[(plotData$reference=="all")&(plotData$algorithm!="greedy"),]
) +
  geom_col(
    mapping = aes(x = algorithm, y = diff, group = position, fill = position),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Viterbi", "MNase-seq")) +
  scale_y_continuous(
    name = "Deviation from Expectation",
    limits = c(-0.2, 0.2),
    breaks = seq(-0.2, 0.2, 0.05),
    expand = c(0,0),
    labels = paste0(unlist(lapply(-4:4*5, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("Minor-In", "Minor-Out")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.015 - 0.0275*(position == "minor_in"),
      group = position,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width = 0.9),
    vjust = 0
  ) +
  ggtitle(bquote(underline(bold("All Nucleosomes")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 30)
  )

P

png(
  "PaperFigures/viterbi_rotational_barplot_all.png",
  width = 900,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Viterbi Strong ----

P <- ggplot(
  data=plotData[(plotData$reference=="strong")&(plotData$algorithm!="greedy"),]
) +
  geom_col(
    mapping = aes(x = algorithm, y = diff, group = position, fill = position),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Viterbi", "MNase-seq")) +
  scale_y_continuous(
    name = "Deviation from Expectation",
    limits = c(-0.3, 0.3),
    breaks = seq(-0.3, 0.3, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(-3:3*10, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("Minor-In", "Minor-Out")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.025 - 0.035*(position == "minor_in"),
      group = position,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width = 0.9),
    vjust = 0
  ) +
  ggtitle(bquote(underline(bold("Strongly-Positioned")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 30)
  )

P

png(
  "PaperFigures/viterbi_rotational_barplot_strong.png",
  width = 900,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Greedy Minor-Out ----

P <- ggplot(
  data=plotData[
    (plotData$algorithm!="viterbi")&(plotData$position == "minor_out"),
  ]
) +
  geom_col(
    mapping = aes(x=algorithm, y=diff, group=reference, fill=reference),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Greedy", "MNase-seq")) +
  scale_y_continuous(
    name = "Proportion of Dyad Positions",
    limits = c(0.1-E_out, 0.7-E_out),
    breaks = seq(0.1-E_out, 0.7-E_out, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(1:7*10, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("All Nucleosomes", "Strongly-Positioned\nNucleosomes")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.015,
      group = reference,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width=0.9),
    vjust = 0
  ) +
  annotate(
    geom = "text",
    label = "🡄",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.4,
    size = 5
  ) +
  annotate(
    geom = "text",
    label = "    Random",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.5
  ) +
  coord_cartesian(clip = "off") +
  ggtitle(bquote(underline(bold("Minor-Out")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 90),
    axis.text.x = element_text(margin = margin(t = -75, l = 0, r = 0, b = 75))
  ) +
  geom_errorbarh(
    data = data.frame(x = c("greedy", "Weiner")),
    mapping = aes(
      xmin = stage(start = x, after_scale = xmin - 0.4),
      xmax = stage(start = x, after_scale = xmax + 0.4),
      y = -0.04
    ),
    width = 0,
    linewidth = 1.25
  )

P

png(
  "PaperFigures/greedy_rotational_barplot_minor_out.png",
  width = 1250,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Greedy Minor-In ----

P <- ggplot(
  data=plotData[
    (plotData$algorithm!="viterbi")&(plotData$position == "minor_in"),
  ]
) +
  geom_col(
    mapping = aes(x=algorithm, y=diff, group=reference, fill=reference),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Greedy", "MNase-seq")) +
  scale_y_continuous(
    name = "Proportion of Dyad Positions",
    limits = c(0.1-E_in, 0.7-E_in),
    breaks = seq(0.1-E_in, 0.7-E_in, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(1:7*10, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("All Nucleosomes", "Strongly-Positioned\nNucleosomes")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.015,
      group = reference,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width=0.9),
    vjust = 1
  ) +
  annotate(
    geom = "text",
    label = "🡄",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.4,
    size = 5
  ) +
  annotate(
    geom = "text",
    label = "    Random",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.5
  ) +
  coord_cartesian(clip = "off") +
  ggtitle(bquote(underline(bold("Minor-In")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 90)
  ) +
  geom_errorbarh(
    data = data.frame(x = c("greedy", "Weiner")),
    mapping = aes(
      xmin = stage(start = x, after_scale = xmin - 0.4),
      xmax = stage(start = x, after_scale = xmax + 0.4),
      y = -Inf
    ),
    width = 0,
    linewidth = 1.25
  )

P

png(
  "PaperFigures/greedy_rotational_barplot_minor_in.png",
  width = 1250,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Viterbi Minor-Out ----

P <- ggplot(
  data=plotData[
    (plotData$algorithm!="greedy")&(plotData$position == "minor_out"),
  ]
) +
  geom_col(
    mapping = aes(x=algorithm, y=diff, group=reference, fill=reference),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Viterbi", "MNase-seq")) +
  scale_y_continuous(
    name = "Proportion of Dyad Positions",
    limits = c(0.1-E_out, 0.7-E_out),
    breaks = seq(0.1-E_out, 0.7-E_out, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(1:7*10, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("All Nucleosomes", "Strongly-Positioned\nNucleosomes")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.015,
      group = reference,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width=0.9),
    vjust = 0
  ) +
  annotate(
    geom = "text",
    label = "🡄",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.4,
    size = 5
  ) +
  annotate(
    geom = "text",
    label = "    Random",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.5
  ) +
  coord_cartesian(clip = "off") +
  ggtitle(bquote(underline(bold("Minor-Out")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 90),
    axis.text.x = element_text(margin = margin(t = -75, l = 0, r = 0, b = 75))
  ) +
  geom_errorbarh(
    data = data.frame(x = c("viterbi", "Weiner")),
    mapping = aes(
      xmin = stage(start = x, after_scale = xmin - 0.4),
      xmax = stage(start = x, after_scale = xmax + 0.4),
      y = -0.04
    ),
    width = 0,
    linewidth = 1.25
  )

P

png(
  "PaperFigures/viterbi_rotational_barplot_minor_out.png",
  width = 1250,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Viterbi Minor-In ----

P <- ggplot(
  data=plotData[
    (plotData$algorithm!="greedy")&(plotData$position == "minor_in"),
  ]
) +
  geom_col(
    mapping = aes(x=algorithm, y=diff, group=reference, fill=reference),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name = NULL, labels = c("Viterbi", "MNase-seq")) +
  scale_y_continuous(
    name = "Proportion of Dyad Positions",
    limits = c(0.1-E_in, 0.7-E_in),
    breaks = seq(0.1-E_in, 0.7-E_in, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(1:7*10, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("All Nucleosomes", "Strongly-Positioned\nNucleosomes")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  geom_text(
    mapping = aes(
      x = algorithm,
      y = diff - 0.015,
      group = reference,
      label = signif
    ),
    size = 8,
    position = position_dodge2(width=0.9),
    vjust = 1
  ) +
  annotate(
    geom = "text",
    label = "🡄",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.4,
    size = 5
  ) +
  annotate(
    geom = "text",
    label = "    Random",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.5
  ) +
  coord_cartesian(clip = "off") +
  ggtitle(bquote(underline(bold("Minor-In")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.65,1),
    plot.margin = margin(l = 10, t = 10, b = 10, r = 90)
  ) +
  geom_errorbarh(
    data = data.frame(x = c("viterbi", "Weiner")),
    mapping = aes(
      xmin = stage(start = x, after_scale = xmin - 0.4),
      xmax = stage(start = x, after_scale = xmax + 0.4),
      y = -Inf
    ),
    width = 0,
    linewidth = 1.25
  )

P

png(
  "PaperFigures/viterbi_rotational_barplot_minor_in.png",
  width = 1250,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/rotational_barplots.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(data, P, inData, outData, plotData, E_in, E_out)

##### Greedy vs. Reference vs. Scores IGV Track #####

# Nothing to generate here...

##### Explanation of Viterbi Placement #####
X <- -500:500 / 100
N <- length(X)
W <- c(10, pi, sqrt(2), 1, 3, 8, 1.0258)
A <- c(1, 0.9, 0.95, 1.25, 1.1, 0.6, 1.23)
B <- c(0.1, 2, 4, sqrt(3), exp(1), 27, 0)
M <- 7 # length(W)
f <- function(x, W, A, B) {
  return(sum(A*sin(W*x+B)))
}
plotData <- data.frame(
  x = X,
  y = unlist(lapply(X, f, W, A, B))
)

ggplot(data = plotData) +
  geom_line(mapping = aes(x = x, y = y), linewidth = 1) +
  annotate(
    geom = "line",
    x = c(0, 0),
    y = c(-Inf, Inf),
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(-0.2, -0.2),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(-4, -4),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 15, face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.border = element_blank()
  ) +
  guides(linetype = guide_none()) +
  scale_y_continuous(name = "Score") +
  scale_x_continuous(
    name = "Genome Position (bp)",
    limits = c(-5,5),
    expand = c(0,0)
  )

ggsave(
  "PaperFigures/ViterbiExplanation.png",
  width = 2000,
  height = 1000,
  units = "px"
)

rm(X, plotData, A, B, M, N, W, f)

##### Viterbi Nucleosome Positions Relative to Reference Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = rep(-80:80, 2),
  frequency = c(
    data$frequency[data$position %in% -80:80],
    weiner_dist$frequency[weiner_dist$position %in% -80:80]
  ),
  dataset = rep(c("viterbi", "weiner"), each = 161)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("#990099", "black"),
    labels = c("Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_dyad_positions.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/viterbi_dyad_positions.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_dyad_positions.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string)

##### Viterbi Nucleosome Offset Relative to Reference Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = rep(-100:100, 2),
  frequency = c(
    data$frequency[data$position %in% -100:100],
    weiner_offset$frequency[weiner_offset$position %in% -100:100]
  ),
  dataset = rep(c("viterbi", "weiner"), each = 201)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("#990099", "black"),
    labels = c("Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/viterbi_dyad_offset.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string)

##### Viterbi Placement TSS and Scatter Plot #####
viterbi_TSS <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/TSS_evals",
    "/viterbi_frequencies.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(viterbi_TSS) <- c("position", "frequency")

regressionData <- data.frame(
  x = brogaard_TSS$frequency,
  y = viterbi_TSS$frequency
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Reference Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Viterbi Dyad Frequency",
    limits = c(0, 65),
    breaks = 0:6*10,
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 100,
    y = 10,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme

P

png(
  "PaperFigures/viterbi_scatter.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/viterbi_scatter.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  position = viterbi_TSS$position,
  frequency = viterbi_TSS$frequency
)

s <- 5 # Smoothing window radius
smooth_plotData = data.frame(
  position = (R[1]+s):(R[2]-s),
  frequency = rep(0, R[2]-R[1]-(2*s-1))
)

for (i in 1:nrow(smooth_plotData)) {
  smooth_plotData$frequency[i] <- mean(
    plotData$frequency[
      (plotData$position >= smooth_plotData$position[i]-s) &
        (plotData$position <= smooth_plotData$position[i]+s)
    ]
  )
}

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = 0,
      ymax = y * 60/100
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 0.75,
    color = "#FFB0FF"
  ) +
  geom_line(
    data = smooth_plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "#990099"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Viterbi Dyad Frequency",
    limits = c(0, 60),
    breaks = 0:6*10,
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/60,
      name = "Reference Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "#990099")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = -100,
    y = 55,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/viterbi_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

data <- plotData
data$smooth_frequency <- rep(NA, nrow(data))
data$smooth_frequency[
  data$position %in% smooth_plotData$position
] <- smooth_plotData$frequency

write.table(
  data,
  "PaperFigures/FigureData/viterbi_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

rm(C, data, viterbi_TSS, model, P, plotData, regressionData, smooth_plotData, i,
   p, r, s, X)

##### Greedy vs. Viterbi vs. Brogaard vs. Weiner IGV Track #####

# Nothing to generate here...


##### Explanation of Viterbi Placement with 20bp Overlap #####
X <- -500:500 / 100
N <- length(X)
W <- c(10, pi, sqrt(2), 1, 3, 8, 1.0258)
A <- c(1, 0.9, 0.95, 1.25, 1.1, 0.6, 1.23)
B <- c(0.1, 2, 4, sqrt(3), exp(1), 27, 0)
M <- 7 # length(W)
f <- function(x, W, A, B) {
  return(sum(A*sin(W*x+B)))
}
plotData <- data.frame(
  x = X,
  y = unlist(lapply(X, f, W, A, B))
)

ggplot(data = plotData) +
  geom_line(mapping = aes(x = x, y = y), linewidth = 1) +
  annotate(
    geom = "line",
    x = c(0, 0),
    y = c(-Inf, Inf),
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(-0.2, -0.2),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(-3.6, -3.6),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 15, face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.border = element_blank()
  ) +
  guides(linetype = guide_none()) +
  scale_y_continuous(name = "Score") +
  scale_x_continuous(
    name = "Genome Position (bp)",
    limits = c(-5,5),
    expand = c(0,0)
  )

ggsave(
  "PaperFigures/ViterbiWithOverlapExplanation.png",
  width = 2000,
  height = 1000,
  units = "px"
)

rm(X, plotData, A, B, M, N, W, f)

##### Viterbi with 20bp Overlap Position Relative to Reference Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_20b_overlap_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = rep(-80:80, 2),
  frequency = c(
    data$frequency[data$position %in% -80:80],
    weiner_dist$frequency[weiner_dist$position %in% -80:80]
  ),
  dataset = rep(c("viterbi", "weiner"), each = 161)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("#990099", "black"),
    labels = c("Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_20bp_overlap_dyad_positions.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/viterbi_20bp_overlap_dyad_positions.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_20bp_overlap_dyad_positions.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string)

##### Viterbi with 20bp Overlap Offset Relative to Reference Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_20b_overlap_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = rep(-100:100, 2),
  frequency = c(
    data$frequency[data$position %in% -100:100],
    weiner_offset$frequency[weiner_offset$position %in% -100:100]
  ),
  dataset = rep(c("viterbi", "weiner"), each = 201)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("#990099", "black"),
    labels = c("Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_20bp_overlap_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/viterbi_20bp_overlap_dyad_offset.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_20bp_overlap_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string)

##### Viterbi with 20bp Overlap Placement TSS and Scatter Plot #####
viterbi_TSS <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/TSS_evals",
    "/viterbi_20b_overlap_frequencies.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(viterbi_TSS) <- c("position", "frequency")

regressionData <- data.frame(
  x = brogaard_TSS$frequency,
  y = viterbi_TSS$frequency
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Reference Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Viterbi Dyad Frequency",
    limits = c(0, 65),
    breaks = 0:6*10,
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 100,
    y = 10,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme

P

png(
  "PaperFigures/viterbi_20bp_overlap_scatter.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/viterbi_20bp_overlap_scatter.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  position = viterbi_TSS$position,
  frequency = viterbi_TSS$frequency
)

s <- 5 # Smoothing window radius
smooth_plotData = data.frame(
  position = (R[1]+s):(R[2]-s),
  frequency = rep(0, R[2]-R[1]-(2*s-1))
)

for (i in 1:nrow(smooth_plotData)) {
  smooth_plotData$frequency[i] <- mean(
    plotData$frequency[
      (plotData$position >= smooth_plotData$position[i]-s) &
        (plotData$position <= smooth_plotData$position[i]+s)
    ]
  )
}

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = 0,
      ymax = y * 60/100
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 0.75,
    color = "#FFB0FF"
  ) +
  geom_line(
    data = smooth_plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "#990099"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Viterbi Dyad Frequency",
    limits = c(0, 60),
    breaks = 0:6*10,
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/60,
      name = "Reference Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "#990099")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = -100,
    y = 55,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/viterbi_20bp_overlap_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

data <- plotData
data$smooth_frequency <- rep(NA, nrow(data))
data$smooth_frequency[
  data$position %in% smooth_plotData$position
] <- smooth_plotData$frequency

write.table(
  data,
  "PaperFigures/FigureData/viterbi_20bp_overlap_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

rm(C, data, viterbi_TSS, model, P, plotData, regressionData, smooth_plotData, i,
   p, r, s, X)

##### Dynamic Probability Distribution Explanation #####
X <- -500:500 / 100
N <- length(X)
W <- c(10, pi, sqrt(2), 1, 3, 8, 1.0258)
A <- c(1, 0.9, 0.95, 1.25, 1.1, 0.6, 1.23)
B <- c(0.1, 2, 4, sqrt(3), exp(1), 27, 0)
M <- 7 # length(W)
f <- function(x, W, A, B) {
  return(sum(A*sin(W*x+B)))
}
plotData <- data.frame(
  x = X,
  y = unlist(lapply(X, f, W, A, B))
)

ggplot(data = plotData) +
  geom_line(mapping = aes(x = x, y = y), linewidth = 1) +
  annotate(
    geom = "line",
    x = c(0, 0),
    y = c(-Inf, Inf),
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(-0.2, -0.2),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(-4, -4),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 15, face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.border = element_blank()
  ) +
  guides(linetype = guide_none()) +
  scale_y_continuous(name = "Score") +
  scale_x_continuous(
    name = "Genome Position (bp)",
    limits = c(-5,5),
    expand = c(0,0)
  )

ggsave(
  "PaperFigures/DynamicProbExplanation_1.png",
  width = 2000,
  height = 1000,
  units = "px"
)

ggplot(data = plotData) +
  geom_line(mapping = aes(x = x, y = y), linewidth = 1) +
  annotate(
    geom = "line",
    x = c(0, 0),
    y = c(-Inf, Inf),
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(0.2, 0.2),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(4, 4),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 15, face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.border = element_blank()
  ) +
  guides(linetype = guide_none()) +
  scale_y_continuous(name = "Score") +
  scale_x_continuous(
    name = "Genome Position (bp)",
    limits = c(-5,5),
    expand = c(0,0)
  )

ggsave(
  "PaperFigures/DynamicProbExplanation_2.png",
  width = 2000,
  height = 1000,
  units = "px"
)

ggplot(data = plotData[plotData$x <= 0,]) +
  geom_line(mapping = aes(x = 2*x, y = y), linewidth = 1) +
  annotate(
    geom = "line",
    x = c(0, 0),
    y = c(-Inf, Inf),
    linetype = "longdash",
    linewidth = 1
  ) +
  annotate(
    geom = "line",
    x = c(4, 4),
    y = c(-Inf, Inf),
    linetype = "dotted",
    linewidth = 1
  ) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 15, face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF"),
    panel.border = element_blank()
  ) +
  guides(linetype = guide_none()) +
  scale_y_continuous(name = "Probability") +
  scale_x_continuous(
    name = "Genome Position (bp)",
    limits = c(-5,5),
    expand = c(0,0)
  )

ggsave(
  "PaperFigures/DynamicProbExplanation_3.png",
  width = 2000,
  height = 1000,
  units = "px"
)

rm(plotData, A, B, M, N, W, X, f)

##### Probability Relative to Reference Dyad #####
periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

# Load Probabilities
chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "dynamic_probability_distribution/chromosomes.txt"
  ),
  header = FALSE
)$V1
probs <- list()
for (chrom in chroms) {
  probs[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/",
      "dynamic_probability_distribution/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(probs[[chrom]]) <- c("position", "prob")
}

chroms <- read.table(
  "PreprocessedData/all_nucleosomes/chromosomes.txt",
  header = FALSE
)$V1
nucleosomes <- list()
for (chrom in chroms) {
  nucleosomes[[chrom]] <- read.table(
    paste0("PreprocessedData/all_nucleosomes/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(nucleosomes[[chrom]]) <- c("start", "stop", "dyad")
}

W <- c(-100,100) # The window width to graph around the true dyad

total <- rep(0, W[2]-W[1] + 1)
n <- 0
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  for (i in 1:nrow(nucleosomes[[chrom]])) {
    dyad <- nucleosomes[[chrom]]$dyad[i]
    theseProbs <- probs[[chrom]]$prob[
      (probs[[chrom]]$position >= dyad + W[1]) &
        (probs[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseProbs) == W[2] - W[1] + 1) {
      total <- total + theseProbs
      n <- n + 1
    }
  }
}

plotData <- data.frame(
  relative_position = W[1]:W[2],
  average_prob = total / n
)

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

P <- ggplot(
  data = plotData,
  mapping = aes(x = relative_position, y = average_prob)
) +
  geom_line(
    data = annotationData[annotationData$x != 0,],
    mapping = aes(x = x, y = y, group = group),
    linetype = "dotted",
    linewidth = 0.75
  ) +
  annotate(
    geom = "line",
    linewidth = 1.25,
    linetype = "dashed",
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  geom_line(linewidth = 1.25) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Probability",
    limits = c(0.015,0.04),
    breaks = seq(0.015,0.04,0.005),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/probability_dyad_positions.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))

write(string, "PaperFigures/FigureData/probability_dyad_positions.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/probability_dyad_positions.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, G, nucleosomes, P, periodicity, plotData, probs, chrom,
   chroms, dyad, i, n, string, theseProbs, total, W)

##### Dynamic Probability Distribution TSS and Scatter Plot #####
prob_TSS <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/TSS_evals",
    "/mean_prob.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(prob_TSS) <- c("position", "mean_probability")

regressionData <- data.frame(
  x = brogaard_TSS$frequency,
  y = prob_TSS$mean_probability
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Probability",
    limits = c(0.018, 0.027),
    breaks = seq(0.018, 0.027, 0.003),
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = bquote(bold("r = ")*bold(.(toString(signif(r, 3)))*"*")),
    x = 100,
    y = 0.02,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme

P

png(
  "PaperFigures/probability_scatter.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/probability_scatter.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  position = prob_TSS$position,
  mean_probability = prob_TSS$mean_probability
)

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = 0.018,
      ymax = y * 0.009/100 + 0.018
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = mean_probability),
    linewidth = 0.75,
    color = "#990099"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Probability",
    limits = c(0.018, 0.027),
    breaks = seq(0.018, 0.027, 0.003),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/0.009 - 0.018*100/0.009,
      name = "Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "#990099")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  annotate(
    geom = "text",
    label = bquote(bold("r = ")*bold(.(toString(signif(r, 3)))*"*")),
    x = -100,
    y = 0.025,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/probability_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  "PaperFigures/FigureData/probability_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

rm(C, model, P, plotData, prob_TSS, regressionData, p, r, X)

##### Dynamic Probability IGV Track #####

# Nothing to generate here...

##### Strongly-Positioned Nucleosome log2 Ratio Pattern #####
pattern <- read.delim(
  "NucleosomePattern/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

periodicity <- read.delim(
  "NucleosomePattern/cl_ice_log2/periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,])[ncol(pattern):1], unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.45, 0.2),
    breaks = c(-0.4, -0.2, 0, 0.2),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/strong_log2_pattern.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#-: ",
  toString(periodicity$period[periodicity$strand == "-"][1]),
  "\n#+: ",
  toString(periodicity$period[periodicity$strand == "+"][1])
)

write(string, "PaperFigures/FigureData/strong_log2_pattern.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/strong_log2_pattern.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(pattern, periodicity, plotData, annotationData, string, P, G)

##### Strongly-Positioned Nucleosome Periodogram log2 Ratios #####
periodogram <- read.table(
  "NucleosomePattern/cl_ice_log2/periodogram.tsv",
  header = TRUE,
  sep = "\t"
)

periodicity <- read.delim(
  "NucleosomePattern/cl_ice_log2/periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

annotationData <- data.frame(
  x = c(
    rep(periodicity$period[periodicity$strand == "-"][1], 2),
    rep(periodicity$period[periodicity$strand == "+"][1], 2)
  ),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("-","+"), each = 2)
)

labs <- rep("", 14)
labs[c(4,9,14)] = c(5,10,15)

periodicity_annotations <- c(
  as.character(round(annotationData$x[annotationData$group == "-"][1], 2)),
  as.character(round(annotationData$x[annotationData$group == "+"][1], 2))
)
for (i in 1:length(periodicity_annotations)) {
  item <- periodicity_annotations[i]
  if (str_length(unlist(strsplit(item, "[.]"))[2]) < 2) {
    periodicity_annotations[i] <- paste0(periodicity_annotations[i], "0")
  }
}

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group, linetype = group),
    linewidth = 0.75
  ) +
  geom_line(
    data = periodogram,
    mapping = aes(x = period, y = normalized_power, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (nt)",
    limits = c(2,15),
    breaks = 2:15,
    labels = labs,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,0.65),
    breaks = (0:6)/10,
    expand = c(0,0)
  ) +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand", "Plus Strand")
  ) +
  annotate(
    "text",
    x = 10.5,
    y = 0.5,
    label = paste0(periodicity_annotations[1], "nt"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "bottom",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.5,
    y = 0.5,
    label = paste0(periodicity_annotations[2], "nt"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "blue"
  ) +
  guides(linetype = guide_none())

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/strong_log2_periododgram.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#-: ",
  toString(periodicity$period[periodicity$strand == "-"][1]),
  "\n#+: ",
  toString(periodicity$period[periodicity$strand == "+"][1])
)

write(string, "PaperFigures/FigureData/strong_log2_periodogram.tsv")

write.table(
  periodogram,
  file = "PaperFigures/FigureData/strong_log2_periodogram.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(annotationData, G, P, periodicity, periodogram, labs, string)

##### Strongly-Positioned Nucleosome Absolute Difference Pattern #####
pattern <- read.delim(
  "NucleosomePattern/cl_ice_abs_diff/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,])[ncol(pattern):1], unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(30, 55),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/strong_abs_diff_pattern.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/strong_abs_diff_pattern.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(pattern, G, P, plotData)

##### Greedy Nucleosome Position Relative to Strongly-Positioned Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/strong_greedy_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

weiner_data <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/strong_Weiner_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(weiner_data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = rep(-80:80, 2),
  frequency = c(
    data$frequency[data$position %in% -80:80],
    weiner_data$frequency[weiner_data$position %in% -80:80]
  ),
  dataset = rep(c("greedy", "weiner"), each = 161)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "black"),
    labels = c("Greedy", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_strong_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/greedy_strong_dyad_position.tsv")
write.table(
  cbind(
    rbind(data, weiner_data),
    source = c(rep("greedy", nrow(data)), rep("weiner", nrow(weiner_data)))
  ),
  file = "PaperFigures/FigureData/greedy_strong_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string, weiner_data)

##### Greedy Nucleosome Offset Relative to Strongly-Positioned Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/strong_greedy_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

weiner_data <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/strong_Weiner_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(weiner_data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = rep(-100:100, 2),
  frequency = c(
    data$frequency[data$position %in% -100:100],
    weiner_data$frequency[weiner_data$position %in% -100:100]
  ),
  dataset = rep(c("greedy", "weiner"), each = 201)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "black"),
    labels = c("Greedy", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest to Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_strong_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/greedy_strong_dyad_offset.tsv")
write.table(
  cbind(
    rbind(data, weiner_data),
    source = c(rep("greedy", nrow(data)), rep("weiner", nrow(weiner_data)))
  ),
  file = "PaperFigures/FigureData/greedy_strong_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string, weiner_data)

##### Alternate Score vs. Greedy vs. Brogaard IGV Track #####

# Nothing to generate here...

##### Viterbi Nucleosome Position Relative to Strongly-Positioned Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/strong_viterbi_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

weiner_data <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/strong_Weiner_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(weiner_data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = rep(-80:80, 2),
  frequency = c(
    data$frequency[data$position %in% -80:80],
    weiner_data$frequency[weiner_data$position %in% -80:80]
  ),
  dataset = rep(c("viterbi", "weiner"), each = 161)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("#990099", "black"),
    labels = c("Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_strong_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/viterbi_strong_dyad_position.tsv")
write.table(
  cbind(
    rbind(data, weiner_data),
    source = c(rep("greedy", nrow(data)), rep("weiner", nrow(weiner_data)))
  ),
  file = "PaperFigures/FigureData/viterbi_strong_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string, weiner_data)

##### Viterbi Nucleosome Offset Relative to Strongly-Positioned Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/strong_viterbi_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

weiner_data <- read.delim(
  paste0(
    "NucleosomePattern/weiner_brogaard_distance_evals",
    "/strong_Weiner_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(weiner_data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = rep(-100:100, 2),
  frequency = c(
    data$frequency[data$position %in% -100:100],
    weiner_data$frequency[weiner_data$position %in% -100:100]
  ),
  dataset = rep(c("viterbi", "weiner"), each = 201)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("#990099", "black"),
    labels = c("Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_strong_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/viterbi_strong_dyad_offset.tsv")
write.table(
  cbind(
    rbind(data, weiner_data),
    source = c(rep("greedy", nrow(data)), rep("weiner", nrow(weiner_data)))
  ),
  file = "PaperFigures/FigureData/viterbi_strong_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string, weiner_data)

##### Alternate Viterbi vs. Greedy vs. Brogaard vs. Weiner IGV #####

# Nothing to generate here...

##### Random Log2 Pattern #####
pattern <- read.delim(
  "RandomPattern/cl_ice_log2/0/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,])[ncol(pattern):1], unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.45, 0.2),
    breaks = c(-0.4, -0.2, 0, 0.2),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/random_log2_pattern.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/random_log2_pattern.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(pattern, plotData, P, G)

##### Random Abs Diff Pattern #####
pattern <- read.delim(
  "RandomPattern/cl_ice_abs_diff/0/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,])[ncol(pattern):1], unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(30, 55),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/random_abs_diff_pattern.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/FigS15.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(pattern, G, P, plotData)

##### Nucleosome Log2 Std Relation #####
pattern <- read.delim(
  "NucleosomePattern/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  mean_val = c(unlist(pattern[4,]), unlist(pattern[5,])),
  std_val = c(sqrt(unlist(pattern[6,])), sqrt(unlist(pattern[7,]))),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val, color = strand),
    size = 1.5
  ) +
  scale_x_continuous(
    name = bquote(bold("Mean log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5,0.25),
    breaks = seq(-0.5, 0.25, 0.1),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(
      bold("Std Dev of log"[bold("2")])*bold("(Cellular/Naked DNA)")
    ),
    limits = c(0.8, 1.1),
    breaks = seq(0.8, 1.1, 0.1),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand", "Plus Strand")
  ) +
  plotTheme +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 5))
  ) +
  theme(legend.position.inside = c(0.01,1))

P

png(
  "PaperFigures/nucleosome_log2_std_relation.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/nucleosome_log2_std_relation.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  mean_val = unlist(pattern[2,]),
  std_val = sqrt(unlist(pattern[3,]))
)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val),
    size = 1.5
  ) +
  scale_x_continuous(
    name = bquote(bold("Mean log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5,0.25),
    breaks = seq(-0.5, 0.25, 0.1),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(
      bold("Std Dev of log"[bold("2")])*bold("(Cellular/Naked DNA)")
    ),
    limits = c(0.8, 1.1),
    breaks = seq(0.8, 1.1, 0.1),
    expand = c(0,0)
  ) +
  plotTheme

P

png(
  "PaperFigures/nucleosome_log2_std_relation_overall.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/nucleosome_log2_std_relation_overall.tsv",
  quote = FALSE,
  row.names = FALSE
)

rm(P, pattern, plotData)

##### Nucleosome Abs Diff Std Relation #####
pattern <- read.delim(
  "NucleosomePattern/cl_ice_abs_diff/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  mean_val = c(unlist(pattern[4,]), unlist(pattern[5,])),
  std_val = c(sqrt(unlist(pattern[6,])), sqrt(unlist(pattern[7,]))),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val, color = strand),
    size = 1.5
  ) +
  scale_x_continuous(
    name = "Mean |Cellular - Naked DNA|",
    limits = c(30,55),
    breaks = seq(30, 55, 5),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Std Dev of |Cellular - Naked DNA|",
    limits = c(25, 75),
    breaks = seq(25, 75, 10),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand", "Plus Strand")
  ) +
  plotTheme +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 5))
  ) +
  theme(legend.position.inside = c(0.01,1))

P

png(
  "PaperFigures/nucleosome_abs_diff_std_relation.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/nucleosome_abs_diff_std_relation.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  mean_val = unlist(pattern[2,]),
  std_val = sqrt(unlist(pattern[3,]))
)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val),
    size = 1.5
  ) +
  scale_x_continuous(
    name = "Mean |Cellular - Naked DNA|",
    limits = c(30,55),
    breaks = seq(30, 55, 5),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Std Dev of |Cellular - Naked DNA|",
    limits = c(25, 75),
    breaks = seq(25, 75, 10),
    expand = c(0,0)
  ) +
  plotTheme

P

png(
  "PaperFigures/nucleosome_abs_diff_std_relation_overall.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/nucleosome_abs_diff_std_relation_overall.tsv",
  quote = FALSE,
  row.names = FALSE
)

rm(P, pattern, plotData)

##### Random Log2 Std Relation #####
pattern <- read.delim(
  "RandomPattern/cl_ice_log2/0/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  mean_val = c(unlist(pattern[4,]), unlist(pattern[5,])),
  std_val = c(sqrt(unlist(pattern[6,])), sqrt(unlist(pattern[7,]))),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val, color = strand),
    size = 1.5
  ) +
  scale_x_continuous(
    name = bquote(bold("Mean log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.28,-0.15),
    breaks = seq(-0.28, -0.15, 0.05),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(
      bold("Std Dev of log"[bold("2")])*bold("(Cellular/Naked DNA)")
    ),
    limits = c(0.9, 1),
    breaks = seq(0.9, 1, 0.02),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand", "Plus Strand")
  ) +
  plotTheme +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 5))
  ) +
  theme(legend.position.inside = c(0.75,1))

P

png(
  "PaperFigures/random_log2_std_relation.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/random_log2_std_relation.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  mean_val = unlist(pattern[2,]),
  std_val = sqrt(unlist(pattern[3,]))
)

L <- lm(std_val ~ mean_val, plotData)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val),
    size = 1.5
  ) +
  scale_x_continuous(
    name = bquote(bold("Mean log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.26,-0.2),
    breaks = seq(-0.26, -0.2, 0.05),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(
      bold("Std Dev of log"[bold("2")])*bold("(Cellular/Naked DNA)")
    ),
    limits = c(0.915, 0.975),
    breaks = seq(0.915, 0.975, 0.015),
    expand = c(0,0)
  ) +
  plotTheme

P

png(
  "PaperFigures/random_log2_std_relation_overall.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/random_log2_std_relation_overall.tsv",
  quote = FALSE,
  row.names = FALSE
)

rm(P, pattern, plotData)

##### Random Abs Diff Std Relation #####
pattern <- read.delim(
  "RandomPattern/cl_ice_abs_diff/0/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  mean_val = c(unlist(pattern[4,]), unlist(pattern[5,])),
  std_val = c(sqrt(unlist(pattern[6,])), sqrt(unlist(pattern[7,]))),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val, color = strand),
    size = 1.5
  ) +
  scale_x_continuous(
    name = "Mean |Cellular - Naked DNA|",
    limits = c(32,38),
    breaks = seq(32, 38, 2),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Std Dev of |Cellular - Naked DNA|",
    limits = c(34, 44),
    breaks = seq(34, 44, 2),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand", "Plus Strand")
  ) +
  plotTheme +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 5))
  )

P

png(
  "PaperFigures/random_abs_diff_std_relation.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/random_abs_diff_std_relation.tsv",
  quote = FALSE,
  row.names = FALSE
)

plotData <- data.frame(
  mean_val = unlist(pattern[2,]),
  std_val = sqrt(unlist(pattern[3,]))
)

P <- ggplot(data = plotData) +
  geom_point(
    mapping = aes(x = mean_val, y = std_val),
    size = 1.5
  ) +
  scale_x_continuous(
    name = "Mean |Cellular - Naked DNA|",
    limits = c(33,37),
    breaks = seq(33, 37, 2),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Std Dev of |Cellular - Naked DNA|",
    limits = c(34, 42),
    breaks = seq(34, 42, 2),
    expand = c(0,0)
  ) +
  plotTheme

P

png(
  "PaperFigures/random_abs_diff_std_relation_overall.png",
  width = 2000,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/random_abs_diff_std_relation_overall.tsv",
  quote = FALSE,
  row.names = FALSE
)

rm(P, pattern, plotData)

##### TSS and Scatter Plot of Differences #####

# Load Data
chroms <- unlist(read.table("PreprocessedData/TSS/chromosomes.txt"))

TSS_data <- list()
N <- 0
diff_minus <- list()
diff_plus <- list()
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  TSS_data[[chrom]] <- read.table(
    paste0("PreprocessedData/TSS/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(TSS_data[[chrom]]) <- c(
    "position",
    "gene",
    "internal",
    "flag",
    "strand"
  )
  N <- N + nrow(TSS_data[[chrom]])
  
  diff_minus[[chrom]] <- read.table(
    paste0("PreprocessedData/cl_ice_diff_minus/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(diff_minus[[chrom]]) <- c("position", "value")
  
  diff_plus[[chrom]] <- read.table(
    paste0("PreprocessedData/cl_ice_diff_plus/", chrom, ".tsv"),
    header = FALSE,
    sep = "\t"
  )
  colnames(diff_plus[[chrom]]) <- c("position", "value")
}

# Get TSS Spots
values <- matrix(nrow = N, ncol = R[2]-R[1])
inbetween_positions <- R[1]:(R[2]-1) + 0.5
I <- 1
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  for (i in 1:nrow(TSS_data[[chrom]])) {
    position <- TSS_data[[chrom]]$position[i]
    if (TSS_data[[chrom]]$strand[i] == "+") {
      theseDiffs <- rbind(
        diff_minus[[chrom]][
          (diff_minus[[chrom]]$position >= position + R[1]) &
            (diff_minus[[chrom]]$position <= position + R[2]),
        ],
        plus_vals <- diff_plus[[chrom]][
          (diff_plus[[chrom]]$position >= position + R[1]) &
            (diff_plus[[chrom]]$position <= position + R[2]),
        ]
      )
      theseDiffs <- rbind(
        theseDiffs,
        data.frame(
          position = (position + inbetween_positions)[
            !((position+inbetween_positions) %in% theseDiffs$position)
          ],
          value = NA
        )
      )
      
      theseDiffs <- theseDiffs[order(theseDiffs$position),]
      
      values[I,] <- theseDiffs$value
    } else {
      theseDiffs <- rbind(
        diff_minus[[chrom]][
          (diff_minus[[chrom]]$position >= position - R[2]) &
            (diff_minus[[chrom]]$position <= position - R[1]),
        ],
        plus_vals <- diff_plus[[chrom]][
          (diff_plus[[chrom]]$position >= position - R[2]) &
            (diff_plus[[chrom]]$position <= position - R[1]),
        ]
      )
      theseDiffs <- rbind(
        theseDiffs,
        data.frame(
          position = (position - inbetween_positions)[
            !((position-inbetween_positions) %in% theseDiffs$position)
          ],
          value = NA
        )
      )
      
      theseDiffs <- theseDiffs[
        order(theseDiffs$position, decreasing = TRUE),
      ]
      
      values[I,] <- theseDiffs$value
    }
    
    I <- I + 1
  }
}

plotData <- data.frame(
  position = inbetween_positions,
  diff = colMeans(values, na.rm = TRUE)
)

regressionData <- data.frame(
  x = (
    brogaard_TSS$frequency[1:(nrow(brogaard_TSS)-1)] +
      brogaard_TSS$frequency[2:nrow(brogaard_TSS)]
  ) / 2 ,
  y = plotData$diff
)

model <- lm(y ~ x, data = regressionData)

C <- cor.test(regressionData$x, regressionData$y)
r <- C$estimate
p <- C$p.value

X <- seq(0, max(brogaard_TSS$frequency), 0.1)
P <- ggplot(data = regressionData) +
  geom_point(
    mapping = aes(x = x, y = y),
    size = 0.75
  ) +
  annotate(
    geom = "line",
    x = X,
    y = model$coefficients[1] + model$coefficients[2]*X,
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Dyad Frequency",
    limits = c(0,120),
    breaks = 0:6*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Cellular - Naked DNA",
    limits = c(-26, 10),
    breaks = ((-5):2) * 5,
    expand = c(0,0)
  ) +
  annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = 100,
    y = -20,
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  plotTheme

P

png(
  "PaperFigures/diff_scatter.png",
  width = 1700,
  height = 1300,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  regressionData,
  file = "PaperFigures/FigureData/diff_scatter.tsv",
  quote = FALSE,
  row.names = FALSE
)

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = -26,
      ymax = y * 36/100 - 26
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(
      x = position,
      y = diff
    ),
    linewidth = 0.75,
    color = "blue"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Cellular - Naked DNA",
    limits = c(-26, 10),
    breaks = ((-5):2)*5,
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/36 + 26*100/36,
      name = "Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "blue")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  ) + annotate(
    geom = "text",
    label = paste0("r = ", toString(signif(r, 3)), "*"),
    fontface = "bold",
    x = -100,
    y = 7,
    size = smallTxtSize,
    size.unit = "pt"
  )

P

png(
  "PaperFigures/diff_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

data <- data.frame(
  diff_position = c(inbetween_positions, ""),
  diff_mean = c(plotData$diff, "")
)

write.table(
  data,
  "PaperFigures/FigureData/diff_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(diff_minus, diff_plus, data, P, plotData, plus_vals, theseDiffs, TSS_data,
   values, chrom, chroms, i, I, inbetween_positions, N, position, C, model,
   regressionData, r, p, X)

##### All Nucleosome log2 Ratio Pattern Non-Aligned #####
pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,]), unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5, 0),
    breaks = c(-5:0 / 10),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (3′→5′)", "Plus Strand")
  ) +
  theme(legend.direction="horizontal", legend.position.inside=c(0.005,0.125))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_pattern_unaligned.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#-: ",
  toString(periodicity$period[periodicity$strand == "-"][1]),
  "\n#+: ",
  toString(periodicity$period[periodicity$strand == "+"][1])
)

write(string, "PaperFigures/FigureData/all_log2_pattern_unaligned.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_log2_pattern_unaligned.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(pattern, periodicity, plotData, string, P, G)

##### All Nucleosome Absolute Difference Pattern Non-Aligned #####
pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_abs_diff/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(unlist(pattern[4,]), unlist(pattern[5,])),
  strand = rep(c("-", "+"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(30, 46),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (3′→5′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_pattern_unaligned.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_abs_diff_pattern_unaligned.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(pattern, G, P, plotData)

##### Not Strong Nucleosome log2 Ratio Pattern and Periodogram #####
all_pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)
strong_pattern <- read.delim(
  "NucleosomePattern/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

minus_all_mean = unlist(all_pattern[4,])
minus_all_n = unlist(all_pattern[8,])
plus_all_mean = unlist(all_pattern[5,])
plus_all_n = unlist(all_pattern[9,])
minus_strong_mean = unlist(strong_pattern[4,])
minus_strong_n = unlist(strong_pattern[8,])
plus_strong_mean = unlist(strong_pattern[5,])
plus_strong_n = unlist(strong_pattern[9,])

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = c(
    (
      (
        minus_all_mean*minus_all_n - minus_strong_mean*minus_strong_n
      ) / (minus_all_n - minus_strong_n)
    )[ncol(all_pattern):1],
    (
      (
        plus_all_mean*plus_all_n - plus_strong_mean*plus_strong_n
      ) / (plus_all_n - plus_strong_n)
    )
  ),
  strand = rep(c("-", "+"), each = 146)
)

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5, 0),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/notstrong_log2_pattern.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#-: ",
  toString(periodicity$period[periodicity$strand == "-"][1]),
  "\n#+: ",
  toString(periodicity$period[periodicity$strand == "+"][1])
)

write(string, "PaperFigures/FigureData/notstrong_log2_pattern.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/notstrong_log2_pattern.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

L_minus <- lsp(
  plotData$values[plotData$strand == "-"],
  times = -73:72 + 0.5,
  type = "period",
  to = 15,
  ofac = 100
)
L_plus <- lsp(
  plotData$values[plotData$strand == "+"],
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

annotationData <- data.frame(
  x = c(
    rep(L_minus$peak.at[1], 2),
    rep(L_plus$peak.at[1], 2)
  ),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("-","+"), each = 2)
)

labs <- rep("", 14)
labs[c(4,9,14)] = c(5,10,15)

periodicity_annotations <- c(
  as.character(round(annotationData$x[annotationData$group == "-"][1], 2)),
  as.character(round(annotationData$x[annotationData$group == "+"][1], 2))
)
for (i in 1:length(periodicity_annotations)) {
  item <- periodicity_annotations[i]
  if (str_length(unlist(strsplit(item, "[.]"))[2]) < 2) {
    periodicity_annotations[i] <- paste0(periodicity_annotations[i], "0")
  }
}

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group, linetype = group),
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (nt)",
    limits = c(2,15),
    breaks = 2:15,
    labels = labs,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,0.6),
    breaks = (0:6)/10,
    expand = c(0,0)
  ) +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Minus Strand", "Plus Strand")
  ) +
  annotate(
    "text",
    x = 10.5,
    y = 0.5,
    label = paste0(periodicity_annotations[1], "nt"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "bottom",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.5,
    y = 0.5,
    label = paste0(periodicity_annotations[2], "nt"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "blue"
  ) +
  guides(linetype = guide_none())

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/notstrong_log2_periododgram.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#-: ",
  toString(plotData$period[plotData$strand == "-"][1]),
  "\n#+: ",
  toString(plotData$period[plotData$strand == "+"][1])
)

write(string, "PaperFigures/FigureData/notstrong_log2_periodogram.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/notstrong_log2_periodogram.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

##### Log2 Pattern with Shift #####

shifted_pattern <- read.delim(
  "NucleosomePattern_shifted/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

unshifted_pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = rep(-73:72 + 0.5, 6),
  values = c(
    unlist(shifted_pattern[4,])[ncol(shifted_pattern):1],
    unlist(shifted_pattern[5,]),
    unlist(shifted_pattern[2,]),
    unlist(unshifted_pattern[4,])[ncol(unshifted_pattern):1],
    unlist(unshifted_pattern[5,]),
    unlist(unshifted_pattern[2,])
  ),
  strand = rep(rep(c("-", "+", "0"), each = 146), 2),
  shifted = rep(c(TRUE, FALSE), each = 146*3)
)

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$strand == "+",],
    mapping = aes(x = position, y = values, color = shifted),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5, 0),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    labels = c("Unshifted", "10bp Shift"),
    type = c("blue", "red")
  ) +
  theme(legend.direction = "horizontal", legend.position.inside = c(0.7,0.3)) +
  ggtitle("Plus Strand")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_plus_shifted.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_log2_plus_shifted.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$strand == "-",],
    mapping = aes(x = position, y = values, color = shifted),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5, 0),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    labels = c("Unshifted", "10bp Shift"),
    type = c("blue", "red")
  ) +
  theme(legend.direction = "horizontal", legend.position.inside = c(0.7,0.3)) +
  ggtitle("Minus Strand (5′→3′)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_minus_shifted.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_log2_minus_shifted.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$strand == "0",],
    mapping = aes(x = position, y = values, color = shifted),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5, 0),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    labels = c("Unshifted", "10bp Shift"),
    type = c("blue", "red")
  ) +
  theme(legend.direction = "horizontal", legend.position.inside = c(0.7,0.3)) +
  ggtitle("Average Pattern (Both Strands)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_pattern_shifted.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_log2_pattern_shifted.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

diffData <- data.frame(
  position = rep(-73:72 + 0.5, 3),
  strand = rep(c("-", "+", "0"), each = 146),
  diff = c(
    plotData$values[(plotData$strand == "-") & plotData$shifted] - 
      plotData$values[(plotData$strand == "-") & (!plotData$shifted)],
    plotData$values[(plotData$strand == "+") & plotData$shifted] - 
      plotData$values[(plotData$strand == "+") & (!plotData$shifted)],
    plotData$values[(plotData$strand == "0") & plotData$shifted] - 
      plotData$values[(plotData$strand == "0") & (!plotData$shifted)]
  )
)

P <- ggplot() +
  geom_line(
    data = diffData,
    mapping = aes(x = position, y = abs(diff), color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("Deviation from log"[bold("2")])*bold(" Pattern")),
    limits = c(0, 0.5),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue", "black"),
    labels = c(
      "Minus Strand (5′→3′)",
      "Plus Strand",
      "Average Pattern (Both Strands)"
    )
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_log2_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "+",],
    mapping = aes(x = position, y = abs(diff)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("Deviation from log"[bold("2")])*bold(" Pattern")),
    limits = c(0, 0.5),
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Plus Strand")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_plus_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_log2_plus_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "-",],
    mapping = aes(x = position, y = abs(diff)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("Deviation from log"[bold("2")])*bold(" Pattern")),
    limits = c(0, 0.5),
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Minus Strand (5′→3′)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_minus_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_log2_minus_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "0",],
    mapping = aes(x = position, y = abs(diff)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("Deviation from log"[bold("2")])*bold(" Pattern")),
    limits = c(0, 0.5),
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Average Pattern (Both Strands)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_pattern_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_log2_pattern_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "0",],
    mapping = aes(x = position, y = diff),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("Deviation from log"[bold("2")])*bold(" Pattern")),
    limits = c(-0.5, 0.5),
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Average Pattern (Both Strands)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_pattern_shifted_signed_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_log2_pattern_shifted_signed_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(shifted_pattern, unshifted_pattern, plotData, annotationData, diffData, P, G)

##### Absolute Difference Pattern with Shift #####

shifted_pattern <- read.delim(
  "NucleosomePattern_shifted/cl_ice_abs_diff/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

unshifted_pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_abs_diff/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = rep(-73:72 + 0.5, 6),
  values = c(
    unlist(shifted_pattern[4,])[ncol(shifted_pattern):1],
    unlist(shifted_pattern[5,]),
    unlist(shifted_pattern[2,]),
    unlist(unshifted_pattern[4,])[ncol(unshifted_pattern):1],
    unlist(unshifted_pattern[5,]),
    unlist(unshifted_pattern[2,])
  ),
  strand = rep(rep(c("-", "+", "0"), each = 146), 2),
  shifted = rep(c(TRUE, FALSE), each = 146*3)
)

P <- ggplot() +
  geom_line(
    data = plotData[plotData$strand == "+",],
    mapping = aes(x = position, y = values, color = shifted),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = c(-60 + 20*(0:6)),
    labels = c(-60 + 20*(0:6)),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(30, 46),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    labels = c("Unshifted","10bp Shift"),
    type = c("blue", "red")
  ) +
  plotTheme +
  ggtitle("Plus Strand")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_plus_shifted.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_abs_diff_plus_shifted.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = plotData[plotData$strand == "-",],
    mapping = aes(x = position, y = values, color = shifted),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = c(-60 + 20*(0:6)),
    labels = c(-60 + 20*(0:6)),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(30, 46),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    labels = c("Unshifted","10bp Shift"),
    type = c("blue", "red")
  ) +
  plotTheme +
  ggtitle("Minus Strand (5′→3′)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_minus_shifted.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_abs_diff_minus_shifted.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = plotData[plotData$strand == "0",],
    mapping = aes(x = position, y = values, color = shifted),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = c(-60 + 20*(0:6)),
    labels = c(-60 + 20*(0:6)),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "|Cellular - Naked DNA|",
    limits = c(30, 46),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    labels = c("Unshifted","10bp Shift"),
    type = c("blue", "red")
  ) +
  plotTheme +
  ggtitle("Average Pattern (Both Strands)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_pattern_shifted.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_abs_diff_pattern_shifted.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

diffData <- data.frame(
  position = rep(-73:72 + 0.5, 3),
  strand = rep(c("-", "+", "0"), each = 146),
  diff = c(
    plotData$values[(plotData$strand == "-") & plotData$shifted] - 
      plotData$values[(plotData$strand == "-") & (!plotData$shifted)],
    plotData$values[(plotData$strand == "+") & plotData$shifted] - 
      plotData$values[(plotData$strand == "+") & (!plotData$shifted)],
    plotData$values[(plotData$strand == "0") & plotData$shifted] - 
      plotData$values[(plotData$strand == "0") & (!plotData$shifted)]
  )
)

P <- ggplot() +
  geom_line(
    data = diffData,
    mapping = aes(x = position, y = abs(diff), color = strand),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Deviation from\nAbsolute Difference Pattern",
    limits = c(0, 16),
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue", "black"),
    labels = c("Minus Strand (5′→3′)", "Plus Strand", "Pattern")
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_abs_diff_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "+",],
    mapping = aes(x = position, y = abs(diff)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Deviation from\nAbsolute Difference Pattern",
    limits = c(0, 16),
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Plus Strand")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_plus_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_abs_diff_plus_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "-",],
    mapping = aes(x = position, y = abs(diff)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Deviation from\nAbsolute Difference Pattern",
    limits = c(0, 16),
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Minus Strand (5′→3′)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_minus_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_abs_diff_minus_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "0",],
    mapping = aes(x = position, y = abs(diff)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Deviation from\nAbsolute Difference Pattern",
    limits = c(0, 16),
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Average Pattern (Both Strands)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_pattern_shifted_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_abs_diff_pattern_shifted_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_line(
    data = diffData[diffData$strand == "0",],
    mapping = aes(x = position, y = diff),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Deviation from\nAbsolute Difference Pattern",
    limits = c(-16, 16),
    breaks = -3:3 * 5,
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle("Average Pattern (Both Strands)")

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_abs_diff_pattern_shifted_signed_diff.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  diffData,
  file = "PaperFigures/FigureData/all_abs_diff_pattern_shifted_signed_diff.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(shifted_pattern, unshifted_pattern, plotData, maxes, diffData, P, G)

##### Greedy Nucleosome Position Relative to Reference Dyad vs. Viterbi #####
greedy <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_greedy_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(greedy) <- c("position", "frequency")

viterbi <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(viterbi) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = rep(-80:80, 3),
  frequency = c(
    greedy$frequency[greedy$position %in% -80:80],
    viterbi$frequency[viterbi$position %in% -80:80],
    weiner_dist$frequency[weiner_dist$position %in% -80:80]
  ),
  dataset = factor(
    rep(c("greedy", "viterbi", "weiner"), each = 161),
    levels = c("greedy", "viterbi", "weiner")
  )
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData[annotationData$x != 0,],
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  annotate(
    geom = "line",
    x = c(0,0),
    y = c(-Inf, Inf),
    linewidth = 1.25,
    linetype = "dashed"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "#990099", "black"),
    labels = c("Greedy", "Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_viterbi_weiner_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/greedy_viterbi_weiner_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot(data = plotData[plotData$dataset != "weiner",]) +
  geom_line(
    data = annotationData[annotationData$x != 0,],
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  annotate(
    geom = "line",
    x = c(0,0),
    y = c(-Inf, Inf),
    linewidth = 1.25,
    linetype = "dashed"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "#990099"),
    labels = c("Greedy", "Viterbi")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_viterbi_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/greedy_viterbi_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(annotationData, greedy, viterbi, G, P, plotData)

##### Greedy Nucleosome Offset Relative to Reference Dyad vs. Viterbi #####
greedy <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_greedy_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(greedy) <- c("position", "frequency")

viterbi <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(viterbi) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = rep(-100:100, 3),
  frequency = c(
    greedy$frequency[greedy$position %in% -100:100],
    viterbi$frequency[viterbi$position %in% -100:100],
    weiner_offset$frequency[weiner_offset$position %in% -100:100]
  ),
  dataset = factor(
    rep(c("greedy", "viterbi", "weiner"), each = 201),
    levels = c("greedy", "viterbi", "weiner")
  )
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData[annotationData$x != 0,],
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  annotate(
    geom = "line",
    x = c(0,0),
    y = c(-Inf, Inf),
    linewidth = 1.25,
    linetype = "dashed"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "#990099", "black"),
    labels = c("Greedy", "Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_viterbi_weiner_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/greedy_viterbi_weiner_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot(data = plotData[plotData$dataset != "weiner",]) +
  geom_line(
    data = annotationData[annotationData$x != 0,],
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  annotate(
    geom = "line",
    x = c(0,0),
    y = c(-Inf, Inf),
    linewidth = 1.25,
    linetype = "dashed"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "#990099"),
    labels = c("Greedy", "Viterbi")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 2000),
    breaks = seq(0, 2000, 400),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_viterbi_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/greedy_viterbi_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(annotationData, greedy, viterbi, G, P, plotData)

##### Nucleosome Dipyrimidine Distribution #####

data <- read_json("NucleosomePattern/all_dinucs_stats.json")
plotData <- data.frame(values = unlist(data$data), groups = rep("all", data$n))

P <- ggplot(data = plotData) +
  geom_boxplot(
    mapping = aes(y = values),
    outliers = FALSE,
    linewidth = 1,
    color = "black",
    lineend = "square",
    width = 0.1
  ) +
  geom_point(
    data = plotData[
      (plotData$values > data$Q3 + 1.5*(data$Q3-data$Q1)) |
        (plotData$values < data$Q1 - 1.5*(data$Q3-data$Q1)),
    ],
    mapping = aes(x = 0, y = values),
    position = position_jitter(width = 0.05, height = 0),
    size = 0.75
  ) +
  plotTheme +
  scale_x_continuous(breaks = NULL, name = NULL, expand = c(0.025,0.025)) +
  scale_y_continuous(
    limits=c(0,146),
    expand=c(0,0),
    breaks=seq(0,146,20),
    name = "Number of Dipyrimidines"
  )

P

png(
  "PaperFigures/dinucleotide_boxplot.png",
  width = 600,
  height = 2000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData$values,
  file = "PaperFigures/FigureData/dinucleotide_boxplot.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

##### Viterbi Nucleosome Positions Relative to Q1 Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_Q1_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = -80:80,
  frequency = data$frequency[data$position %in% -80:80]
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad with ≤ 72 Dipyrimidines (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_Q1_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_Q1_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(annotationData, data, G, P, plotData)

##### Viterbi Nucleosome Offset Relative to Q1 Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_Q1_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = -100:100,
  frequency = data$frequency[data$position %in% -100:100]
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad with ≤ 72 Dipyrimidines (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_Q1_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_Q1_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(annotationData, data, G, P, plotData)

##### Viterbi Nucleosome Positions Relative to Q3 Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_Q3_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

plotData <- data.frame(
  position = -80:80,
  frequency = data$frequency[data$position %in% -80:80]
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad with ≥ 81 Dipyrimidines (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_Q3_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_Q3_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(annotationData, data, G, P, plotData)

##### Viterbi Nucleosome Offset Relative to Q3 Dyad #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/all_viterbi_Q3_offset_freq_table.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(data) <- c("position", "frequency")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

plotData <- data.frame(
  position = -100:100,
  frequency = data$frequency[data$position %in% -100:100]
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad with ≥ 81 Dipyrimidines (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 550),
    breaks = seq(0, 550, 100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_Q3_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_Q3_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(annotationData, data, G, P, plotData)

##### Acetylation Barplots #####
load("NucleosomeClassDataTable.rda")
D <- D[D$dataset == "Cellular vs. Ice Combined",2:ncol(D)]

acetylationData <- cbind(
  score = D$score,
  D[16:ncol(D)][grep("ac", colnames(D)[16:ncol(D)])]
)
acetylationData <- acetylationData[!is.na(acetylationData[2]),]

plotData <- data.frame(
  mod = character(0),
  presence = logical(0),
  mean = numeric(0),
  std = numeric(0),
  signif = numeric(0)
)
signifData <- data.frame(
  mod = character(0),
  star = character(0)
)
for (j in 2:ncol(acetylationData)) {
  low <- acetylationData$score[
    acetylationData[[j]] < quantile(acetylationData[[j]], 0.25)
  ]
  high <- acetylationData$score[
    acetylationData[[j]] > quantile(acetylationData[[j]], 0.75)
  ]
  t <- t.test(low, high)
  p <- t$p.value * (ncol(acetylationData) - 1)
  char = "n.s."
  if (p < 0.05) { char = "*" }
  if (p < 0.01) { char = "**" }
  if (p < 0.001) { char = "***" }
  plotData <- rbind(
    plotData,
    data.frame(
      mod = rep(colnames(acetylationData)[j], 2),
      presence = c(TRUE, FALSE),
      mean = c(mean(high), mean(low)),
      std = c(sd(high) / sqrt(length(high)), sd(low) / sqrt(length(low))),
      signif = rep(p, 2)
    )
  )
  if (p < 0.05) {
    signifData = rbind(
      signifData,
      data.frame(
        mod = colnames(acetylationData)[j],
        mean = max(mean(high), mean(low)),
        star = char
      )
    )
  }
}
plotData$presence <- factor(plotData$presence, levels = c(TRUE, FALSE))

P <- ggplot(data = plotData) +
  geom_col(
    mapping=aes(x=mod,y=mean,fill=signif<0.05,group=presence,color=signif<0.05),
    position = position_dodge(),
  ) +
  plotTheme +
  scale_fill_discrete(type = c("grey", "white")) +
  scale_color_discrete(type = c("#6D6D6D", "black")) +
  guides(fill = "none", color = "none") +
  scale_y_continuous(
    name = "Mean Rotational Score",
    expand = c(0, 0),
    limits = c(-10, 135),
    breaks = c(0,45,90,135)
  ) +
  scale_x_discrete(name = "Posttranslational Modification") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 15)),
    axis.ticks.x = element_blank()
  ) +
  geom_text(
    mapping = aes(
      x=stage(mod,after_scale=x+rep(c(-0.2,0.2),length(unique(plotData$mod)))),
      y = -7,
      label = rep(c("🡅","🡇"), length(unique(plotData$mod)))
    ),
    fontface = "bold",
    size = 8
  ) +
  coord_cartesian(clip = FALSE, ylim = c(0,135)) +
  geom_errorbar(
    mapping = aes(
      x = mod,
      ymin = mean - std,
      ymax = mean + std,
      group = presence,
      color = signif < 0.05
    ),
    position = position_dodge()
  ) +
  geom_errorbarh(
    data = signifData,
    mapping = aes(
      xmin = stage(mod, after_scale = xmin - 0.4),
      xmax = stage(mod, after_scale = xmax + 0.4),
      y = max(mean) + 5,
      group = mod
    ),
    lineend = "square",
    height = 0
  ) +
  geom_text(
    data = signifData,
    mapping = aes(x = mod, y = max(mean) + 6, label = star),
    fontface = "bold"
  )

P

png(
  "PaperFigures/acetylation_barplots.png",
  width = 2000,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/acetylation_barplots.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(D,acetylationData,signifData,P,plotData,t,j,high,low,char,p)

##### Methylation Barplots #####
load("NucleosomeClassDataTable.rda")
D <- D[D$dataset == "Cellular vs. Ice Combined",2:ncol(D)]

methylationData <- cbind(
  score = D$score,
  D[16:ncol(D)][grep("me", colnames(D)[16:ncol(D)])]
)
methylationData <- methylationData[!is.na(methylationData[2]),]
methylationData <- methylationData[
  !(1:ncol(methylationData) %in% grep("H4", colnames(methylationData)))
]

plotData <- data.frame(
  mod = character(0),
  presence = logical(0),
  mean = numeric(0),
  std = numeric(0),
  signif = numeric(0)
)
signifData <- data.frame(
  mod = character(0),
  star = character(0)
)
for (j in 2:ncol(methylationData)) {
  low <- methylationData$score[
    methylationData[[j]] < quantile(methylationData[[j]], 0.25)
  ]
  high <- methylationData$score[
    methylationData[[j]] > quantile(methylationData[[j]], 0.75)
  ]
  t <- t.test(low, high)
  p <- t$p.value * (ncol(methylationData) - 1)
  char = "n.s."
  if (p < 0.05) { char = "*" }
  if (p < 0.01) { char = "**" }
  if (p < 0.001) { char = "***" }
  plotData <- rbind(
    plotData,
    data.frame(
      mod = rep(colnames(methylationData)[j], 2),
      presence = c(TRUE, FALSE),
      mean = c(mean(high), mean(low)),
      std = c(sd(high) / sqrt(length(high)), sd(low) / sqrt(length(low))),
      signif = rep(p, 2)
    )
  )
  if (p < 0.05) {
    signifData = rbind(
      signifData,
      data.frame(
        mod = colnames(methylationData)[j],
        mean = max(mean(high), mean(low)),
        star = char
      )
    )
  }
}
plotData$presence <- factor(plotData$presence, levels = c(TRUE, FALSE))

P <- ggplot(data = plotData) +
  geom_col(
    mapping=aes(x=mod,y=mean,fill=signif<0.05,group=presence,color=signif<0.05),
    position = position_dodge(),
  ) +
  plotTheme +
  scale_fill_discrete(type = c("grey", "white")) +
  scale_color_discrete(type = c("#6D6D6D", "black")) +
  guides(fill = "none", color = "none") +
  scale_y_continuous(
    name = "Mean Rotational Score",
    expand = c(0, 0),
    limits = c(-10, 135),
    breaks = c(0,45,90,135)
  ) +
  scale_x_discrete(name = "Posttranslational Modification") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 15)),
    axis.ticks.x = element_blank()
  ) +
  geom_text(
    mapping = aes(
      x=stage(mod,after_scale=x+rep(c(-0.2,0.2),length(unique(plotData$mod)))),
      y = -7,
      label = rep(c("🡅","🡇"), length(unique(plotData$mod)))
    ),
    fontface = "bold",
    size = 8
  ) +
  coord_cartesian(clip = FALSE, ylim = c(0,135)) +
  geom_errorbar(
    mapping = aes(
      x = mod,
      ymin = mean - std,
      ymax = mean + std,
      group = presence,
      color = signif < 0.05
    ),
    position = position_dodge()
  ) +
  geom_errorbarh(
    data = signifData,
    mapping = aes(
      xmin = stage(mod, after_scale = xmin - 0.4),
      xmax = stage(mod, after_scale = xmax + 0.4),
      y = max(mean) + 5,
      group = mod
    ),
    lineend = "square",
    height = 0
  ) +
  geom_text(
    data = signifData,
    mapping = aes(x = mod, y = max(mean) + 6, label = star),
    fontface = "bold"
  )

P

png(
  "PaperFigures/methylation_barplots.png",
  width = 1500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/methylation_barplots.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(D,methylationData,signifData,P,plotData,t,j,high,low,char,p)

##### Phosphorylation Barplots #####
load("NucleosomeClassDataTable.rda")
D <- D[D$dataset == "Cellular vs. Ice Combined",2:ncol(D)]

phosphorylationData <- cbind(
  score = D$score,
  D[16:ncol(D)][grep("ph", colnames(D)[16:ncol(D)])]
)
phosphorylationData <- phosphorylationData[!is.na(phosphorylationData[2]),]

plotData <- data.frame(
  mod = character(0),
  presence = logical(0),
  mean = numeric(0),
  std = numeric(0),
  signif = numeric(0)
)
signifData <- data.frame(
  mod = character(0),
  star = character(0)
)
for (j in 2:ncol(phosphorylationData)) {
  low <- phosphorylationData$score[
    phosphorylationData[[j]] < quantile(phosphorylationData[[j]], 0.25)
  ]
  high <- phosphorylationData$score[
    phosphorylationData[[j]] > quantile(phosphorylationData[[j]], 0.75)
  ]
  t <- t.test(low, high)
  p <- t$p.value * (ncol(phosphorylationData) - 1)
  char = "n.s."
  if (p < 0.05) { char = "*" }
  if (p < 0.01) { char = "**" }
  if (p < 0.001) { char = "***" }
  plotData <- rbind(
    plotData,
    data.frame(
      mod = rep(colnames(phosphorylationData)[j], 2),
      presence = c(TRUE, FALSE),
      mean = c(mean(high), mean(low)),
      std = c(sd(high) / sqrt(length(high)), sd(low) / sqrt(length(low))),
      signif = rep(p, 2)
    )
  )
  if (p < 0.05) {
    signifData = rbind(
      signifData,
      data.frame(
        mod = colnames(phosphorylationData)[j],
        mean = max(mean(high), mean(low)),
        star = char
      )
    )
  }
}
plotData$presence <- factor(plotData$presence, levels = c(TRUE, FALSE))

P <- ggplot(data = plotData) +
  geom_col(
    mapping=aes(x=mod,y=mean,fill=signif<0.05,group=presence,color=signif<0.05),
    position = position_dodge(),
  ) +
  plotTheme +
  scale_fill_discrete(type = c("grey", "white")) +
  scale_color_discrete(type = c("#6D6D6D", "black")) +
  guides(fill = "none", color = "none") +
  scale_y_continuous(
    name = "Mean Rotational Score",
    expand = c(0, 0),
    limits = c(-10, 135),
    breaks = c(0,45,90,135)
  ) +
  scale_x_discrete(name = "Posttranslational Modification") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 15)),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(hjust = 1.5),
    plot.margin = margin(t = 10, b = 10, r = 100, l = 100)
  ) +
  geom_text(
    mapping = aes(
      x=stage(mod,after_scale=x+rep(c(-0.2,0.2),length(unique(plotData$mod)))),
      y = -8,
      label = rep(c("🡅","🡇"), length(unique(plotData$mod)))
    ),
    fontface = "bold",
    size = 8
  ) +
  coord_cartesian(clip = FALSE, ylim = c(0,135)) +
  geom_errorbar(
    mapping = aes(
      x = mod,
      ymin = mean - std,
      ymax = mean + std,
      group = presence,
      color = signif < 0.05
    ),
    position = position_dodge()
  ) +
  geom_errorbarh(
    data = signifData,
    mapping = aes(
      xmin = stage(mod, after_scale = xmin - 0.4),
      xmax = stage(mod, after_scale = xmax + 0.4),
      y = max(mean) + 5,
      group = mod
    ),
    lineend = "square",
    height = 0
  ) +
  geom_text(
    data = signifData,
    mapping = aes(x = mod, y = max(mean) + 6, label = star),
    fontface = "bold"
  )

P

png(
  "PaperFigures/phosphorylation_barplots.png",
  width = 1500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/phosphorylation_barplots.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(D,phosphorylationData,signifData,P,plotData,t,j,high,low,char,p)

##### Nucleosome Class Barplots #####
load("NucleosomeClassDataTable.rda")
D <- D[D$dataset == "Cellular vs. Ice Combined",2:ncol(D)]

reformat_data <- function(scores, truths, class) {
  true <- scores[truths]
  false <- scores[!truths]
  
  t <- t.test(true, false)
  
  p <- data.frame(
    class = c(class, class),
    truth = c(TRUE, FALSE),
    mean = c(mean(true), mean(false)),
    std = c(sd(true) / sqrt(length(true)), sd(false) / sqrt(length(false))),
    signif = rep(t$p.value, 2)
  )
  
  return(p)
}

# Genic
L <- reformat_data(D$score, D$in_gene, "in_gene")
plotData <- L

# Plus One
L <- reformat_data(D$score[D$in_gene], D$X.1[D$in_gene], "plus_one")
plotData <- rbind(plotData, L)

# ARS
L <- reformat_data(D$score, D$ARS, "ARS")
plotData <- rbind(plotData, L)

# tRNA
L <- reformat_data(D$score, D$tRNA, "tRNA")
plotData <- rbind(plotData, L)

# snoRNA
L <- reformat_data(D$score, D$snoRNA, "snoRNA")
plotData <- rbind(plotData, L)

# High vs. Low Expression Genes
scores <- c(D$score[D$in_high_expr_gene], D$score[D$in_low_expr_gene])
truths<-c(rep(TRUE,sum(D$in_high_expr_gene)),rep(FALSE,sum(D$in_low_expr_gene)))
L <- reformat_data(scores, truths, "gene_expr")
plotData <- rbind(plotData, L)

plotData$signif <- plotData$signif * length(unique(plotData$class))

signifData<-data.frame(class=unique(plotData$class),mean=unique(plotData$mean))
chars = c()
for (class in signifData$class) {
  char = "n.s."
  if (plotData$signif[plotData$class == class][1] < 0.05) { char = "*" }
  if (plotData$signif[plotData$class == class][1] < 0.01) { char = "**" }
  if (plotData$signif[plotData$class == class][1] < 0.001) { char = "***" }
  chars = c(chars, char)
}
signifData$star <- chars
signifData <- signifData[signifData$star != "n.s.",]

L <- c("tRNA", "snoRNA", "in_gene", "plus_one", "gene_expr", "ARS")
TrueNames <- c(
  "In tRNA Gene",
  "In snoRNA Gene",
  "In Protein-Coding Gene",
  "+1 Nucleosome",
  "In High-Expression Gene",
  "In ARS"
)
FalseNames <- c(
  "Not in tRNA Gene",
  "Not in snoRNA Gene",
  "Intergenic",
  "Distal Nucleosome",
  "In Low-Expression Gene",
  "Not in ARS"
)
plotData$class <- factor(plotData$class, levels = L)
plotData$truth <- factor(plotData$truth, levels = c(TRUE, FALSE))

P <- ggplot(data = plotData) +
  geom_col(
    mapping = aes(
      x = class,
      y = mean,
      fill = signif < 0.05,
      group = truth,
      color = signif < 0.05
    ),
    position = position_dodge(),
  ) +
  plotTheme +
  scale_fill_discrete(type = c("grey", "white")) +
  scale_color_discrete(type = c("#6D6D6D", "black")) +
  guides(fill = "none", color = "none") +
  scale_y_continuous(
    name = "Mean Rotational Score",
    expand = c(0, 0),
    limits = c(-10, 120),
    breaks = c(0,40,80,120)
  ) +
  scale_x_discrete(name = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(margin = margin(t = 80)),
    axis.title.y = element_text(hjust = 0.6)
  ) +
  geom_errorbar(
    mapping = aes(
      x = class,
      ymin = mean - std,
      ymax = mean + std,
      group = truth,
      color = signif < 0.05
    ),
    position = position_dodge()
  ) +
  geom_errorbarh(
    data = signifData,
    mapping = aes(
      xmin = stage(class, after_scale = xmin - 0.4),
      xmax = stage(class, after_scale = xmax + 0.4),
      y = max(mean) + 3
    ),
    lineend = "square",
    height = 0
  ) +
  geom_text(
    data = signifData,
    mapping = aes(x = class, y = max(mean) + 4, label = star),
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off", ylim = c(0, 120)) +
  geom_text(
    data = data.frame(
      class = rep(L, 2),
      label = c(TrueNames, FalseNames)
    ),
    mapping = aes(
      x = stage(class, after_scale = x + rep(c(-0.2, 0.2), each = length(L))),
      label = label,
      y = -3
    ),
    fontface = "bold",
    angle = 45,
    hjust = 1
  )

P

png(
  "PaperFigures/class_barplots.png",
  width = 1700,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/class_barplots.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(D,L,signifData,P,plotData,t,j,high,low,char,p,TrueNames,FalseNames,
   reformat_data,truths,scores,chars,class)

##### Nucleosome Rotational Score Histogram #####
load("NucleosomeClassDataTable.rda")
D <- D[D$dataset == "Cellular vs. Ice Combined",2:ncol(D)]

P <- ggplot(data = D) +
  geom_histogram(mapping = aes(x = score), bins = 200) +
  plotTheme +
  scale_x_continuous(
    name = "Rotational Score",
    expand = c(0,0),
    breaks = seq(-100, 550, 100)
  ) +
  scale_y_continuous(name = "Count", expand = c(0,0), limits = c(0,1500))

P

png(
  "PaperFigures/rotational_score_histogram.png",
  width = 2000,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  D$score,
  file = "PaperFigures/FigureData/rotational_score_histogram.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(P, D)

##### H3K4me3 Mean Score Plot #####
load("NucleosomeClassDataTable.rda")
D <- D[D$dataset == "Cellular vs. Ice Combined",2:ncol(D)]

chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "genome_scores/chromosomes.txt"
  ),
  header = FALSE
)$V1
scores <- list()
for (chrom in chroms) {
  scores[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/genome_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[chrom]]) <- c("position", "score")
}

W <- c(-100,100) # The window width to graph around the nucleosome dyad

get_mean_scores <- function(scores_list, dyads, window, chroms) {
  total <- rep(0, window[2]-window[1] + 1)
  n <- 0
  for (chrom in chroms) {
    print(paste0("> ", chrom, "..."))
    theseDyads <- dyads$dyad[dyads$chrom == chrom]
    for (i in 1:length(theseDyads)) {
      dyad <- theseDyads[i]
      theseScores <- scores_list[[chrom]]$score[
        (scores_list[[chrom]]$position >= dyad + window[1]) &
          (scores_list[[chrom]]$position <= dyad + window[2])
      ]
      
      if (length(theseScores) == window[2] - window[1] + 1) {
        total <- total + theseScores
        n <- n + 1
      }
    }
  }
  
  return(
    data.frame(
      relative_position = window[1]:window[2],
      average_score = total / n
    )
  )
}

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

modData <- D[!is.na(D[16]),]

mod <- "H3K4me3"

highThresh = quantile(modData[[mod]], 0.75)
highDyads <- modData[modData[[mod]] > highThresh,c(1,2)]
highData <- get_mean_scores(scores, highDyads, W, chroms)

lowThresh <- quantile(modData[[mod]], 0.25)
lowDyads <- modData[modData[[mod]] < lowThresh,c(1,2)]
lowData <- get_mean_scores(scores, lowDyads, W, chroms)

plotData <- cbind(
  rbind(highData, lowData),
  presence = c(rep("+", nrow(highData)), rep("-", nrow(lowData)))
)

P <- ggplot() +
  geom_line(
    data = annotationData[annotationData$x != 0,],
    mapping = aes(x = x, y = y, group = group),
    linetype = "dotted",
    linewidth = 0.75
  ) +
  annotate(
    geom = "line",
    linewidth = 1.25,
    linetype = "dashed",
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = relative_position, y = average_score, color = presence),
    linewidth = 1
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.6,2.9),
    breaks = seq(-2.5,2.9,0.5),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = mod,
    type = c("blue", "red"),
    labels = c("🡇", "🡅")
  ) +
  plotTheme +
  theme(legend.direction="horizontal", legend.title=element_text(face="bold"))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  paste0("PaperFigures/", mod, "_score.png"),
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = paste0("PaperFigures/FigureData/", mod, "_score.tsv"),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(P,G,D,scores,highData,highDyads,lowData,lowDyads,plotData,modData,
   chrom,chroms,mod,highThresh,lowThresh,get_mean_scores,annotationData,W)

##### ARS Mean Score Plot #####
load("NucleosomeClassDataTable.rda")
D <- D[D$dataset == "Cellular vs. Ice Combined",2:ncol(D)]

chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "genome_scores/chromosomes.txt"
  ),
  header = FALSE
)$V1
scores <- list()
for (chrom in chroms) {
  scores[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/genome_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[chrom]]) <- c("position", "score")
}

W <- c(-100,100) # The window width to graph around the nucleosome dyad

get_mean_scores <- function(scores_list, dyads, window, chroms) {
  total <- rep(0, window[2]-window[1] + 1)
  n <- 0
  for (chrom in chroms) {
    print(paste0("> ", chrom, "..."))
    theseDyads <- dyads$dyad[dyads$chrom == chrom]
    for (i in 1:length(theseDyads)) {
      dyad <- theseDyads[i]
      theseScores <- scores_list[[chrom]]$score[
        (scores_list[[chrom]]$position >= dyad + window[1]) &
          (scores_list[[chrom]]$position <= dyad + window[2])
      ]
      
      if (length(theseScores) == window[2] - window[1] + 1) {
        total <- total + theseScores
        n <- n + 1
      }
    }
  }
  
  return(
    data.frame(
      relative_position = window[1]:window[2],
      average_score = total / n
    )
  )
}

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

class <- "ARS"

trueDyads <- D[D[[class]], c(1,2)]
trueData <- get_mean_scores(scores, highDyads, W, chroms)

falseDyads <- D[!D[[class]], c(1,2)]
falseData <- get_mean_scores(scores, lowDyads, W, chroms)

plotData <- cbind(
  rbind(trueData, falseData),
  presence = c(rep("+", nrow(trueData)), rep("-", nrow(falseData)))
)

P <- ggplot() +
  geom_line(
    data = annotationData[annotationData$x != 0,],
    mapping = aes(x = x, y = y, group = group),
    linetype = "dotted",
    linewidth = 0.75
  ) +
  annotate(
    geom = "line",
    linewidth = 1.25,
    linetype = "dashed",
    x = c(0,0),
    y = c(-Inf, Inf)
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = relative_position, y = average_score, color = presence),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.4,2.7),
    breaks = seq(-2.4,2.7,0.6),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = class,
    type = c("blue", "red"),
    labels = c("Not in ARS", "In ARS")
  ) +
  plotTheme +
  theme(legend.direction="horizontal", legend.title=element_blank())

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  paste0("PaperFigures/", class, "_score.png"),
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = paste0("PaperFigures/FigureData/", class, "_score.tsv"),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(P,G,D,scores,trueData,trueDyads,falseData,falseDyads,plotData,modData,
   chrom, chroms, class, get_mean_scores, annotationData, W)

##### Venn Diagrams #####
data <- read.table("VennDiagrams.txt",sep="\t")
colnames(data) <- c("A", "B", "C", "A&B", "A&C", "B&C", "A&B&C")
data$A <- as.integer(str_split_i(data$A, ": ", 2))

E <- euler(unlist(as.vector(data[1,])), shape = "ellipse")

P <- plot(
  E,
  fills = list(fill=c("red", "green", "blue"), alpha=c(0.25,0.25,0.25)),
  edges = list(lex = 2),
  labels = NULL
);P

png(
  paste0("PaperFigures/venn_greedy_weiner_brogaard.png"),
  width = 1000,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

E <- euler(unlist(as.vector(data[2,])), shape = "ellipse")

P <- plot(
  E,
  fills = list(fill=c("red", "green", "blue"), alpha=c(0.25,0.25,0.25)),
  edges = list(lex = 2),
  labels = NULL
);P

png(
  paste0("PaperFigures/venn_viterbi_weiner_brogaard.png"),
  width = 1000,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

E <- euler(unlist(as.vector(data[3,])), shape = "ellipse")

P <- plot(
  E,
  fills = list(fill=c("red", "green", "blue"), alpha=c(0.25,0.25,0.25)),
  edges = list(lex = 2),
  labels = NULL
);P

png(
  paste0("PaperFigures/venn_viterbi20bp_weiner_brogaard.png"),
  width = 1000,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

rownames(data) <- c(
  "Greedy_Weiner_Brogaard",
  "Viterbi_Weiner_Brogaard",
  "Viterbi20bp_Weiner_Brogaard"
)
write.table(
  data,
  file = paste0("PaperFigures/FigureData/venn.tsv"),
  col.names = TRUE,
  row.names = TRUE,
  sep = "\t",
  quote = FALSE
)

##### Positions Relative to Redundant Dyads #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/redundant_freqTable.tsv"
  ),
  sep = "\t"
)

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/detrended_periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

# > Greedy ----
plotData <- data.frame(
  position = rep(-100:100, 2),
  frequency = c(data$greedyVsRedund, data$weinerVsRedund),
  dataset = rep(c("greedy", "weiner"), each = 201)
)
plotData <- plotData[(plotData$position >= -80) & (plotData$position <= 80),]

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "black"),
    labels = c("Greedy", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 6300),
    breaks = seq(0, 6300, 900),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_redundant_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/greedy_redundant_dyad_position.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/greedy_redundant_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

# > Viterbi ----
plotData <- data.frame(
  position = rep(-100:100, 2),
  frequency = c(data$viterbiVsRedund, data$weinerVsRedund),
  dataset = rep(c("viterbi", "weiner"), each = 201)
)
plotData <- plotData[(plotData$position >= -80) & (plotData$position <= 80),]

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = position, y = frequency, color = dataset),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("#990099", "black"),
    labels = c("Viterbi", "MNase-seq")
  ) +
  scale_x_continuous(
    name = "Position Relative to Reference Dyad (bp)",
    limits = c(-80, 80),
    breaks = -4:4*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 6600),
    breaks = seq(0, 6600, 1100),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_redundant_dyad_position.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0("#Period: ", toString(periodicity$period[1]))
write(string, "PaperFigures/FigureData/viterbi_redundant_dyad_position.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_redundant_dyad_position.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE,
)

rm(annotationData, data, G, P, periodicity, plotData, string)

##### Offsets Relative to Redundant Dyads #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/redundant_offsetFreqTable.tsv"
  ),
  sep = "\t"
)[1:4]
colnames(data) <- c("offset", "greedy", "viterbi", "weiner")

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

write.table(
  data,
  file = "PaperFigures/FigureData/redundant_dyad_offsets.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Greedy ----
P <- ggplot(data = data) +
  geom_line(
    mapping = aes(x = offset, y = greedy),
    linewidth = 1.25
  ) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 6500),
    breaks = seq(0, 6500, 1300),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/greedy_redundant_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

# > Viterbi ----
P <- ggplot(data = data) +
  geom_line(
    mapping = aes(x = offset, y = viterbi),
    linewidth = 1.25
  ) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 7000),
    breaks = seq(0, 7000, 1000),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_redundant_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

# > Weiner ----
P <- ggplot(data = data) +
  geom_line(
    mapping = aes(x = offset, y = weiner),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 4500),
    breaks = seq(0, 4500, 500),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/weiner_redundant_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

rm(data, G, P)

##### Brogaard et al. Redundant Center-to-Center Plot #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/redundant_center_to_centers.tsv"
  ),
  sep = "\t",
  header = FALSE
)
colnames(data) <- c("distance", "frequency")
data <- data[2:nrow(data),]

P <- ggplot(data = data) +
  geom_line(
    mapping = aes(x = distance, y = frequency),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Center to Center Distance (bp)",
    limits = c(0, 100),
    breaks = 0:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 40000),
    breaks = seq(0, 40000, 10000),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/redundant_center_to_center.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  data,
  file = "PaperFigures/FigureData/redundant_center_to_center.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(data, P, G)

##### Viterbi-Redundant Center-to-Center Plot #####

data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/redundant_freqTable.tsv"
  ),
  sep = "\t"
)

positive <- data$viterbiVsRedund[data$distance > 0]
negative <- data$viterbiVsRedund[data$distance < 0]
plotData <- data.frame(
  distance = 1:100,
  frequency = positive + negative[100:1]
)

P <- ggplot(data = plotData) +
  geom_line(
    mapping = aes(x = distance, y = frequency),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Center to Center Distance (bp)",
    limits = c(0, 100),
    breaks = 0:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 10000),
    breaks = seq(0, 10000, 2500),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/redundant_viterbi_center_to_center.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/redundant_viterbi_center_to_center.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(data, P, G, plotData, positive, negative)

##### Dinucleotide Frequency Plots #####
directory <- "NucleosomePattern/composite_models/composite_ice"

N <- read.table(paste0(directory, "/all_sequence_analysis/n.txt"))$V1
all_WW = read.table(
  paste0(directory, "/all_sequence_analysis/WW.tsv"),
  sep = "\t"
)$V1
all_SS = read.table(
  paste0(directory, "/all_sequence_analysis/SS.tsv"),
  sep = "\t"
)$V1

P_N <- read.table(paste0(directory, "/P1_sequence_analysis/n.txt"))$V1
P1_WW = read.table(
  paste0(directory, "/P1_sequence_analysis/WW.tsv"),
  sep = "\t"
)$V1
P1_SS = read.table(
  paste0(directory, "/P1_sequence_analysis/SS.tsv"),
  sep = "\t"
)$V1
P3_WW = read.table(
  paste0(directory, "/P3_sequence_analysis/WW.tsv"),
  sep = "\t"
)$V1
P3_SS = read.table(
  paste0(directory, "/P3_sequence_analysis/SS.tsv"),
  sep = "\t"
)$V1

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

# > WW Frequency ----

plotData <- data.frame(
  frequency = c(P1_WW, P3_WW),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("P1", "P9"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(1700, 2900),
    breaks = c(seq(1700, 2900, 200)),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL,
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  theme(legend.direction="horizontal", legend.position.inside=c(0.025, 1.05))

P

png(
  "PaperFigures/WW_frequency.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW Frequency (Weak) ----

plotData <- data.frame(
  frequency = P1_WW,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(1700, 2900),
    breaks = c(seq(1700, 2900, 200)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/WW_frequency_weak.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW Frequency (Strong) ----

plotData <- data.frame(
  frequency = P3_WW,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(1700, 2900),
    breaks = c(seq(1700, 2900, 200)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/WW_frequency_strong.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency ----

plotData <- data.frame(
  frequency = c(P1_SS, P3_SS),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("P1", "P9"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(700, 1200),
    breaks = c(seq(700, 1200, 100)),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL,
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  theme(legend.direction="horizontal", legend.position.inside=c(0.025, 1.05))

P

png(
  "PaperFigures/SS_frequency.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Weak) ----

plotData <- data.frame(
  frequency = P1_SS,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(700, 1200),
    breaks = c(seq(700, 1200, 100)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/SS_frequency_weak.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Strong) ----

plotData <- data.frame(
  frequency = P3_SS,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(700, 1200),
    breaks = c(seq(700, 1200, 100)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/SS_frequency_strong.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW vs. SS Frequency (Weak) ----

plotData <- data.frame(
  prop = c(P1_WW, (P1_SS - 700)*12/5 + 1700),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("WW", "SS"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = prop, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(1700, 2900),
    breaks = c(seq(1700, 2900, 200)),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 5/12 + 700 - 1700*5/12,
      name = "GG/GC/CG/CC Frequency",
      breaks = c(700.00000000001, #GGPLOT WHY???? 
                 seq(800, 1200, 100))
    )
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL
  ) +
  guides(color = guide_none()) +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.y.left = element_text(colour = "blue"),
    axis.title.y.right = element_text(colour = "red")
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/WW_SS_frequency_weak.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW vs. SS Frequency (Strong) ----

plotData <- data.frame(
  prop = c(P3_WW, (P3_SS - 700)*12/5 + 1700),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("WW", "SS"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = prop, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(1700, 2900),
    breaks = c(seq(1700, 2900, 200)),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 5/12 + 700 - 1700*5/12,
      name = "GG/GC/CG/CC Frequency",
      breaks = c(700.00000000001, #GGPLOT WHY???? 
                 seq(800, 1200, 100))
    )
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL
  ) +
  guides(color = guide_none()) +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.y.left = element_text(colour = "blue"),
    axis.title.y.right = element_text(colour = "red")
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/WW_SS_frequency_strong.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Periodograms ----
L_P1_WW <- lsp(P1_WW, times=-73:72+0.5, type="period", to=15, ofac=100)
L_P3_WW <- lsp(P3_WW, times=-73:72+0.5, type="period", to=15, ofac=100)
L_P1_SS <- lsp(P1_SS, times=-73:72+0.5, type="period", to=15, ofac=100)
L_P3_SS <- lsp(P3_SS, times=-73:72+0.5, type="period", to=15, ofac=100)

# >> WW Periodogram ----
plotData <- data.frame(
  period = c(L_P1_WW$scanned, L_P3_WW$scanned),
  normalized_power = c(L_P1_WW$power, L_P3_WW$power),
  group=c(rep("P1",length(L_P1_WW$scanned)),rep("P9",length(L_P3_WW$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_P1_WW$peak.at[1], L_P3_WW$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("P1", "P9"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group),
    linetype = "dashed",
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  annotate(
    "text",
    x = 4,
    y = 0.7,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[3],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("AA/AT/TA/TT Frequency"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_frequency_periodogram.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#P1: ",
  toString(annotationData$x[1]),
  "\n#P9: ",
  toString(annotationData$x[3])
)

write(string, "PaperFigures/FigureData/WW_frequency_periodogram.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_periodogram.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

# >> SS Periodogram ----
plotData <- data.frame(
  period = c(L_P1_SS$scanned, L_P3_SS$scanned),
  normalized_power = c(L_P1_SS$power, L_P3_SS$power),
  group=c(rep("P1",length(L_P1_SS$scanned)),rep("P9",length(L_P3_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_P1_SS$peak.at[1], L_P3_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("P1", "P9"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group),
    linetype = "dashed",
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  annotate(
    "text",
    x = 13.4,
    y = 0.7,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "right",
    vjust = "top",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10,
    y = 0.7,
    label = paste0(round(annotationData$x[3],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "right",
    vjust = "top",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("GG/GC/CG/CC Frequency"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/SS_frequency_periodogram.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#P1: ",
  toString(annotationData$x[1]),
  "\n#P9: ",
  toString(annotationData$x[3])
)

write(string, "PaperFigures/FigureData/SS_frequency_periodogram.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_periodogram.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

# >> Low Rotational Score Periodogram ----
plotData <- data.frame(
  period = c(L_P1_WW$scanned, L_P1_SS$scanned),
  normalized_power = c(L_P1_WW$power, L_P1_SS$power),
  group=c(rep("WW",length(L_P1_WW$scanned)),rep("SS",length(L_P1_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_P1_WW$peak.at[1], L_P1_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("WW", "SS"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_SS_frequency_periodogram_weak.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_periodogram_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# >> High Rotational Score Periodogram ----
plotData <- data.frame(
  period = c(L_P3_WW$scanned, L_P3_SS$scanned),
  normalized_power = c(L_P3_WW$power, L_P3_SS$power),
  group=c(rep("WW",length(L_P3_WW$scanned)),rep("SS",length(L_P3_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_P3_WW$peak.at[1], L_P3_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("WW", "SS"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group, linetype = group),
    linewidth = 0.75
  ) +
  scale_linetype_manual(
    values = c("dotted", "dashed")
  ) +
  guides(linetype = guide_none()) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[3],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "bottom",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_SS_frequency_periodogram_strong.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#WW: ",
  toString(annotationData$x[1]),
  "\n#SS: ",
  toString(annotationData$x[3])
)

write(string, "PaperFigures/FigureData/WW_SS_frequency_periodogram_strong.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_periodogram_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(annotationData,plotData,P,all_WW,all_SS,N,P_N,Q,Q1_SS,Q1_WW,P1_SS,P1_WW,G,
   Q3_WW,Q3_SS,P3_WW,P3_SS,directory,string,L_P1_SS,L_P1_WW,L_P3_WW,L_P3_SS)

##### Dinucleotide Frequency Plots (Quartiles) #####
directory <- "NucleosomePattern/composite_models/composite_ice"

N <- read.table(paste0(directory, "/all_sequence_analysis/n.txt"))$V1
all_WW = read.table(
  paste0(directory, "/all_sequence_analysis/WW.tsv"),
  sep = "\t"
)$V1
all_SS = read.table(
  paste0(directory, "/all_sequence_analysis/SS.tsv"),
  sep = "\t"
)$V1

Q <- read.table(paste0(directory, "/Q1_sequence_analysis/n.txt"))$V1
Q1_WW = read.table(
  paste0(directory, "/Q1_sequence_analysis/WW.tsv"),
  sep = "\t"
)$V1
Q1_SS = read.table(
  paste0(directory, "/Q1_sequence_analysis/SS.tsv"),
  sep = "\t"
)$V1
Q3_WW = read.table(
  paste0(directory, "/Q3_sequence_analysis/WW.tsv"),
  sep = "\t"
)$V1
Q3_SS = read.table(
  paste0(directory, "/Q3_sequence_analysis/SS.tsv"),
  sep = "\t"
)$V1

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

# > WW Frequency ----

plotData <- data.frame(
  frequency = c(Q1_WW, Q3_WW),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("Q1", "Q3"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(4400, 7100),
    breaks = c(seq(4400, 7100, 900)),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL,
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  theme(legend.direction="horizontal", legend.position.inside=c(0.025, 1.05))

P

png(
  "PaperFigures/WW_frequency_quartiles.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_quartiles.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW Frequency (Weak) ----

plotData <- data.frame(
  frequency = Q1_WW,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(4400, 7100),
    breaks = c(seq(4400, 7100, 900)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/WW_frequency_Q1.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_Q1.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW Frequency (Strong) ----

plotData <- data.frame(
  frequency = Q3_WW,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(4400, 7100),
    breaks = c(seq(4400, 7100, 900)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/WW_frequency_Q3.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_Q3.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency ----

plotData <- data.frame(
  frequency = c(Q1_SS, Q3_SS),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("Q1", "Q3"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(1800, 2800),
    breaks = c(seq(1800, 2800, 200)),
    expand = c(0,0)
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL,
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  theme(legend.direction="horizontal", legend.position.inside=c(0.025, 1.05))

P

png(
  "PaperFigures/SS_frequency_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Weak) ----

plotData <- data.frame(
  frequency = Q1_SS,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(1800, 2800),
    breaks = c(seq(1800, 2800, 200)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/SS_frequency_Q1.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_Q1.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Strong) ----

plotData <- data.frame(
  frequency = Q3_SS,
  position = seq(-72.5, 72.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(1800, 2800),
    breaks = c(seq(1800, 2800, 200)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/SS_frequency_Q3.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_Q3.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW vs. SS Frequency (Weak) ----

plotData <- data.frame(
  prop = c(Q1_WW, (Q1_SS - 1800)*27/10 + 4400),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("WW", "SS"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = prop, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(4400, 7100),
    breaks = c(seq(4400, 7100, 900)),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 10/27 + 1800 - 4400*10/27,
      name = "GG/GC/CG/CC Frequency",
      breaks = seq(1800, 2800, 200)
    )
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL
  ) +
  guides(color = guide_none()) +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.y.left = element_text(colour = "blue"),
    axis.title.y.right = element_text(colour = "red")
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/WW_SS_frequency_Q1.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_Q1.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW vs. SS Frequency (Strong) ----

plotData <- data.frame(
  prop = c(Q3_WW, (Q3_SS - 1800)*27/10 + 4400),
  position = rep(seq(-72.5, 72.5, 1), 2),
  group = rep(c("WW", "SS"), each = 146)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = prop, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = -7:7 * 10,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(4400, 7100),
    breaks = c(seq(4400, 7100, 900)),
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 10/27 + 1800 - 4400*10/27,
      name = "GG/GC/CG/CC Frequency",
      breaks = seq(1800, 2800, 200)
    )
  ) +
  scale_color_discrete(
    type = c("red", "blue"),
    name = NULL
  ) +
  guides(color = guide_none()) +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.y.left = element_text(colour = "blue"),
    axis.title.y.right = element_text(colour = "red")
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/WW_SS_frequency_Q3.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_Q3.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Periodograms ----
L_Q1_WW <- lsp(Q1_WW, times=-73:72+0.5, type="period", to=15, ofac=100)
L_Q3_WW <- lsp(Q3_WW, times=-73:72+0.5, type="period", to=15, ofac=100)
L_Q1_SS <- lsp(Q1_SS, times=-73:72+0.5, type="period", to=15, ofac=100)
L_Q3_SS <- lsp(Q3_SS, times=-73:72+0.5, type="period", to=15, ofac=100)

# >> WW Periodogram ----
plotData <- data.frame(
  period = c(L_Q1_WW$scanned, L_Q3_WW$scanned),
  normalized_power = c(L_Q1_WW$power, L_Q3_WW$power),
  group=c(rep("Q1",length(L_Q1_WW$scanned)),rep("Q3",length(L_Q3_WW$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_Q1_WW$peak.at[1], L_Q3_WW$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("Q1", "Q3"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group),
    linetype = "dashed",
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  annotate(
    "text",
    x = 4,
    y = 0.7,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[3],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("AA/AT/TA/TT Frequency"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_frequency_periodogram_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#Q1: ",
  toString(annotationData$x[1]),
  "\n#Q3: ",
  toString(annotationData$x[3])
)

write(string, "PaperFigures/FigureData/WW_frequency_periodogram_quartile.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_periodogram_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

# >> SS Periodogram ----
plotData <- data.frame(
  period = c(L_Q1_SS$scanned, L_Q3_SS$scanned),
  normalized_power = c(L_Q1_SS$power, L_Q3_SS$power),
  group=c(rep("Q1",length(L_Q1_SS$scanned)),rep("Q3",length(L_Q3_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_Q1_SS$peak.at[1], L_Q3_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("Q1", "Q3"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group),
    linetype = "dashed",
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  annotate(
    "text",
    x = 13.4,
    y = 0.7,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "right",
    vjust = "top",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10,
    y = 0.7,
    label = paste0(round(annotationData$x[3],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "right",
    vjust = "top",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("GG/GC/CG/CC Frequency"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/SS_frequency_periodogram_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#Q1: ",
  toString(annotationData$x[1]),
  "\n#Q3: ",
  toString(annotationData$x[3])
)

write(string, "PaperFigures/FigureData/SS_frequency_periodogram_quartile.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_periodogram_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

# >> Low Rotational Score Periodogram ----
plotData <- data.frame(
  period = c(L_Q1_WW$scanned, L_Q1_SS$scanned),
  normalized_power = c(L_Q1_WW$power, L_Q1_SS$power),
  group=c(rep("WW",length(L_Q1_WW$scanned)),rep("SS",length(L_Q1_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_Q1_WW$peak.at[1], L_Q1_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("WW", "SS"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_SS_frequency_periodogram_Q1.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_periodogram_Q1.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# >> High Rotational Score Periodogram ----
plotData <- data.frame(
  period = c(L_Q3_WW$scanned, L_Q3_SS$scanned),
  normalized_power = c(L_Q3_WW$power, L_Q3_SS$power),
  group=c(rep("WW",length(L_Q3_WW$scanned)),rep("SS",length(L_Q3_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_Q3_WW$peak.at[1], L_Q3_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("WW", "SS"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group, linetype = group),
    linewidth = 0.75
  ) +
  scale_linetype_manual(
    values = c("dotted", "dashed")
  ) +
  guides(linetype = guide_none()) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[3],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "bottom",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_SS_frequency_periodogram_Q3.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#WW: ",
  toString(annotationData$x[1]),
  "\n#SS: ",
  toString(annotationData$x[3])
)

write(string, "PaperFigures/FigureData/WW_SS_frequency_periodogram_Q3.tsv")

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_periodogram_Q3.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(annotationData,plotData,P,all_WW,all_SS,N,P_N,Q,Q1_SS,Q1_WW,P1_SS,P1_WW,G,
   Q3_WW,Q3_SS,P3_WW,P3_SS,directory,string,L_P1_SS,L_P1_WW,L_P3_WW,L_P3_SS)


##### A/T Bias at -3/+3 #####
load_seqs <- function(directory) {
  chroms <- read.table(paste0(directory, "/chromosomes.txt"))$V1
  data <- data.frame(chrom = character(0), loc = integer(0))
  Cnms <- c("chrom", "loc")
  for (i in -73:73) {
    data <- cbind(data, V = character(0))
    colnames(data) <- c(Cnms, toString(i))
    Cnms <- colnames(data)
  }
  for (chrom in chroms) {
    seqData <- read.table(paste0(directory, "/", chrom, ".tsv"))
    #data <- c(data, apply(seqData, 1, paste0, collapse = ""))
    seqData <- cbind(chrom = rep(chrom, nrow(seqData)), seqData)
    colnames(seqData) <- Cnms
    data <- rbind(data, seqData)
  }
  
  return(data)
}

viterbi <- load_seqs(
  "NucleosomePattern/composite_models/composite_ice/all_sequence_analysis/seqs"
)
brogaard <- load_seqs(
  "NucleosomePattern/ReferenceSequenceAnalysis/brogaard_all/seqs"
)
weiner <- load_seqs("NucleosomePattern/ReferenceSequenceAnalysis/weiner/seqs")
redund <- load_seqs("NucleosomePattern/ReferenceSequenceAnalysis/redund/seqs")
pombe_unique <- load_seqs("pombe_sequence_analysis/unique/seqs")
pombe_redund <- load_seqs("pombe_sequence_analysis/redundant/seqs")

redund_noUnique <- redund[0,]
for (chrom in unique(redund$chrom)) {
  subset <- redund[redund$chrom == chrom,]
  bad <- brogaard[brogaard$chrom == chrom,]
  redund_noUnique <- rbind(
    redund_noUnique,
    subset[!(subset$loc %in% bad$loc),]
  )
}

mouse_stats <- read.delim("mouse_AT3_data.txt", header = FALSE)

# > A/T Bias Bar Graphs ----

E <- 0.62 - 0.31^2
E_pombe <- 0.64 - 0.32^2
AT_mouse <- as.integer(unlist(str_split(mouse_stats$V1[1], ": "))[2])
N_mouse <- as.numeric(unlist(str_split(mouse_stats$V1[2], ": "))[2])
E_mouse <- AT_mouse / N_mouse - (AT_mouse/(2*N_mouse))^2
plotData <- data.frame(
  prop = c(
    sum((viterbi$`-3` == "A") | (viterbi$`3` == "T")) / nrow(viterbi) - E,
    sum((brogaard$`-3` == "A") | (brogaard$`3` == "T")) / nrow(brogaard) - E,
    sum((weiner$`-3` == "A") | (weiner$`3` == "T")) / nrow(weiner) - E
  ),
  source = c("viterbi", "brogaard", "weiner")
)

labels <- c(
  "**Chemical**<br>**Cleavage**",
  "**CPD-seq**",
  "**MNase-seq**"
)

P <- ggplot(data = plotData) +
  geom_col(
    mapping = aes(
      x = factor(source, levels = c("brogaard", "viterbi", "weiner")),
      y = prop
    ),
    color = "black",
    fill = "white",
    linewidth = 1.25,
    lineend = "square"
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0,0),
    linetype = "solid",
    color = "black",
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_discrete(
    name = NULL,
    labels = labels
  ) +
  scale_y_continuous(
    name = "Proportion of Nucleosomes\nwith A/T Nucleotide at -3/+3",
    limits = c(0.5 - E, 1 - E),
    breaks = seq(0.5-E, 1-E, 0.1),
    labels = paste0(unlist(lapply(5:10*10, toString)), "%"),
    expand = c(0,0)
  ) +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(hjust = 1),
    plot.margin = margin(t = 10, b = 10, l = 10, r = 60)
  ) +
  coord_cartesian(clip = "off") +
  annotate(
    geom = "text",
    label = "🡄",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.4,
    size = 5
  ) +
  annotate(
    geom = "text",
    label = "    Random",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.5
  )


P

png(
  "PaperFigures/AT3_bars.png",
  width = 1250,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

string <- paste0("# Expected Proportion: ", toString(E), "\n")

write(string, "PaperFigures/FigureData/AT3_bars.tsv")
write.table(
  plotData,
  file = "PaperFigures/FigureData/AT3_bars.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

# > Viterbi Logo Near Dyad ----

seqs <- character(nrow(viterbi))
for (i in -5:5) {
  seqs <- paste0(seqs, viterbi[[toString(i)]])
}

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg
) +
  scale_y_continuous(
    limits = c(0, 0.4),
    expand = c(0,0),
    name = "Relative Entropy (bits)"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("CPD-seq Dyads"))))

P

png(
  "PaperFigures/viterbi_seqLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

write.table(
  seqs,
  file = "PaperFigures/FigureData/viterbi_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-4000,6000),
    breaks = seq(-4000, 6000, 2000)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("CPD-seq Dyads"))))

P

png(
  "PaperFigures/viterbi_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Brogaard Logo Near Dyad ----

seqs <- character(nrow(brogaard))
for (i in -5:5) {
  seqs <- paste0(seqs, brogaard[[toString(i)]])
}

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg
) +
  scale_y_continuous(
    limits = c(0, 0.4),
    expand = c(0,0),
    name = "Relative Entropy (bits)"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(underline(bold("Chemical Cleavage Dyads")))
  )

P

png(
  "PaperFigures/brogaard_seqLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

write.table(
  seqs,
  file = "PaperFigures/FigureData/brogaard_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-4000,6000)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(underline(bold("Chemical Cleavage Dyads")))
  )

P

png(
  "PaperFigures/brogaard_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > MNase Logo Near Dyad ----

seqs <- character(nrow(weiner))
for (i in -5:5) {
  seqs <- paste0(seqs, weiner[[toString(i)]])
}

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg
) +
  scale_y_continuous(
    limits = c(0, 0.4),
    expand = c(0,0),
    name = "Relative Entropy (bits)"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("MNase-Seq"))))

P

png(
  "PaperFigures/weiner_seqLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

write.table(
  seqs,
  file = "PaperFigures/FigureData/weiner_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-4000,6000),
    breaks = seq(-4000, 6000, 2000)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("MNase-Seq"))))

P

png(
  "PaperFigures/weiner_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Redundant Logo Near Dyad ----

seqs <- character(nrow(redund))
for (i in -5:5) {
  seqs <- paste0(seqs, redund[[toString(i)]])
}

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg
) +
  scale_y_continuous(
    limits = c(0, 0.4),
    expand = c(0,0),
    name = "Relative Entropy (bits)"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(underline(bold("Chemical Cleavage Redundant Dyads")))
  )

P

png(
  "PaperFigures/redundant_seqLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

write.table(
  seqs,
  file = "PaperFigures/FigureData/redundant_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-10000, 20000),
    breaks = seq(-10000, 20000, 5000)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(underline(bold("Brogaard"~bolditalic("et al.")~"Redundant Dyads")))
  )

P

png(
  "PaperFigures/redundant_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Pombe Logo Near Dyad ----

seqs <- character(nrow(pombe_unique))
for (i in -5:5) {
  seqs <- paste0(seqs, pombe_unique[[toString(i)]])
}
seqs <- seqs[!grepl("N", seqs)]

P <- plot_sequence_logo(
  seqs,
  bg_freqs = c("A"=0.32, "T"=0.32, "C"=0.18, "G"=0.18)
) +
  scale_y_continuous(
    limits = c(0, 0.4),
    expand = c(0,0),
    name = "Relative Entropy (bits)"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(
      underline(bold(bolditalic("Schizosaccharomyces pombe")~"Unique Dyads"))
    )
  )

P

png(
  "PaperFigures/pombe_seqLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

write.table(
  seqs,
  file = "PaperFigures/FigureData/pombe_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = c("A"=0.32, "T"=0.32, "C"=0.18, "G"=0.18),
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-4000, 5000),
    breaks = seq(-4000, 5000, 1000)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(
      underline(bold(bolditalic("Schizosaccharomyces pombe")~"Unique Dyads"))
    )
  )

P

png(
  "PaperFigures/pombe_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Pombe Redundant Logo Near Dyad ----

seqs <- character(nrow(pombe_redund))
for (i in -5:5) {
  seqs <- paste0(seqs, pombe_redund[[toString(i)]])
}
seqs <- seqs[!grepl("N", seqs)]

P <- plot_sequence_logo(
  seqs,
  bg_freqs = c("A"=0.32, "T"=0.32, "C"=0.18, "G"=0.18)
) +
  scale_y_continuous(
    limits = c(0, 0.4),
    expand = c(0,0),
    name = "Relative Entropy (bits)"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(
      underline(bold(bolditalic("Schizosaccharomyces pombe")~"Redundant Dyads"))
    )
  )

P

png(
  "PaperFigures/pombe_redundant_seqLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

write.table(
  seqs,
  file = "PaperFigures/FigureData/pombe_redundant_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = c("A"=0.32, "T"=0.32, "C"=0.18, "G"=0.18),
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-8000, 12000),
    breaks = seq(-8000, 12000, 4000)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(
    bquote(
      underline(bold(bolditalic("Schizosaccharomyces pombe")~"Redundant Dyads"))
    )
  )

P

png(
  "PaperFigures/pombe_redundant_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Mouse Logo Near Dyad ----

seqs <- matrix(0, nrow=4,ncol=11)
rownames(seqs) <- c("A", "C", "G", "T")
colnames(seqs) <- unlist(lapply(-5:5, toString))
for (chrom in list.files("mouse_sequence_analysis", ".tsv")) {
  print(paste0("> ", chrom, "..."))
  theseSeqs <- read.table(
    paste0("mouse_sequence_analysis/", chrom),
    sep = "\t"
  )
  colnames(theseSeqs) <- c("loc", unlist(lapply(-73:73, toString)))
  for (i in -5:5) {
    tab <- table(theseSeqs[[toString(i)]])
    for (N in rownames(seqs)) {
      seqs[N, toString(i)] <- seqs[N, toString(i)] + tab[N]
    }
  }
}

AT <- as.numeric(unlist(str_split(mouse_stats$V1[1],": "))[2])
N <- as.numeric(unlist(str_split(mouse_stats$V1[2],": "))[2])
A <- (AT/N)/2
P <- plot_sequence_logo(
  seqs,
  bg_freqs = c("A"=A, "T"=A, "C"=(1-2*A)/2, "G"=(1-2*A)/2)
) +
  scale_y_continuous(
    name = "Relative Entropy (bits)",
    limits = c(0,0.4),
    expand = c(0,0)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Mouse Nucleosomes"))))

P

png(
  "PaperFigures/mouse_seqLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

write.table(
  seqs,
  file = "PaperFigures/FigureData/mouse_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = c("A"=A, "T"=A, "C"=(1-2*A)/2, "G"=(1-2*A)/2),
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-200000, 300000),
    
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:11,
    labels = -5:5,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Mouse Nucleosomes"))))

P

png(
  "PaperFigures/mouse_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Offset Plots ----

data <- data.frame(chrom = character(0), loc = integer(0), dist = integer(0))
chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice",
    "/viterbi_dyad_calls/chromosomes.txt"
  )
)$V1
for (chrom in chroms) {
  theseDyads <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/viterbi_dyad_calls/",
      chrom,
      ".tsv"
    )
  )$V1
  data <- rbind(
    data,
    data.frame(
      chrom = rep(chrom, length(theseDyads)),
      loc = theseDyads,
      dist = read.table(
        paste0(
          "NucleosomePattern/composite_models/composite_ice/distance_evals",
          "/all_viterbi_distances/",
          chrom,
          ".tsv"
        )
      )$V1
    )
  )
}

FALSE %in% (viterbi$chrom == data$chrom) # Should be FALSE
FALSE %in% (viterbi$loc == data$loc) # Should be FALSE
data <- cbind(data, viterbi[3:ncol(viterbi)])

data <- data[abs(data$dist) <= 100,]

annotationData <- data.frame(
  x = rep(10.1*(-9:9), 2),
  y = rep(c(-Inf, Inf), each = 19),
  group = rep(-9:9, 2)
)

A.3 <- data[data$`-3` == "A",]
T3 <- data[data$`3` == "T",]
AorT <- data[(data$`-3` == "A") | (data$`3` == "T"),]

NoAnorT <- data[(data$`-3` != "A") & (data$`3` != "T"),]

A_table <- table(A.3$dist)
T_table <- table(T3$dist)
AorT_table <- table(AorT$dist)

NoAnorT_table <- table(NoAnorT$dist)

N_AorT <- nrow(AorT)
N <- nrow(data)

# >> A at -3 Offsets ----
plotData <- data.frame(
  dist = as.integer(names(A_table)),
  frequency = as.integer(A_table)
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = dist, y = frequency),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1500),
    breaks = seq(0, 1500, 500),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_Aminus3_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_Aminus3_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# >> T at +3 Offsets ----
plotData <- data.frame(
  dist = as.integer(names(T_table)),
  frequency = as.integer(T_table)
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = dist, y = frequency),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1500),
    breaks = seq(0, 1500, 500),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_Tplus3_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_Tplus3_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# >> A/T at -3/+3 Offsets ----
plotData <- data.frame(
  dist = as.integer(names(AorT_table)),
  frequency = as.integer(AorT_table)
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = dist, y = frequency),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1800),
    breaks = seq(0, 1800, 600),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_AorT_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_AorT_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# >> No A/T at -3/+3 Offsets ----
plotData <- data.frame(
  dist = as.integer(names(NoAnorT_table)),
  frequency = as.integer(NoAnorT_table)
)

P <- ggplot(data = plotData) +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    mapping = aes(x = dist, y = frequency),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Offset from Nearest Reference Dyad (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1800),
    breaks = seq(0, 1800, 600),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/viterbi_NoAnorT_dyad_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_NoAnorT_dyad_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Reversed Brogaard Unique Offset Plot ----
D<-"NucleosomePattern/composite_models/composite_ice/reverse_AT3_distance_evals"
data <- read.table(paste0(D, "/brogaard_AT3.tsv"))
colnames(data) <- c("distance", "frequency")
data$ref <- "AT3"

data2 <- read.table(paste0(D, "/brogaard_notAT3.tsv"))
colnames(data2) <- c("distance", "frequency")
data2$ref <- "no_AT3"

plotData <- data.frame(
  distance = rep(-100:100, 2),
  frequency = c(
    data[abs(data$distance) <= 100,]$frequency / N_AorT,
    data2[abs(data2$distance) <= 100,]$frequency / (N - N_AorT)
  ),
  ref = rep(c("AT3", "no_AT3"), each = 201)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = distance, y = frequency, color = ref),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("with A/T at -3/+3", "without A/T at -3/+3")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest CPD-seq Dyad with(out) A/T at -3/+3 (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Frequency",
    limits = c(0, 0.06),
    breaks = seq(0, 0.06, 0.01),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/brogaard_viterbi_AT3_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/brogaard_viterbi_AT3_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Reversed Brogaard Redundant Offset Plots ----
D<-"NucleosomePattern/composite_models/composite_ice/reverse_AT3_distance_evals"
data <- read.table(paste0(D, "/redund_AT3.tsv"))
colnames(data) <- c("distance", "frequency")
data$ref <- "AT3"

data2 <- read.table(paste0(D, "/redund_notAT3.tsv"))
colnames(data2) <- c("distance", "frequency")
data2$ref <- "no_AT3"

plotData <- data.frame(
  distance = rep(-100:100, 2),
  frequency = c(
    data[abs(data$distance) <= 100,]$frequency / N_AorT,
    data2[abs(data2$distance) <= 100,]$frequency / (N - N_AorT)
  ),
  ref = rep(c("AT3", "no_AT3"), each = 201)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = distance, y = frequency, color = ref),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("with A/T at -3/+3", "without A/T at -3/+3")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest CPD-seq Dyad with(out) A/T at -3/+3 (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Frequency",
    limits = c(0, 0.2),
    breaks = seq(0, 0.2, 0.05),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/redundant_viterbi_AT3_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/redundant_viterbi_AT3_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Reversed Weiner Redundant Offset Plots ----
D<-"NucleosomePattern/composite_models/composite_ice/reverse_AT3_distance_evals"
data <- read.table(paste0(D, "/weiner_AT3.tsv"))
colnames(data) <- c("distance", "frequency")
data$ref <- "AT3"

data2 <- read.table(paste0(D, "/weiner_notAT3.tsv"))
colnames(data2) <- c("distance", "frequency")
data2$ref <- "no_AT3"

plotData <- data.frame(
  distance = rep(-100:100, 2),
  frequency = c(
    data[abs(data$distance) <= 100,]$frequency / N_AorT,
    data2[abs(data2$distance) <= 100,]$frequency / (N - N_AorT)
  ),
  ref = rep(c("AT3", "no_AT3"), each = 201)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = distance, y = frequency, color = ref),
    linewidth = 1.25
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("with A/T at -3/+3", "without A/T at -3/+3")
  ) +
  scale_x_continuous(
    name = "Offset from Nearest CPD-seq Dyad with(out) A/T at -3/+3 (bp)",
    limits = c(-100, 100),
    breaks = -5:5*20,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Frequency",
    limits = c(0, 0.012),
    breaks = seq(0, 0.012, 0.002),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/weiner_viterbi_AT3_offset.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/weiner_viterbi_AT3_offset.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(data, plotData, P, G, A_table, AorT_table, chrom, chroms, D, E, labels,
   A.3, T3, AorT, annotationData, brogaard, P1, P3, viterbi, weiner, load_seqs,
   data_AT3, data_noAT3, data2, redund, all, NoAnorT, NoAnorT_table, N, N_AorT,
   T_table, theseDyads, mouse_stats, pombe_redund, pombe_unique,redund_noUnique,
   E_mouse, E_pombe, AT_mouse, N_mouse, string)

##### Average Nucleosome Score with Huge Window #####

# Load Scores
chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "genome_scores/chromosomes.txt"
  ),
  header = FALSE
)$V1
scores <- list()
for (chrom in chroms) {
  scores[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/genome_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[chrom]]) <- c("position", "score")
}

chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/viterbi_dyad_calls",
    "/chromosomes.txt"
  ),
  header = FALSE
)$V1
nucleosomes <- list()
for (chrom in chroms) {
  nucleosomes[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/viterbi_dyad_calls/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(nucleosomes[[chrom]]) <- c("dyad")
}

W <- c(-500,500) # The window width to graph around the true dyad

total <- rep(0, W[2]-W[1] + 1)
n <- 0
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  for (i in 1:nrow(nucleosomes[[chrom]])) {
    dyad <- nucleosomes[[chrom]]$dyad[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      total <- total + theseScores
      n <- n + 1
    }
  }
};plotData <- data.frame(
  relative_position = W[1]:W[2],
  average_score = total / n
);write.table(
  plotData,
  file = "PaperFigures/FigureData/score_dyad_position_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

annotationData <- data.frame(
  x = rep(10.1*-24:24, each = 2),
  y = rep(c(-Inf, Inf), 49),
  group = rep(1:49, each = 2)
)

P <- ggplot(
  data = plotData[
    (plotData$relative_position <= 250) & (plotData$relative_position >= -250),
  ],
  mapping = aes(x = relative_position, y = average_score)
) +
  geom_line(
    data = annotationData,
    mapping = aes(x=x, y=y, group=group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(linewidth = 0.75) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-250, 250),
    breaks = -5:5*50,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2.5,2.5),
    breaks = seq(-2.5,2.5,0.5),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/score_dyad_position_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

rm(chrom,scores,chroms,nucleosomes,W,total,n,plotData,P,G,theseScores,dyad,i)

##### Linker Length Distributions #####
rot <- read.table(
  "NucleosomePattern/composite_models/composite_ice/rotational_scores.tsv",
  col.names = c("chrom", "loc", "score")
)
Rot <- data.frame(chrom = character(0), loc = integer(0), score = numeric(0))
for (chrom in unique(rot$chrom)) {
  rotSubset <- rot[rot$chrom == chrom,]
  Rot <- rbind(Rot,rotSubset[sort(rotSubset$loc, index.return=TRUE)$ix,])
}

all_linker_lens <- c()
for (chrom in unique(Rot$chrom)) {
  subset <- Rot[Rot$chrom == chrom,]
  all_linker_lens <- c(
    all_linker_lens,
    subset$loc[2:nrow(subset)] - subset$loc[1:(nrow(subset)-1)] - 146 - 1
  )
}
AT <- table(all_linker_lens)

Prot <- rot[rev(sort(rot$score, index.return = TRUE)$ix),]

N <- nrow(rot)
P <- floor(N/10)
P9 <- Prot[1:P,]
P1 <- Prot[(N-P):N,]
Q_N <- floor(N/4)
Q3 <- Prot[1:Q_N,]
Q1 <- Prot[(N-Q_N):N,]

get_linker_lens <- function(S, R) {
  result <- c()
  for (chrom in unique(S$chrom)) {
    subset <- S[S$chrom == chrom,]
    rotSubset <- R[R$chrom == chrom,]
    fronts <- c()
    backs <- c()
    for (i in 1:nrow(subset)) {
      j <- which(rotSubset$loc == subset$loc[i])
      L <- c()
      if (j > 1) {
        f <- rotSubset$loc[j]
        b <- rotSubset$loc[j - 1]
        
        if (!(f %in% fronts) || !(TRUE %in% (b %in% backs[fronts == f]))) {
          L <- c(L, f - b)
          fronts <- c(fronts, f)
          backs <- c(backs, b)
        }
      }
      if (j < nrow(rotSubset)) {
        f <- rotSubset$loc[j + 1]
        b <- rotSubset$loc[j]
        
        if (!(f %in% fronts) || !(TRUE %in% (b %in% backs[fronts == f]))) {
          L <- c(L, f - b)
          fronts <- c(fronts, f)
          backs <- c(backs, b)
        }
      }
      result <- c(result, L - 146 - 1 )
    }
  }
  return(result)
}

P9_linker_lens <- get_linker_lens(P9, Rot)
P1_linker_lens <- get_linker_lens(P1, Rot)
Q1_linker_lens <- get_linker_lens(Q1, Rot)
Q3_linker_lens <- get_linker_lens(Q3, Rot)

P1T <- table(P1_linker_lens)
P9T <- table(P9_linker_lens)
Q1T <- table(Q1_linker_lens)
Q3T <- table(Q3_linker_lens)

# > All ----
plotData <- data.frame(
  len = 0:100,
  frequency = as.integer(AT[1:101])
)

annotationData <- data.frame(
  x = rep(10.1*0:100 - 146 - 1, 2),
  y = rep(c(-Inf, Inf), each = 101),
  group = rep(0:100, 2)
)
annotationData <- annotationData[
  (annotationData$x >= 0) & (annotationData$x <= 100),
]

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = len, y = frequency),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Linker Length (bp)",
    limits = c(0, 100),
    breaks = c(0:10*10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    breaks = c(0:3*900),
    expand = c(0, 0)
  ) +
  plotTheme +
  theme(
    legend.direction = "horizontal",
    plot.margin = margin(t = 10, r = 10, l = 10, b = 10),
    legend.position.inside = c(0.075, 0.95),
    plot.title = element_textbox(
      size = tinyTxtSize,
      fill = "white",
      margin = margin(b = 9, l = 5, r = 5),
      padding = margin(l = 5, r = 5, b = 0, t = 0),
      vjust = 0.5,
      hjust = 0.5
    )
  ) +
  geom_text(
    data = annotationData,
    mapping = aes(x = x),
    y = 2750,
    label = "🡇",
    vjust = 0
  ) +
  annotate(
    geom = "line",
    x = c(
      annotationData$x[1] - 1,
      annotationData$x[1] - 1,
      annotationData$x[nrow(annotationData)] + 1,
      annotationData$x[nrow(annotationData)] + 1
    ),
    y = c(2800, 2900, 2900, 2800),
    linewidth = 1.25,
    lineend = "square"
  ) +
  coord_cartesian(clip = "off", ylim = c(0, 2700)) +
  ggtitle("Minor Groove In-Phase Positions")

P

png(
  "PaperFigures/linker_len_distribution_all.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/linker_len_distribution_all.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Deciles ----
plotData <- data.frame(
  len = 0:100,
  decile = rep(c(1,9), each = 101),
  frequency = c(as.integer(P1T[1:101]), as.integer(P9T[1:101]))
)

annotationData <- data.frame(
  x = rep(10.1*0:100 - 146 - 1, 2),
  y = rep(c(-Inf, Inf), each = 101),
  group = rep(0:100, 2)
)
annotationData <- annotationData[
  (annotationData$x >= 0) & (annotationData$x <= 100),
]

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = len, y = frequency, color = factor(decile)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Linker Length (bp)",
    limits = c(0, 100),
    breaks = c(0:10*10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    breaks = c(0:3*250),
    expand = c(0, 0)
  ) +
  scale_colour_discrete(
    name = NULL,
    type = c("#CC79A7", "#009E73"),
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  plotTheme +
  theme(
    legend.direction = "horizontal",
    plot.margin = margin(t = 10, r = 10, l = 10, b = 10),
    legend.position.inside = c(0.075, 0.95),
    plot.title = element_textbox(
      size = tinyTxtSize,
      fill = "white",
      margin = margin(b = 9, l = 5, r = 5),
      padding = margin(l = 5, r = 5, b = 0, t = 0),
      vjust = 0.5,
      hjust = 0.5
    )
  ) +
  geom_text(
    data = annotationData,
    mapping = aes(x = x),
    y = 760,
    label = "🡇",
    vjust = 0
  ) +
  annotate(
    geom = "line",
    x = c(
      annotationData$x[1] - 1,
      annotationData$x[1] - 1,
      annotationData$x[nrow(annotationData)] + 1,
      annotationData$x[nrow(annotationData)] + 1
    ),
    y = c(775, 800, 800, 775),
    linewidth = 1.25,
    lineend = "square"
  ) +
  coord_cartesian(clip = "off", ylim = c(0, 750)) +
  ggtitle("Minor Groove In-Phase Positions")

P

png(
  "PaperFigures/linker_len_distribution.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/linker_len_distribution.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Quartiles ----
plotData <- data.frame(
  len = 0:100,
  quartile = rep(c(1,3), each = 101),
  frequency = c(as.integer(Q1T[1:101]), as.integer(Q3T[1:101]))
)

annotationData <- data.frame(
  x = rep(10.1*0:100 - 146 - 1, 2),
  y = rep(c(-Inf, Inf), each = 101),
  group = rep(0:100, 2)
)
annotationData <- annotationData[
  (annotationData$x >= 0) & (annotationData$x <= 100),
]

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = len, y = frequency, color = factor(quartile)),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Linker Length (bp)",
    limits = c(0, 100),
    breaks = c(0:10*10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    breaks = c(0:4*400),
    expand = c(0, 0)
  ) +
  scale_colour_discrete(
    name = NULL,
    type = c("#CC79A7", "#009E73"),
    labels = c("Low Rotational Score", "High Rotational Score")
  ) +
  plotTheme +
  theme(
    legend.direction = "horizontal",
    plot.margin = margin(t = 10, r = 10, l = 10, b = 10),
    legend.position.inside = c(0.075, 0.95),
    plot.title = element_textbox(
      size = tinyTxtSize,
      fill = "white",
      margin = margin(b = 12, l = 5, r = 5),
      padding = margin(l = 5, r = 5, b = 0, t = 0),
      vjust = 0.5,
      hjust = 0.5
    )
  ) +
  geom_text(
    data = annotationData,
    mapping = aes(x = x),
    y = 1650,
    label = "🡇",
    vjust = 0
  ) +
  annotate(
    geom = "line",
    x = c(
      annotationData$x[1] - 1,
      annotationData$x[1] - 1,
      annotationData$x[nrow(annotationData)] + 1,
      annotationData$x[nrow(annotationData)] + 1
    ),
    y = c(1700, 1750, 1750, 1700),
    linewidth = 1.25,
    lineend = "square"
  ) +
  coord_cartesian(clip = "off", ylim = c(0, 1600)) +
  ggtitle("Minor Groove In-Phase Positions")

P

png(
  "PaperFigures/linker_len_distribution_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/linker_len_distribution_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

rm(P,Rot,rot,rotSubset,subset,i,j,N,L,chrom,plotData,annotationData,P9T,P1T,P1,
   P9, P1_linker_lens, P9_linker_lens, Prot)

##### Minor-Out vs. Minor-In Positions Bar Graphs with Redundant Set #####
data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/periodic_proportions.txt"
  ),
  sep = "\t",
  header = FALSE
)
E_out <- as.numeric(unlist(strsplit(data$V1[1], "="))[2])

data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/redundant_periodic.tsv"
  ),
  sep = "\t",
  header = FALSE
)
colnames(data) <- c("comparison", "successes", "total")

outData <- data.frame(
  comparison = data$comparison,
  diff = data$successes / data$total - E_out
)
outData <- cbind(outData, position = rep("minor_out", nrow(outData)))

data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/anti_periodic_proportions.txt"
  ),
  sep = "\t",
  header = FALSE
)
E_in <- as.numeric(unlist(strsplit(data$V1[1], "="))[2])

data <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/distance_evals",
    "/redundant_anti_periodic.tsv"
  ),
  sep = "\t",
  header = FALSE
)
colnames(data) <- c("comparison", "successes", "total")

inData <- data.frame(
  comparison = data$comparison,
  diff = data$successes / data$total - E_in
)
inData <- cbind(inData, position = rep("minor_in", nrow(inData)))

plotData <- rbind(outData, inData)
plotData <- plotData[
  !(
    plotData$comparison %in% c(
      "greedyVsRedund",
      "redundVsBrogaard",
      "brogaardVsRedund"
    )
  ),
]

# > Altogether ----

P <- ggplot(data=plotData) +
  geom_col(
    mapping = aes(x = comparison, y = diff, group = position, fill = position),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square"
  ) +
  scale_x_discrete(name=NULL, labels=c("Redundant","Viterbi","MNase-seq")) +
  scale_y_continuous(
    name = "Deviation from Expectation",
    limits = c(-0.2, 0.2),
    breaks = seq(-0.2, 0.2, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(-2:2*10, toString)), "%")
  ) +
  scale_fill_discrete(
    name = NULL,
    type = c("black", "white"),
    labels = c("Minor-In", "Minor-Out")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  ggtitle(bquote(underline(bold("Redundant")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position.inside = c(0.7,1),
    plot.margin = margin(t = 10, l = 10, b = 10, r = 20)
  )

P

png(
  "PaperFigures/viterbi_redundant_rotational_barplot.png",
  width = 1200,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/viterbi_redundant_rotational_barplot.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Minor-Out ----

P <- ggplot(data=plotData[plotData$position == "minor_out",]) +
  geom_col(
    mapping = aes(x = comparison, y = diff),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square",
    fill = "white"
  ) +
  scale_x_discrete(
    name = NULL,
    labels = c("Chemical\nCleavage", "CPD-seq", "MNase-seq")
  ) +
  scale_y_continuous(
    name = "Proportion of Dyad Positions",
    limits = c(0.2 - E_out, 0.6 - E_out),
    breaks = seq(0.2 - E_out, 0.6 - E_out, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(2:6*10, toString)), "%")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(t = 10, b = 10, l = 10, r = 60),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      margin = margin(t = -60, b = 60, r = 0, l = 0)
    )
  ) +
  annotate(
    geom = "text",
    label = "🡄",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.4,
    size = 5
  ) +
  annotate(
    geom = "text",
    label = "    Random",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.5
  ) +
  coord_cartesian(clip = "off") +
  ggtitle(bquote(underline(bold("Minor-Out"))))

P

png(
  "PaperFigures/viterbi_redundant_rotational_barplot_minor_out.png",
  width = 1500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Minor-In ----

P <- ggplot(data=plotData[plotData$position == "minor_in",]) +
  geom_col(
    mapping = aes(x = comparison, y = diff),
    position = position_dodge(),
    color = "black",
    linewidth = 1,
    lineend = "square",
    fill = "white"
  ) +
  scale_x_discrete(
    name = NULL,
    labels = c("Chemical\nCleavage","CPD-seq","MNase-seq")
  ) +
  scale_y_continuous(
    name = "Proportion of Dyad Positions",
    limits = c(0.2 - E_in, 0.6 - E_in),
    breaks = seq(0.2 - E_in, 0.6 - E_in, 0.1),
    expand = c(0,0),
    labels = paste0(unlist(lapply(2:6*10, toString)), "%")
  ) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0, 0),
    linewidth = 1.25,
    linetype = "solid"
  ) +
  ggtitle(bquote(underline(bold("Redundant")))) +
  plotTheme +
  theme(
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(t = 10, b = 10, l = 10, r = 60),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  annotate(
    geom = "text",
    label = "🡄",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.4,
    size = 5
  ) +
  annotate(
    geom = "text",
    label = "    Random",
    x = Inf,
    y = 0,
    hjust = "left",
    fontface = "bold",
    vjust = 0.5
  ) +
  coord_cartesian(clip = "off") +
  ggtitle(bquote(underline(bold("Minor-In"))))

P

png(
  "PaperFigures/viterbi_redundant_rotational_barplot_minor_in.png",
  width = 1500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

rm(P, plotData, inData, outData, E_out, E_in, data)

##### Explanation of Rotationally In-Phase Nucleosomes #####

X <- seq(0, 10, 0.01)
Y <- cos(2*pi*X)
nearPeak <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearTrough <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearTrough] <- "trough"
colours[nearPeak] <- "peak"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X <= 3.25) | (5.75 <= X)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

boundaries <- plotData[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

plotData$c[plotData$x <= 3.25] <- paste0("red", plotData$c[plotData$x <= 3.25])
plotData$c[plotData$x >= 5.75] <- paste0("blue", plotData$c[plotData$x >= 5.75])

ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_ellipse(
    data = data.frame(
      x = c(-2,11),
      y = c(3.1, 3.1),
      c = c("redinbetween", "blueinbetween")
    ),
    mapping = aes(
      x0 = x,
      y0 = y,
      a = 5.25,
      b = 1,
      angle = 0,
      fill = c,
      color = c
    )
  ) +
  scale_x_continuous(name = NULL, expand = c(0,0)) +
  coord_cartesian(clip = "on", xlim = c(0,10)) +
  guides(linetype = guide_none(), color = guide_none(), fill = guide_none()) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  ) +
  scale_color_discrete(
    type = c(
      "blue",
      "#BBBBFF",
      "#000033",
      "purple",
      "#FFAAFF",
      "red",
      "#FFAAAA",
      "#440000", 
      "#440044"
    )
  ) +
  scale_fill_discrete(type = c("blue", "red"))

ggsave(
  "PaperFigures/RotationallyInPhaseExplanation.png",
  width = 2000,
  height = 500,
  units = "px"
)

rm(X, Y, colours, linetype, groups, plotData, boundaries, nearPeak)

##### Explanation of Rotationally Out-of-Phase Nucleosomes #####

X <- seq(0, 6.25, 0.01)
Y <- cos(2*pi*X)
nearPeak <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearTrough <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearTrough] <- "trough"
colours[nearPeak] <- "peak"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X <= 3.25)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

boundaries <- plotData[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

plotData$c <- paste0("red", plotData$c)

X <- seq(3.25, 10, 0.01)
Y <- cos(2*pi*X + pi)
nearTrough <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearPeak <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearPeak] <- "peak"
colours[nearTrough] <- "trough"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X >= 6.25)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData2 <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData2[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x + pi)
plotData2 <- rbind(
  plotData2,
  boundaries
)

boundaries <- plotData2[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x + pi)
plotData2 <- rbind(
  plotData2,
  boundaries
)

plotData2$c <- paste0("blue", plotData2$c)

ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_line(
    data = plotData2,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_ellipse(
    data = data.frame(
      x = c(-2,11.5),
      y = c(3.1, 3.1),
      c = c("redinbetween", "blueinbetween")
    ),
    mapping = aes(
      x0 = x,
      y0 = y,
      a = 5.25,
      b = 1,
      angle = 0,
      fill = c,
      color = c
    )
  ) +
  scale_x_continuous(name = NULL, expand = c(0,0)) +
  coord_cartesian(clip = "on", xlim = c(0,10)) +
  guides(linetype = guide_none(), color = guide_none(), fill = guide_none()) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  ) +
  scale_color_discrete(
    type = c(
      "blue",
      "#BBBBFF",
      "#000033",
      "red",
      "#FFAAAA",
      "#440000"
    )
  ) +
  scale_fill_discrete(type = c("blue", "red"))

ggsave(
  "PaperFigures/RotationallyOutOfPhaseExplanation.png",
  width = 2000,
  height = 500,
  units = "px"
)

rm(X, Y, colours, linetype, groups, plotData, boundaries, nearPeak, plotData2)

##### PDB Structure Pie Chart #####
prop <- 44/181
data <- data.frame(
  xmin = c(0, prop),
  xmax = c(prop, 1),
  c = c("red", "blue")
)

P <- ggplot(data = data) +
  geom_rect(mapping = aes(xmin=xmin, xmax=xmax, ymin=0, ymax=1, fill=c)) +
  coord_polar() +
  scale_fill_discrete(
    name = NULL,
    type = c("blue", "red"),
    labels = c("No A/T at -3/+3", "A/T at -3/+3")
  ) +
  plotTheme +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.box.margin = margin(0,0,0,0),
    legend.margin = margin(0,0,0,0),
    legend.spacing = unit(0,"pt"),
    legend.box.spacing = unit(-40,"pt"),
    plot.margin=margin(t = 10, b = 10, r = 10, l = 30)
  ) +
  annotate(
    geom = "text",
    label = paste0(toString(round(prop*100,2)), "%"),
    x = prop/2,
    y = 0.5,
    fontface = "bold",
    size = 20,
    size.unit = "pt",
    hjust = 0.35,
    color = "white"
  ) +
  annotate(
    geom = "text",
    label = paste0(toString(100-round(prop*100,2)), "%"),
    x = (1-prop)/2 + prop,
    y = 0.5,
    fontface = "bold",
    size = 20,
    size.unit = "pt",
    hjust = 0.25,
    color = "white"
  ) +
  ggtitle(
    bquote(bold(underline(bolditalic("in vitro")~"Nucleosome Structures")))
  )

P

png(
  "PaperFigures/PDB_structures_pie.png",
  width = 1200,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

rm(P, data, prop)

##### Transcription Aligned with +1 Nucleosome Plots #####
chroms = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII",
           "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV",
           "chrXV", "chrXVI")
directory <- "NucleosomePattern/composite_models/composite_ice"

D <- data.frame()
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  data <- read.table(
    paste0("PreprocessedData/TSS/", chrom, ".tsv"),
    sep = "\t"
  )[c(1,2,5)]
  colnames(data) <- c("TSS", "gene", "strand")
  
  data$plus_one <- read.table(
    paste0(directory, "/plusOnes/", chrom, ".tsv"),
    sep = "\t"
  )$V1
  
  data$rotational <- read.table(
    paste0(directory, "/plusOnes_rotationalScores/", chrom, ".tsv"),
    sep = "\t"
  )$V1
  
  scoreData <- read.table(
    paste0(directory, "/plusOnes_alignedScores/", chrom, ".tsv"),
    sep = "\t"
  )
  colnames(scoreData) <- paste0("score", unlist(lapply(-500:650, toString)))
  
  nucData <- read.table(
    paste0(directory, "/plusOnes_alignedNucs/", chrom, ".tsv"),
    sep = "\t"
  )
  colnames(nucData) <- paste0("nuc", unlist(lapply(-500:650, toString)))
  
  dinucData <- read.table(
    paste0(directory, "/plusOnes_alignedDinucs/", chrom, ".tsv"),
    sep = "\t"
  )
  colnames(dinucData) <- paste0("dinuc", unlist(lapply(-499.5:649.5,toString)))
  
  thisD <- cbind(
    chrom = rep(chrom, nrow(data)),
    data,
    scoreData,
    nucData,
    dinucData
  )
  
  if (ncol(D) < 1) {
    D <- thisD
  } else {
    D <- rbind(D, thisD)
  }
}

valid <- D$plus_one >= 0
D <- D[valid,]

G <- read.table("../Data/GeneExpression_holstege98.tsv", sep="\t", header=TRUE)
geneTable <- data.frame(
  ExpressionLevel = numeric(0),
  HalfLife = numeric(0),
  TranscriptionalFreq = numeric(0)
)
for (i in 1:nrow(D)) {
  j <- which(G$ORF == D$gene[i])
  if (length(j) == 1) {
    geneTable <- rbind(geneTable, G[j, 2:ncol(G)])
  } else {
    geneTable <- rbind(
      geneTable,
      data.frame(
        ExpressionLevel = NA,
        HalfLife = NA,
        TranscriptionalFreq = NA
      )
    )
  }
}
D <- cbind(D, geneTable)

# > Overall ----

S <- colMeans(D[paste0("score", unlist(lapply(0:650, toString)))])

plotData <- data.frame(
  x = -500:650,
  y = unlist(D[paste0("score", unlist(lapply(-500:650, toString)))][which.max(D$rotational),])
)
ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = x, y = y)
  )

Ps <- c()
Pwrs <- c()
W <- 73
for (i in (W+1):(nrow(plotData)-W)) {
  print(paste0("> ", toString(plotData$x[i]), "..."))
  L <- lsp(
    plotData$y[(i-W):(i+W)],
    times = (-W):W,
    type = "period",
    from = 8,
    to = 12,
    ofac = 100,
    plot = FALSE
  )
  Ps <- c(Ps,L$peak.at[1])
  Pwrs <- c(Pwrs, L$peak[1])
}
plotData$p <- c(rep(NA, W), Ps, rep(NA, W))
plotData$pwr <- c(rep(NA, W), Pwrs, rep(NA, W))

ggplot() +
  geom_line(
    data = data.frame(
      x = c(-Inf, Inf, -73, -73, 73, 73),
      y = c(10.1, 10.1, -Inf, Inf, -Inf, Inf),
      group = c(1, 1, 2, 2, 3, 3)
    ),
    mapping = aes(x = x, y = y, group = group),
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = x, y = p)
  )

L <- lsp(S, times=0:650, type="period", to=15, ofac=100)

Dns <- D[paste0("dinuc",unlist(lapply(0.5:649.5,toString)))]
Dns <- (Dns == "AA") | (Dns == "TT") | (Dns == "AT") | (Dns == "TA")
Dns <- colSums(Dns)

L_D <- lsp(Dns, times=0.5:649.5, type="period", to=15, ofac=100)

N <- D[paste0("nuc", unlist(lapply(0:650, toString)))]
plus2s <- apply(N[2:ncol(N)], 1, which.max)
linker1_2 <- plus2s - 147
linker1_2 <- linker1_2[(linker1_2 >= 0) & (linker1_2 <= 100)]

Q1 <- D[D$rotational < quantile(D$rotational, 0.25),]
S_Q1 <- colMeans(Q1[paste0("score", unlist(lapply(0:650, toString)))])

L_Q1 <- lsp(S_Q1, times=0:650, type="period", to=15, ofac=100)

Dns_Q1 <- Q1[paste0("dinuc",unlist(lapply(0.5:649.5,toString)))]
Dns_Q1 <- (Dns_Q1=="AA") | (Dns_Q1=="TT") | (Dns_Q1=="AT") | (Dns_Q1=="TA")
Dns_Q1 <- colSums(Dns_Q1)

L_D_Q1 <- lsp(Dns_Q1, times=0.5:649.5, type="period", to=15, ofac=100)

Q3 <- D[D$rotational > quantile(D$rotational, 0.75),]
S_Q3 <- colMeans(Q3[paste0("score", unlist(lapply(0:650, toString)))])

L_Q3 <- lsp(S_Q3, times=0:650, type="period", to=15, ofac=100)

Dns_Q3 <- Q3[paste0("dinuc",unlist(lapply(0.5:649.5,toString)))]
Dns_Q3 <- (Dns_Q3=="AA") | (Dns_Q3=="TT") | (Dns_Q3=="AT") | (Dns_Q3=="TA")
Dns_Q3 <- colSums(Dns_Q3)

L_D_Q3 <- lsp(Dns_Q3, times=0.5:649.5, type="period", to=15, ofac=100)

# > Stratify by Gene Expression Level ----

D2 <- D[!is.na(D$ExpressionLevel),]
summary(D2$ExpressionLevel)

low <- D2[D2$ExpressionLevel < quantile(D2$ExpressionLevel, 0.25),]
high <- D2[D2$ExpressionLevel > quantile(D2$ExpressionLevel, 0.75),]

S_low <- colMeans(low[paste0("score", unlist(lapply(0:650, toString)))])
S_high <- colMeans(high[paste0("score", unlist(lapply(0:650, toString)))])

L_low <- lsp(S_low, times=0:650, type="period", to=15, ofac=100)
L_high <- lsp(S_high, times=0:650, type="period", to=15, ofac=100)

Dns_low <- low[paste0("dinuc",unlist(lapply(0.5:649.5,toString)))]
Dns_low <- (Dns_low=="AA")|(Dns_low=="TT")|(Dns_low=="AT")|(Dns_low=="TA")
Dns_low <- colSums(Dns_low)

L_D_low <- lsp(Dns_low, times=0.5:649.5, type="period", to=15, ofac=100)

Dns_high <- high[paste0("dinuc",unlist(lapply(0.5:649.5,toString)))]
Dns_high <- (Dns_high=="AA")|(Dns_high=="TT")|(Dns_high=="AT")|(Dns_high=="TA")
Dns_high <- colSums(Dns_high)

L_D_high <- lsp(Dns_high, times=0.5:649.5, type="period", to=15, ofac=100)

N_low <- low[paste0("nuc", unlist(lapply(0:650, toString)))]
plus2s_low <- apply(N_low[2:ncol(N)], 1, which.max)
linker1_2_low <- plus2s_low - 147
linker1_2_low <- linker1_2_low[(linker1_2_low >= 0) & (linker1_2_low <= 100)]

N_high <- high[paste0("nuc", unlist(lapply(0:650, toString)))]
plus2s_high <- apply(N_high[2:ncol(N)], 1, which.max)
linker1_2_high <- plus2s_high - 147
linker1_2_high <- linker1_2_high[(linker1_2_high >= 0)&(linker1_2_high <= 100)]

##### Dinucleotide Frequency Plots with Wide Windows #####
directory <- "NucleosomePattern/composite_models/composite_ice"

P_N <- read.table(paste0(directory, "/P1_sequence_analysis_wide/n.txt"))$V1
P1_WW = read.table(
  paste0(directory, "/P1_sequence_analysis_wide/WW.tsv"),
  sep = "\t"
)$V1
P1_SS = read.table(
  paste0(directory, "/P1_sequence_analysis_wide/SS.tsv"),
  sep = "\t"
)$V1
P3_WW = read.table(
  paste0(directory, "/P3_sequence_analysis_wide/WW.tsv"),
  sep = "\t"
)$V1
P3_SS = read.table(
  paste0(directory, "/P3_sequence_analysis_wide/SS.tsv"),
  sep = "\t"
)$V1
Q_N <- read.table(paste0(directory, "/Q1_sequence_analysis_wide/n.txt"))$V1
Q1_WW = read.table(
  paste0(directory, "/Q1_sequence_analysis_wide/WW.tsv"),
  sep = "\t"
)$V1
Q1_SS = read.table(
  paste0(directory, "/Q1_sequence_analysis_wide/SS.tsv"),
  sep = "\t"
)$V1
Q3_WW = read.table(
  paste0(directory, "/Q3_sequence_analysis_wide/WW.tsv"),
  sep = "\t"
)$V1
Q3_SS = read.table(
  paste0(directory, "/Q3_sequence_analysis_wide/SS.tsv"),
  sep = "\t"
)$V1

annotationData <- data.frame(
  x = rep(10.1*(-14:14), 2),
  y = rep(c(-Inf, Inf), each = 29),
  group = rep(-14:14, 2)
)

ins <- round(5.05*-99:99)
ins <- c(ins, ins - 1)
outs <- round(10.1*-49:49)
outs <- c(outs, outs - 1)

X <- -499.5:499.5
colors <- c()
for (x in X) {
  if (min(abs(x - outs)) < 1) {
    colors <- c(colors, 1)
  } else if (min(abs(x - ins)) < 1) {
    colors <- c(colors, -1)
  } else {
    colors <- c(colors, 0)
  }
}

# > WW Frequency (Weak) ----

plotData <- data.frame(
  frequency = P1_WW,
  position = seq(-499.5, 499.5, 1),
  color = colors
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -150 & plotData$position <= 150,],
    mapping = aes(x = position, y = frequency, color = color),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-150,150,50),
    limits = c(-150, 150),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(1700, 2900),
    breaks = c(seq(1700, 2900, 200)),
    expand = c(0,0)
  ) +
  guides(color = guide_none()) +
  scale_color_continuous(
    palette = c("#1a8cff", "black", "red")
  )

P

png(
  "PaperFigures/WW_frequency_weak_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_weak_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW Frequency (Strong) ----

plotData <- data.frame(
  frequency = P3_WW,
  position = seq(-499.5, 499.5, 1),
  color = colors
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -150 & plotData$position <= 150,],
    mapping = aes(x = position, y = frequency, color = color),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-150,150,50),
    limits = c(-150, 150),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(1700, 2900),
    breaks = c(seq(1700, 2900, 200)),
    expand = c(0,0)
  ) +
  guides(color = guide_none()) +
  scale_color_continuous(
    palette = c("#1a8cff", "black", "red")
  )

P

png(
  "PaperFigures/WW_frequency_strong_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_strong_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Weak) ----

plotData <- data.frame(
  frequency = P1_SS,
  position = seq(-499.5, 499.5, 1),
  color = colors
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -150 & plotData$position <= 150,],
    mapping = aes(x = position, y = frequency, color = color),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-150,150,50),
    limits = c(-150, 150),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(700, 1200),
    breaks = c(seq(700, 1200, 100)),
    expand = c(0,0)
  ) +
  guides(color = guide_none()) +
  scale_color_continuous(
    palette = c("#1a8cff", "black", "red")
  )

P

png(
  "PaperFigures/SS_frequency_weak_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_weak_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Strong) ----

plotData <- data.frame(
  frequency = P3_SS,
  position = seq(-499.5, 499.5, 1),
  color = colors
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -150 & plotData$position <= 150,],
    mapping = aes(x = position, y = frequency, color = color),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-150,150,50),
    limits = c(-150, 150),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(700, 1200),
    breaks = c(seq(700, 1200, 100)),
    expand = c(0,0)
  ) +
  guides(color = guide_none()) +
  scale_color_continuous(
    palette = c("#1a8cff", "black", "red")
  )

P

png(
  "PaperFigures/SS_frequency_strong_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_strong_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW Frequency (Weak Quartile) ----

plotData <- data.frame(
  frequency = Q1_WW,
  position = seq(-499.5, 499.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -250 & plotData$position <= 250,],
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "blue"
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-250,250,50),
    limits = c(-250, 250),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(4500, 7000),
    breaks = c(seq(4500, 7000, 500)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/WW_frequency_weak_wide_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_weak_wide_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > WW Frequency (Strong Quartile) ----

plotData <- data.frame(
  frequency = Q3_WW,
  position = seq(-499.5, 499.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -250 & plotData$position <= 250,],
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "blue"
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-250,250,50),
    limits = c(-250, 250),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "AA/AT/TA/TT Frequency",
    limits = c(4700, 7200),
    breaks = c(seq(4700, 7200, 500)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/WW_frequency_strong_wide_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_frequency_strong_wide_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Weak Quartile) ----

plotData <- data.frame(
  frequency = Q1_SS,
  position = seq(-499.5, 499.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -250 & plotData$position <= 250,],
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "red"
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-250,250,50),
    limits = c(-250, 250),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(1800, 2800),
    breaks = c(seq(1800, 2800, 200)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/SS_frequency_weak_wide_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_weak_wide_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > SS Frequency (Strong Quartile) ----

plotData <- data.frame(
  frequency = Q3_SS,
  position = seq(-499.5, 499.5, 1)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[plotData$position >= -250 & plotData$position <= 250,],
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "red"
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-250,250,50),
    limits = c(-250, 250),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "GG/GC/CG/CC Frequency",
    limits = c(1800, 2800),
    breaks = c(seq(1800, 2800, 200)),
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/SS_frequency_strong_wide_quartile.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/SS_frequency_strong_wide_quartile.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# > Periodograms ----
selected <- (-500:499 + 0.5 >= -150) & (-500:499 + 0.5 <= 150)
L_P1_WW <- lsp(P1_WW[selected],times=-150:149+0.5,type="period",to=15,ofac=100)
L_P3_WW <- lsp(P3_WW[selected],times=-150:149+0.5,type="period",to=15,ofac=100)
L_P1_SS <- lsp(P1_SS[selected],times=-150:149+0.5,type="period",to=15,ofac=100)
L_P3_SS <- lsp(P3_SS[selected],times=-150:149+0.5,type="period",to=15,ofac=100)

# >> Low Rotational Score Periodogram ----
plotData <- data.frame(
  period = c(L_P1_WW$scanned, L_P1_SS$scanned),
  normalized_power = c(L_P1_WW$power, L_P1_SS$power),
  group=c(rep("WW",length(L_P1_WW$scanned)),rep("SS",length(L_P1_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_P1_WW$peak.at[1], L_P1_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("WW", "SS"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_SS_frequency_periodogram_weak_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_periodogram_weak_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# >> High Rotational Score Periodogram ----
plotData <- data.frame(
  period = c(L_P3_WW$scanned, L_P3_SS$scanned),
  normalized_power = c(L_P3_WW$power, L_P3_SS$power),
  group=c(rep("WW",length(L_P3_WW$scanned)),rep("SS",length(L_P3_SS$scanned)))
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(c(L_P3_WW$peak.at[1], L_P3_SS$peak.at[1]), each = 2),
  y = rep(c(-Inf, Inf), 2),
  group = rep(c("WW", "SS"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = group, linetype = group),
    linewidth = 0.75
  ) +
  scale_linetype_manual(
    values = c("dotted", "dashed")
  ) +
  guides(linetype = guide_none()) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power, color = group),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("red", "blue"),
    labels = c("GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[3],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top",
    color = "red"
  ) +
  annotate(
    "text",
    x = 10.2,
    y = 0.7,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "bottom",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/WW_SS_frequency_periodogram_strong_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- paste0(
  "#WW: ",
  toString(annotationData$x[1]),
  "\n#SS: ",
  toString(annotationData$x[3])
)

write(
  string,
  "PaperFigures/FigureData/WW_SS_frequency_periodogram_strong_wide.tsv"
)

write.table(
  plotData,
  file = "PaperFigures/FigureData/WW_SS_frequency_periodogram_strong_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(annotationData,plotData,P,all_WW,all_SS,N,P_N,Q,Q1_SS,Q1_WW,P1_SS,P1_WW,G,
   Q3_WW,Q3_SS,P3_WW,P3_SS,directory,string,L_P1_SS,L_P1_WW,L_P3_WW,L_P3_SS)

##### Average Nucleosome Score with Huge Window - By Rotational Score #####

rot <- read.table(
  "NucleosomePattern/composite_models/composite_ice/rotational_scores.tsv",
  col.names = c("chrom", "loc", "score")
)
rot <- rot[rev(sort(rot$score, index.return = TRUE)$ix),]

N <- nrow(rot)
P <- floor(N/10)
P9 <- rot[1:P,]
P1 <- rot[(N-P):N,]

# Load Scores
chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "genome_scores/chromosomes.txt"
  ),
  header = FALSE
)$V1
scores <- list()
for (chrom in chroms) {
  scores[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/genome_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[chrom]]) <- c("position", "score")
}

W <- c(-500,500) # The window width to graph around the true dyad

# > Weak ----
total <- rep(0, W[2]-W[1] + 1)
n <- 0
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  nucleosomes <- P1[P1$chrom == chrom,]
  for (i in 1:nrow(nucleosomes)) {
    dyad <- nucleosomes$loc[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      total <- total + theseScores
      n <- n + 1
    }
  }
};plotData <- data.frame(
  relative_position = W[1]:W[2],
  average_score = total / n
);write.table(
  plotData,
  file = "PaperFigures/FigureData/score_dyad_position_wide_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

annotationData <- data.frame(
  x = rep(10.1*-14:14, each = 2),
  y = rep(c(-Inf, Inf), 29),
  group = rep(1:29, each = 2)
)

colors <- integer(1001)
ins <- round(5.05*-99:99)
ins <- c(ins, ins - 1)
colors[-500:500 %in% ins] <- -1
outs <- round(10.1*-49:49)
outs <- c(outs, outs - 1)
colors[-500:500 %in% outs] <- 1
plotData$color <- colors

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x=x, y=y, group=group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[
      (plotData$relative_position<=150) & (plotData$relative_position>=-150),
    ],
    mapping = aes(x=relative_position, y=average_score, color=color),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-150, 150),
    breaks = -3:3*50,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-2,5),
    breaks = seq(-2,2,1),
    expand = c(0,0)
  ) +
  plotTheme +
  guides(color = guide_none()) +
  scale_color_continuous(
    palette = c("#1a8cff", "black", "red")
  ) +
  coord_cartesian(
    clip = "off",
    xlim = c(-150, 150),
    ylim = c(-2,2)
  ) +
  geom_ellipse(
    data = data.frame(
      centers = 1
    ),
    mapping = aes(
      x0 = centers,
      y0 = 2.35,
      a = 73,
      b = 0.4,
      angle = 0
    ),
    fill = "black"
  ) +
  theme(
    plot.title = element_text(colour = "#FFFFFF")
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "title", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/score_dyad_position_wide_weak.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

# > Strong ----
total <- rep(0, W[2]-W[1] + 1)
n <- 0
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  nucleosomes <- P9[P9$chrom == chrom,]
  for (i in 1:nrow(nucleosomes)) {
    dyad <- nucleosomes$loc[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      total <- total + theseScores
      n <- n + 1
    }
  }
};plotData <- data.frame(
  relative_position = W[1]:W[2],
  average_score = total / n
);write.table(
  plotData,
  file = "PaperFigures/FigureData/score_dyad_position_wide_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

annotationData <- data.frame(
  x = rep(10.1*-14:14, each = 2),
  y = rep(c(-Inf, Inf), 29),
  group = rep(1:29, each = 2)
)

colors <- integer(1001)
ins <- round(5.05*-99:99)
ins <- c(ins, ins - 1, ins[ins < -5] - 2, ins[ins > 5] + 1)
colors[-500:500 %in% ins] <- -1
outs <- round(10.1*-49:49)
outs <- c(outs, outs - 1, outs[outs < -5] - 2, outs[outs > 5] + 1)
colors[-500:500 %in% outs] <- 1
plotData$color <- colors

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x=x, y=y, group=group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData[
      (plotData$relative_position<=150) & (plotData$relative_position>=-150),
    ],
    mapping = aes(x = relative_position, y = average_score, color = color),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-300, 300),
    breaks = -3:3*50,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Mean Score",
    limits = c(-5,10),
    breaks = seq(-5,5,1),
    expand = c(0,0)
  ) +
  plotTheme +
  guides(color = guide_none(), alpha = guide_none()) +
  scale_color_continuous(
    palette = c("#1a8cff", "black", "red")
  ) +
  coord_cartesian(
    clip = "off",
    xlim = c(-150, 150),
    ylim = c(-5,5)
  ) +
  geom_ellipse(
    data = data.frame(
      centers = 1
    ),
    mapping = aes(
      x0 = centers,
      y0 = 6,
      a = 73,
      b = 1,
      angle = 0
    ),
    fill = "black"
  ) +
  theme(
    plot.title = element_text(colour = "#FFFFFF")
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "title", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/score_dyad_position_wide_strong.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

rm(chrom,scores,chroms,nucleosomes,W,total,n,plotData,P,G,theseScores,dyad,i,
   colors, ins, outs)

##### All Nucleosome log2 Ratio Pattern Non-Aligned Combined #####
pattern <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/damagePattern.tsv",
  header = FALSE,
  sep = "\t"
)

periodicity <- read.delim(
  "NucleosomePattern_all/cl_ice_log2/periodogram_peaks.tsv",
  header = TRUE,
  sep = "\t"
)

plotData <- data.frame(
  position = -73:72 + 0.5,
  values = (
    unlist(pattern[4,])*unlist(pattern[8,]) +
      unlist(pattern[5,])*unlist(pattern[9,])
  ) / (unlist(pattern[8,]) + unlist(pattern[9,]))
)

annotationData <- data.frame(
  x = rep(10.1*(-7:7), 2),
  y = rep(c(-Inf, Inf), each = 15),
  group = rep(-7:7, 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = values),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (nt)",
    limits = c(-73,73),
    breaks = -60 + 20*(0:6),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = bquote(bold("log"[bold("2")])*bold("(Cellular/Naked DNA)")),
    limits = c(-0.5, 0),
    expand = c(0,0)
  ) +
  plotTheme

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_pattern_nonaligned_combined.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_log2_pattern_nonaligned_combined.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

L <- lsp(plotData$values, times=-73:72+0.5, type="period", to=15, ofac=100)

plotData <- data.frame(
  period = L$scanned,
  normalized_power = L$power
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = c(L$peak.at[1], L$peak.at[1]),
  y = c(-Inf, Inf)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y),
    linewidth = 0.75,
    linetype = "dashed"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = normalized_power),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  guides(linetype = guide_none()) +
  annotate(
    "text",
    x = 10.2,
    y = 0.9,
    label = paste0(round(annotationData$x[1],2), "bp"),
    fontface = "bold",
    size = smallTxtSize,
    size.unit = "pt",
    hjust = "left",
    vjust = "top"
  )

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/all_log2_periodogram_nonaligned_combined.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

string <- toString(L$peak.at[1])
write(
  string,
  "PaperFigures/FigureData/all_log2_periodogram_nonaligned_combined.tsv"
)

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_log2_periodogram_nonaligned_combined.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  append = TRUE
)

rm(annotationData, G, P, periodicity, periodogram, labs, string)

##### Redundant vs. Viterbi Dyad Positions with High/Low Rotational Score #####
data <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice",
    "/rotational_center_to_center/redundant.tsv"
  ),
  sep = "\t"
)
colnames(data) <- c("x", "P1", "P9")

annotationData <- data.frame(
  x = rep(10.1*(-14:14), 2),
  y = rep(c(-Inf, Inf), each = 29),
  group = rep(-14:14, 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = data[data$x >= -150 & data$x <= 150,],
    mapping = aes(x = x, y = P1),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-150,150,50),
    limits = c(-150, 150),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 400),
    breaks = 0:4*100,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/redundant_viterbi_frequency_weak_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    data = data[data$x >= -150 & data$x <= 150,],
    mapping = aes(x = x, y = P9),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = seq(-150,150,50),
    limits = c(-150, 150),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Frequency",
    limits = c(0, 1000),
    breaks = 0:4*250,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/redundant_viterbi_frequency_strong_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  data,
  file = "PaperFigures/FigureData/redundant_viterbi_frequency_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

##### TSS Plot by Rotational Score #####
TSS <- read.delim(
  paste0(
    "NucleosomePattern/composite_models/composite_ice",
    "/rotational_center_to_center/TSS.tsv"
  ),
  header = FALSE,
  sep = "\t"
)
colnames(TSS) <- c("x", "P1", "P9")

plotData <- data.frame(
  position = TSS$x,
  frequency = TSS$P1
)

s <- 5 # Smoothing window radius
smooth_plotData = data.frame(
  position = (R[1]+s):(R[2]-s),
  frequency = rep(0, R[2]-R[1]-(2*s-1))
)

for (i in 1:nrow(smooth_plotData)) {
  smooth_plotData$frequency[i] <- mean(
    plotData$frequency[
      (plotData$position >= smooth_plotData$position[i]-s) &
        (plotData$position <= smooth_plotData$position[i]+s)
    ]
  )
}

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = 0,
      ymax = y * 12/100
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 0.75,
    color = "#555555"
  ) +
  geom_line(
    data = smooth_plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Low Rotational Score\nDyad Frequency",
    limits = c(0, 12),
    breaks = 0:6*2,
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/12,
      name = "Reference Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "black")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  )

P

png(
  "PaperFigures/weak_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

data <- plotData
data$smooth_frequency <- rep(NA, nrow(data))
data$smooth_frequency[
  data$position %in% smooth_plotData$position
] <- smooth_plotData$frequency

write.table(
  data,
  "PaperFigures/FigureData/weak_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

plotData <- data.frame(
  position = TSS$x,
  frequency = TSS$P9
)

s <- 5 # Smoothing window radius
smooth_plotData = data.frame(
  position = (R[1]+s):(R[2]-s),
  frequency = rep(0, R[2]-R[1]-(2*s-1))
)

for (i in 1:nrow(smooth_plotData)) {
  smooth_plotData$frequency[i] <- mean(
    plotData$frequency[
      (plotData$position >= smooth_plotData$position[i]-s) &
        (plotData$position <= smooth_plotData$position[i]+s)
    ]
  )
}

P <- ggplot() +
  geom_rect(
    data = brogaard_plotData,
    mapping = aes(
      xmin = x,
      xmax = x + w,
      ymin = 0,
      ymax = y * 15/100
    ),
    fill = "#AAAAAA"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 0.75,
    color = "#555555"
  ) +
  geom_line(
    data = smooth_plotData,
    mapping = aes(x = position, y = frequency),
    linewidth = 1.25,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Distance from TSS (bp)",
    limits = R,
    breaks = seq(R[1], R[2], 100),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "High Rotational Score\nDyad Frequency",
    limits = c(0, 15),
    breaks = 0:5*3,
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~. * 100/15,
      name = "Reference Dyad Frequency"
    )
  ) +
  plotTheme +
  theme(axis.title.y.left = element_text(color = "black")) +
  annotate(
    geom = "line",
    linetype = "dashed",
    linewidth = 1.25,
    x = c(0,0),
    y = c(-Inf, Inf)
  )

P

png(
  "PaperFigures/strong_TSS.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

data <- plotData
data$smooth_frequency <- rep(NA, nrow(data))
data$smooth_frequency[
  data$position %in% smooth_plotData$position
] <- smooth_plotData$frequency

write.table(
  data,
  "PaperFigures/FigureData/strong_TSS.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

rm(data, TSS, model, P, plotData, regressionData, smooth_plotData, i)

###### Wide Window Score, A/T, and G/C Periodograms #####
score_weak <- read.table(
  file = "PaperFigures/FigureData/score_dyad_position_wide_weak.tsv",
  header = TRUE,
  sep = "\t"
)
score_strong <- read.table(
  file = "PaperFigures/FigureData/score_dyad_position_wide_strong.tsv",
  header = TRUE,
  sep = "\t"
)
WW_weak <- read.table(
  file = "PaperFigures/FigureData/WW_frequency_weak_wide.tsv",
  header = TRUE,
  sep = "\t"
)
WW_strong <- read.table(
  file = "PaperFigures/FigureData/WW_frequency_strong_wide.tsv",
  header = TRUE,
  sep = "\t"
)
SS_weak <- read.table(
  file = "PaperFigures/FigureData/SS_frequency_weak_wide.tsv",
  header = TRUE,
  sep = "\t"
)
SS_strong <- read.table(
  file = "PaperFigures/FigureData/SS_frequency_strong_wide.tsv",
  header = TRUE,
  sep = "\t"
)

score_weak <- score_weak$average_score[
  (score_weak$relative_position >= -150)&(score_weak$relative_position <= 150)
]
score_strong <- score_strong$average_score[
  (score_strong$relative_position>=-150)&(score_strong$relative_position<=150)
]
WW_weak <- WW_weak$frequency[
  (WW_weak$position >= -150) & (WW_weak$position <= 150)
]
WW_strong <- WW_strong$frequency[
  (WW_strong$position >= -150) & (WW_strong$position <= 150)
]
SS_weak <- SS_weak$frequency[
  (SS_weak$position >= -150) & (SS_weak$position <= 150)
]
SS_strong <- SS_strong$frequency[
  (SS_strong$position >= -150) & (SS_strong$position <= 150)
]

L_weak <- lsp(score_weak, times=-150:150, type="period", to=15, ofac=100)
L_strong <- lsp(score_strong, times=-150:150, type="period", to=15, ofac=100)
L_WW_weak <- lsp(WW_weak, times=-150:149+0.5, type="period", to=15, ofac=100)
L_WW_strong <- lsp(WW_strong,times=-150:149+0.5,type="period",to=15,ofac=100)
L_SS_weak <- lsp(SS_weak, times=-150:149+0.5, type="period", to=15, ofac=100)
L_SS_strong <- lsp(SS_strong,times=-150:149+0.5,type="period",to=15,ofac=100)

plotData <- data.frame(
  period = c(L_weak$scanned, L_WW_weak$scanned, L_SS_weak$scanned),
  power = c(L_weak$power, L_WW_weak$power, L_SS_weak$power),
  dataset = c(
    rep("score", length(L_weak$scanned)),
    rep("WW", length(L_WW_weak$scanned)),
    rep("SS", length(L_SS_weak$scanned))
  )
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = c(L_weak$peak.at[1], L_weak$peak.at[1]),
  y = c(-Inf, Inf)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y),
    linetype = "dashed",
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = power, color = dataset),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("black", "red", "blue"),
    labels = c("Mean Score", "GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  annotate(
    geom = "text",
    x = 10,
    y = 0.6,
    label = toString(round(annotationData$x[1], 2)),
    fontface = "bold",
    hjust = "left",
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  guides(linetype = guide_none()) +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/score_WW_SS_periodogram_weak_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/score_WW_SS_periodogram_weak_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

plotData <- data.frame(
  period = c(L_strong$scanned, L_WW_strong$scanned, L_SS_strong$scanned),
  power = c(L_strong$power, L_WW_strong$power, L_SS_strong$power),
  dataset = c(
    rep("score", length(L_strong$scanned)),
    rep("WW", length(L_WW_strong$scanned)),
    rep("SS", length(L_SS_strong$scanned))
  )
)
plotData <- plotData[plotData$period >= 2,]

annotationData <- data.frame(
  x = rep(
    c(L_strong$peak.at[1], L_WW_strong$peak.at[1], L_SS_strong$peak.at[1]),
    each = 2
  ),
  y = rep(c(-Inf, Inf), 3),
  dataset = rep(c("score", "WW", "SS"), each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, color = dataset, linetype = dataset),
    linewidth = 0.75
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = period, y = power, color = dataset),
    linewidth = 1.25
  ) +
  scale_x_continuous(
    name = "Period (bp)",
    limits = c(2,15),
    breaks = 2:15,
    labels = 2:15,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized Power",
    limits = c(0,1),
    breaks = (0:5)/5,
    expand = c(0,0)
  ) +
  plotTheme +
  scale_color_discrete(
    name = NULL,
    type = c("black", "red", "blue"),
    labels = c("Mean Score", "GG/GC/CG/CC", "AA/AT/TA/TT")
  ) +
  annotate(
    geom = "text",
    x = 10.2,
    y = 0.925,
    label = toString(round(annotationData$x[1], 2)),
    fontface = "bold",
    hjust = "left",
    vjust = "bottom",
    size = smallTxtSize,
    size.unit = "pt"
  ) +
  annotate(
    geom = "text",
    x = 10.2,
    y = 0.875,
    label = toString(round(annotationData$x[5], 2)),
    fontface = "bold",
    hjust = "left",
    size = smallTxtSize,
    size.unit = "pt",
    color = "red"
  ) +
  annotate(
    geom = "text",
    x = 10.2,
    y = 0.825,
    label = toString(round(annotationData$x[3], 2)),
    fontface = "bold",
    hjust = "left",
    vjust = "top",
    size = smallTxtSize,
    size.unit = "pt",
    color = "blue"
  ) +
  guides(linetype = guide_none()) +
  scale_linetype_manual(values = c("dashed", "longdash", "dotdash")) +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

G <- ggplotGrob(P)
G$layout[G$layout$name == "panel", "z"] <- max(G$layout$z) + 1
G$layout[G$layout$name == "guide-box-inside", "z"] <- max(G$layout$z) + 1
grid.draw(G)

png(
  "PaperFigures/score_WW_SS_periodogram_strong_wide.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

grid.draw(G)

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/score_WW_SS_periodogram_strong_wide.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

##### Non-Stratified Dinucleotide Frequency Plots #####
directory <- "NucleosomePattern/composite_models/composite_ice"

N <- read.table(paste0(directory, "/all_wide_sequence_analysis/n.txt"))$V1
all_WW = read.table(
  paste0(directory, "/all_wide_sequence_analysis/WW.tsv"),
  sep = "\t"
)$V1
all_SS = read.table(
  paste0(directory, "/all_wide_sequence_analysis/SS.tsv"),
  sep = "\t"
)$V1

chroms <- read.table(
  paste0(
    "NucleosomePattern/WideReferenceSequenceAnalysis",
    "/brogaard_all/seqs/chromosomes.txt"
  )
)$V1
brogaard_N <- 0
for (chrom in chroms) {
  brogaard_N <- brogaard_N + nrow(
    read.table(
      paste0(
        "NucleosomePattern/WideReferenceSequenceAnalysis",
        "/brogaard_all/seqs/",
        chrom,
        ".tsv"
      )
    )
  )
}
brogaard_WW = read.table(
  "NucleosomePattern/WideReferenceSequenceAnalysis/brogaard_all/WW.tsv",
  sep = "\t"
)$V1
brogaard_SS = read.table(
  "NucleosomePattern/WideReferenceSequenceAnalysis/brogaard_all/SS.tsv",
  sep = "\t"
)$V1

chroms <- read.table(
  paste0(
    "NucleosomePattern/WideReferenceSequenceAnalysis",
    "/weiner/seqs/chromosomes.txt"
  )
)$V1
weiner_N <- 0
for (chrom in chroms) {
  weiner_N <- weiner_N + nrow(
    read.table(
      paste0(
        "NucleosomePattern/WideReferenceSequenceAnalysis",
        "/weiner/seqs/",
        chrom,
        ".tsv"
      )
    )
  )
}
weiner_WW = read.table(
  "NucleosomePattern/WideReferenceSequenceAnalysis/weiner/WW.tsv",
  sep = "\t"
)$V1
weiner_SS = read.table(
  "NucleosomePattern/WideReferenceSequenceAnalysis/weiner/SS.tsv",
  sep = "\t"
)$V1

annotationData <- data.frame(
  x = rep(10.1*(-14:14), 2),
  y = rep(c(-Inf, Inf), each = 29),
  group = rep(-14:14, 2)
)

# > WW Frequency ----
plotData <- data.frame(
  frequency = c(
    all_WW[351:650] / N,
    brogaard_WW[351:650] / brogaard_N,
    weiner_WW[351:650] / weiner_N
  ),
  position = rep(seq(-149.5, 149.5, 1), 3),
  group = rep(c("CPD", "brogaard", "weiner"), each = 300)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x = x, y = y, group = group),
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  geom_line(
    data = plotData,
    mapping = aes(x = position, y = frequency, color = group),
    linewidth = 1.25
  ) +
  plotTheme +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-150, 150),
    breaks = seq(-150, 150, 50),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Normalized\nAA/AT/TA/TT Frequency",
    limits = c(0.25, 0.5),
    breaks = 5:10 * 0.05,
    expand = c(0,0)
  ) +
  scale_color_discrete(
    name = NULL,
    type = c("blue", "red", "black"),
    labels = c("Chemical Cleavage", "CPD-seq", "MNase-seq")
  ) +
  theme(legend.direction="horizontal", legend.position.inside=c(0.025, 0.1))

P

png(
  "PaperFigures/all_sets_WW_frequency.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

write.table(
  plotData,
  file = "PaperFigures/FigureData/all_sets_WW_frequency.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

##### pLogos Near Nucleosome Edge #####
load_seqs <- function(directory) {
  chroms <- read.table(paste0(directory, "/chromosomes.txt"))$V1
  data <- data.frame(chrom = character(0), loc = integer(0))
  Cnms <- c("chrom", "loc")
  for (i in -500:500) {
    data <- cbind(data, V = character(0))
    colnames(data) <- c(Cnms, toString(i))
    Cnms <- colnames(data)
  }
  for (chrom in chroms) {
    seqData <- read.table(paste0(directory, "/", chrom, ".tsv"))
    #data <- c(data, apply(seqData, 1, paste0, collapse = ""))
    seqData <- cbind(chrom = rep(chrom, nrow(seqData)), seqData)
    colnames(seqData) <- Cnms
    data <- rbind(data, seqData)
  }
  
  return(data)
}

viterbi <- load_seqs(
  paste0(
    "NucleosomePattern/composite_models/composite_ice",
    "/all_wide_sequence_analysis/seqs"
  )
)
brogaard <- load_seqs(
  "NucleosomePattern/WideReferenceSequenceAnalysis/brogaard_all/seqs"
)
weiner<-load_seqs("NucleosomePattern/WideReferenceSequenceAnalysis/weiner/seqs")
redund<-load_seqs("NucleosomePattern/WideReferenceSequenceAnalysis/redund/seqs")

# > Weiner Left ----
seqs <- character(nrow(weiner))
for (i in -88:-69) {
  seqs <- paste0(seqs, weiner[[toString(i)]])
}

write.table(
  seqs,
  file = "PaperFigures/FigureData/weiner_edge_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-80,120),
    breaks = seq(-80, 120, 40)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:20,
    labels = -88:-69,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("MNase-seq Nucleosome Edge"))))

P

png(
  "PaperFigures/weiner_edge_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Viterbi Left ----
seqs <- character(nrow(viterbi))
for (i in -88:-69) {
  seqs <- paste0(seqs, viterbi[[toString(i)]])
}

write.table(
  seqs,
  file = "PaperFigures/FigureData/viterbi_edge_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-80,120),
    breaks = seq(-80, 120, 40)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:20,
    labels = -88:-69,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("CPD-seq Nucleosome Edge"))))

P

png(
  "PaperFigures/viterbi_edge_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Brogaard Left ----
seqs <- character(nrow(brogaard))
for (i in -88:-69) {
  seqs <- paste0(seqs, brogaard[[toString(i)]])
}

write.table(
  seqs,
  file = "PaperFigures/FigureData/brogaard_edge_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-80,120),
    breaks = seq(-80, 120, 40)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:20,
    labels = -88:-69,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Chemical Cleavage Nucleosome Edge"))))

P

png(
  "PaperFigures/brogaard_edge_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Weiner Right ----
seqs <- character(nrow(weiner))
for (i in 69:88) {
  seqs <- paste0(seqs, weiner[[toString(i)]])
}

write.table(
  seqs,
  file = "PaperFigures/FigureData/weiner_edge2_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-80,120),
    breaks = seq(-80, 120, 40),
    position = "right"
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:20,
    labels = 69:88,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("MNase-seq Nucleosome Edge"))))

P

png(
  "PaperFigures/weiner_edge2_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Viterbi Right ----
seqs <- character(nrow(viterbi))
for (i in 69:88) {
  seqs <- paste0(seqs, viterbi[[toString(i)]])
}

write.table(
  seqs,
  file = "PaperFigures/FigureData/viterbi_edge2_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-80,120),
    breaks = seq(-80, 120, 40)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:20,
    labels = 69:88,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("CPD-seq Nucleosome Edge"))))

P

png(
  "PaperFigures/viterbi_edge2_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

# > Brogaard Left ----
seqs <- character(nrow(brogaard))
for (i in 69:88) {
  seqs <- paste0(seqs, brogaard[[toString(i)]])
}

write.table(
  seqs,
  file = "PaperFigures/FigureData/brogaard_edge2_seqLogo.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- plot_sequence_logo(
  seqs,
  bg_freqs = yeast_bg,
  method = "pLogo",
  doAnnotate = FALSE
) +
  scale_y_continuous(
    name = "Log Odds",
    expand = c(0,0),
    limits = c(-80,120),
    breaks = seq(-80, 120, 40)
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    breaks = 1:20,
    labels = 69:88,
    expand = c(0,0)
  ) +
  ggtitle(bquote(underline(bold("Chemical Cleavage Nucleosome Edge"))))

P

png(
  "PaperFigures/brogaard_edge2_seqpLogo.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)
P
dev.off()

##### Allowed Nucleosomes By Rotational Score #####

rot <- read.table(
  "NucleosomePattern/composite_models/composite_ice/rotational_scores.tsv",
  col.names = c("chrom", "loc", "score")
)
rot <- rot[rev(sort(rot$score, index.return = TRUE)$ix),]

N <- nrow(rot)
P <- floor(N/10)
P9 <- rot[1:P,]
P1 <- rot[(N-P):N,]

# Load Scores
chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "genome_scores/chromosomes.txt"
  ),
  header = FALSE
)$V1
scores <- list()
for (chrom in chroms) {
  scores[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/genome_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[chrom]]) <- c("position", "score")
}

W <- c(-73,73) # The window width to graph around the true dyad

# > Weak ----
allowed <- rep(0, W[2]-W[1] + 1)
n <- 0
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  nucleosomes <- P1[P1$chrom == chrom,]
  for (i in 1:nrow(nucleosomes)) {
    dyad <- nucleosomes$loc[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      allowed[theseScores > 0] <- allowed[theseScores > 0] + 1
      n <- n + 1
    }
  }
};plotData <- data.frame(
  relative_position = W[1]:W[2],
  frequency = allowed / n
);write.table(
  plotData,
  file = "PaperFigures/FigureData/allowed_dyads_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

annotationData <- data.frame(
  x = rep(10.1*-7:7, each = 2),
  y = rep(c(-Inf, Inf), 15),
  group = rep(1:15, each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x=x, y=y, group=group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_rect(
    data = plotData,
    mapping = aes(
      xmin = relative_position - 0.5,
      xmax = relative_position + 0.5,
      ymin = 0,
      ymax = frequency,
      fill = frequency > 0.5
    )
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-73.5, 73.5),
    breaks = -1:1*73,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Proportion with Nucleosome-\nFavoring Score Value",
    limits = c(0,1),
    breaks = 0:4/4,
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle(bquote(underline(bold("Low Rotational Score")))) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0.5, 0.5),
    linetype = "solid",
    linewidth = 1.25
  ) +
  scale_fill_discrete(
    type = c("#777777", "black")
  ) +
  guides(fill = guide_none())

P

png(
  "PaperFigures/allowed_dyads_weak.png",
  width = 1250,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Strong ----
allowed <- rep(0, W[2]-W[1] + 1)
n <- 0
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  nucleosomes <- P9[P9$chrom == chrom,]
  for (i in 1:nrow(nucleosomes)) {
    dyad <- nucleosomes$loc[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      allowed[theseScores > 0] <- allowed[theseScores > 0] + 1
      n <- n + 1
    }
  }
};plotData <- data.frame(
  relative_position = W[1]:W[2],
  frequency = allowed / n
);write.table(
  plotData,
  file = "PaperFigures/FigureData/allowed_dyads_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

annotationData <- data.frame(
  x = rep(10.1*-7:7, each = 2),
  y = rep(c(-Inf, Inf), 15),
  group = rep(1:15, each = 2)
)

P <- ggplot() +
  geom_line(
    data = annotationData,
    mapping = aes(x=x, y=y, group=group),
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_rect(
    data = plotData,
    mapping = aes(
      xmin = relative_position - 0.5,
      xmax = relative_position + 0.5,
      ymin = 0,
      ymax = frequency,
      fill = frequency > 0.5
    )
  ) +
  scale_x_continuous(
    name = "Position Relative to Dyad (bp)",
    limits = c(-73.5, 73.5),
    breaks = -1:1*73,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Proportion with Nucleosome-\nFavoring Score Value",
    limits = c(0,1),
    breaks = 0:4/4,
    expand = c(0,0)
  ) +
  plotTheme +
  ggtitle(bquote(underline(bold("High Rotational Score")))) +
  annotate(
    geom = "line",
    x = c(-Inf, Inf),
    y = c(0.5, 0.5),
    linetype = "solid",
    linewidth = 1.25
  ) +
  scale_fill_discrete(
    type = c("#777777", "black")
  ) +
  guides(fill = guide_none())

P

png(
  "PaperFigures/allowed_dyads_strong.png",
  width = 1250,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

rm(chrom,scores,chroms,nucleosomes,W,allowed,n,plotData,P,theseScores,dyad,i,
   disallowed, P1, P9, rot)

##### Nucleosomes Score Matrices By Rotational Score #####

rot <- read.table(
  "NucleosomePattern/composite_models/composite_ice/rotational_scores.tsv",
  col.names = c("chrom", "loc", "score")
)
rot <- rot[rev(sort(rot$score, index.return = TRUE)$ix),]

N <- nrow(rot)
P <- floor(N/10)
P9 <- rot[1:P,]
P1 <- rot[(N-P):N,]

# Load Scores
chroms <- read.table(
  paste0(
    "NucleosomePattern/composite_models/composite_ice/",
    "genome_scores/chromosomes.txt"
  ),
  header = FALSE
)$V1
scores <- list()
for (chrom in chroms) {
  scores[[chrom]] <- read.table(
    paste0(
      "NucleosomePattern/composite_models/composite_ice/genome_scores/",
      chrom,
      ".tsv"
    ),
    header = FALSE,
    sep = "\t"
  )
  colnames(scores[[chrom]]) <- c("position", "score")
}

W <- c(-73,73) # The window width to graph around the true dyad

# > Weak ----
data <- matrix(nrow = 0, ncol = W[2] - W[1] + 1)
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  nucleosomes <- P1[P1$chrom == chrom,]
  for (i in 1:nrow(nucleosomes)) {
    dyad <- nucleosomes$loc[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      data <- rbind(data, theseScores)
      rownames(data) <- c()
    }
  }
};write.table(
  data,
  file = "PaperFigures/FigureData/score_matrix_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

S <- 0:100/100*10
full <- data.frame(
  y = rep(0:100, each = W[2] - W[1] + 1 + 1),
  score = rep(S, each = W[2] - W[1] + 1 + 1),
  positions = rep(0:(W[2] - W[1] + 1), 101),
  count = 0
)
for (s in S) {
  these <- rowSums(data >= s)
  counts <- table(these)
  for (j in names(counts)) {
    full$count[
      (full$score == s) & (full$positions == as.integer(j))
    ] <- counts[j] / nrow(data)
  }
}

write.table(
  full,
  file = "PaperFigures/FigureData/score_heatmap_weak.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_rect(
    data = full,
    mapping = aes(
      xmin = positions - 0.5,
      xmax = positions + 0.5,
      ymin = y - 0.5,
      ymax = y + 0.5,
      fill = log10(count)
    )
  ) +
  scale_x_continuous(
    name = "Number of positions within 73 bp of dyad",
    limits = c(-0.5, 147.5),
    breaks = seq(0, 147, 49),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Nucleosome Score",
    limits = c(-0.5,101.5),
    expand = c(0,0),
    breaks = c(0,25,50,75,100),
    labels = c(0,2.5,5,7.5,10)
  ) +
  scale_fill_continuous(
    type = "viridis",
    na.value = "#440154",
    name = bquote(bold("Proportion ≥ Score\n(log10)")),
    breaks = seq(-4,0,1),
    limits = c(-4,0),
    labels = c("<-4", "-3", "-2", "-1", "0")
  ) +
  plotTheme +
  theme(legend.position = "right") +
  ggtitle(bquote(underline(bold("Low Rotational Score"))))

P

png(
  "PaperFigures/dyad_heatmap_weak.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

# > Strong ----
data <- matrix(nrow = 0, ncol = W[2] - W[1] + 1)
for (chrom in chroms) {
  print(paste0("> ", chrom, "..."))
  nucleosomes <- P9[P9$chrom == chrom,]
  for (i in 1:nrow(nucleosomes)) {
    dyad <- nucleosomes$loc[i]
    theseScores <- scores[[chrom]]$score[
      (scores[[chrom]]$position >= dyad + W[1]) &
        (scores[[chrom]]$position <= dyad + W[2])
    ]
    
    if (length(theseScores) == W[2] - W[1] + 1) {
      data <- rbind(data, theseScores)
      rownames(data) <- c()
    }
  }
};write.table(
  data,
  file = "PaperFigures/FigureData/score_matrix_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

S <- 0:100/100*10
full <- data.frame(
  y = rep(0:100, each = W[2] - W[1] + 1 + 1),
  score = rep(S, each = W[2] - W[1] + 1 + 1),
  positions = rep(0:(W[2] - W[1] + 1), 101),
  count = 0
)
for (s in S) {
  these <- rowSums(data >= s)
  counts <- table(these)
  for (j in names(counts)) {
    full$count[
      (full$score == s) & (full$positions == as.integer(j))
    ] <- counts[j] / nrow(data)
  }
}

write.table(
  full,
  file = "PaperFigures/FigureData/score_heatmap_strong.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

P <- ggplot() +
  geom_rect(
    data = full,
    mapping = aes(
      xmin = positions - 0.5,
      xmax = positions + 0.5,
      ymin = y - 0.5,
      ymax = y + 0.5,
      fill = log10(count)
    )
  ) +
  scale_x_continuous(
    name = "Number of positions within 73 bp of dyad",
    limits = c(-0.5, 147.5),
    breaks = seq(0, 147, 49),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = "Nucleosome Score",
    limits = c(-0.5,101.5),
    expand = c(0,0),
    breaks = c(0,25,50,75,100),
    labels = c(0,2.5,5,7.5,10)
  ) +
  scale_fill_continuous(
    type = "viridis",
    na.value = "#440154",
    name = bquote(bold("Proportion ≥ Score\n(log10)")),
    breaks = seq(-4,0,1),
    limits = c(-4,0),
    labels = c("<-4", "-3", "-2", "-1", "0")
  ) +
  plotTheme +
  theme(legend.position = "right") +
  ggtitle(bquote(underline(bold("High Rotational Score"))))

P

png(
  "PaperFigures/dyad_heatmap_strong.png",
  width = 2500,
  height = 1000,
  units = "px",
  res = 300
)

P

dev.off()

rm(chrom,scores,chroms,nucleosomes,W,allowed,n,plotData,P,theseScores,dyad,i,
   disallowed)

##### Model of Rotationally In-Phase Nucleosomes Pre-Shift #####

X <- seq(0, 10, 0.01)
Y <- cos(2*pi*X)
nearPeak <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearTrough <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearTrough] <- "trough"
colours[nearPeak] <- "peak"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X <= 3.25) | (5.75 <= X)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

boundaries <- plotData[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

plotData$c[plotData$x <= 3.25] <- paste0("red", plotData$c[plotData$x <= 3.25])
plotData$c[plotData$x >= 5.75] <- paste0("blue", plotData$c[plotData$x >= 5.75])

m <- 0
t <- 1.5
tt <- 1.75
teeth <- data.frame(
  x = c(0.25, 0.5, 0.75,
        1.25, 1.5, 1.75,
        2.25, 2.5, 2.75,
        6.25, 6.5, 6.75,
        7.25, 7.5, 7.75,
        8.25, 8.5, 8.75,
        9.25, 9.5, 9.75),
  y = c(t, m, t,
        t, m, tt,
        tt, m, tt,
        tt, m, tt,
        tt, m, t,
        t, m, t,
        t, m, t
  ),
  id = rep(1:7, each = 3),
  c = rep(c(rep("redinbetween", 3), rep("blueinbetween", 4)), each = 3)
)

ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_ellipse(
    data = data.frame(
      x = c(-2,11),
      y = c(2.1, 2.1),
      c = c("redinbetween", "blueinbetween")
    ),
    mapping = aes(
      x0 = x,
      y0 = y,
      a = 5.25,
      b = 1,
      angle = 0,
      fill = c,
      color = c
    )
  ) +
  scale_x_continuous(name = NULL, expand = c(0,0)) +
  coord_cartesian(clip = "on", xlim = c(0,10)) +
  guides(linetype = guide_none(), color = guide_none(), fill = guide_none()) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  ) +
  scale_color_discrete(
    type = c(
      "blue",
      "#BBBBFF",
      "#000033",
      "purple",
      "#FFAAFF",
      "red",
      "#FFAAAA",
      "#440000", 
      "#440044"
    )
  ) +
  scale_fill_discrete(type = c("blue", "red")) +
  geom_polygon(
    data = teeth,
    mapping = aes(
      x = x,
      y = y,
      group = id,
      fill = c
    )
  )

ggsave(
  "PaperFigures/RotationallyInPhaseModel_1.png",
  width = 2000,
  height = 500,
  units = "px"
)

rm(X, Y, colours, linetype, groups, plotData, boundaries, nearPeak)

##### Model of Rotationally Out-of-Phase Nucleosomes Pre-Shift #####

X <- seq(0, 6.25, 0.01)
Y <- cos(2*pi*X)
nearPeak <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearTrough <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearTrough] <- "trough"
colours[nearPeak] <- "peak"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X <= 3.25)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

boundaries <- plotData[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

plotData$c <- paste0("red", plotData$c)

X <- seq(3.25, 10, 0.01)
Y <- cos(2*pi*X + pi)
nearTrough <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearPeak <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearPeak] <- "peak"
colours[nearTrough] <- "trough"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X >= 6.25)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData2 <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData2[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x + pi)
plotData2 <- rbind(
  plotData2,
  boundaries
)

boundaries <- plotData2[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x + pi)
plotData2 <- rbind(
  plotData2,
  boundaries
)

plotData2$c <- paste0("blue", plotData2$c)

m <- 0
t <- 1.5
tt <- 1.75
teeth <- data.frame(
  x = c(0.25, 0.5, 0.75,
        1.25, 1.5, 1.75,
        2.25, 2.5, 2.75,
        6.75, 7, 7.25,
        7.75, 8, 8.25,
        8.75, 9, 9.25,
        9.75, 10, 10.25),
  y = c(t, m, t,
        t, m, tt,
        tt, m, tt,
        tt, m, tt,
        tt, m, t,
        t, m, t,
        t, m, t
  ),
  id = rep(1:7, each = 3),
  c = rep(c(rep("redinbetween", 3), rep("blueinbetween", 4)), each = 3)
)

ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_line(
    data = plotData2,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_ellipse(
    data = data.frame(
      x = c(-2,11.5),
      y = c(2.1, 2.1),
      c = c("redinbetween", "blueinbetween")
    ),
    mapping = aes(
      x0 = x,
      y0 = y,
      a = 5.25,
      b = 1,
      angle = 0,
      fill = c,
      color = c
    )
  ) +
  scale_x_continuous(name = NULL, expand = c(0,0)) +
  coord_cartesian(clip = "on", xlim = c(0,10)) +
  guides(linetype = guide_none(), color = guide_none(), fill = guide_none()) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  ) +
  scale_color_discrete(
    type = c(
      "blue",
      "#BBBBFF",
      "#000033",
      "red",
      "#FFAAAA",
      "#440000"
    )
  ) +
  scale_fill_discrete(type = c("blue", "red")) +
  geom_polygon(
    data = teeth,
    mapping = aes(
      x = x,
      y = y,
      group = id,
      fill = c
    )
  )

ggsave(
  "PaperFigures/RotationallyOutOfPhaseModel_1.png",
  width = 2000,
  height = 500,
  units = "px"
)

rm(X, Y, colours, linetype, groups, plotData, boundaries, nearPeak, plotData2)


##### Model for Rotationally In-Phase Nucleosomes #####

X <- seq(0, 10, 0.01)
Y <- cos(2*pi*X)
nearPeak <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearTrough <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearTrough] <- "trough"
colours[nearPeak] <- "peak"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X <= 3.25) | (5.75 <= X)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

boundaries <- plotData[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

plotData$c[plotData$x <= 3.25] <- paste0("red", plotData$c[plotData$x <= 3.25])
plotData$c[plotData$x >= 5.75] <- paste0("blue", plotData$c[plotData$x >= 5.75])

m <- 0
t <- 1.5
tt <- 1.75
teeth <- data.frame(
  x = c(0.25, 0.5, 0.75,
        1.25, 1.5, 1.75,
        2.25, 2.5, 2.75,
        3.25, 3.5, 3.75,
        4.25, 4.5, 4.75,
        5.25, 5.5, 5.75,
        6.25, 6.5, 6.75,
        7.25, 7.5, 7.75,
        9.25, 9.5, 9.75),
  y = c(t, m, t,
        t, m, t,
        t, m, t,
        t, m, t,
        t, m, t,
        t, m, t,
        t, m, tt,
        tt, m, tt,
        tt, m, tt
  ),
  id = rep(1:9, each = 3),
  c = rep(c(rep("redinbetween", 8), "blueinbetween"), each = 3)
)

ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_ellipse(
    data = data.frame(
      x = c(3, 14),
      y = c(2.1, 2.1),
      c = c("redinbetween", "blueinbetween")
    ),
    mapping = aes(
      x0 = x,
      y0 = y,
      a = 5.25,
      b = 1,
      angle = 0,
      color = c,
      fill = c
    ),
  ) +
  scale_x_continuous(name = NULL, expand = c(0,0)) +
  coord_cartesian(clip = "on", xlim = c(0,10)) +
  guides(linetype = guide_none(), color = guide_none(), fill = guide_none()) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  ) +
  scale_color_discrete(
    type = c(
      "blue",
      "#BBBBFF",
      "#000033",
      "purple",
      "#FFAAFF",
      "red",
      "#FFAAAA",
      "#440000", 
      "#440044"
    )
  ) +
  scale_fill_discrete(type = c("blue", "red")) +
  geom_polygon(
    data = teeth,
    mapping = aes(
      x = x,
      y = y,
      group = id,
      fill = c
    )
  )

ggsave(
  "PaperFigures/RotationallyInPhaseModel.png",
  width = 2000,
  height = 500,
  units = "px"
)

rm(X, Y, colours, linetype, groups, plotData, boundaries, nearPeak)

##### Model for Rotationally Out-of-Phase Nucleosomes #####

X <- seq(0, 6.25, 0.01)
Y <- cos(2*pi*X)
nearPeak <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearTrough <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearTrough] <- "trough"
colours[nearPeak] <- "peak"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X <= 3.25)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

boundaries <- plotData[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x)
plotData <- rbind(
  plotData,
  boundaries
)

plotData$c <- paste0("red", plotData$c)

X <- seq(3.25, 10, 0.01)
Y <- cos(2*pi*X + pi)
nearTrough <- (X %% 1 < 0.101) | (X %% 1 > 0.899)
nearPeak <- (X %% 0.5 < 0.101) | (X %% 0.5 > 0.399)
colours <- rep("inbetween", length(X))
colours[nearPeak] <- "peak"
colours[nearTrough] <- "trough"
groups <- cumsum(
  c(FALSE, colours[2:length(colours)] != colours[1:(length(colours)-1)])
)
linetype <- rep("outside", length(X))
linetype[(X >= 6.25)] <- "inside"

groups[linetype == "inside"] <- groups[linetype == "inside"] + max(groups) + 1

plotData2 <- data.frame(
  x = X,
  y = Y,
  c = colours,
  g = groups,
  l = linetype
)

boundaries <- plotData2[
  (colours != "inbetween") & (
    c(colours[2:length(colours)] == "inbetween", FALSE) |
      c(FALSE, colours[1:(length(colours)-1)] == "inbetween")
  ),
]
boundaries$x <- boundaries$x - 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x + pi)
plotData2 <- rbind(
  plotData2,
  boundaries
)

boundaries <- plotData2[
  (linetype != "outside") & (
    c(linetype[2:length(colours)] == "outside", FALSE) |
      c(FALSE, linetype[1:(length(colours)-1)] == "outside")
  ),
]
boundaries$x <- boundaries$x + 0.01*(2*(1:nrow(boundaries) %% 2) - 1)
boundaries$y <- cos(2*pi*boundaries$x + pi)
plotData2 <- rbind(
  plotData2,
  boundaries
)

plotData2$c <- paste0("blue", plotData2$c)

m <- 0
t <- 1.5
tt <- 1.75
teeth <- data.frame(
  x = c(0.25, 0.5, 0.75,
        1.25, 1.5, 1.75,
        2.25, 2.5, 2.75,
        3.25, 3.5, 3.75,
        4.25, 4.5, 4.75,
        5.25, 5.5, 5.75,
        6.25, 6.5, 6.75,
        7.25, 7.5, 7.75,
        9.75, 10, 10.25),
  y = c(t, m, t,
        t, m, t,
        t, m, t,
        t, m, t,
        t, m, t,
        t, m, t,
        t, m, tt,
        tt, m, tt,
        tt, m, tt
  ),
  id = rep(1:9, each = 3),
  c = rep(c(rep("redinbetween", 8), "blueinbetween"), each = 3)
)

ggplot() +
  geom_line(
    data = plotData,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_line(
    data = plotData2,
    mapping = aes(
      x = x,
      y = y,
      colour = c,
      linetype = l,
      group = g
    ),
    linewidth = 1.25
  ) +
  geom_ellipse(
    data = data.frame(
      x = c(3, 14.5),
      y = c(2.1, 2.1),
      c = c("redinbetween", "blueinbetween")
    ),
    mapping = aes(
      x0 = x,
      y0 = y,
      a = 5.25,
      b = 1,
      angle = 0,
      color = c,
      fill = c
    ),
  ) +
  scale_x_continuous(name = NULL, expand = c(0,0)) +
  coord_cartesian(clip = "on", xlim = c(0,10)) +
  guides(linetype = guide_none(), color = guide_none(), fill = guide_none()) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  ) +
  scale_color_discrete(
    type = c(
      "blue",
      "#BBBBFF",
      "#000033",
      "red",
      "#FFAAAA",
      "#440000"
    )
  ) +
  scale_fill_discrete(type = c("blue", "red")) +
  geom_polygon(
    data = teeth,
    mapping = aes(
      x = x,
      y = y,
      group = id,
      fill = c
    )
  )

ggsave(
  "PaperFigures/RotationallyOutOfPhaseModel.png",
  width = 2000,
  height = 500,
  units = "px"
)

rm(X, Y, colours, linetype, groups, plotData, boundaries, nearPeak, plotData2)

