##########################################################################
# This R script examines the periodicity of scores at +1 nucleosomes.    #
##########################################################################

##### SETUP #####
require(ggplot2)
require(jsonlite)
require(stringr)
require(rstudioapi)
require(lomb)

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
directory <- "NucleosomePattern/composite_models/composite_ice"
chroms<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII",
          "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
D <- data.frame(
  dataset = character(0),
  chrom = character(0),
  plus_one_dyad = integer(0),
  TSS_offset = integer(0)
)
for (i in -500:650) {
  D[[paste0("score_", toString(i))]] <- numeric(0)
}

for (chrom in chroms) {
  thisD <- read.table(paste0(directory, "/plusOnes/", chrom, ".tsv"), sep="\t")
  colnames(thisD) <- c("plus_one_dyad", "TSS_offset")
  thisD$dataset <- rep("plus_one", nrow(thisD))
  thisD$chrom <- rep(chrom, nrow(thisD))
  
  thisScoreD <- read.table(
    paste0(directory, "/plusOnes_alignedScores/", chrom, ".tsv"),
    sep = "\t"
  )
  colnames(thisScoreD) <- paste0("score_", unlist(lapply(-500:650, toString)))
  
  D <- rbind(D, cbind(thisD, thisScoreD))
  
  thisD<-read.table(paste0(directory,"/mock_plusOnes/",chrom,".tsv"),sep="\t")
  colnames(thisD) <- c("plus_one_dyad")
  thisD$TSS_offset <- rep(NA, nrow(thisD))
  thisD$dataset <- rep("mock", nrow(thisD))
  thisD$chrom <- rep(chrom, nrow(thisD))
  
  thisScoreD <- read.table(
    paste0(directory, "/mock_plusOnes_alignedScores/", chrom, ".tsv"),
    sep = "\t"
  )
  colnames(thisScoreD) <- paste0("score_", unlist(lapply(-500:650, toString)))
  
  D <- rbind(D, cbind(thisD, thisScoreD))
}

##### PERIODICITY #####

P <- data.frame(
  dataset = character(1005*nrow(D)),
  chrom = character(1005*nrow(D)),
  dyad = character(1005*nrow(D)),
  position = integer(1005*nrow(D)),
  peak_at = numeric(1005*nrow(D)),
  peak = numeric(1005*nrow(D))
)
p <- 1
for (i in 1:nrow(D)) {
  print(paste0("> ", toString(i), "/", toString(nrow(D)), "..."))
  for (j in -427:577) {
    print(paste0(">> ", toString(j), "..."))
    x <- unlist(D[i,paste0("score_", unlist(lapply((j-73):(j+73), toString)))])
    names(x) <- c()
    L <- lsp(x, -73:73, type = "period", to = 15, ofac = 100, plot = FALSE)
    
    P$dataset[p] = D$dataset[i]
    P$chrom[p] = D$chrom[i]
    P$dyad[p] = D$plus_one_dyad[i]
    P$position[p] = j
    P$peak_at[p] = L$peak.at[1]
    P$peak[p] = L$peak[1]
    
    p <- p + 1
  }
};write.table(
  P,
  paste0(directory, "/plusOne_periodicity.tsv"),
  sep = "\t",
  row.names = FALSE
)
