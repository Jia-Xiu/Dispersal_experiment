# Using relativ abundance =< 0.1% cutoffs to define the rare and common biopsheres
# Author: Xiu Jia
# Date: 21-02-2020

# load directory --------------------------------------------------------------------------------------------
setwd()

cutoff = 0.1/100

# source the trucate function
source("TruncateTable.r") 

years <- c("0", "70")

for (year in years) {
  
  cat("# load feature table from", year, "years soil")
  
  com <- read.csv(paste("feature_table", year, "year.csv", sep = "_"), sep=",",  header=1, row.names=1, check.names = FALSE)
  com <- t(com)
  str(com)
  
  # The truncated datasets can be stored as follows: 
  common <-TruncateTable(com, cutoff, typem="dominant") 
  write.csv(t(common), paste("feature_table", year, "year", cutoff, "cutoff_common.csv", sep = "_"))
  
  rare <- TruncateTable(com, cutoff, typem="rare") 
  write.csv(t(rare), paste("feature_table", year, "year", cutoff, "cutoff_rare.csv", sep = "_"))
}


# Using sample-specific rarity cutoffs (rank abundance curve based approach) to define the rare and common biopsheres
# Author: Xiu Jia
# Date: 24-12-2019

# load library
library(vegan)

years <- c("0", "70")

cutoff <- "rac"

for (year in years) {
  
  cat("# load feature table from", year, "years soil")
  
  com <- read.csv(paste("feature_table", year, "year.csv", sep = "_"), sep=",",  header=1, row.names=1, check.names = FALSE)
  com <- t(com)
  str(com)
  
  # calculate Chao1 for estimation of the sequencing depth
  # to be notice, Chao1 could only represent the lower bundary of richness estimation
  Chao <- as.data.frame(t(estimateR(com)))
  Chao$slope <- Chao$S.obs/Chao$S.chao1
  head(Chao)
  
  
  # built a empty matrix to store the sample-specific rarity cutoffs
  cutoffs <- matrix(NA, nrow(com), 3)
  row.names(cutoffs) <- row.names(com)
  cutoffs[,2] <- Chao$slope
  cutoffs[,3] <- Chao$S.obs
  colnames(cutoffs) <- c("Rarity.cutoffs", "Slopes", "S.obs")
  
  # using method calculate H-index to calculate rarity cutoff per sample
  for (j in 1:nrow(com)) {
    com_j = sort(as.numeric(com[j,]), decreasing = TRUE)
    com_j <- com_j[com_j!=0]
    slope <- cutoffs[j,2]
    for (i in 1:length(com_j)){
      if (com_j[i]>=i*slope){
        H=i
      }
    }
    cutoffs[j, 1] <- H
  }
  cutoffs <- as.data.frame(cutoffs)
  head(cutoffs)
  
  
  # generate the dataset of the rare biosphere
  df <- com
  for (j in 1:nrow(df)) {
    for (i in 1:ncol(df)) {
      if (df[j, i] > cutoffs[j,1]) {
        df[j, i] <- NA
      }
    }
  }
  
  df[is.na(df)] <- 0
  df <- df[, colSums(df)!=0] 
  df <- df[, rowSums(df)!=0] 
  write.csv(t(df), paste("feature_table", year, "year_rare.csv", sep = "_"))
  
  # generate the dataset of the common biosphere
  df <- com
  for (j in 1:nrow(df)) {
    for (i in 1:ncol(df))
      if (df[j, i] <= cutoffs[j,1]) {
        df[j,i] <- NA
      }
  }
  
  df[is.na(df)] <- 0
  df <- df[, colSums(df)!=0] 
  df <- df[, rowSums(df)!=0] 
  write.csv(t(df),  paste("feature_table", year, "year", cutoff, "cutoff_common.csv", sep = "_"))
}
