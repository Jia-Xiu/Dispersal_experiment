# Multivariate Homogeneity Of Groups Dispersions (Variances) analysis for dispersal experiment
# Author: Xiu Jia
# Date: 05-02-2020

rm(list=ls())

# Load libraries
library(vegan) # for multivariable analysis
library(ggplot2) # for graphing
library(ggpubr) # combine figures
library(RColorBrewer) # for color bar
display.brewer.all()

# change directory

mytheme<-theme_bw()+
  theme(text = element_text(size=12),
        plot.title = element_text(size=12, face="bold"),
        strip.background = element_blank(),
        strip.text=element_text(color="#666666", face="bold", size = 12),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  


# load the rarefied otu table -----------------------------------------------------------------------
com <- read.csv("feature-table-rarified-nontax.csv", sep=";", header=1, row.names=1, check.names = FALSE)
colnames(com) <- gsub("\\-", "\\_", colnames(com))
com <- t(com)
dim(com)
com[1:5, 1:2]
cat("the range of total number of sequences of each species is:", range(apply(com,2,sum)))
cat("the range of total sequences per sample is:", range(apply(com,1,sum)))

# ========================================================================================================
# sea & strl water
# ========================================================================================================

disp.year.list <- list()
  disp.freq.list <- list()
permanova.year.list <- list()
  permanova.freq.list <- list()


years <- c("0", "70")
frequencies <- c("f1", "f4")

for (yr in years) {
  
  comyr <- subset(com, !grepl("___0", row.names(com)))
  comyr <- subset(comyr, grepl(paste("^", yr, "_", sep = ""), row.names(comyr)))
  
  for (fr in frequencies) {
    #subcom <- subset(subcom, grepl("_f1_", row.names(subcom)) | grepl("_f4_", row.names(subcom)))
    subcom <- subset(comyr,  grepl(paste("_", fr, "_", sep = ""), row.names(comyr)))
    dim(subcom)
    
    # bray curtis distance
    dist <- as.data.frame(as.matrix(vegdist(subcom, method="bray", diag=1)))
    dist[1:5, 1:5]
    
    df <- data.frame(row.names=rownames(dist), t(as.data.frame(strsplit(rownames(dist),"_"))))
    
    df$Soil <- factor(df$X1, levels = c("0", "70"), labels = c("0 years", "70 years"))
    df$Frequency <- factor(df$X2, levels = c("f1", "f2", "f3", "f4"), 
                           labels = c("1x / week", "1x / 3 days", "1x / day", "2x / day"))
    df$Water <- factor(df$X3, levels = c("sea", "strl"), 
                       labels = c("Sea Water", "Sterile Water"))
    df$Day <- factor(df$X4, levels=c("2", "4", "8", "12", "16", "20"))
    head(df)
    
    cat("\nDispersion between frequency for year", yr, "& frequency", fr) 
    dist <- as.dist(as.matrix(dist))
    (dispersion <- betadisper(dist, type = c("centroid"), group = df$Water))
    disp.freq.list[[fr]] <- permutest(dispersion, permutation=9999)$tab
    #str(dispersion)
    
  }
  disp.year.list[[yr]] <- disp.freq.list
  
  comyrwater <- subset(comyr, grepl("_f1_", row.names(comyr)) | grepl("_f4_", row.names(comyr)))
  
  dist <- as.data.frame(as.matrix(vegdist(comyrwater, method="bray", diag=1)))
  dist[1:5, 1:5]
  
  df <- data.frame(row.names=rownames(dist), t(as.data.frame(strsplit(rownames(dist),"_"))))
  
  df$Soil <- factor(df$X1, levels = c("0", "70"), labels = c("0 years", "70 years"))
  df$Frequency <- factor(df$X2, levels = c("f1", "f2", "f3", "f4"), 
                         labels = c("1x / week", "1x / 3 days", "1x / day", "2x / day"))
  df$Water <- factor(df$X3, levels = c("sea", "strl"), 
                     labels = c("Sea Water", "Sterile Water"))
  df$Day <- factor(df$X4, levels=c("2", "4", "8", "12", "16", "20"))
  head(df)
 
  dist <- as.dist(as.matrix(dist))
  
  cat("\npermanova (Frequency) for year", yr)
  set.seed(123)
  result_whole <- adonis(dist ~ Water*Frequency, data=df, permutation=9999) 
  result_whole$aov.tab
  permanova.freq.list[[yr]] <- as.data.frame(result_whole$aov.tab)
  
}

#lapply(disp.year.list, function(x) write.table( data.frame(x), 'permutation_test_for_dispersion.csv'  , append= T, sep=',' ))

#lapply(permanova.freq.list, function(x) write.table( data.frame(x), 'PERMANOVA_soil_sea_strl.csv'  , append= T, sep=',' ))



# all -------------------------------------------
subcom <- subset(com, !grepl("___0", row.names(com)))
subcom <- subset(subcom, grepl("_f1_", row.names(subcom)) | grepl("_f4_", row.names(subcom)))
dim(subcom)

# bray curtis distance
dist <- as.data.frame(as.matrix(vegdist(subcom, method="bray", diag=1)))
dist[1:5, 1:5]

df <- data.frame(row.names=rownames(dist), t(as.data.frame(strsplit(rownames(dist),"_"))))

df$Soil <- factor(df$X1, levels = c("0", "70"), labels = c("0 years", "70 years"))
df$Frequency <- factor(df$X2, levels = c("f1", "f4"), 
                       labels = c("1x / week", "2x / day"))
df$Water <- factor(df$X3, levels = c("sea", "strl"), 
                   labels = c("Sea Water", "Sterile Water"))
df$Day <- factor(df$X4, levels=c("2", "4", "8", "12", "16", "20"))
df$Groups <- factor(paste(df$Soil, df$Water, df$Frequency))
levels(df$Groups)
head(df)
str(df)

dist <- as.dist(as.matrix(dist))

cat("\nDispersion between frequency for year", yr, "& frequency", fr) 
(dispersion <- betadisper(dist, type = c("centroid"), group = df$Groups))
plot(dispersion, hull=FALSE, ellipse=TRUE, main = "Dispersion")
boxplot(dispersion, xlab="Species", notch=TRUE)
permutest(dispersion, permutation=9999)$tab
str(dispersion)
head(dispersion)


df <- as.data.frame(dispersion$distances)

groups <- data.frame(row.names=rownames(df), t(as.data.frame(strsplit(rownames(df),"_"))))

df$Soil <- factor(groups$X1, levels = c("0", "70"), labels = c("0 years", "70 years"))
df$Frequency <- factor(groups$X2, levels = c("f1", "f4"), 
                       labels = c("1x / week", "2x / day"))
df$Water <- factor(groups$X3, levels = c("sea", "strl"), 
                   labels = c("Sea Water", "Sterile Water"))
df$Day <- factor(groups$X4, levels=c("2", "4", "8", "12", "16", "20"))
df$Groups <- factor(paste(df$Soil, df$Water, df$Frequency))
levels(df$Groups)
head(df)
str(df)

(my_palette1 = c("#999999", "#000000")) # gray vs. black
(my_palette2 = c("#999999", "#FFFFFF")) # gray vs. white

(p <- ggplot(df, aes(x = Frequency, y = dispersion$distances, fill = Water, color = Water, shape = Water))+ 
    geom_point(size = 3, alpha = 0.3, position = position_jitterdodge()) +
    geom_boxplot(alpha = 0.6, outlier.size=-1) +
    scale_shape_manual(values=c(21, 1)) +
    scale_colour_manual(values = my_palette1) +
    scale_fill_manual(values = my_palette2, guide=guide_legend(override.aes = list(shape=21))) +
    stat_compare_means(aes(group = Water), label = "p.format") + #, label = "p.signif"
    labs(x=" ", y="Distance to centroid") + 
    facet_grid(. ~ Soil) +
    mytheme)

ggsave("Dispersion_distance_to_center_sea_strl.pdf", width = 8, height = 3.5, p, scale = 1, dpi = 600)
  
  
  
