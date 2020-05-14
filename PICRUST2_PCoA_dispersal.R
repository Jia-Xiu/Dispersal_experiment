# analyze results from functional annotation by PICRUST2  
# Author: Jia Xiu
# Date: 2020-05-08

rm(list=ls())

# change directory   
directory = '~/Dropbox/'
subfolder = 'Dispersal/community_analysis'

setwd(paste(directory, subfolder, sep="/"))
getwd()

library(vegan) 
library(ape) 
library(gridExtra)
library(ggplot2) 
library(ggpubr)
library(KEGGREST)
listDatabases()


mytheme <- theme_bw()+
  theme(text = element_text(size=12),
        plot.title = element_text(size=12, face="bold"),
        strip.background = element_blank(),
        strip.text=element_text(face="bold", color="#666666", size = 12),
        legend.box.background = element_rect(),
        legend.title = element_text(face = "bold"),
        legend.box.margin = margin(1, 1, 1, 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  

### PCoA plot start from here


# All
meg_0 <- read.table("pred_metagenome_unstrat_0_year.tsv", sep = '\t', header = 1, row.names = 1, check.names = FALSE)
meg_70 <- read.table("pred_metagenome_unstrat_70_year.tsv", sep = '\t', header = 1, row.names = 1, check.names = FALSE)
meg <- transform(merge(meg_0, meg_70, all = TRUE, by="row.names"), row.names=Row.names, Row.names=NULL)
meg[is.na(meg)] <- 0
colnames(meg) <- gsub("^X", "", colnames(meg))
meg <- t(meg)
meg <- as.data.frame(meg)
# warning("remove day 0")
# meg <- subset(meg, !grepl('___', row.names(meg)))
meg <- meg[, !colSums(meg)==0]
dim(meg)
meg[1:5, 1:3]

# PCoA analysis based on bray curtis or jaccard distance
dist <- vegdist(meg, method = "bray", binary=FALSE, diag=1) 
re <- pcoa(dist, correction="none", rn=NULL)
str(re)

group_info <- data.frame(row.names = rownames(re$vectors), t(as.data.frame(strsplit(rownames(re$vectors),"_"))))
head(group_info)

df <- data.frame(pc1 = re$vectors[,1], 
                 pc2 = re$vectors[,2],
                 pc3 = re$vectors[,3],
                 Soil = as.factor(group_info[,1]),
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]))

df$Soil <- factor(df$Soil, levels = c("0", "70"), labels = c("0 years", "70 years"))
df$Day <- factor(df$Day, levels=c("0", "2", "4", "8", "12", "16", "20"))
head(df); str(df)

my_palette <- c("#999999", "#E69F00", "#56B4E9","#F0E442",  "#009E73", "#CC79A7", "#0072B2")

(p0 <- ggplot(df, aes(pc1, pc2, shape = Soil, fill=Day))+
    geom_point(size=5, alpha=0.7)+
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 22))+ 
    stat_ellipse(geom = "path", type="norm", linetype = 2, alpha=0.4, aes(group=Soil)) +
    annotate("text", x = -0.08, y = -0.08, label = "70 years") +
    annotate("text", x = 0.05, y = -0.1, label = "0 years") +
    labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1] * 100, 2), "%)", sep=""), 
         y=paste("PCoA2 (", round(re$values$Relative_eig[2] * 100, 2), "%)", sep=""),
         title = "")+
    mytheme)

# ggsave("PCoA_metagenome_bray_Year_day.jpg", width = 8.5, height = 6.5, units = "cm", device = "jpeg", p0, scale = 1.5, dpi = 300)

# PERMANOVA
# (result_whole <- adonis(dist ~ Soil*Day*Water*Frequency, data=df, permutation=9999))
# result_whole$aov.tab
# write.csv(result_whole$aov.tab, "permanova_metagenome_all.csv")



# each soil seperately
p <- list()
permanova.list <- list()
years <- c("0", "70")

dist_methods <- readline("\ndistance methods is? e.g. bray, jaccard etc.")

for (yr in years) {
  if (!dist_methods == "bray" | dist_methods == "jaccard") 
    stop("select a distance method for PCoA analysis!")
  
  cat("\nanalysis for year", yr)
  
  # load predicted metagenome 
  meg <- read.table(paste("pred_metagenome_unstrat", yr, "year.tsv", sep = "_"), sep = '\t', header = 1, row.names = 1, check.names = FALSE)
  meg <- t(meg)
  # warning("remove day 0")
  # meg <- subset(meg, !grepl('___', rownames(meg)))
  # meg <- meg[, !colSums(meg)==0]
  str(meg)
  
  # PCoA analysis based on bray curtis or jaccard distance
  dist <- vegdist(meg, method = dist_methods, binary=FALSE, diag=1) 
  re <- pcoa(dist, correction="none", rn=NULL)
  str(re)
  
  group_info <- data.frame(row.names = rownames(re$vectors), t(as.data.frame(strsplit(rownames(re$vectors),"_"))))
  head(group_info)
  
  df <- data.frame(pc1 = re$vectors[,1], 
                   pc2 = re$vectors[,2],
                   pc3 = re$vectors[,3],
                   Frequency = as.factor(group_info[,2]),
                   Water = as.factor(group_info[,3]),
                   Day = as.factor(group_info[,4]),
                   replicates = as.factor(group_info[,5]))
  
  df$Day <- factor(df$Day, levels=c("0", "2", "4", "8", "12", "16", "20"))
  df$Frequency <- factor(df$Frequency, levels = c("", "f1", "f2", "f3", "f4"), 
                         labels = c("", "1x / week", "1x / 3 days", "1x / day", "2x / day"))
  df$Water <- factor(df$Water, levels = c("", "sea", "strl"), 
                     labels = c("", "Sea water", "Sterile water"))
  head(df); str(df)
  
  # calculate average for each treatment
  dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se = sd(x)/ sqrt(length(x))))
  df1 <- summaryBy(pc1 ~ Soil + Frequency + Day + Water, data=df, FUN=dstats)
  df2 <- summaryBy(pc2 ~ Soil + Frequency + Day + Water, data=df, FUN=dstats)
  df0 <- cbind(df1, df2[, c(5, 7)])
  df0$Groups <- factor(paste(df0$Water, df0$Frequency),
                       levels = c(" ", "Sea water 1x / week", "Sea water 1x / 3 days",
                                  "Sea water 1x / day", "Sea water 2x / day", "Sterile water 1x / week", 
                                  "Sterile water 2x / day" ),
                       labels = c("Day 0 untreated soil", "Sea water 1x / week", "Sea water 1x / 3 days",
                                  "Sea water 1x / day", "Sea water 2x / day", "Sterile water 1x / week", 
                                  "Sterile water 2x / day" ))
  levels(df0$Groups)
  head(df0); str(df0)
  
  my_palette <- c("#999999", "#E69F00", "#56B4E9","#F0E442",  "#009E73", "#CC79A7", "#0072B2")
  
  (p[[yr]] <- ggplot(df0, aes(x = pc1.mean, y = pc2.mean, shape = Groups, color = Day, fill = Day))+
      geom_errorbar(aes(ymin = pc2.mean - pc2.se,ymax = pc2.mean + pc2.se), color = "#999999") + 
      geom_errorbarh(aes(xmin = pc1.mean - pc1.se, xmax = pc1.mean + pc1.se), color = "#999999") +
      geom_point(size = 5, alpha = .9)+
      scale_colour_manual(values = my_palette) +
      scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
      scale_shape_manual(values=c(24, 22, 23, 25, 21, 0, 1))+ 
      labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1]*100, 2), "%)", sep=""), 
           y=paste("PCoA2 (", round(re$values$Relative_eig[2]*100, 2), "%)", sep=""),
           title = paste(yr, 'years soil', sep=" "))+
      mytheme)

  # set.seed(123)
  # (result_whole <- adonis(dist ~ Water*Frequency*Day, data=df, permutation=9999))
  # result_whole$aov.tab
  # permanova.list[[yr]] <- as.data.frame(result_whole$aov.tab)
}

#p0 <- do.call(grid.arrange, p, ncol = 2)
(p1 <- (ggarrange(p[[1]], p[[2]], labels = c( "B", "C"), 
                  common.legend = TRUE, legend = "right", ncol = 2, nrow = 1)))
(p2 <- ggarrange(p0, NULL, labels = c( "A", ""), widths = c(1, 0.9), ncol = 2, nrow = 1))
(p2 <- ggarrange(p2, p1, labels = c( "", ""), ncol = 1, nrow = 2))

ggsave(paste("PCoA_metagenome", dist_methods, "all_raw.pdf", sep = "_"), width = 19.5, height = 15, 
       units = "cm", p2, scale = 1.5, dpi = 300)

# str(permanova.list)
# df <- do.call(rbind.data.frame, permanova.list)
# df$Year <- gsub("\\..*", "", row.names(df))
# df$groups <- sub('.*\\.', '', row.names(df))
# write.csv(df, paste("permanova_metagenome_", dist_methods, ".csv", sep = ""))

