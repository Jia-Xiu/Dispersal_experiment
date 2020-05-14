# beta diversity analysis for dispersal experiment
# Author: Xiu Jia
# Date: 29-01-2020

rm(list=ls())

# Load libraries
library(vegan) # for multivariable analysis
library(ape) # for pcoa 
library(ggplot2) # for graphing
library(ggpubr) # combine figures
library(RColorBrewer) # for color bar
library(plyr) # for rename
library(doBy)
display.brewer.all()

# change directory
setwd()


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


# load the rarefied otu table -----------------------------------------------------------------------
com <- read.csv("feature-table-rarified-nontax.csv", sep=";", header=1, row.names=1, check.names = FALSE)
colnames(com) <- gsub("\\-", "\\_", colnames(com))
write.csv(com, "feature-table-rarified.csv")
com <- t(com)
dim(com)
com[1:5, 1:2]
cat("the range of total number of sequences of each species is:", range(apply(com,2,sum)))
cat("the range of total sequences per sample is:", range(apply(com,1,sum)))

# ========================================================================================================
# for all samples 
# ========================================================================================================

# PCoA analysis  bray curtis distance
dist <- vegdist(com, method="bray", binary=FALSE, diag=1) 
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
df$Water <- factor(df$Water, levels = c("sea", "strl"), labels = c("Sea water", "Sterile water"))
df$Day <- factor(df$Day, levels=c("0", "2", "4", "8", "12", "16", "20"))
df$Frequency <- factor(df$Frequency, levels = c("f1", "f2", "f3", "f4"), 
                       labels = c("1x per week", "1x per 3 days", "1x per day", "2x per day"))
head(df); str(df)

my_palette <- c("#999999", "#E69F00", "#56B4E9","#F0E442",  "#009E73", "#CC79A7", "#0072B2")

(p0 <- ggplot(df, aes(pc1, pc2, shape=Soil, fill=Day))+
    geom_point(size=5, alpha=0.7)+
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 22))+ 
    stat_ellipse(geom = "path", type="norm", linetype = 2, alpha=0.4, aes(group=Soil)) +
    annotate("text", x = -0.4, y = .25, label = "0 years") +
    annotate("text", x = 0.4, y = .1, label = "70 years") + 
    labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1] * 100, 2), "%)", sep=""), 
         y=paste("PCoA2 (", round(re$values$Relative_eig[2] * 100, 2), "%)", sep=""),
         title = "")+
    mytheme)

#ggsave("PCoA_bray_Year_day.jpg", width = 8.5, height = 6.5, units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)

# ========================================================================================================
# for each soil types
# ========================================================================================================

p <- list()
permanova.list <- list()
years <- c("0", "70")

for (yr in years) {
  
cat("analysis for year", yr)
subcom <- subset(com, grepl(paste('^', yr, '_', sep=""), rownames(com)))
subcom <- subcom[, !colSums(subcom) == 0]
#write.csv(t(subcom), paste("feature_table", yr, "year.csv", sep = "_"))

# warning("remove day 0")
# subcom <- subset(subcom, !grepl('___', rownames(subcom)))
# subcom <- subcom[, !colSums(subcom) == 0]
dim(subcom)
subcom[1:5, 1:2]

# PCoA analysis  bray curtis distance
dist <- vegdist(subcom, method="bray", binary=FALSE, diag=1) 
re <- pcoa(dist, correction="none", rn=NULL)
str(re)

group_info <- data.frame(row.names = rownames(re$vectors), t(as.data.frame(strsplit(rownames(re$vectors),"_"))))
head(group_info)

df <- data.frame(pc1 = re$vectors[,1], 
                 pc2 = re$vectors[,2],
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]))

df$Day <- factor(df$Day, levels=c("0", "2", "4", "8", "12", "16", "20"))
df$Frequency <- factor(df$Frequency, levels = c("", "f1", "f2", "f3", "f4"), 
                       labels = c("", "1x / week", "1x / 3 days", "1x / day", "2x / day"))
df$Water <- factor(df$Water, levels = c("", "sea", "strl"), 
                   labels = c("", "Sea water", "Sterile water"))

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

( p[[yr]] <- ggplot(df0, aes(x = pc1.mean, y = pc2.mean, shape = Groups, color = Day, fill = Day))+
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

}


(p1 <- (ggarrange(p[[1]], p[[2]], labels = c( "B", "C"), 
                                  common.legend = TRUE, legend = "right", ncol = 2, nrow = 1)))
(p2 <- ggarrange(p0, NULL, labels = c( "A", ""), widths = c(1, 0.9), ncol = 2, nrow = 1))
(p2 <- ggarrange(p2, p1, labels = c( "", ""), ncol = 1, nrow = 2))

ggsave("PCoA_bray_all_raw.pdf", width = 19, height = 15, 
       units = "cm", p2, scale = 1.5, dpi = 300)


# ========================================================================================================
# Day in Facet
# ========================================================================================================

# daya in facet
subcom <- subset(com, !grepl('___', rownames(com)))
dim(subcom)
subcom[1:5, 1:2]

group_info <- data.frame(row.names=rownames(subcom),t(as.data.frame(strsplit(rownames(subcom),"_"))))
head(group_info)

# PCoA analysis  bray curtis distance
dist <- vegdist(subcom, method="bray", binary=FALSE, diag=1) 
re <- pcoa(dist, correction="none", rn=NULL)
str(re)

df <- data.frame(pc1 = re$vectors[,1], 
                 pc2 = re$vectors[,2],
                 Soil = as.factor(group_info[,1]),
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]))

df$Day <- factor(df$Day, levels=c("2", "4", "8", "12", "16", "20"))
df$Soil <- factor(df$Soil, levels = c("0", "70"), labels = c("0 years", "70 years"))
df$Frequency <- factor(df$Frequency, levels = c("f1", "f2", "f3", "f4", "f0"), 
                       labels = c("1x / week", "1x / 3 days", "1x / day", "2x / day", "untreated soil"))
df$Water <- factor(df$Water, levels = c("sea", "strl"), 
                   labels = c("Sea Water", "Sterile water"))

df <- df[df$Water == "Sea Water", ]

str(df)

(my_palette = c("#00AFBB", "#E7B800"))
(my_palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))
(p <- ggplot(df, aes(x = pc1, y = pc2, shape = Frequency, color = Frequency, fill = Frequency))+
   facet_grid(Soil ~ Day, scales="free") +
   geom_point(size = 3, alpha = .9)+
   scale_colour_manual(values = my_palette) +
   scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
   scale_shape_manual(values=c(24, 22, 23, 25, 21, 0, 1))+ 
   labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1]*100, 2), "%)", sep=""), 
        y=paste("PCoA2 (", round(re$values$Relative_eig[2]*100, 2), "%)", sep=""))+
   mytheme)

ggsave("PCoA_bray_all_seawater_raw.pdf", width = 28, height = 10, units = "cm", p, scale = 1)


# ========================================================================================================
# compare sea water and strl water 
# ========================================================================================================
(my_palette = c("#00AFBB", "#E7B800"))

cat("change the year for analysis")
year <- "70"
subcom <- subset(com, grepl(paste('^', year, '_', sep=""), rownames(com)))
subcom <- subset(subcom, grepl('f1', rownames(subcom)) | grepl('f4', rownames(subcom)))
dim(subcom)
group_info <- data.frame(row.names=rownames(subcom),t(as.data.frame(strsplit(rownames(subcom),"_"))))
head(group_info)

# PCoA analysis  bray curtis distance
dist <- vegdist(subcom, method="bray", binary=FALSE, diag=1) 
re <- pcoa(dist, correction="none", rn=NULL)
str(re)

df <- data.frame(pc1 = re$vectors[,1], 
                 pc2 = re$vectors[,2],
                 Soil = as.factor(group_info[,1]),
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]))

df$Day <- factor(df$Day, levels = c("2", "4", "8", "12", "16", "20"), 
                 labels = c("2 days", "4 days", "8 days", "12 days", "16 days", "20 days"))
df$Frequency <- factor(df$Frequency, levels = c("f1", "f4"), labels = c("1x / week", "2x / day"))
df$Water <- factor(df$Water, levels = c("sea", "strl"), labels = c("Sea water", "Sterile water"))
head(df);str(df)

(p4 <- ggplot(df, aes(pc1, pc2, shape=Frequency, fill=Water))+
    geom_point(size=3, alpha=0.7)+ 
    stat_ellipse(type = "norm", aes(colour = Water, linetype = Frequency)) + 
    facet_grid( ~ Day)+
    labs(x = paste("PCoA1 (", round(re$values$Relative_eig[1] * 100, 2), "%)", sep=""), 
         y = paste("PCoA2 (", round(re$values$Relative_eig[2] * 100, 2), "%)", sep=""),
         title = paste(year, 'years soil', sep=" "))+
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21)))+
    scale_color_manual(values = my_palette)+
    scale_shape_manual(values=c(23, 21))+ 
    mytheme)

(p0 <- annotate_figure((ggarrange(p1, p2, labels = c("A", "B"), 
                                 common.legend = TRUE, legend = "right", ncol = 2, nrow = 1)), right = " "))
(p <- annotate_figure((ggarrange(p0, p3, p4, labels = c( "", "C", "D"), legend = "none",
                                 heights = c(1, 0.7, 0.7), ncol = 1, nrow = 3)), right = " "))

ggsave("PCoA_bray_Sea_Strl.jpg", width = 16, height = 7, units = "cm", device = "jpeg", p0, scale = 1.5, dpi = 300)
ggsave("PCoA_bray_Sea_Strl_facet.jpg", width = 16, height = 16, units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)



# ========================================================================================================
# PERMANOVA
# ========================================================================================================
# remove sample from day 0
com_treatments <- com[!grepl("___0", rownames(com)), ]
cat("the range of total number of sequences of each species is:", range(apply(com_treatments,2,sum)), 
    "\nthe range of total sequences per sample is:", range(apply(com_treatments,1,sum)))

# remove species only appears in day 0
com_treatments <- com_treatments[ , colSums(com_treatments) != 0]
str(com_treatments)
com_treatments[1:5, 1:2]

group_info <- data.frame(row.names=rownames(com_treatments),t(as.data.frame(strsplit(rownames(com_treatments),"_"))))
head(group_info)

# bray curtis distance
dist <- vegdist(com_treatments, method="bray", binary=FALSE, diag=1) 

df <- data.frame(Soil = as.factor(group_info[,1]),
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]),
                 row.names = row.names(group_info))

df$Soil <- factor(df$Soil, levels = c("0", "70"), labels = c("0 years", "70 years"))
df$Water <- factor(df$Water, levels = c("sea", "strl"), labels = c("Sea water", "Sterile water"))
df$Day <- factor(df$Day, levels=c("2", "4", "8", "12", "16", "20"))
df$Frequency <- factor(df$Frequency, levels = c("f1", "f2", "f3", "f4"), 
                       labels = c("1x / week", "1x / 3 days", "1x / day", "2x / day"))
head(df); str(df)


# change as distance matirx
dist <- as.dist(as.matrix(dist))
str(dist)

# four way permanova
set.seed(123)
(result_whole <- adonis(dist ~ Soil+Day+Water+Frequency+Soil:Day+Soil:Water, data=df, permutation=9999))
result_whole
#write.csv(result_whole$aov.tab, "PERMANOVA_bray_whole_without_day0_updated.csv")

# ====================================================================================================
cat("change the year for analysis")
year <- "70"
subcom <- subset(com_treatments, grepl(paste('^', year, '_', sep=""), rownames(com_treatments)))
subcom <- subcom[ , colSums(subcom) != 0]
dim(subcom)
subcom[1:5, 1:2]

group_info <- data.frame(row.names=rownames(subcom),t(as.data.frame(strsplit(rownames(subcom),"_"))))

df <- data.frame(Soil = as.factor(group_info[,1]),
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]),
                 row.names = row.names(group_info))

df$Water <- factor(df$Water, levels = c("sea", "strl"), labels = c("Sea water", "Sterile water"))
df$Day <- factor(df$Day, levels=c("2", "4", "8", "12", "16", "20"))
df$Frequency <- factor(df$Frequency, levels = c("f1", "f2", "f3", "f4"), 
                       labels = c("1x / week", "1x / 3 days", "1x / day", "2x / day"))
head(df); str(df)

set.seed(123)
# bray curtis distance
dist <- vegdist(subcom, method="bray", binary=FALSE, diag=1) 
dist <- as.dist(as.matrix(dist))
(result_subcom <- adonis(dist ~ Day*Water*Frequency, data=df, permutation=9999))
result_subcom
write.csv(result_subcom$aov.tab, paste("PERMANOVA_bray_year", year, "without_day0.csv", sep = "_"))


