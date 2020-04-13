# analyze  rare biosphere and common biosphere for the dispersal project
# Date: 31-12-2019
# Author: Jia Xiu


# load packages --------------------------------------------------------------------------------------------
rm(list=ls())

library(ggplot2)
library(RColorBrewer); display.brewer.all()
library(ggpubr) # combine figures
library(vegan)
library(ape) # pcoa 
library(reshape2) # melt
library(doBy) # se function
library(plyr) # rbind.fill
library(VennDiagram)
library(cowplot)
library(picante) # tree
library(GUniFrac) # for UniFrac
#library(rmarkdown)

# load directory --------------------------------------------------------------------------------------------
#directory = 'C:/Users/P278113/Dropbox'
directory = '~/Dropbox/' 
subfolder = 'Dispersal/community_analysis'

setwd(paste(directory, subfolder, sep="/"))
getwd()

mytheme<-theme_bw()+
  theme(text = element_text(size=12),
        plot.title = element_text(size=12, face="bold"),
        strip.background = element_blank(),
        strip.text=element_text(color="#666666", size = 12),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  

dstats<-function(x)(c(n=length(x), mean=mean(x), se=sd(x)/sqrt(length(x))))

# load the rare and common biospheres ----------------------------------------------------------------------------
cutoff <- readline("Which cutoff was used to define the rare biosphere\n
          if it is relative abundance based approach, type: eg. 0.001 \n
          if it is rank abundance curve based approach, type: rac")

years <- c('0', '70')

# a loop to calculte the fraction of the rare and common biosphere for each dataset
for (year in years) {
  rare <- read.csv(paste("feature_table", year, "year", cutoff, "cutoff_rare.csv", sep = "_"), sep=",", header=1, row.names=1, check.names=FALSE)
  rare <- t(rare)
  rare[is.na(rare)] <- 0
  rare <- rare[, colSums(rare != 0) > 0]
  str(rare)
  cat("\nthe number of samples is:", nrow(rare), "\nthe number of species/ASVs is:", ncol(rare),
      "\nthe range of sequence number among samples is:", range(rowSums(rare)))
  
  common <- read.csv(paste("feature_table", year, "year", cutoff, "cutoff_common.csv", sep = "_"), sep=",", header=1, row.names=1, check.names=FALSE) 
  common <- t(common)
  common[is.na(common)] <- 0
  common <- common[, colSums(common != 0) > 0]
  str(common)
  cat("\nthe number of samples is:", nrow(common), "\nthe number of species/ASVs is:", ncol(common),
      "\nthe range of sequence number among samples is:", range(rowSums(common)))
  
  # Venn digram for the overlap between rare and common --------------------------------------------------------
  venn.diagram(list(Rare = colnames(rare),
                    Common = colnames(common)), 
               filename = paste("venn_diagram", cutoff, "cutoff", year, "year.png", sep="_"), 
               height = 1200, 
               width = 1200,
               resolution =400, 
               imagetype="png", col="white", lwd=0.6, 
               fill=c("#66c2a5", "#fc8d62"),
               alpha = 0.70, 
               cex=0.9, 
               cat.cex=0.9,
               compression = 'lzw',
               units = 'px',
               lty = 'blank',
               fontface = "bold",
               fontfamily = "sans",
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.dist = c(-0.06, -0.035),
               cat.fontfamily = "sans")
}



# ============================================================================================================
# number of ASVs in the rare and common biosphere 
# ============================================================================================================

binary_or_not <- readline("do you want to sum up number of ASVs (binary? type: y, or yes\n
                          or number of indiciduals in each biosphere?, type: n, or no\n")

years <- c('0', '70') 

datalist <- list()

# a loop to calculte the fraction of the rare and common biosphere for each dataset
for (year in years) {
  
  cat("\n\nFor dataset:", year, "\n")
  
  rare <- read.csv(paste("feature_table", year, "year", cutoff, "cutoff_rare.csv", sep = "_"), sep=",", header=1, row.names=1, check.names=FALSE)
  rare <- t(rare)
  if (substring(binary_or_not, 1) == "y") {
    rare[rare > 0] <- 1
  }
  str(rare)
  cat("\nin rare biophere: \nthe number of samples is:", nrow(rare), "\nthe number of species/ASVs is:", ncol(rare),
      "\nthe range of  rowSums among samples is:", range(rowSums(rare)))
  
  common <- read.csv(paste("feature_table", year, "year", cutoff, "cutoff_common.csv", sep = "_"), sep=",", header=1, row.names=1, check.names=FALSE) 
  common <- t(common)
  if (substring(binary_or_not, 1) == "y") {
    common[common > 0] <- 1
  }
  str(common)
  cat("\nin common biophere: \nthe number of samples is:", nrow(common), "\nthe number of species/ASVs is:", ncol(common),
      "\nthe range of rowSums among samples is:", range(rowSums(common)))
  
  df <- data.frame(Common = rowSums(common), Rare = rowSums(rare))
  group_info <- data.frame(row.names=row.names(df), 
                           t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))
  
  df <- data.frame(Soil = as.factor(group_info[,1]),
                   Frequency = as.factor(group_info[,2]),
                   Water = as.factor(group_info[,3]),
                   Day = as.factor(group_info[,4]),
                   replicates= as.factor(group_info[,5]),
                   df, row.names = row.names(df))
  
  df$Soil <- factor(df$Soil, levels = c("0", "70"), labels = c("0 years", "70 years"))
  df$Frequency <- factor(df$Frequency, levels = c("", "f1", "f2", "f3", "f4"), 
                         labels = c("0d", "1x / week", "1x / 3 days", "1x / day", "2x / day"))
  df$Water <- factor(df$Water, levels = c("", "sea", "strl"), 
                     labels = c("0d", "Sea Water", "Sterile Water"))
  #df$Day <- factor(df$Day, levels = c("0", "2", "4", "8", "12", "16", "20"))
  
  df <- melt(df, id.vars = c("Soil", "Frequency", "Water", "Day", "replicates"))
  df$variable <- factor(df$variable, levels = c("Common", "Rare"), 
                        labels = c("Common biosphere", "Rare biosphere"))
  head(df)
  str(df)
  datalist[[year]] <- df
  
  dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
  data <- summaryBy(value ~ Soil+Frequency+Water+Day+variable, data=df, FUN=dstats)
  head(data)
  
  #datalist[[year]] <- data
}

str(datalist)
df <- do.call(rbind.data.frame, datalist)
head(df)
str(df)

# average of the richeness in each biosphere for two types of soil
(data <- summaryBy(value ~ Soil + variable, data=df, FUN=dstats))

# average of the richeness in each biosphere for two types of soil
df_wide <- dcast(df, Soil + Frequency + Water + Day + replicates ~ variable, value.var="value")
df_wide$sum <- df_wide$`Common biosphere` + df_wide$`Rare biosphere`
df_wide$asv_common <- round(df_wide$`Common biosphere`*100/df_wide$sum, 2)
df_wide$asv_rare <- round(df_wide$`Rare biosphere`*100/df_wide$sum, 2)
head(df_wide)
(data <- summaryBy(asv_common ~ Soil, data=df_wide, FUN=dstats))
(data <- summaryBy(asv_rare ~ Soil, data=df_wide, FUN=dstats))

# -----------------------------------------------------------
# correlation analysis
df <- df[df$Water != "0d", ]
df$Day <- as.numeric(as.character(df$Day))
head(df)
str(df)

pd <- position_dodge(0.1) 
# (my_palette = c("#00AFBB", "#E7B800")) # blue vs. yellow
(my_palette = c("#999999", "#000000")) # gray vs. black

df1 <- df[df$variable == "Common biosphere", ]
(p1 <- ggplot(df1, aes(x = Day, y = value, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 2, alpha = 0.6) + #, shape = 21, colour = "#666666"
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Number of ASVs", title = "Common biosphere") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 32, size = 3) + 
    mytheme)

df2 <- df[df$variable != "Common biosphere", ]
nrow(df2)
head(df2)
(p2 <- ggplot(df2, aes(x = Day, y = value, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 2, alpha = 0.6) + #, shape = 21, colour = "#666666"
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Number of ASVs", title = "Rare biosphere") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 180, size = 3) + 
    mytheme)

(p <- annotate_figure((ggarrange(p1, p2, labels = c("A", "B"), ncol = 1, nrow = 2, 
                                 align = "v", common.legend = TRUE, legend = "right")), right = " "))

ggsave(paste("rare_common_ASVs", cutoff, "cutoff_raw.pdf", sep = "_"), width = 17, height = 16, units = "cm", p, scale = 1.5)


# ============================================================================================================
# Relative abundance of rare and common biosphere 
# ============================================================================================================
years = c('0', '70') 

datalist <- list()

# a loop to calculte the fraction of the rare and common biosphere for each dataset
for (year in years) {
  
  cat("\n\nFor dataset:", year, "\n")
  rare <- read.csv(paste("feature_table", year, "year", cutoff, "cutoff_rare.csv", sep = "_"), sep=",", header=1, row.names=1, check.names=FALSE)
  rare <- t(rare)
  str(rare)
  cat("\nthe number of samples is:", nrow(rare), "\nthe number of species/ASVs is:", ncol(rare),
      "\nthe range of sequence number among samples is:", range(rowSums(rare)))
  
  rare <- read.csv(paste("feature_table", year, "year", cutoff, "cutoff_common.csv", sep = "_"), sep=",", header=1, row.names=1, check.names=FALSE)
  common <- t(common)
  str(common)
  cat("\nthe number of samples is:", nrow(common), "\nthe number of species/ASVs is:", ncol(common),
      "\nthe range of sequence number among samples is:", range(rowSums(common)))
  
  df <- data.frame(Common = rowSums(common), Rare = rowSums(rare))
  group_info <- data.frame(row.names=row.names(df), 
                           t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))
  
  df <- data.frame(Soil = as.factor(group_info[,1]),
                   Frequency = as.factor(group_info[,2]),
                   Water = as.factor(group_info[,3]),
                   Day = as.factor(group_info[,4]),
                   replicates= as.factor(group_info[,5]),
                   df, row.names = row.names(df))
  
  df$Sum <- df$Common+df$Rare
  df$Common <- df$Common*100/df$Sum
  df$Rare <- df$Rare*100/df$Sum
  
  df$Soil <- factor(df$Soil, levels = c("0", "70"), labels = c("0 years", "70 years"))
  df$Frequency <- factor(df$Frequency, levels = c("", "f1", "f2", "f3", "f4"), 
                         labels = c("0d", "1x / week", "1x / 3 days", "1x / day", "2x / day"))
  df$Water <- factor(df$Water, levels = c("", "sea", "strl"), 
                     labels = c("0d", "Sea Water", "Sterile Water"))
  df$Day <- factor(df$Day, levels = c("0", "2", "4", "8", "12", "16", "20"))
  
  
  df <- melt(df[, -8], id.vars = c("Soil", "Frequency", "Water", "Day", "replicates"))
  df$variable <- factor(df$variable, levels = c("Common", "Rare"), 
                        labels = c("Common biosphere", "Rare biosphere"))
  head(df)
  str(df)
  datalist[[year]] <- df
  
  dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
  data <- summaryBy(value ~ Soil+Frequency+Water+Day+variable, data=df, FUN=dstats)
  head(data)
  
  #datalist[[year]] <- data
}

str(datalist)
df <- do.call(rbind.data.frame, datalist)

# -----------------------------------------------------------------------------------
# correlation analysis
df <- df[df$Water != "0d", ]
df$Day <- as.numeric(as.character(df$Day))
head(df)
str(df)

pd <- position_dodge(0.1) 
(my_palette = c("#999999", "#000000")) 

df1 <- df[df$variable == "Common biosphere", ]
str(df1)
(p1 <- ggplot(df1, aes(x = Day, y = value, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 2, alpha = 0.6) + #, shape = 21, colour = "#666666"
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Total relative abundance (%)", title = "Common biosphere") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 50, size = 3) + 
    mytheme)

df2 <- df[df$variable != "Common biosphere", ]
(p2 <- ggplot(df2, aes(x = Day, y = value, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 2, alpha = 0.6) + #, shape = 21, colour = "#666666"
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Total relative abundance (%)", title = "Rare biosphere") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 45, size = 3) + 
    mytheme)

(p <- annotate_figure((ggarrange(p1, p2, labels = c("A", "B"), ncol = 1, nrow = 2, 
                                 align = "v", common.legend = TRUE, legend = "right")), right = " "))

ggsave("rare_common_ASVs_abundance_raw.pdf", width = 17, height = 16, units = "cm", p, scale = 1.5)



# stacked-bar plot  
df0 <- df[df$Soil == "0", ]
str(df0)

(f1 <- ggplot(df0, aes(x=Day, y=value.mean, fill=variable)) + 
    geom_bar(stat="identity", width=0.8, colour="black") +
    #scale_fill_manual(values=my_palette) + # my_palette = c(brewer.pal(6, "Set2")[c(1:2)]) 
    scale_fill_brewer(palette="Set2") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
    facet_grid(Water~Frequency, scales = "free_x", space = "free_x") +
    labs(x = "Day", y = "Relative abundance (%)", title = "0 year soil")+
    mytheme)

# bar plot with error bar
(p1 <- ggplot(df0, aes(x=Day, y=value.mean, fill=variable))+
    facet_grid(Water~Frequency, scales = "free_x", space = "free_x") +
    geom_bar(position = position_dodge(width = 0.7, preserve = "single"), stat="identity",
             colour="black", size=.3) + 
    geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se), colour="black", width=.3, size=.3,    # Thinner lines
                  position = position_dodge(width = 0.7, preserve = "single")) +
    scale_fill_brewer(palette="Set2") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
    labs(x = "Day", y = "Relative abundance (%)", title = "0 year soil")+
    mytheme)

df70 <- df[df$Soil == "70", ]
str(df70); head(df70)

# stacked-bar plot  
(f2 <- ggplot(df70, aes(x=Day, y=value.mean, fill=variable)) + 
    geom_bar(stat="identity", width=0.8, colour="black") +
    #scale_fill_manual(values=my_palette) + # my_palette = c(brewer.pal(6, "Set2")[c(1:2)]) 
    scale_fill_brewer(palette="Set2") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
    facet_grid(Water~Frequency, scales = "free_x", space = "free_x") +
    labs(x = "Day", y = "Relative abundance (%)", title = "70 year soil")+
    mytheme)

# bar plot with error bar
(p2 <- ggplot(df70, aes(x=Day, y=value.mean, fill=variable))+
    facet_grid(Water~Frequency, scales = "free_x", space = "free_x") +
    geom_bar(position = position_dodge(width = 0.7, preserve = "single"), stat="identity",
             colour="black", size=.3) + 
    geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se), colour="black", width=.3, size=.3,    # Thinner lines
                  position = position_dodge(width = 0.7, preserve = "single")) +
    scale_fill_brewer(palette="Set2") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
    labs(x = "Day", y = "Relative abundance (%)", title = "70 year soil")+
    mytheme)

(f <- ggarrange(f1, f2, labels = c("A", "B"), common.legend = TRUE, legend = "right", ncol = 1))
ggsave("Abundance_rare_common_stacked_barplot_raw.pdf", width = 10, height = 18, units = "cm", f, scale = 2)
ggsave("Abundance_rare_common_stacked_barplot_raw.jpg", width = 10, height = 18, units = "cm", f, scale = 2, dpi = 300)
(p <- ggarrange(p1, p2, labels = c("A", "B"), common.legend = TRUE, legend = "right", ncol = 1))
ggsave("Abundance_rare_common_barplot_raw.pdf", width = 12, height = 16, units = "cm", p, scale = 2)
ggsave("Abundance_rare_common_barplot_raw.jpg", width = 12, height = 16, units = "cm", p, scale = 2, dpi = 300)

# scatter plot
df_bio <- dcast(df[, 1:6], Year+Month+replicates+variable ~ Dataset)
df_bio$ratio <- df_bio$cDNA/df_bio$DNA
df_bio$Month <- factor(df_bio$Month, levels = c("M", "J", "S", "N"), labels=c("May", "Jun", "Sep", "Nov"))
head(df_bio)

(f <- ggplot(df_bio, aes(Year, ratio, shape = Month, fill = variable)) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1.2, alpha=0.7) +
    #geom_point(size=4, alpha=0.7)+
    scale_shape_manual(values=c(24, 22, 21, 23))+
    geom_hline(yintercept = c(1), linetype = "dashed") + 
    scale_fill_brewer(palette="Set2", guide=guide_legend(override.aes = list(shape=21)))+
    labs(x = "Stage of succession (Years)", y = "cDNA/DNA ratio", fill = "Biosphere", title = " ")+
    theme_bw()+
    theme(text = element_text(size=12),
          legend.box.background = element_rect(),
          legend.box.margin = margin(1, 1, 1, 1),
          legend.title = element_text(face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) )

ggsave("Abundance_rare_common_ratio.png", width = 8.5, height = 6, units = "cm", f, scale = 1.5, dpi = 300)
ggsave("Abundance_rare_common_ratio.pdf", width = 8.5, height = 6, units = "cm", f, scale = 1.5)
