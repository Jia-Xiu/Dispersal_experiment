# alpha diversity analysis for dispersal experiment
# Author: Xiu Jia
# Date: 24-01-2020 

rm(list=ls())

# load the directory
#directory = 'C:/Users/P278113/Dropbox'
directory = '~/Dropbox' 
subfolder = 'Dispersal/community_analysis'

setwd(paste(directory, subfolder, sep="/"))
getwd()

# load packages
library(vegan)
library(ggplot2)
library(RColorBrewer) 
display.brewer.all()
library(ggpubr) # combine figures & add significant values
library(reshape2)
library(doBy) # for stat
library(picante)
library(ggpmisc) # for stat_poly_eq regresssion labeling
library(ggpubr) # for stat_cor correlation labeling
library(plyr) # for ddply

mytheme<-theme_bw()+
  theme(text = element_text(size=12),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  


# load the rarefied otu table -----------------------------------------------------------------------
com <- read.csv("feature-table-rarified-nontax.csv", sep=";", header=1, row.names=1, check.names = FALSE)
colnames(com) <- gsub("\\-", "\\_", colnames(com))
com <- t(com)
dim(com)
com[1:5, 1:2]

shannon <- diversity(com)
simpson <- diversity(com, "simpson")
richness <- specnumber(com)
Chao <- estimateR(com)
Chao1 <- as.data.frame(t(Chao))
# Pielou's J = H'/ln(S) where H' is Shannon Weiner diversity and S is the total number of species in a sample,
pielous.evenness <- shannon/log(richness)

# PD ---
phylo <- read.tree("tree.nwk")
match.phylo.com <- match.phylo.data(phylo, t(com)); # species as rows, samples as columns for com table
str(match.phylo.com$phy)
#write.tree(match.phylo.com$phy, "pruned_tree.tre")

pd <- pd(t(match.phylo.com$data), match.phylo.com$phy,include.root = F)
str(pd)

group_info <- data.frame(row.names=rownames(com),t(as.data.frame(strsplit(rownames(com),"_"))))
head(group_info)

df <- data.frame(cbind(richness, Chao1, pd, shannon, simpson, pielous.evenness),
                 Soil = as.factor(group_info[,1]),
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]))

df$Soil <- factor(df$Soil, levels = c("0", "70"), labels = c("0 years", "70 years"))
df$Water <- factor(df$Water, levels = c("", "sea", "strl"), 
                   labels = c("0d", "Sea Water", "Sterile Water"))
#df$Day <- factor(df$Day, levels=c("0", "2", "4", "8", "12", "16", "20"))
df$Day <- as.numeric(as.character(df$Day))
df$Frequency <- factor(df$Frequency, levels = c("", "f1", "f2", "f3", "f4"), 
                       labels = c("0d", "1x / week", "1x / 3 days", "1x / day", "2x / day"))
df <- df[, -c(2, 4, 6, 8)]
head(df)


# compare richness of 0 and 70-year soil in seawater treatments
df1 <- df
df1$Day <- factor(df1$Day)
df1 <- df1[df1$Water != "Sterile Water", ]
str(df1)

(my_palette = c("#F0E442", "#0072B2")) # blue vs. yellow
(p4 <- ggplot(df1, aes(x = Day, y = PD, fill = Soil))+ 
    geom_boxplot(position=position_dodge(0.8), alpha = 0.8, outlier.size=-1) +
    geom_point(size = 2, shape = 16, alpha = 0.3, position = position_jitterdodge()) +
    #facet_grid(. ~ Water) +
    scale_fill_manual(values = my_palette) + 
    labs(x = "Day", y = "Phylogenetic diversity") +
    stat_compare_means(aes(group = Soil), label = "p.signif") +  
    #By default method = ???wilcox.test??? (non-parametric test)
    mytheme)

(p <- annotate_figure((ggarrange(p1, p2, p3, labels = c("A", "B", "C"), 
                                 common.legend = TRUE, legend = "right", ncol = 3, nrow = 1)), right = " "))

ggsave("richness_0_70yr_soil.jpg", width = 24, height = 6, units = "cm", device = "jpeg", p, scale = 1.5)


# correlation for all treatments except day 0 ---------------------------
dfs <- df[df$Water != "0d", ]
str(dfs)

# compare two types of soil
(my_palette = c("#56B4E9", "#009E73", "#D55E00", "#CC79A7")) 
(p0 <- ggplot(dfs, aes(x = Day, y = richness, fill = Frequency, color = Frequency, shape = Soil)) + 
    geom_smooth(aes(color = Frequency, linetype=Soil), method=lm, se=FALSE, size = .8) +
    geom_point(size = 3, alpha = 0.7) +
    facet_grid(. ~ Water) +
    labs(x = "Day", y = "Richness") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Frequency), method = "spearman", label.x = 8, label.y = 500, size = 3) + 
    mytheme)


custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352",
                "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pd <- position_dodge(0.1) 

(my_palette = c("#999999", "#000000")) # gray vs. black
(p1 <- ggplot(dfs, aes(x = Day, y = richness, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 3, alpha = 0.7) + #, shape = 21, colour = "#666666"
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Richness") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 500, size = 3) + 
    mytheme)

(p2 <- ggplot(dfs, aes(x = Day, y = PD, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 3, alpha = 0.7) +
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Phylogenetic diversity") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 42, size = 3) + 
    mytheme)

(p3 <- ggplot(dfs, aes(x = Day, y = shannon, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 3, alpha = 0.7) + 
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Shannon") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 5.5, size = 3) + 
    mytheme)

(p4 <- ggplot(dfs, aes(x = Day, y = pielous.evenness, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 3, alpha = 0.7) + 
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Pielou's evenness") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, label.y = 0.81, size = 3) + 
    mytheme)

(p <- annotate_figure((ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), 
                                 common.legend = TRUE, legend = "right", ncol = 1, nrow = 4)), right = " "))

ggsave("alpha_div_day_correlation_sea_strl_water_raw.pdf", width = 17, height = 32, units = "cm", p, scale = 1.5)




# regression scatter plot - sea water -----------------------------------

dfs <- df[df$Water == "Sea Water", ]
str(dfs)

pd <- position_dodge(0.1) 
(my_palette = c("#E7B800", "#FC4E07", "#C3D7A4", "#52854C"))


formula <- y ~ x

(p <- ggplot(dfs, aes(x = Day, y = shannon, fill = Frequency)) + #shape = Frequency, 
    geom_point(size = 2, alpha = 0.7, shape = 21, colour = "#666666") + 
    facet_grid(Soil ~ Frequency) +
    geom_smooth(aes(color = Frequency), method=lm, se=FALSE) +
    labs(x = "Day", y = "Shannon") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = 0.15,
                 formula = formula, parse = TRUE, size = 3) +
   mytheme)

ggsave("regression_shannon_day_seawater.jpg", width = 17, height = 8, units = "cm", p, scale = 1.5, device = "jpeg",dpi = 600)

# correlation ------
(my_palette = c("#E7B800", "#FC4E07", "#C3D7A4", "#52854C"))
(p <- ggplot(dfs, aes(x = Day, y = pielous.evenness, fill = Frequency)) + 
    geom_smooth(aes(color = Frequency), method=lm, se=FALSE) +
    geom_point(size = 2, alpha = 0.7, shape = 21, colour = "#666666") + 
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Pielou's evenness") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Frequency), method = "spearman", label.x = 8, label.y = 0.81, size = 3) + 
    # richness label.y = 500, shannon label.y = 5.5
    mytheme)

ggsave("correlation_pielous.evenness_day_seawater_raw.pdf", width = 17, height = 8, units = "cm", p, scale = 1.5)



# -----------------------------------------------------------------------------------------------------------
#  seawater & strl water

dfs <- df[df$Frequency == "1x / week" | df$Frequency == "2x / day", ]
str(dfs)

pd <- position_dodge(0.1) 
(my_palette = c("#00AFBB", "#E7B800"))


(p <- ggplot(dfs, aes(x = Day, y = shannon, fill = Water)) + 
    geom_smooth(aes(color = Water), method=lm, se=FALSE) +
    geom_point(size = 2, alpha = 0.7, shape = 21, colour = "#666666") + 
    facet_grid(Soil ~ Frequency) +
    labs(x = "Day", y = "Shannon") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 8, size = 3) + 
    mytheme)

ggsave("correlation_shannon_day_sea_strl_water_raw.pdf", width = 11, height = 8, units = "cm", p, scale = 1.5)

# -----------------------------------------------------------------------------------------------------------
#  seawater 
dfs <- df[df$Water == "Sea Water", ]
dfs <- dfs[, c(1, 7, 9, 11, 12, 13, 15)]
dfs <- melt(dfs, id.vars = c("Soil", "Frequency", "Day"))
dfs$variable <- factor(dfs$variable, levels = c("richness", "PD", "shannon", "pielous.evenness"), 
                       labels = c("Richness", "Phylogenetic diversity", "Shannon", "Pielou's evenness"))
dfs <- summaryBy(value ~ Soil + Frequency + Day + variable, data=dfs, FUN=dstats)
head(dfs)
str(dfs)

pd <- position_dodge(0.1)  ##pd <- position_dodge(width = 0.2, preserve = "single")
(my_palette = c("#E7B800", "#FC4E07", "#C3D7A4", "#52854C"))
(f <- ggplot(dfs, aes(x = Day, y = value.mean, group = Frequency, color = Frequency)) + 
    geom_errorbar(aes(group = Frequency, 
                      ymin = value.mean - value.se, 
                      ymax = value.mean + value.se), 
                  width=0.6, position=pd, alpha=0.9) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2, alpha=0.9) +
    scale_colour_manual(values = my_palette) +
    facet_grid(variable ~ Soil, scales = "free_y", switch="y") +
    labs(x="Days",y="", title=" ")+
    mytheme +
    theme(strip.placement = "outside"))

ggsave("alpha_div_four_seawater.jpg", width = 15, height = 17, units = "cm", f, scale = 1, device = "jpeg", dpi = 600)


# -----------------------------------------------------------------------------------------------------------------------
# scatter plot
(p <- ggplot(df, aes(x=Day, y=shannon, group=Frequency, color=Frequency))+ 
   geom_point()+
   geom_smooth(method = "loess", formula = y ~ x, se = FALSE)+
   scale_color_brewer(palette="Dark2")+
   facet_grid(Soil ~ Water, scales = "free_x", space = "free_x") +
   labs(x="Days",y="Shannon", title=" ")+
   mytheme)

ggsave("Shannon_lineplot.jpg", width = 17, height = 10, units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)

# box plot
(f <- ggplot(df, aes(x=Day, y=shannon, fill=Frequency))+
    geom_boxplot(position=position_dodge(0.8))+
    geom_point(size = 2, shape = 16, alpha = 0.8, position = position_jitterdodge()) +
    scale_fill_brewer(palette="Set3")+
    facet_grid(Soil ~ Water, scales = "free_x", space = "free_x") +
    labs(x="Days", y="Shannon", title=" ")+
    mytheme)

ggsave("Shannon_boxplot.jpg", width = 17, height = 10, units = "cm", device = "jpeg", f, scale = 1.5, dpi = 300)



# scatter plot
(p <- ggplot(df, aes(x=Day, y=pielous.evenness, group=Frequency, color=Frequency))+ 
    geom_point()+
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE)+
    scale_color_brewer(palette="Dark2")+
    facet_grid(Soil ~ Water, scales = "free_x", space = "free_x") +
    labs(x="Days",y="Pielou's evenness", title=" ")+
    mytheme)

ggsave("pielous.evenness_lineplot.jpg", width = 17, height = 10, units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)

# box plot
(f <- ggplot(df, aes(x=Day, y=pielous.evenness, fill=Frequency))+
    geom_point(size = 2, shape = 16, alpha = 0.8, position = position_jitterdodge()) +
    geom_boxplot(position=position_dodge(0.8))+
    scale_fill_brewer(palette="Set3")+
    facet_grid(Soil ~ Water, scales = "free_x", space = "free_x") +
    labs(x="Days", y="Pielou's evenness", title=" ")+
    mytheme)

ggsave("pielous.evenness_boxplot.jpg", width = 17, height = 10, units = "cm", device = "jpeg", f, scale = 1.5, dpi = 300)

















