# beta diversity (turnover) analysis for dispersal experiment
# Author: Xiu Jia
# Date: 29-01-2020

rm(list=ls())

# Load libraries
library(vegan) # for multivariable analysis
library(ape) # for pcoa 
library(ggplot2) # for graphing
library(ggpubr) # combine figures
library(RColorBrewer) # for color bar
library(plyr) # for ddply
display.brewer.all()

# change directory

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


# load the rarefied otu table -----------------------------------------------------------------------
com <- read.csv("feature-table-rarified-nontax.csv", sep=";", header=1, row.names=1, check.names = FALSE)
colnames(com) <- gsub("\\-", "\\_", colnames(com))
com <- t(com)
dim(com)
com[1:5, 1:2]
cat("the range of total number of sequences of each species is:", range(apply(com,2,sum)))
cat("the range of total sequences per sample is:", range(apply(com,1,sum)))

# ========================================================================================================
# Temporal changes for all samples 
# ========================================================================================================

# bray curtis distance
bray <- vegdist(com, method="bray", diag=1) 
df <- data.frame(as.table(as.matrix(bray)))[lower.tri(bray, diag = FALSE), ]
row.names(df) <- paste(df$Var1, df$Var2, sep = "_")

group <- data.frame(row.names=rownames(df), t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))
head(group)

df <- transform(merge(df, group[, -c(5, 10)], by="row.names"), row.names=Row.names, Row.names=NULL, check.names=FALSE)
df <- df[which(df$X1==df$X6),]
df <- df[which(df$X2==df$X7),]
df <- df[which(df$X3==df$X8),]

df$Soil <- factor(df$X1, levels = c("0", "70"), labels = c("0 years", "70 years"))
df$Frequency <- factor(df$X2, levels = c("f1", "f2", "f3", "f4"), 
                       labels = c("1x / week", "1x / 3 days", "1x / day", "2x / day"))
df$Water <- factor(df$X3, levels = c("sea", "strl"), labels = c("Sea Water", "Sterile Water"))
df$Day1 <- as.numeric(as.character(df$X4))
df$Day2 <- as.numeric(as.character(df$X9))

df <- df[, c(3, 12:16)]
df <- na.omit(df)
df$daydif <- abs(df$Day1-df$Day2)
head(df)
str(df)

# correlation
pd <- position_dodge(0.1) 
# (my_palette = c("#00AFBB", "#E7B800")) # blue vs. yellow
(my_palette = c("#999999", "#000000")) # gray vs. black
(p <- ggplot(df, aes(x = daydif, y = Freq, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 2, alpha = 0.6) + #, shape = 21, colour = "#666666"
    facet_grid(Soil ~ Frequency) +
    labs(x = expression(Delta~"Time"~"(days)"),    
         y = "Community dissimilarity between samples (Bray-Curtis)") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    stat_cor(aes(color = Water), method = "spearman", label.x = 4, label.y = 0.2, size = 3) + 
    mytheme)

ggsave("Temporal_turnover_correlation_raw.pdf", width = 17, height = 8, units = "cm", p, scale = 1.5)


# regression ---------------------------------------------------------------------------------
df$y <- df$Freq
df$x <- df$daydif

lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

eq <- ddply(df, .(Soil, Frequency, Water), lm_eqn)

(my_palette = c("#00AFBB", "#E7B800"))

(p <- ggplot(df, aes(daydif, Freq, fill = Water, color = Water, shape = Water)) + 
    geom_smooth(aes(color = Water, linetype=Water), method=lm, se=FALSE, size = .8) +
    geom_point(size = 2, alpha = 0.6) + #, shape = 21, colour = "#666666"
    facet_grid(Soil ~ Frequency) +
    labs(x = expression(Delta~"Time"~"(days)"),    
         y = "Community dissimilarity between samples (Bray-Curtis)") +
    scale_colour_manual(values = my_palette) +
    scale_fill_manual(values = my_palette, guide=guide_legend(override.aes = list(shape=21))) +
    scale_shape_manual(values=c(21, 1))+ 
    scale_x_continuous(breaks=c(2, 4, 8, 12, 16, 20)) +
    geom_text(data=eq, aes(color = Water, x = 8, y = 0.2, label=V1), size = 3, parse = TRUE, inherit.aes=FALSE) +
    mytheme)

ggsave("Temporal_turnover_regression_raw.pdf", width = 17, height = 8, units = "cm", p, scale = 1.5)


