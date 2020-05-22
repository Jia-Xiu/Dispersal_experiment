# analyze results from functional annotation by PICRUST2  ??? heatmap
# Author: Jia Xiu
# Date: 2020-05-08

rm(list=ls())

# change directory   
getwd()

library(vegan) 
library(reshape2)
library(doBy)
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


# load ko reference
ko <- keggList("ko")
names(ko) <- gsub("ko:", "", names(ko))  


# All
meg_0 <- read.table("pred_metagenome_unstrat_0_year.tsv", sep = '\t', header = 1, row.names = 1, check.names = FALSE)
meg_70 <- read.table("pred_metagenome_unstrat_70_year.tsv", sep = '\t', header = 1, row.names = 1, check.names = FALSE)
meg <- transform(merge(meg_0, meg_70, all = TRUE, by="row.names"), row.names=Row.names, Row.names=NULL)
meg[is.na(meg)] <- 0
colnames(meg) <- gsub("^X", "", colnames(meg))
meg <- t(meg)
meg <- as.data.frame(meg)
warning("remove day 0")
meg <- subset(meg, !grepl('___', row.names(meg)))
meg <- meg[, !colSums(meg)==0]
range(colSums(meg))
dim(meg)
meg[1:5, 1:3]

# keep ko functions in my dataset
ko_clean <- ko[names(ko) %in% colnames(meg)]
str(ko_clean)
class(ko_clean)
cat("These ko were not found in the dowload dataset:\n", colnames(meg)[!colnames(meg) %in% names(ko)])

# Remove all before and up to "; "
ko_clean <- gsub(".*; ", "", ko_clean)
# Remove square brackets and text within
ko_clean <- gsub("\\[.*\\]", "", ko_clean)
# Remove leading/trailing whitespace
ko_clean <- trimws(ko_clean, which = c("both"))
head(ko_clean)
# write.csv(as.data.frame(ko_clean), "ko_list_dispersal.csv")

data <- read.csv("kegg_orthology_list_dispersal.csv")
data$groups <- ifelse(grepl('*methyltransferase*', data$ko_clean), "methyltransferase", as.character(data$groups))
data$groups <- ifelse(grepl('*MFS*', data$ko_clean), "MFS_transporter", as.character(data$groups))
data$groups <- ifelse(grepl('*aminotransferase*', data$ko_clean), "aminotransferase", as.character(data$groups))
data$groups <- ifelse(grepl('*antibiotic*', data$ko_clean), "ARGs", as.character(data$groups))
levels(factor(data$groups))
data[55:70,]
# write.csv(data, "kegg_orthology_list_dispersal_ARGs.csv")

# heatmap ---------------------------------------------------------------------------------------------------------
functional_group <- "ARGs" #"chemotaxis_flagellar_assembly"
df_nr <- read.csv(paste("ko", functional_group, "nr.csv", sep = "_"))
# nr_chemotaxis <- c("K03406", "K03407", "K03408", "K05874", "K05875", "K03412", "K03413", "K00575")
ko_sub <- ko_clean[names(ko_clean) %in% df_nr$KEGG_orthology]
ko_sub <- as.data.frame(ko_sub)
ko_sub$Var2 <- row.names(ko_sub)
head(ko_sub)

df <- meg[, colnames(meg) %in% row.names(ko_sub)]
dim(df)
df[1:5, 1:4]

# Z-score transformtion for all functions
cols <- colnames(df)
df[cols] <- scale(df[cols])

# heatmap
df1 <- melt(as.matrix(df))

group_info <- data.frame(t(as.data.frame(strsplit(as.character(df1$Var1), "_"))))
head(group_info)

df1 <- data.frame(df1,
                 Soil = as.factor(group_info[,1]),
                 Frequency = as.factor(group_info[,2]),
                 Water = as.factor(group_info[,3]),
                 Day = as.factor(group_info[,4]),
                 replicates = as.factor(group_info[,5]))

dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se = sd(x)/ sqrt(length(x))))
df1 <- summaryBy(value ~ Var2 + Soil + Water + Frequency + Day, data=df1, FUN=dstats)

df1 <- merge(df1, ko_sub, by="Var2", all=TRUE)
df1 <- df1[with(df1, order(Day, Frequency, Water, Soil)), ]
df1$Treaments <- factor(paste(df1$Water, df1$Frequency, df1$Day, sep = "_"), 
                        levels = c("sea_f1_2", "sea_f1_4",  "sea_f1_8", "sea_f1_12", "sea_f1_16", "sea_f1_20", 
                                   "sea_f2_2", "sea_f2_4",  "sea_f2_8", "sea_f2_12", "sea_f2_16", "sea_f2_20",
                                   "sea_f3_2", "sea_f3_4",  "sea_f3_8", "sea_f3_12", "sea_f3_16", "sea_f3_20", 
                                   "sea_f4_2", "sea_f4_4",  "sea_f4_8", "sea_f4_12", "sea_f4_16", "sea_f4_20",
                                   "strl_f1_2", "strl_f1_4", "strl_f1_8", "strl_f1_12", "strl_f1_16", "strl_f1_20",
                                   "strl_f4_2", "strl_f4_4", "strl_f4_8", "strl_f4_12", "strl_f4_16", "strl_f4_20"))
df1$Soil <- factor(df1$Soil, levels = c("0", "70"), labels = c("0 years", "70 years"))
head(df1)
str(df1)

base_size = 9
(p <- ggplot(df1, aes(Treaments, ko_sub)) + 
    geom_tile(aes(fill = value.mean), colour = "lightgray") + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", space = "Lab", 
                         name="z-score(abundance)", midpoint = 0, limit = c(-3, 10)) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    facet_grid(~ Soil, scales = "free_x", space = "free_x") +
    labs(x = "Treatment", y = "", title = "") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
          legend.position = "top", 
          legend.justification='left',
          legend.direction='horizontal'))

ggsave(paste("Heatmap_", functional_group, ".jpg", sep = ""), width = 18, height = 8, 
       units = "cm", device = "jpeg", p, scale = 1.5, dpi = 300)

