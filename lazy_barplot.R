library("phyloseq")
library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("reshape2")
library("indicspecies")
library("ggnetwork")
library("ape")
library("microbiome")
library("ggthemes")
library("cowplot")
library("ggsignif")
library("hrbrthemes")
library(gcookbook)
library(tidyverse)
rm(list=ls())
setwd("~/EPA_LK")

# load the rarefied, renamed, relative abundance data table with all species
load(file = "ps_rel.RData") #renamed, rarefied, relative abundance transformed ps object

sample_sums(ps_rel) 

ps_rel <- subset_samples(ps_rel, Site.Name == "Lindsays.Patch" & Objective == "cor" & Sample.Type == "Lesion")

ps_rel <- subset_samples(ps_rel, Site.Name == "Lindsays.Patch" & Objective == "water")

pseq <- aggregate_rare(ps_rel, level = "Genus", detection = 2/100, prevalence = 25/100)

pseq@sam_data$Site.Status <- factor(pseq@sam_data$Site.Status, levels = c("Vulnerable", "Epidemic", "Endemic")) 

p <- plot_composition(pseq,
                      taxonomic.level = "Genus",
                      otu.sort = "abundance",
                      sample.sort = "Site.Status",
                      x.label = "Site.Status") +
  scale_fill_manual(values=c("#dedede", "#006BFF", "#ffd700", "#277B00", "#ff0000", "#97DE00", "#00ffff", "#6495ed" , '#843BFF', "#BF5B9A", "#003E94", "#FF9200", "#60f542", '#F96D67', '#9067F9', '#F967C1', "#0072B2", "#CC3300","#660066", "#E69F00", "#56B4E9","#2168FF", "#006BFF", "#ffd700", "#277B00", "#ff0000", "#97DE00", "#00ffff")) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Site Status", y = "Relative abundance (%)",
       title = "MCAV Lesion Samples at Lindsay's Patch over Time") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
print(p)  
ggsave(filename="MCAV_lesion_lindsay.png", plot=p, device="png",  width = 7, height = 6, dpi=500)
