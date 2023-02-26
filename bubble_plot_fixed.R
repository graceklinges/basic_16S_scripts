# libraries
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

load(file = "~/Mote_nutrient_experiment/data/ps_rel.RData")
setwd("~/Mote_nutrient_experiment/data/")

mapfile = "~/Mote_nutrient_experiment/data/Mapping-file-full-renamed.txt" #confirm correct mapping file associated
ps_rel = subset_samples(ps_rel, SampleID != "HP-7-2-T2")
ps_rel = subset_samples(ps_rel, SampleID != "MP-7-3-T2")
map = import_qiime_sample_data(mapfile)
sample_data(ps_rel) <- map

ps_genus <- tax_glom(ps_rel, "Genus")
ps_melt <- psmelt(ps_genus)

# summarize by Genus, sort, get top taxa
# Get families with mean relative abundance >0.01 across all samples 
genus <- ps_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
top <- genus[which(genus$Aver > 0.015),]
names <- top$Genus
names
ps_melt_top <- ps_melt[ps_melt$Genus %in% names,]
names

# plot
#note that ASV 2 in this dataset is the same as ASV1 in geno 7 only dataset- they are named by ranked abundance. BLASTed to confirm

ps_melt_top$Genus <- factor(ps_melt_top$Genus, levels = c("Unclassified_ASV_2", "Aquarickettsia", "Family_Helicobacteraceae", "[Caedibacter] taeniospiralis group",
                                                          "Spirochaeta 2"))

colorblind_pallette = c("#0072B2", "#58B8EC", "#CC3300", "#660066", "#39CC0A")

#basic format
nozero<- subset(ps_melt_top, Abundance > 0)
p <- ggplot(nozero, aes(x = Sample, y = reorder(OTU, Abundance))) +
  geom_point(aes(size=Abundance, fill = Genus), shape = 21) + #shape = 21 tells it to do colored circles
  facet_grid(. ~ Phenotype, scales = "free", space = "free", switch = "y") + #splits into resistant v susceptible categories
  theme_facet() +
  scale_fill_manual(values = colorblind_pallette, 
                    breaks = c("Aquarickettsia", "Unclassified_ASV_2", "Family_Helicobacteraceae", "[Caedibacter] taeniospiralis group", "Spirochaeta 2"),
                    labels = c("Aquarickettsia", "Unclassified ASV 1", "Family Helicobacteraceae", "Cysteiniphilum", "Spirochaeta")) +  
  guides(fill = guide_legend(override.aes = list(size=4))) + #makes legend circles bigger
  labs(fill = "Genus", size = "Relative \nAbundance")
p

ggsave(filename="~/Mote_nutrient_experiment/genotype_7/plots/geno_compare_bubble.svg", plot=p, device="svg",  width = 7, height = 4, dpi=500)



