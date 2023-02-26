setwd("~/Desktop/Lugol KMD figures") #wherever your data is

library(ggplot2)
library(dunn.test)
library(ggpubr)
library(dplyr)
library(agricolae)

rm(list=ls()) #handy command to clear your R environment

sa <- read.csv("lugol_growth.csv") #read in your file of raw data
#sa <- sa[sa$species == "OFAV",]
sa <- sa[sa$species == "APAL",]

sa$geno_trt <-  paste(sa$geno,sa$Trt) #making a new column that is genotype and treatment so I can do all my stats together

sa$geno <- factor(sa$geno, levels = c("AP16", "AP20")) 
#sa$geno <- factor(sa$geno, levels = c("OF3", "OF61")) 

#take a look at the data
p <- ggplot(data=sa, aes(x=geno_trt, y=delta)) +
  geom_bar(stat = "summary", fun = "mean")
p

#SA is approx normally distributed
ggdensity(sa$delta)
shapiro.test(sa$delta) #passes shapiro test, high p value
model <- lm(delta ~ Trt, data = sa)

plot(model, which = 2) #q-q plot
plot(model, which = 3) #variance
 
#modeling relationship of tle to treatment
summary(model)
anova(model)
#omnibus test shows main effect of treatment is significant if run on both species together (p value less than 0.05)
#treatment is NS if run on OFAV alone
#sig for apal alone

#within species
aov.model <- aov(delta ~ geno_trt, data = sa)
summary(aov.model)
tukey_stats <- TukeyHSD(aov.model)
tukey_groups <- HSD.test(aov.model, "geno_trt", group=TRUE)
tukey_groups <- tukey_groups$groups #helps with letters for plots
tukey_groups

sign <- aov.model[aov.model$p.value < 0.05]

tukey_stats <- as.data.frame(tukey_stats$geno_trt)
write.table(tukey_stats, file = "ofav_SA_lugol.txt", sep = "\t")

colors <- c( '#0045EE','#EE0000', '#FBD92E', '#E26060', '#BE0032', '#74A4D1', '#106EC8', '#70BD64','#1A9906')

boxplot <- ggplot(sa, aes(x=Trt, y=delta, fill=Trt)) +
  geom_boxplot() + #I personally think the boxplot looks better and is more informative
  #geom_bar(stat = "summary", fun = "mean") + #but you could also do a barplot 
  facet_grid(~geno) +
  theme_bw() +
  scale_fill_manual(values = c("red", "blue"))  +
  guides(fill=guide_legend(title="Treatment Type")) +
  labs(x="Treatment Type", y = "Change in Surface Area (mm^2): T0 to T4") +
  stat_summary(geom = 'text', label = c("b", "c", "a", "ab"), fun = max, vjust = -1, size = 3.5) 
boxplot

ggsave(filename="sa_boxplot_ofav_lugol.png", plot=boxplot, device="png",  width = 6, height = 5, dpi=500)
