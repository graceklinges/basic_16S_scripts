library("dada2")
library("seqinr")
library("biomformat")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("cowplot")
library("phyloseq")

setwd("~/Mote_nutrient_experiment/data")

# functions
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

# phyloseq object output from Dada2
ps_full <- readRDS("ps_object.rds") #two different ways to load in data, change to your phyloseq object
load(file = "ps_full.RData") 
ps_full <- ps #this just saves your full ps object so you can compare reads lost later on

# import metadata and merge into phyloseq object
mapfile = "Mapping-file-full-renamed.txt" #change to your mapping file
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

# summary of data
ps
summary(sample_data(ps))

ntaxa(ps)
nsamples(ps)
rank_names(ps)
sample_names(ps)[1:5]
sample_variables(ps)

# remove mitochondria and chloroplasts, is.na important becuase if not included
# this command will also remove all Family = NA or Order = NA
ps_with_mito = subset_taxa(ps, (Order!="Chloroplast") | is.na(Order))
ps_no_mito = subset_taxa(ps_with_mito, (Family!="Mitochondria") | is.na(Family))
ps_no_Eukaryota = subset_taxa(ps_no_mito, (Kingdom!="Eukaryota") | is.na(Kingdom)) #lost 0 taxa

ntaxa(ps) - ntaxa(ps_no_Eukaryota)
ps = ps_no_Eukaryota

summary(taxa_sums(ps)) # 5 was first quantile, change 5 down below to whatever your first quantile is
# Filter on prevalence or total counts
pst = fast_melt(ps)
prevdt = pst[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = taxaID]
keepTaxa = prevdt[(Prevalence >=0 & TotalCounts >5), taxaID]
ps_pruned = prune_taxa(keepTaxa,ps)
ps_pruned
sum(sample_sums(ps_pruned)) #5,149,563
summary(sample_sums(ps_pruned))
summary(taxa_sums(ps_pruned))
ps <- ps_pruned
save(ps, file = "ps_pruned.RData")

clr <- microbiome::transform(ps, 'clr') #centered log ratio transformation, I use for beta diversity

min(sample_sums(ps))

#plot otus by sequencing depth
observed <- estimate_richness(ps, measures = c('Observed'))
explore.df <- cbind(observed, sample_sums(ps), sample_data(ps)$Region) #change Region to some variable you have
colnames(explore.df) <- c('Observed', 'Sample_Sums', 'Region')
observed_mean <- mean(explore.df$Observed)
sample_sum_mean <- mean(explore.df$Sample_Sums)
observed_plot <- ggplot(data = explore.df, aes(x = Sample_Sums, y = Observed, color = Region)) + 
  geom_point() +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95, 
              inherit.aes = F, mapping = aes(Sample_Sums, Observed),
              data = explore.df) +
  ylab("Observed OTUs") +
  scale_colour_colorblind()
ggsave(filename="observed_by_sums.png", plot=observed_plot, device="png", dpi=500)

# Let's use ampvis2 again so we can easily make a rarefaction curve

# Need to convert from phyloseq to ampvis
av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(ps)@.Data)),
                           t(phyloseq::otu_table(ps)@.Data),
                           phyloseq::tax_table(ps)@.Data,
                           check.names = F
)

tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
av2_tax<-as.data.frame(tax)

#NOTE THAT IF YOU HAVE A VARIABLE CALLED SPECIES ITS GONNA CAUSE ISSUES
#Extract metadata from the phyloseq object:
av2_metadata <- data.frame(phyloseq::sample_data(ps), 
                           check.names = F
)

av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)

#Load the data with amp_load:
av2_obj <- amp_load(av2_otutable, metadata = av2_metadata,taxonomy = av2_tax)

# RARE CURVE
rare_plot_amp <- amp_rarecurve(data = av2_obj, color_by = "Region")
rare_curve_plot <- rare_plot_amp + ylab('Observed ASVs (count)') + 
  geom_vline(xintercept=min(sample_sums(ps)), linetype='dashed') +
  scale_colour_colorblind() +
  xlim(c(0, 35000))
plot(rare_plot_amp)
rare_curve_plot
ggsave(filename="unrarefied.png", plot=rare_plot_amp, device="png", dpi=500)
ggsave(filename="rarefy.png", plot=rare_curve_plot, device="png", dpi=500)


#set sample size (5098) to be your minimum sequencing depth from the above command
ps_rarefied <- rarefy_even_depth(ps, sample.size = 5098, rngseed = 999) 
ps_rarefied 
sum(sample_sums(ps_rarefied)) 
sample_sums(ps_rarefied)
#relative abundance transform
ps_rel <- transform(ps_rarefied, "compositional")

# save ps objects
save(ps, file = "ps_rename_g50.RData") #unrarefied, pruned
save(ps_rel, file = "ps_rel_g50.RData") #rel abundance, rarefied, pruned
save(clr, file = "ps_clr_g50.RData") #clr transformed, pruned
