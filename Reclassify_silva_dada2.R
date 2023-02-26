library(phyloseq)
library(seqinr)
library(dada2)

#starting with rep seq file from MR DNA DADA2 pipeline
#as qza format from qiime
#in bash, qiime tools export --input-path --output-path
#exported rep seq file as fasta
#katie sent me feature table as text file
#converted biom txt file to json formatted biom using biom_convert
#biom convert -i 021920EM515F-table-dada2-feature-table.biom.txt -o json.biom --table-type="OTU table" --to-json
ps <- import_biom("json.biom", refseqfilename="dna-sequences.fasta")
seqtab<-t(as(otu_table(ps),"matrix"))

otu_names <- colnames(seqtab)

#imported req seps fasta file to excel, copied col 1 to col 2, filtered col 1 to 
#choose one, begins with >, copied to new sheet --> tab separated fasta, export as .txt
sequences <- as.matrix(read.table("dna-sequences.txt", sep="\t", header=FALSE))
colnames(seqtab) <- sequences[,2] #changes column names of seqtab to exact ASV instead of qiime ID based on rep seqs file

# Assign taxonomy based on silva reference database at genus level, you must have the appropriate Silva database downloaded
tax_silva <- assignTaxonomy(seqtab, "/Volumes/Grace\ External\ 2/Mote_Nutrient_Exp/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Assign taxonomy based on silva reference database at species (100%) level
silva_sp <- addSpecies(tax_silva, "/Volumes/Grace\ External\ 2/Mote_Nutrient_Exp/silva_species_assignment_v132.fa.gz", returnMultiple=TRUE)

# Export sequence table with genus and species assignments as phyloseq objects
ps_sctld <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax_silva))
ps_sctld_sp <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(silva_sp))

#remove chloroplasts, mitochondria, eukaryotes
ps_with_mito = subset_taxa(ps_sctld, (Order!="Chloroplast") | is.na(Order)) #removed 223 taxa
ps_no_mito = subset_taxa(ps_with_mito, (Family!="Mitochondria") | is.na(Family)) #removed 104 taxa
ps_no_Eukaryota = subset_taxa(ps_no_mito, (Kingdom!="Eukaryota") | is.na(Kingdom)) #removed 46 taxa
ps_sctld <- ps_no_Eukaryota

ps_with_mito = subset_taxa(ps_sctld_sp, (Order!="Chloroplast") | is.na(Order)) #removed 223 taxa
ps_no_mito = subset_taxa(ps_with_mito, (Family!="Mitochondria") | is.na(Family)) #removed 104 taxa
ps_no_Eukaryota = subset_taxa(ps_no_mito, (Kingdom!="Eukaryota") | is.na(Kingdom)) #removed 46 taxa
ps_sctld_sp <- ps_no_Eukaryota

# Save as RDS objects
saveRDS(ps_sctld, file = "ps_sctld.rds")
saveRDS(ps_sctld_sp, file = "ps_sctld_sp.rds")
saveRDS(seqtab, file = "seqtab.rds")

#unique sequences
uniqueSeqs <- as.list(colnames(seqtab))
write.fasta(uniqueSeqs, uniqueSeqs, "uniqueSeqs.fasta")

#taxonomy text file
tax <-as(tax_table(ps_sctld_sp),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax <-as.data.frame(tax)
tax $taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "taxonomy_sctld.txt", quote=FALSE, col.names=FALSE, sep="\t")

#confidence calls, merged these with taxonomy in excel for ease of reading
write.table(bootstraps, "taxonomy_confidence.txt", quote=FALSE, col.names=TRUE, sep="\t")

#write to csv
OTU1 = as(otu_table(ps_sctld_sp), "matrix")
OTU1 <- t(OTU1)
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
write.table(OTUdf, "feature-table.txt", quote=FALSE, col.names=TRUE, sep="\t")

