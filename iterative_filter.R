# Place filtered files in filtered/ subdirectory
# filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(pathR, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,200),
#                      maxN=0, maxEE=c(2,2), truncQ=10, rm.phix=TRUE,
#                      compress=TRUE, multithread=TRUE)
# write.table(out, file = "filter.txt", sep="\t")


# Assessing read loss with iterative filtration steps. How much are we losing with each filtration/trimming step
filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(pathR, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     compress=TRUE, multithread=TRUE, matchIDs = TRUE)
write.table(out, file = "truncLen_out.txt", sep="\t")

filtFs_2 <- file.path(pathF, "filtered2", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_2 <- file.path(pathR, "filtered2", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs, filtFs_2, filtRs, filtRs_2,
                     maxN=0, compress=TRUE, multithread=TRUE, matchIDs = TRUE)
write.table(out, file = "maxN.txt", sep="\t")

filtFs_3 <- file.path(pathF, "filtered3", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_3 <- file.path(pathR, "filtered3", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs_2, filtFs_3, filtRs_2, filtRs_3,
                     maxEE=c(2,2), compress=TRUE, multithread=TRUE, matchIDs = TRUE)
write.table(out, file = "maxEE.txt", sep="\t")

filtFs_4 <- file.path(pathF, "filtered4", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_4 <- file.path(pathR, "filtered4", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs_3, filtFs_4, filtRs_3, filtRs_4, truncQ=10, compress=TRUE, multithread=TRUE, matchIDs = TRUE)
write.table(out, file = "truncQ_out.txt", sep="\t")

filtFs_5 <- file.path(pathF, "filtered5", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_5 <- file.path(pathR, "filtered5", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(filtFs_4, filtFs_5, filtRs_4, filtRs_5, rm.phix=TRUE, compress=TRUE, multithread=TRUE, matchIDs = TRUE)
write.table(out, file = "phiX_out.txt", sep="\t")
