####
#### DADA2
#### Re-analyze O'Donnell et al. (2016)
#### 2022.05.12
#### R 4.1.2
####

# Load library and functions
library(dada2); packageVersion("dada2") # 1.22.0, 2022.5.12
library(ShortRead); packageVersion("ShortRead") # 1.52.0, 2022.5.12
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2022.5.12
source("funcs/F01_HelperFunctions.R")

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
dir.create("00_SessionInfo")

#------------------------------------------------#
# Quality filtering
#------------------------------------------------#
# Generate output folder
#wdir <- basename(rstudioapi::getSourceEditorContext()$path)
#(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
output_folder <- "02_DADA2Out"
dir.create(output_folder)

# Load sequence reads
path <- "/Users/ushio/Desktop/ODonnell_2016/Analysis/seqdata_seq_reversed"
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = T)) # Forward read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
(sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1))

# Visualize quality
#plotQualityProfile(fnFs[2])
#plotQualityProfile(fnRs[2])

# Performing filtering and trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     #truncLen=c(110, 110), # Truncate the end of reads
                     trimRight = c(15, 15),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = T,
                     compress = TRUE, multithread = TRUE) # On Windows set multithread=FALSE
#head(out)
out
#plotQualityProfile(filtFs[1:4])
#plotQualityProfile(filtRs[1:4])

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample_names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(filt_path, paste0(sample_names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample_names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
min_nbases <- 100 * sum(out[,2]) # Use a small number of bases to speed up the analysis
#min_nbases <- 1e+09 # Use the small number of bases to speed up the analysis
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)

# Visualize errors
#plotErrors(errF, nominalQ = T)
#plotErrors(errR, nominalQ = T)
ggsave(sprintf("%s/errF.pdf", output_folder), plotErrors(errF, nominalQ = T), width = 10, height = 10)
ggsave(sprintf("%s/errR.pdf", output_folder), plotErrors(errR, nominalQ = T), width = 10, height = 10)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_names[out[,2]>0 & out[,1]>0]
names(derepRs) <- sample_names[out[,2]>0 & out[,1]>0]

# Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")
#dadaFs[[1]]

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, trimOverhang = TRUE)
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, trimOverhang = TRUE, minOverlap = 20, maxMismatch = 1)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Cutting unexpected length sequences
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(60,160)]
#seqtab2 <- seqtab
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtab2)

# Track reads thourhg the pipeline
out2 <- out[out[,2]>0 & out[,1]>0,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab_nochim),  rowSums(seqtab_nochim)/out2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample_names[out[,2]>0 & out[,1]>0]
head(track)

# Taxa output for claident tax assignment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">ASV%05d", 1:length(seqs)), seqs)), ncol = 1)
seqtab_nochim_for_csv <- seqtab_nochim
colnames(seqtab_nochim_for_csv) <- sprintf("ASV%05d", 1:length(seqs))

# Save outputs
write.csv(seqs, paste0(output_folder, "/seq_only.csv"), row.names = colnames(seqtab_nochim))
write.table(seqs_out, paste0(output_folder, "/ASV_seqs.fa"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(seqtab_nochim_for_csv, paste0(output_folder, "/seqtab_nochim.csv"))
write.csv(track, paste0(output_folder, "/track.csv"))

# Save workspace
#load("02_DADA2_ITSOut/02_DADA2_ITSOut.RData")
rm(derepFs)
rm(derepRs)
rm(dadaFs)
rm(dadaRs)
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

