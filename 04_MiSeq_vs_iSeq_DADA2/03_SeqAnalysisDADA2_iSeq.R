####
#### DADA2 analysis of iSeq-style fastq files
#### 2020.7.21 Ushio
#### R 4.0.2
####

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
dir.create("03_SessionInfo")
output_folder <- "03_SeqAnalysisDADA2_iSeqOut"
dir.create(output_folder)

# Load library and functions
library(dada2); packageVersion("dada2") #1.16.0, 2020.6.25
library(ShortRead); packageVersion("ShortRead") #1.46.0, 2020.6.25
library(tidyverse); packageVersion("tidyverse") #1.3.0, 2020.6.25
source("00_Helpers/DADA2Helpers.R")

# Load sequence reads
path <- paste0(getwd(), "/02_DemultiplexedOut_iSeqStyle")
fnFs <- sort(list.files(path, pattern=".forward.fastq", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern=".reverse.fastq", full.names = T)) # Reverse read files

# Identify primers
FWD <- as.character(read.table("00_sample_info/F_primer.txt")[2,])
REV <- as.character(read.table("00_sample_info/R_primer.txt")[2,])
FWD_orients <- AllOrients(FWD)#; AllOrients(FWD_rc)
REV_orients <- AllOrients(REV)#; AllOrients(REV_rc)

# Pre-filtering to remove Ns (filter and compress)
fnFs_filtN <- str_sub(file.path(path, "01_filtN", basename(fnFs)), 1, -4) # Put N-filterd files in filtN/ subdirectory
fnRs_filtN <- str_sub(file.path(path, "01_filtN", basename(fnRs)), 1, -4)

filterAndTrim(fnFs, fnFs_filtN, fnRs, fnRs_filtN, maxN = 0, multithread = TRUE, compress = FALSE)
system(sprintf("pigz %s/*.fastq", file.path(path, "01_filtN")))
fnFs_filtN <- file.path(path, "01_filtN", basename(fnFs)) # Rename file path
fnRs_filtN <- file.path(path, "01_filtN", basename(fnRs)) # Rename file path

# Identify primers
seq_id <- 3
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs_filtN[[seq_id]]),
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs_filtN[[seq_id]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs_filtN[[seq_id]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs_filtN[[seq_id]]))

# Tell the path to cutadapt command
cutadapt <- "/usr/local/bin/cutadapt"
system2(cutadapt, args = "--version")

# Remove primers
path_cut <- file.path(path, "02_cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)
fnFs_cut <- file.path(path_cut, basename(fnFs))
fnRs_cut <- file.path(path_cut, basename(fnRs))

FWD_RC <- dada2:::rc(FWD)
REV_RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1_flags <- paste("-g", FWD, "-a", REV_RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2_flags <- paste("-G", REV, "-A", FWD_RC)

# Run Cutadapt
# Cutadapt commands may crash from R console. In that case, excute it from terminal.
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-j 72", # Multithred option
                             R1_flags, R2_flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs_cut[i], "-p", fnRs_cut[i], # output files
                             fnFs_filtN[i], fnRs_filtN[i])) # input files
}

# Sanity check
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs_cut[[seq_id]]), 
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs_cut[[seq_id]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs_cut[[seq_id]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs_cut[[seq_id]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path_cut, pattern = ".forward.fastq", full.names = TRUE))
cutRs <- sort(list.files(path_cut, pattern = ".reverse.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
sample_names <- unname(sapply(cutFs, function(x) strsplit(basename(x), "_")[[1]][2]))
head(sample_names)
# Visualize quality# Visualize quality
#plotQualityProfile(fnFs[1:3])
#plotQualityProfile(fnRs[1:3])

# Perform quality filtering
filtFs <- file.path(path, "03_filtered", paste0(sample_names, "_F_filt.fastq"))
filtRs <- file.path(path, "03_filtered", paste0(sample_names, "_R_filt.fastq"))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     # Output sequences of Claident already trimmed Ns and primers
                     #trimLeft = c(0,0),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=F, minLen = 20,
                     compress=FALSE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
system(sprintf("pigz %s/*.fastq", file.path(path, "03_filtered")))

# Exclude 0 seq samples, rename filtFs and filtRs
valid_sample_id <- out[,2]>0 & out[,1]>0
filtFs <- file.path(path, "03_filtered", paste0(sample_names[valid_sample_id], "_F_filt.fastq.gz"))
filtRs <- file.path(path, "03_filtered", paste0(sample_names[valid_sample_id], "_R_filt.fastq.gz"))

# Learn the error rates
min_nbases <- 5e+09
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)

# Visualize errors
ggsave("03_SeqAnalysisDADA2_iSeqOut/errF_iSeq.pdf", plotErrors(errF, nominalQ = T), width = 10, height = 10)
ggsave("03_SeqAnalysisDADA2_iSeqOut/errR_iSeq.pdf", plotErrors(errR, nominalQ = T), width = 10, height = 10)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_names[valid_sample_id]
names(derepRs) <- sample_names[valid_sample_id]

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dadaFs[[1]]

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, trimOverhang = TRUE, minOverlap = 5, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Cutting unexpected length sequences
seqtab2 <- seqtab
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtab)

# Track reads thourhg the pipeline
out2 <- out[valid_sample_id,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab_nochim),  rowSums(seqtab_nochim)/out2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample_names[out[,2] > 0 & out[,1] > 0]
head(track)

# Taxa output for claident tax assginment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">Taxa%05d", 1:length(seqs)), seqs)), ncol = 1)
write.table(seqs_out, paste0(output_folder, "/ASV_seqs.fa"), col.names = FALSE, row.names = FALSE, quote = FALSE)

# Save workspace
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("03_SessionInfo/SessionInfo_iSeq", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))
write.csv(seqtab_nochim, "seqtab_nochim_iseq.csv", row.names = F)

system("rm -r 02_DemultiplexedOut_iSeqStyle/01_filtN")
system("rm -r 02_DemultiplexedOut_iSeqStyle/02_cutadapt")
