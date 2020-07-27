####
#### Import as phyloseq objects (MiSeq and iSeq test)
#### 2020.7.27 Ushio
#### R 4.0.2
####

# Set random seeds (for reproduction)
output_folder <- "05_PhyloseqObjOut"
dir.create(output_folder)

# Load library and functions
library(phyloseq); packageVersion("phyloseq") #1.32.0, 2020.7.22

# Load workspace of MiSeq-version fastq
load("03_SeqAnalysisDADA2_MiSeqOut/03_SeqAnalysisDADA2_MiSeqOut.RData")

# Load sample sheet and tax information
sample_sheet_all <- read.csv("00_sample_info/SampleSheet.csv")
sample_sheet_all$fastq_type <- "MiSeq"
tax_claident <- read.delim("03_SeqAnalysisDADA2_MiSeqOut/ASV_merge_classigntax")

# Compile metadata
tax_claident$seq <- colnames(seqtab_nochim)
tax_claident$seqlen <- nchar(colnames(seqtab_nochim))
rownames(tax_claident) <- colnames(seqtab_nochim)
sample_sheet <- sample_sheet_all[match(rownames(seqtab_nochim), sample_sheet_all$Sample_Name2),]
rownames(sample_sheet) <- paste0("MiSeq_", sample_sheet$Sample_Name2)
rownames(seqtab_nochim) <- rownames(sample_sheet)

# Preparetion to import to phyloseq
dim(sample_sheet); dim(tax_claident); dim(seqtab_nochim)
all(rownames(sample_sheet) == rownames(seqtab_nochim)) # sample name check
all(rownames(tax_claident) == colnames(seqtab_nochim)) # taxa name check

# Import data to phyloseq
ps_miseq <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
                     sample_data(sample_sheet),
                     tax_table(as.matrix(tax_claident)))
saveRDS(ps_miseq, file = "05_PhyloseqObjOut/ps_miseq.obj")

# Delete all object
rm(list = ls())
# Load workspace of iSeq-version fastq
load("03_SeqAnalysisDADA2_iSeqOut/03_SeqAnalysisDADA2_iSeqOut.RData")

# Load sample sheet and tax information
sample_sheet_all <- read.csv("00_sample_info/SampleSheet.csv")
sample_sheet_all$fastq_type <- "iSeq"
tax_claident <- read.delim("03_SeqAnalysisDADA2_iSeqOut/ASV_merge_classigntax")

# Compile metadata
tax_claident$seq <- colnames(seqtab_nochim)
tax_claident$seqlen <- nchar(colnames(seqtab_nochim))
rownames(tax_claident) <- colnames(seqtab_nochim)
sample_sheet <- sample_sheet_all[match(rownames(seqtab_nochim), sample_sheet_all$Sample_Name2),]
rownames(sample_sheet) <- paste0("iSeq_", sample_sheet$Sample_Name2)
rownames(seqtab_nochim) <- rownames(sample_sheet)

# Preparetion to import to phyloseq
dim(sample_sheet); dim(tax_claident); dim(seqtab_nochim)
all(rownames(sample_sheet) == rownames(seqtab_nochim)) # sample name check
all(rownames(tax_claident) == colnames(seqtab_nochim)) # taxa name check

# Import data to phyloseq
ps_iseq <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
                     sample_data(sample_sheet),
                     tax_table(as.matrix(tax_claident)))
saveRDS(ps_iseq, file = "05_PhyloseqObjOut/ps_iseq.obj")


