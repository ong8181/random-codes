####
#### OTU Clustering
#### Re-analyze O'Donnell et al. (2016)
#### 2022.05.12
#### R 4.1.2
####

# Load DADA2 results
load("02_DADA2Out/02_DADA2Out.RData")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #1.3.1, 2022.2.15
library(phyloseq); packageVersion("phyloseq") #1.38.0, 2022.2.15
library(ShortRead); packageVersion("ShortRead") #1.52.0, 2022.2.15
library(dada2); packageVersion("dada2") #1.22.0, 2022.2.15
# Packages that are required but not loaded:
# library(DECIPHER)
# library(Biostrings)

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# <-----------------------------------------------------> #
#                       Load data
# <-----------------------------------------------------> #
# Load sample data
sample_sheet <- read.csv("sampledata/SampleSheet_All.csv")
#seqs
dim(seqtab_nochim); sum(seqtab_nochim)


# <-----------------------------------------------------> #
#                  Clustering by DECIPHER
# <-----------------------------------------------------> #
# Preparation
n_cores <- parallel::detectCores() # set to number of cpus/processors to use for the clustering
dna <- Biostrings::DNAStringSet(seqs) # "seqs" is an object from DADA2 output
## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = n_cores)
aln_dist <- DECIPHER::DistanceMatrix(aln, processors = n_cores)
# use `cutoff = 0.03` for a 97% OTU 
clusters <- DECIPHER::IdClusters(aln_dist, method = "complete",
                                 processors = n_cores, type = "clusters",
                                 cutoff = 0.03, showPlot = F)
colnames(clusters) <- "cluster_id"
clusters$cluster_id <- factor(clusters$cluster_id, levels = unique(clusters$cluster_id))
length(unique(clusters$cluster_id))
#clusters <- clusters %>% add_column(sequence = seqs)
## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
merged_seqtab <- seqtab_nochim %>% t %>%
  rowsum(clusters$cluster_id) %>% t
dim(merged_seqtab)


# <-----------------------------------------------------> #
#             Make OTU-based phyloseq objects
# <-----------------------------------------------------> #
# Import DADA2 ASV output to phyloseq
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE))
# Merge taxa in "ps" using cluster information
ps_otu <- speedyseq::merge_taxa_vec(ps, group = clusters$cluster_id)
otu_seqs <- colnames(otu_table(ps_otu))

# Quick taxa assignment
#taxa <- assignTaxonomy(otu_table(ps_otu), "~/DADA2_DB/rdp_train_set_18.fa.gz")
#taxa <- assignSpecies(otu_table(ps_otu), "~/DADA2_DB/rdp_species_assignment_18.fa.gz")
#rdp_species_assignment_18.fa
otu_only <- data.frame(taxa_id = sprintf("OTU%05d", 1:length(otu_seqs)),
                       seq = otu_seqs)
write.csv(otu_only, sprintf("%s/otu_only.csv", output_folder), row.names = T)

# Check corresponence between ASV sequences and OTU representative sequences
clusters$seq_sum <- colSums(seqtab_nochim)
clusters$asv_seq <- seqs
asv_1st <- clusters %>% group_by(cluster_id) %>% summarize(otu_seqs = asv_seq[[1]])
clusters$otu_seq <- unlist(asv_1st[match(clusters$cluster_id, asv_1st$cluster_id), "otu_seqs"])
write.csv(clusters, sprintf("%s/cluster_summary.csv", output_folder), row.names = F)

# Save OTU table
otu_out <- as.matrix(c(rbind(sprintf(">OTU%05d", 1:length(otu_seqs)), otu_seqs)), ncol=1)
write.table(otu_out, sprintf("%s/OTU_seqs.fa", output_folder), col.names = FALSE, row.names = FALSE, quote = FALSE)

otu_mat <- as.matrix(otu_table(ps_otu)@.Data)
colnames(otu_mat) <- sprintf("OTU%05d", 1:length(otu_seqs))
write.csv(otu_mat, sprintf("%s/otu_table.csv", output_folder), row.names = T)


# <-----------------------------------------------------> #
#                     Save workspace
# <-----------------------------------------------------> #
# Save workspace
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))


