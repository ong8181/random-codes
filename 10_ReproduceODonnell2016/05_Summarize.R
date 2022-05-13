####
#### Visualize patterns
#### R 4.1.2
#### 2022.5.13 Ushio
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2022.2.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2022.2.16
library(cowplot); packageVersion("cowplot") # 1.1.1, 2022.2.16
library(ggsci); packageVersion("ggsci") # 2.9
#library(cols4all); packageVersion("cols4all") # 0.2
theme_set(theme_cowplot())
source("funcs/F02_HelperFunctions.R") # Helper function for visualization

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# <-----------------------------------------------------> #
#  Load data
# <-----------------------------------------------------> #
# Load sample data
seqtab_data <- read.csv("03_OTUClusteringOut/otu_table.csv", row.names = 1)
sample_sheet <- read.csv("sampledata/SampleSheet_All.csv")
tax_sheet <- read.delim("04_TaxaAssignmentOut/merge_classigntax")
tax_seq <- read.csv("03_OTUClusteringOut/otu_only.csv", row.names = 1)

# Check structure
dim(seqtab_data)
dim(sample_sheet)
dim(tax_sheet); dim(tax_seq)


# <-----------------------------------------------------> #
#  Compile data
# <-----------------------------------------------------> #
# Check whether there is 0 sequences samples
(zero_sample <- sample_sheet$Sample_Name2[is.na(match(sample_sheet$Sample_Name2, rownames(seqtab_data)))])
if(length(zero_sample) > 0){
  # Generate dummy data frame
  seqtab_data_tmp <- matrix(0, ncol = ncol(seqtab_data), nrow = nrow(sample_sheet))
  rownames(seqtab_data_tmp) <- sample_sheet$Sample_Name2
  colnames(seqtab_data_tmp) <- colnames(seqtab_data)
  # Add object
  seqtab_data_tmp[match(rownames(seqtab_data), rownames(seqtab_data_tmp)),] <-
    as.matrix(seqtab_data)
  # Replace seq tables
  seqtab_data <- as.data.frame(seqtab_data_tmp)
}


# <-----------------------------------------------------> #
#  Import to phyloseq
# <-----------------------------------------------------> #
# Adjust row- and col-names
rownames(sample_sheet) <- rownames(seqtab_data) <- sample_sheet$Sample_Name2
tax_sheet$seq <- tax_seq[,2]
tax_sheet$seqlen <- nchar(tax_seq[,2])
colnames(seqtab_data) <- rownames(tax_sheet) <- tax_seq[,1]

# Import to phyloseq
ps_all <- phyloseq(otu_table(seqtab_data, taxa_are_rows = FALSE),
                   sample_data(sample_sheet),
                   tax_table(as.matrix(tax_sheet)))

# Visualize
target_taxa <- "phylum"
ps_sub <- taxa_name_summarize(ps_all, target_taxa, top_taxa_n = 10)
ps_m1 <- speedyseq::psmelt(ps_sub)
ps_m2 <- ps_m1 %>% group_by_at(c("Sample", target_taxa)) %>%
  summarize(sequence_reads = sum(Abundance)) #%>% filter(sequence_reads > 0)
ps_m3 <- ps_m1 %>% group_by(Sample, rep_tax) %>%
  summarize(sequence_reads = sum(Abundance)) #%>% filter(sequence_reads > 0)

# Set colors
#library(cols4all)
#cols4all::c4a_gui()
#c4a_palettes(type = "cat", series = "misc")

# Figures
f1 <- ggplot(ps_m2, aes_(x = as.name("Sample"), y = as.name("sequence_reads"),
                         fill = as.name(target_taxa))) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  scale_fill_igv() +
  #scale_fill_discrete_c4a_cat("cols25") +
  xlab(NULL) + ylab("Sequence reads")

f2 <- ggplot(ps_m3, aes(x = Sample, y = sequence_reads, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_igv(name = target_taxa) +
  #scale_fill_discrete_c4a_cat("cols25", name = target_taxa) +
  #coord_cartesian(ylim = c(0,10000)) +
  NULL

f3 <- plot_richness(ps_all, measures = "Observed") + xlab(NULL) + theme(axis.text.x = element_text(size = 6))

# Extract top taxa > 100 reads
#top_taxa <- taxa_names(ps_all)[taxa_sums(ps_all) > 100]
#ps_top <- prune_taxa(top_taxa, ps_all)
#f4 <- plot_richness(ps_top, measures = "Observed") + xlab(NULL) + theme(axis.text.x = element_text(size = 6))


# --------------------------------------------- #
# Save results
# --------------------------------------------- #
# Output figures
pdf(file = sprintf("%s/Summary.pdf", output_folder), width = 18, height = 6)
plot_grid(f2, f3, ncol = 2, rel_widths = c(1,0.7))
dev.off()

# Re-output data
write.csv(otu_table(ps_all), sprintf("%s/otu_table.csv", output_folder))
write.csv(sample_data(ps_all), sprintf("%s/sample_data.csv", output_folder))
write.csv(as.data.frame(tax_table(ps_all)), sprintf("%s/tax_table.csv", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))
