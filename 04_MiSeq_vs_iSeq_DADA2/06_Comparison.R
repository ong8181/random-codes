####
#### Visualize results
#### 2020.7.27 Ushio
#### R 4.0.2
####

# Set random seeds (for reproduction)
output_folder <- "06_ComparisonOut"
dir.create(output_folder)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #1.3.0, 2020.6.25
library(phyloseq); packageVersion("phyloseq") #1.32.0, 2020.7.22
library(ggsci); packageVersion("ggsci") #2.9, 2020.7.22
library(fantaxtic); packageVersion("fantaxtic") #0.1.0, 7.22
library(cowplot); packageVersion("cowplot") #1.0.0, 2020.7.22
theme_set(theme_cowplot())

# Load workspace
ps_miseq <- readRDS("05_PhyloseqObjOut/ps_miseq.obj")
ps_iseq <- readRDS("05_PhyloseqObjOut/ps_iseq.obj")

# Combine phyloseq objects
ps <- merge_phyloseq(ps_miseq, ps_iseq)

# Visualize patterns (barplot)
ps2 <- get_top_taxa(ps, n = 50, relative = TRUE,
                    discard_other = FALSE, other_label = "Other")
b1 <- fantaxtic_bar(ps2, color_by = "phylum", label_by = "phylum", other_label = "Other", facet_by = "fastq_type")
b1$data$Sample <- factor(b1$data$Sample, levels = rownames(sample_data(ps2)))
b1$data$facet <- factor(b1$data$facet, levels = c("MiSeq", "iSeq"))

# Visualize patterns (scatterplot)
ps_aa_m0 <- psmelt(ps)
ps_aa_m1 <- ps_aa_m0[,c("Sample", "fastq_type", "query", "phylum", "order", "Abundance")]
ps_aa_m1$sample_x_query <- paste0(ps_aa_m1$Sample, "_", ps_aa_m1$query)
ps_aa_miseq <- ps_aa_m1[ps_aa_m1$fastq_type == "MiSeq",]
ps_aa_iseq <- ps_aa_m1[ps_aa_m1$fastq_type == "iSeq",]

# Converting to the relative abundance
ps_ra <- transform_sample_counts(ps, function(x) x / sum(x))
ps_ra_m0 <- psmelt(ps_ra)
ps_ra_m1 <- ps_ra_m0[,c("Sample", "fastq_type", "query", "phylum", "order", "Abundance")]
ps_ra_m1$sample_x_query <- paste0(ps_ra_m1$Sample, "_", ps_ra_m1$query)
ps_ra_miseq <- ps_ra_m1[ps_ra_m1$fastq_type == "MiSeq",]
ps_ra_iseq <- ps_ra_m1[ps_ra_m1$fastq_type == "iSeq",]

# Sort data
ps_aa_iseq <- ps_aa_iseq[match(substr(ps_aa_miseq$sample_x_query,2,20), ps_aa_iseq$sample_x_query),]
ps_ra_iseq <- ps_ra_iseq[match(substr(ps_aa_miseq$sample_x_query,2,20), ps_ra_iseq$sample_x_query),]
ps_ra_miseq <- ps_ra_miseq[match(ps_aa_miseq$sample_x_query, ps_ra_miseq$sample_x_query),]

# Check correspondence
head(ps_aa_iseq); head(ps_aa_miseq); head(ps_ra_iseq); head(ps_ra_miseq)
tail(ps_aa_iseq); tail(ps_aa_iseq); tail(ps_ra_iseq); tail(ps_ra_miseq)
all(ps_aa_iseq$query == ps_aa_miseq$query)
all(ps_aa_iseq$query == ps_ra_miseq$query)
all(ps_aa_iseq$query == ps_ra_iseq$query)

# Compile as a new dataframe
ps_m2 <- data.frame(sample_name = substr(ps_aa_miseq$Sample, 7, 10),
                    query = ps_aa_miseq$query,
                    phylum = ps_aa_miseq$phylum,
                    order = ps_aa_miseq$order,
                    MiSeq_abs_reads = ps_aa_miseq$Abundance,
                    iSeq_abs_reads = ps_aa_iseq$Abundance,
                    MiSeq_rlt_reads = ps_ra_miseq$Abundance,
                    iSeq_rlt_reads = ps_ra_iseq$Abundance)

# Visualize scatter plot
s1 <- ggplot(ps_m2, aes(x = MiSeq_abs_reads, y = iSeq_abs_reads, color = phylum)) +
  geom_point(alpha = 0.7) + scale_color_igv() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, col = "gray50") +
  xlab("Sequence reads (MiSeq fastq)") + ylab("Sequence reads (iSeq fastq)")

s2 <- ggplot(ps_m2, aes(x = MiSeq_rlt_reads, y = iSeq_rlt_reads, color = phylum)) +
  geom_point(alpha = 0.7) + scale_color_igv() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, col = "gray50") +
  xlab("Relative abundance (MiSeq fastq)") + ylab("Relative abundance (iSeq fastq)")


# Save figures
ggsave(sprintf("%s/MiSeq_vs_iSeq_barplot.pdf", output_folder), plot = b1, width = 14, height = 10)
ggsave(sprintf("%s/SequenceReads_scatterplot.pdf", output_folder), plot = s1, width = 8.5, height = 6)
ggsave(sprintf("%s/RelativeAbundance_scatterplot.pdf", output_folder), plot = s2, width = 8.5, height = 6)
