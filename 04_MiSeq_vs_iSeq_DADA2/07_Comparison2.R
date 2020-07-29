####
#### Visualize results
#### 2020.7.27 Ushio
#### R 4.0.2
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #1.3.0, 2020.6.25
library(phyloseq); packageVersion("phyloseq") #1.32.0, 2020.7.22
library(ggsci); packageVersion("ggsci") #2.9, 2020.7.22
library(cowplot); packageVersion("cowplot") #1.0.0, 2020.7.22
theme_set(theme_cowplot())

# Load workspace
load("06_ComparisonOut/06_ComparisonOut.RData")

# Create output directory
output_folder <- "07_Comparison2Out"
dir.create(output_folder)

# Check ASV distributions
tax_sums_df1 <- aggregate(ps_aa_m1$Abundance, by=list(ps_aa_m1$query, ps_aa_m1$fastq_type), sum)
tax_sums_df2 <- aggregate(ps_ra_m1$Abundance, by=list(ps_ra_m1$query, ps_ra_m1$fastq_type), sum)
colnames(tax_sums_df1) <- c("query", "fastq_type", "reads")
colnames(tax_sums_df2) <- c("query", "fastq_type", "relative_abundance")

h1 <- ggplot(tax_sums_df1, aes(x = reads, fill = fastq_type)) +
  geom_histogram(position = "identity", alpha = 0.5) + scale_x_log10() +
  scale_fill_startrek(name = "FASTQ type") +
  ggtitle("ASV reads distribution")

h2 <- ggplot(tax_sums_df2, aes(x = relative_abundance, fill = fastq_type)) +
  geom_histogram(position = "identity", alpha = 0.5) + scale_x_log10() +
  scale_fill_startrek(name = "FASTQ type") +
  ggtitle("ASV relative abundance distribution") + xlab("Relative abundance")

# Save figures
ggsave(sprintf("%s/MiSeq_vs_iSeq_histogram.pdf", output_folder),
       plot = plot_grid(h1, h2, ncol = 1, align = "hv"),
       width = 10, height = 10)

# Save workspace
#save(list = ls(all.names = TRUE),
#     file = sprintf("%s/%s.RData", output_folder, output_folder))


