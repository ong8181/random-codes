####
#### Statistical analyses after DADA2 processing
#### 2021.4.28 Ushio
#### R 4.0.3
####

# Create output directory
dir.create("00_SessionInfo")
dir.create("DADA2PhyloseqOut")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.0
library(phyloseq); packageVersion("phyloseq") # 1.34.0
library(cowplot); packageVersion("cowplot") # 1.1.1
library(ggsci); packageVersion("ggsci") # 2.9
theme_set(theme_cowplot())

# ----------------------------------------------- #
#                 Load demo data
# ----------------------------------------------- #
sample_sheet <- read.csv("data_subset/sample_sheet.csv", row.names = 1)
asv_sheet <- read.csv("data_subset/asv_sheet.csv", row.names = 1)
tax_sheet <- read.csv("data_subset/tax_sheet.csv", row.names = 1)


# ----------------------------------------------- #
#         Check colnames and rownames
# ----------------------------------------------- #
dim(sample_sheet); dim(asv_sheet); dim(tax_sheet)
all(rownames(asv_sheet) == rownames(sample_sheet))
all(colnames(asv_sheet) == rownames(tax_sheet))


# ----------------------------------------------- #
#      Import to phyloseq and data handling
# ----------------------------------------------- #
# Import to phyloseq
ps_all <- phyloseq(otu_table(asv_sheet, taxa_are_rows = FALSE),
                   sample_data(sample_sheet),
                   tax_table(as.matrix(tax_sheet)))

# Data handling
## Select "sea" samples
subset_samples(ps_all, Site == "sea")
## Select samples with total reads > 10^4
prune_samples(sample_sums(ps_all) > 10^4, ps_all)
## Select "Proteobacteria" phylum
subset_taxa(ps_all, phylum == "Proteobacteria")
## Select top 50 taxa
top_taxa <- names(sort(taxa_sums(ps_all), decreasing = TRUE)[1:50])
prune_taxa(top_taxa, ps_all)
## Transform sample counts
transform_sample_counts(ps_all, function(x) 100* x / sum(x))
transform_sample_counts(ps_all, log)


# ----------------------------------------------- #
#               Basic visualization
# ----------------------------------------------- #
# plot_bar
p1 <- plot_bar(ps_all, x = "Sample", fill = "phylum") + scale_fill_igv()
p2 <- plot_bar(ps_all, x = "replicate", fill = "phylum") +
  facet_wrap(.~ Site+Method) + scale_fill_igv()

# plot_richness
p3 <- plot_richness(ps_all, measures = "Observed")

# Visualize patterns (scatterplot)
ps_m1 <- psmelt(ps_all) %>%
  select(OTU, Sample, Abundance, Site, Method) %>%
  filter(OTU == top_taxa[1], Site == "pond")
p4 <- ggplot(ps_m1, aes(y = Abundance, x = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_jitter(height = 0, width = 0.2)


# ----------------------------------------------- #
#              Rarefaction curve
# ----------------------------------------------- #
## From https://stackoverflow.com/questions/47234809/coloring-rarefaction-curve-lines-by-metadata-vegan-package-phyloseq-package
# Rarefy samples
ps_rr <- vegan::rarecurve(as.matrix(otu_table(ps_all)), step = 50, label = TRUE)
rare <- lapply(ps_rr, function(x){
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw_read = rownames(b))
  b$raw_read <- as.numeric(gsub("N", "",  b$raw_read))
  return(b)
})
# Label list
names(rare) <- rownames(sample_data(ps_all))
# Convert to the data frame
rare_df <- map_dfr(rare, function(x) return(data.frame(x)), .id = "sample")
# Visualize with ggplot
p5 <- ggplot(rare_df, aes(x = raw_read, y = ASV, color = sample)) +
  geom_line() + scale_color_igv() +
  xlab("Reads") + ylab("The number of ASV")


# ----------------------------------------------- #
#           Dimensionality reduction
# ----------------------------------------------- #
# Nonmetric dimensional scaling
ps_sea <- subset_samples(ps_all, Site == "sea" & sample_nc == "sample")
set.seed(1234)
ps_bray <- ordinate(ps_sea, "NMDS", "bray")
p6 <- plot_ordination(ps_sea, ps_bray, color = "Method", shape = "Method") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=Method)) +
  geom_point(size = 2) + scale_color_startrek() + scale_fill_startrek() +
  xlab("Axis 1") + ylab("Axis 2") + ggtitle("NMDS")

# t-SNE
#library(tsnemicrobiota)
tsne_res <- tsnemicrobiota::tsne_phyloseq(ps_sea, distance='bray', perplexity = 5, rng_seed = 1234)
p7 <- tsnemicrobiota::plot_tsne_phyloseq(ps_sea, tsne_res, color = "Method", shape = "Method") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=Method)) +
  geom_point(size = 2) + scale_color_startrek() + scale_fill_startrek() +
  xlab("Axis 1") + ylab("Axis 2") + ggtitle("t-SNE")

# UMAP by uwot package
#library(uwot)
set.seed(1234)
umap_res <- uwot::tumap(otu_table(ps_sea)@.Data, n_neighbors = 5, n_components = 2)
umap_df <- cbind(data.frame(sample_data(ps_sea)), umap_res)
colnames(umap_df)[(ncol(umap_df)-1):ncol(umap_df)] <- c("UMAP1", "UMAP2")
p8 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Method, shape = Method)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=Method)) +
  geom_point(size = 2) + scale_color_startrek() + scale_fill_startrek() +
  xlab("Axis 1") + ylab("Axis 2") + ggtitle("UMAP")

# Show all results
p678 <- plot_grid(p6 + theme(legend.position = "none"),
                  p7 + theme(legend.position = "none"),
                  p8 + theme(legend.position = "none"),
                  get_legend(p6),
                  rel_widths = c(1,1,1,0.3),
                  ncol = 4)


# ----------------------------------------------- #
#             Other helpful packages
# ----------------------------------------------- #
# speedyseq (https://github.com/mikemc/speedyseq)
#library(speedyseq)
system.time(psmelt(ps_all))             # 0.103 sec elapsed
system.time(speedyseq::psmelt(ps_all))  # 0.025 sec elapsed

# fantaxtic (https://github.com/gmteunisse/Fantaxtic)
#library(fantaxtic)
ps_top20 <- fantaxtic::get_top_taxa(ps_all, n = 20, relative = TRUE,
                                    discard_other = FALSE, other_label = "Other")
p9 <- fantaxtic::fantaxtic_bar(ps_top20, color_by = "phylum",
                                label_by = "phylum", other_label = "Other")


# ----------------------------------------------- #
#       Save figures and session information
# ----------------------------------------------- #
# Save figures
ggsave(file = "DADA2PhyloseqOut/p1_barplot.pdf", plot = p1, width = 10, height = 6)
ggsave(file = "DADA2PhyloseqOut/p2_barplotfacet.pdf", plot = p2, width = 10, height = 6)
ggsave(file = "DADA2PhyloseqOut/p3_richness.pdf", plot = p3, width = 8, height = 6)
ggsave(file = "DADA2PhyloseqOut/p4_boxplot.pdf", plot = p4, width = 5, height = 4)
ggsave(file = "DADA2PhyloseqOut/p5_rarefaction.pdf", plot = p5, width = 10, height = 6)
ggsave(file = "DADA2PhyloseqOut/p6_dim_reduction.pdf", plot = p678, width = 14, height = 5)
ggsave(file = "DADA2PhyloseqOut/p9_fantaxtic.pdf", plot = p9, width = 10, height = 6)

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
