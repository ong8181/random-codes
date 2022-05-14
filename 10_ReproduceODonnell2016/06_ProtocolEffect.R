####
#### Visualize patterns
#### R 4.1.2
#### 2022.5.13 Ushio
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 1.3.1, 2022.2.16
library(phyloseq); packageVersion("phyloseq") # 1.38.0, 2022.2.16
library(iNEXT); packageVersion("iNEXT") # 2.0.20, 2022.5.13
library(cowplot); packageVersion("cowplot") # 1.1.1, 2022.2.16
library(ggsci); packageVersion("ggsci") # 2.9
library(ggrepel); packageVersion("ggrepel") # 0.9.1, 2022.5.13
library(cols4all); packageVersion("cols4all") # 0.3, 2022.5.13
theme_set(theme_cowplot())
source("funcs/F03_HelperFunctions.R") # Helper function for visualization

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# <-----------------------------------------------------> #
#  Load and compile data
# <-----------------------------------------------------> #
# Load sample data
seqtab_data <- read.csv("05_SummarizeOut/otu_table.csv", row.names = 1)
sample_sheet <- read.csv("05_SummarizeOut/sample_data.csv", row.names = 1)
tax_sheet <- read.csv("05_SummarizeOut/tax_table.csv", row.names = 1)

# Add "double_PCR" and "single_PCR" information
# According to O'Donnell et al. (2016),
# total reads for "double PCR" = 16,635,743 (= library1_S1 + library2_S2)
# total reads for "single PCR" = 13,200,683 (= RK16sTS1_S1 + RK16sTS2_S2 + RK16sTS3_S3)
sample_sheet <- sample_sheet %>% mutate(pcr_method = case_when(
  library == "20150401-lib1" ~ "double_PCR",
  library == "20150401-lib2" ~ "double_PCR",
  library == "20141113-lib1" ~ "single_PCR",
  library == "20141113-lib2" ~ "single_PCR",
  library == "20141113-lib3" ~ "single_PCR"
))


# <-----------------------------------------------------> #
#  (Re-)import to phyloseq and transform data to %
# <-----------------------------------------------------> #
# Check data structure
all(rownames(sample_sheet) == rownames(seqtab_data))
all(colnames(seqtab_data) == rownames(tax_sheet))

# Import to phyloseq
ps_all <- phyloseq(otu_table(seqtab_data, taxa_are_rows = FALSE),
                   sample_data(sample_sheet),
                   tax_table(as.matrix(tax_sheet)))

# Transform read counts
ps_rel_naive <- ps_all %>% 
  subset_samples(sample_type == "environmental") %>%
  prune_samples(sample_sums(.) > 0, .) %>% 
  transform_sample_counts(function(x) x/sum(x))

# Other option: Transform read counts after coverage-based rarefaction
ps_raref0 <- ps_all %>% subset_samples(sample_type == "environmental") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  rarefy_even_coverage(., coverage = 0.9999, include_iNEXT_results = TRUE)
## Visualize rarefaction curve
g0 <- plot_rarefy(ps_raref0, plot_rarefied_point = TRUE)
g1 <- plot_rarefy(ps_raref0, plot_rarefied_point = FALSE)
## Transform read counts
ps_raref <- ps_raref0[[1]] %>% transform_sample_counts(function(x) x/sum(x))


# <-----------------------------------------------------> #
#  Extract dominant OTUs as in O'Donnell et al. (2016)
# <-----------------------------------------------------> #
# Check dominant OTUs
top10_naive1 <- ps_rel_naive %>%
  subset_samples(pcr_method == "single_PCR") %>% 
  subset_taxa(kingdom == "Metazoa") %>% 
  taxa_sums() %>% sort(decreasing = T) %>% names %>% .[1:10]
top10_naive2 <- ps_rel_naive %>%
  subset_samples(pcr_method == "double_PCR") %>% 
  subset_taxa(kingdom == "Metazoa") %>% 
  taxa_sums %>% sort(decreasing = T) %>% names %>% .[1:10]
top10_raref1 <- ps_raref %>%
  subset_samples(pcr_method == "single_PCR") %>% 
  subset_taxa(kingdom == "Metazoa") %>% 
  taxa_sums %>% sort(decreasing = T) %>% names %>% .[1:10]
top10_raref2 <- ps_raref %>%
  subset_samples(pcr_method == "double_PCR") %>% 
  subset_taxa(kingdom == "Metazoa") %>% 
  taxa_sums %>% sort(decreasing = T) %>% names %>% .[1:10]
top10_naive1; top10_naive2; unique(top10_naive1, top10_naive2)
top10_raref1; top10_raref2; unique(top10_raref1, top10_raref2)

# REF:Taxon from O'Donnell et al. (2016)
odonnell_taxa <- c("Elysia pusilla", "Olea hansineensis", "Mytilus trossulus", "Percomorphaceae",
                   "Sessilia", "Cymatogaster aggregata", "Chthamalus", "Percomorphaceae",
                   "Cyphophthalmus", "Teleostei", "Homo sapiens")

# Check taxonomy
ps_top <- ps_raref %>% 
  prune_taxa(unique(top10_raref1, top10_raref2), .) %>% 
  prune_samples(sample_sums(.) > 0, .)
tax_table(ps_top)[,c("phylum", "family", "genus", "species")]


# --------------------------------------------- #
# Visualize patterns
# --------------------------------------------- #
ps_top_melt <- speedyseq::psmelt(ps_top) %>% 
  select(OTU, Sample, Abundance, sequencing_run, env_sample_name,
         index_sequence, library, pcr_method, rarefied_n_taxa,
         phylum, family, genus)

# General pattern
s0 <- ps_top_melt %>% 
  ggplot(aes(x = pcr_method, y = Abundance, color = env_sample_name)) +
  geom_jitter(width = 0.2, height = 0) +
  #scale_color_igv() +
  scale_y_log10() +
  scale_color_manual(values = c4a("wright25")) +
  facet_wrap(~ OTU) + panel_border() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab(NULL) + ylab("Relative abundance (%)") +
  labs(title = "PCR method and sample and the relative abundance of dominant OTUs") +
  NULL

# Three most dominant OTUs
s1 <- ps_top_melt %>% filter(OTU == "OTU00003") %>% 
  ggplot(aes(x = pcr_method, y = Abundance, color = index_sequence)) +
  geom_jitter(width = 0.2, height = 0) +
  scale_y_log10() + 
  scale_color_manual(values = c4a("wright25")) +
  facet_wrap(~ env_sample_name) + panel_border() +
  xlab(NULL) + ylab("Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Effects of tag and PCR method on the relative abundance of OTU00003",
       subtitle = "Each panel indicates each environmental sample") +
  NULL
s2 <- (s1 %+% (ps_top_melt %>% filter(OTU == "OTU00005"))) + 
  ggtitle("Effects of tag and PCR method on the relative abundance of OTU00005")
s3 <- (s1 %+% (ps_top_melt %>% filter(OTU == "OTU00008"))) + 
  ggtitle("Effects of tag and PCR method on the relative abundance of OTU00008")


# --------------------------------------------- #
# Bray-Curtis dissimilarity (NMDS)
# --------------------------------------------- #
# Set random seed for reproducibility
set.seed(1234)

# Perform NMDS
ps_bray <- ordinate(ps_raref, "NMDS", "bray")
n1 <- plot_ordination(ps_raref, ps_bray, color = "pcr_method", shape = "pcr_method") +
  #stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=pcr_method)) +
  geom_point(size = 2) + scale_color_manual(values = c4a("wright25")) +
  geom_text_repel(label = sample_data(ps_raref) %>% pull(index_sequence), size = 3, show.legend = F) +
  #geom_label_repel(label = sample_data(ps_raref) %>% pull(index_sequence), fill = "white", size = 3) +
  guides(color = guide_legend(override.aes =
                                list(label = NULL, shape = c(16,17), size = 3, color = c4a("wright25")[1:2]))) + 
  facet_wrap(~ env_sample_name) + panel_border() +
  xlab("Axis 1") + ylab("Axis 2") +
  labs(title = "NMDS",
       subtitle = "Sample, tag sequence, and PCR method") +
  NULL


# --------------------------------------------- #
# Within- and among primer dissimilarities
# --------------------------------------------- #
# Compile Bray-Curtis dissimilarity
vegan_bray <- vegan::vegdist(as.matrix(otu_table(ps_raref)), method = "bray") %>% as.matrix()
vegan_bray[upper.tri(vegan_bray)] <- NA
vegan_sdf <- sample_data(ps_raref) %>% data.frame %>% select(Sample_Name2, env_sample_name, pcr_method, index_sequence)
vegan_df <- cbind(vegan_sdf, vegan_bray)
vegan_longer <- vegan_df %>%
  pivot_longer(cols = - colnames(vegan_sdf), names_to = "pair_sample", values_to = "bray_value") %>% 
  na.omit()
## Add pair_sample information
vegan_longer$pair_index <- vegan_sdf %>%
  slice(match(vegan_longer$pair_sample, Sample_Name2)) %>% pull(index_sequence)
vegan_longer$pair_pcr_method <- vegan_sdf %>%
  slice(match(vegan_longer$pair_sample, Sample_Name2)) %>% pull(pcr_method)
vegan_longer$pair_env_sample <- vegan_sdf %>%
  slice(match(vegan_longer$pair_sample, Sample_Name2)) %>% pull(env_sample_name)

# Select target rows
## Condition 1: Evaluate dissimilarities within the same sample
vegan_longer_same <- vegan_longer %>%
  filter(env_sample_name == pair_env_sample,
         pcr_method == pair_pcr_method) %>% 
  mutate(index_pair_cat = ifelse(index_sequence == pair_index, "within_index", "among_index"))
vegan_longer_same$pcr_method[vegan_longer_same$pcr_method == "single_PCR"] <- "single PCR"
vegan_longer_same$pcr_method[vegan_longer_same$pcr_method == "double_PCR"] <- "double PCR"
vegan_longer_same$pcr_method <- factor(vegan_longer_same$pcr_method, levels = c("single PCR", "double PCR"))
vegan_longer_same$index_pair_cat[vegan_longer_same$index_pair_cat == "within_index"] <- "within index"
vegan_longer_same$index_pair_cat[vegan_longer_same$index_pair_cat == "among_index"] <- "among index"
vegan_longer_same$index_pair_cat <- factor(vegan_longer_same$index_pair_cat, levels = c("within index", "among index"))

# Visualize dissimilarities
b1 <- vegan_longer_same %>% 
  ggplot(aes(x = index_pair_cat, y = bray_value, facet = pcr_method, color = env_sample_name)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.8) +
  facet_wrap(~ pcr_method) + panel_border() +
  scale_color_manual(values = c4a("wright25"), name = "Sample name") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) + 
  xlab(NULL) +
  ylab("Bray-Curtis dissimilarities") +
  NULL


# --------------------------------------------- #
# Save results
# --------------------------------------------- #
# Save figures
## Relative abundance
ggsave(sprintf("%s/TagMethodEffect_Overall.pdf", output_folder),
       plot = s0, width = 10, height = 8)
ggsave(sprintf("%s/TagMethodEffect_Top1.pdf", output_folder),
       plot = s1, width = 10, height = 8)
ggsave(sprintf("%s/TagMethodEffect_Top2.pdf", output_folder),
       plot = s2, width = 10, height = 8)
ggsave(sprintf("%s/TagMethodEffect_Top3.pdf", output_folder),
       plot = s3, width = 10, height = 8)
## NMDS
ggsave(sprintf("%s/TagMethodEffect_NMDS.pdf", output_folder),
       plot = n1, width = 10, height = 8)
## Reproduce O'Donnell et al. (2016) Fig.2
ggsave(sprintf("%s/Reproduce_ODonnell2016_Fig2.pdf", output_folder),
       plot = b1, width = 10, height = 8)

# Save workspace
save(list = ls(all.names = TRUE),
     file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))
