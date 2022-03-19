####
#### Coverage-based rarefaction
#### 2022.3.8 Ushio
####

#setwd("09_RarefyCoverage/")

# Load libraries and function
library(tidyverse); packageVersion("tidyverse") # 1.3.1
library(phyloseq); packageVersion("phyloseq") # 1.38.0
library(iNEXT); packageVersion("iNEXT") # 2.0.20
library(ggsci); packageVersion("ggsci") # 2.9
source("rarefy_even_coverage.R")


# -------------------------------------------- #
# Load and check demo data
# -------------------------------------------- #
sample_sheet <- read.csv("data_demo/sample_sheet.csv", row.names = 1)
asv_sheet <- read.csv("data_demo/asv_sheet.csv", row.names = 1)
tax_sheet <- read.csv("data_demo/tax_sheet.csv", row.names = 1)
# Check colnames and rownames
dim(sample_sheet); dim(asv_sheet); dim(tax_sheet)
all(rownames(asv_sheet) == rownames(sample_sheet))
all(colnames(asv_sheet) == rownames(tax_sheet))


# ----------------------------------------------- #
# Import to phyloseq and data handling
# ----------------------------------------------- #
# Import to phyloseq
ps_all <- phyloseq(otu_table(asv_sheet, taxa_are_rows = FALSE),
                   sample_data(sample_sheet),
                   tax_table(as.matrix(tax_sheet)))
ps_sample <- ps_all %>%
  subset_samples(sample_nc == "sample" & Site == "sea") %>%
  prune_taxa(taxa_sums(.) > 0, .)
sample_sums(ps_sample)


# ----------------------------------------------- #
# Perform coverage-based rarefaction
# ----------------------------------------------- #
set.seed(1234)
ps_rare_raw <- rarefy_even_coverage(ps_sample,
                                    coverage = 0.97,
                                    knots = 200,
                                    n_rarefy_iter = 100,
                                    rarefy_average_method = "floor",
                                    include_iNEXT_results = TRUE)
ps_rare <- ps_rare_raw[[1]] # Extract phyloseq object
sample_sums(ps_rare)
sample_data(ps_rare)


# ----------------------------------------------- #
# Visualize rarefaction curve
# ----------------------------------------------- #
g1 <- plot_rarefy(ps_rare_raw)
g1 <- g1 + coord_cartesian(xlim = c(0, 4000)) +
  scale_color_igv()
#ggsave(file = "img/rarefy_plot.png",
#       plot = g1, width = 8, height = 6)


# ----------------------------------------------- #
# Appendix
# ----------------------------------------------- #
# ----------------------------------------------- #
# metagMisc-version
# ----------------------------------------------- #
# Should be 'taxa_are_rows = TRUE'
#ps_all2 <- phyloseq(otu_table(t(asv_sheet), taxa_are_rows = TRUE),
#                    sample_data(sample_sheet),
#                    tax_table(as.matrix(tax_sheet)))
#ps_sample2 <- ps_all2 %>% subset_samples(sample_nc == "sample" & Site == "sea")
# Calculate coverage and rarefy
#ps_rare3 <- metagMisc::phyloseq_coverage_raref(ps_sample2, coverage = 0.97, iter = 1)
#sample_sums(ps_rare3)
# Check results
#plot(sample_sums(ps_rare2), sample_sums(ps_sample)); abline(0,1)
#plot(sample_sums(ps_rare2), sample_sums(ps_rare)); abline(0,1)
#plot(sample_sums(ps_rare2), sample_sums(ps_rare3)); abline(0,1)

