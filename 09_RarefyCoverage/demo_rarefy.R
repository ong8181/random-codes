####
#### Coverage-based rarefaction
#### 2022.3.7 Ushio
####

#setwd("09_RarefyCoverage/")

# Load libraries and function
library(tidyverse); packageVersion("tidyverse") # 1.3.1
library(phyloseq); packageVersion("phyloseq") # 1.38.0
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
ps_sample <- ps_all %>% subset_samples(sample_nc == "sample" & Site == "sea")
sample_sums(ps_sample)


# ----------------------------------------------- #
# Perform coverage-based rarefaction
# ----------------------------------------------- #
ps_rare <- rarefy_even_coverage(ps_sample, coverage = 0.99, rarefy_step = 10)
sample_sums(ps_rare)
sample_data(ps_rare)


# ----------------------------------------------- #
# Visualize rarefaction
# ----------------------------------------------- #
(g1 <- plot_rarefy(ps_sample, ps_rare))

