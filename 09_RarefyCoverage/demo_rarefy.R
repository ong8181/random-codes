####
#### Coverage-based rarefaction
#### 2022.3.8 Ushio
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
ps_sample <- ps_all %>%
  subset_samples(sample_nc == "sample" & Site == "sea") %>%
  prune_taxa(taxa_sums(.) > 0, .)
sample_sums(ps_sample)


# ----------------------------------------------- #
# Perform coverage-based rarefaction
# ----------------------------------------------- #
ps_rare <- rarefy_even_coverage(ps_sample, coverage = 0.97, rarefy_step = 10)
sample_sums(ps_rare)
sample_data(ps_rare)


# ----------------------------------------------- #
# Visualize rarefaction curve
# ----------------------------------------------- #
g1 <- plot_rarefy(ps_sample, ps_rare, rarefy_step = 10)
g1 + xlim(0, 4000)
#(g2 <- plot_rarefy(ps_sample, ps_rare, plot_rarefied_point = FALSE, rarefy_step = 10))


# ----------------------------------------------- #
# iNEXT-version
# ----------------------------------------------- #
library(iNEXT)

# iNEXT-version rarefaction
ps_rare3 <- rarefy_even_coverage_inext(ps_sample, coverage = 0.97)
sample_data(ps_rare3)

# iNEXT-version rarefaction
ps_rare_inext <- rarefy_even_coverage_inext(ps_sample, coverage = 0.97,
                                            include_iNEXT_results = TRUE,
                                            knots = 200, nboot = 1)
# Visualize results
g2 <- plot_rarefy_inext(ps_rare_inext)
g2 + xlim(0, 4000)



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
#ps_rare2 <- metagMisc::phyloseq_coverage_raref(ps_sample2, coverage = 0.97, iter = 1)
#sample_sums(ps_rare2)
# Check results
#plot(sample_sums(ps_rare3), sample_sums(ps_sample)); abline(0,1)
#plot(sample_sums(ps_rare3), sample_sums(ps_rare)); abline(0,1)
#plot(sample_sums(ps_rare3), sample_sums(ps_rare2)); abline(0,1)

