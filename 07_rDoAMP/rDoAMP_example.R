####
#### rDoAMP v0.1.0
#### R tool for testing whether your primer set amplifies your target
#### - DOes my primer set AMPlify my target sequences? -
####
#### History:
#### 2021.12.09 < v0.1.0, Ushio: The initial versions, locally stored
#### 2021.12.10   v0.1.0, Ushio: For Github version
####


#------------------ REQIRED packages ------------------#
## seqkit (https://bioinf.shenwei.me/seqkit/)
## "rentres" R package (https://github.com/ropensci/rentrez)
#------------------------------------------------------#


# Load library
library(rentrez); packageVersion("rentrez") # 1.2.3, 2021.12.9
#entrez_db_searchable("nucleotide") # Check query options

# Load functions
source("func/expand_deg_primer.R") # Function to expand degenerate primers
source("func/rDoAMP_v0.1.0.R") # rDoAMP
# Load popular primer sets
popular_primer_set <- readRDS("func/popular_primer_set.obj")
dir.create("0_SessionInfo") # To save session information here


#-----------------------------------------------------------------#
# SET USER-SPECIFIED PARAMETERS
#-----------------------------------------------------------------#
## Fish example
SEARCH_QUERY <- "Trachurus AND mitochondrion AND 1000:20000[SLEN]"
N_RETMAX <- 200 # Maximum number of hits returned by the search
F_PRIMER <- popular_primer_set$MiFish_U$forward # Or, "GTCGGTAAAACTCGTGCCAGC"
R_PRIMER <- popular_primer_set$MiFish_U$reverse # Or, "CATAGTGGGGTATCTAATCCCAGTTTG"
N_MISMATCH <- 0 # Maximum number of primer-template mismatch

## Prokaryote example
#SEARCH_QUERY <- "Erwinia AND 16S AND 500:1000[SLEN]"
#N_RETMAX <- 200 # Maximum number of hits returned by the search
#F_PRIMER <- popular_primer_set$Prok16S_515F_806R$forward # Or, "GTGYCAGCMGCCGCGGTAA"
#R_PRIMER <- popular_primer_set$Prok16S_515F_806R$reverse # Or, "GGACTACNVGGGTWTCTAAT"
#N_MISMATCH <- 0 # Maximum number of primer-template mismatch
#-----------------------------------------------------------------#


#-----------------------------------------------------------------#
# Examples: doamp_auto() and doamp_custom()
#-----------------------------------------------------------------#
## Display primer sequences
F_PRIMER; R_PRIMER

## Automatically download target sequences and extract amplicons
doamp_auto(SEARCH_QUERY,
           F_primer = F_PRIMER,
           R_primer = R_PRIMER,
           n_retmax = N_RETMAX,
           n_mismatch = N_MISMATCH,
           output_dir = "DoAMP1_Out",
           overwrite_output_dir = FALSE)

## Extract amplicons from a user-specified custom FASTA file
## (Use download.fa as a user-specified custom FASTA file)
doamp_custom("DoAMP1_Out/download.fa",
             F_primer = F_PRIMER,
             R_primer = R_PRIMER,
             n_mismatch = N_MISMATCH,
             output_dir = "DoAMP2_Out",
             overwrite_output_dir = FALSE)


#-----------------------------------------------------------------#
# Save session information
#-----------------------------------------------------------------#
writeLines(capture.output(sessionInfo()),
           sprintf("0_SessionInfo/SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

