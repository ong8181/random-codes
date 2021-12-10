####
#### rDoAMP v0.1.0
#### R tool for testing whether your primer set amplifies your target
#### - DOes my primer set AMPlify my target? -
####
#### History:
#### 2021.12.9 < v0.1.0, Ushio: The initial versions, locally stored
#### 2021.12.9   v0.1.0, Ushio: For Github version
####

#-----------------------------------------------------------------#
# doamp_auto(): Automatically download target sequences
#               and extract (possible) amplicon sequences
#-----------------------------------------------------------------#
doamp_auto <- function (search_query,
                        F_primer,
                        R_primer,
                        n_retmax = 20,
                        n_mismatch = 0,
                        output_dir = "rDoAMP_Out",
                        random_sampling = TRUE,
                        random_sampling_seed = 1234,
                        n_retidmax = n_retmax * 10,
                        save_parameter = TRUE,
                        save_stat = TRUE,
                        overwrite_output_dir = FALSE) {
  #-----------------------------------------------------------------#
  # Create output directory
  #-----------------------------------------------------------------#
  time_start <- proc.time() # Measure elapsed time
  # Check and create output directory
  if(overwrite_output_dir) {
    dir.create(output_dir, showWarnings = FALSE)
    if(length(list.files(sprintf("%s/", output_dir))) > 0) system(sprintf("rm %s/*", output_dir)) # Clean files in output_dir
  } else {
    if(dir.exists(output_dir)) {
      stop("Output directory already exists")
    } else {
      dir.create(output_dir)
    }
  }
  
  #-----------------------------------------------------------------#
  # Search data in Entrez
  #-----------------------------------------------------------------#
  if (!random_sampling) n_retidmax <- n_retmax
  rentrez_search <- entrez_search(db = "nucleotide", term = search_query, retmax = n_retidmax)
  # Random sampling from rentrez_search w/ n_retidmax IDs
  if (random_sampling & length(rentrez_search$ids) >= n_retmax) {
    set.seed(random_sampling_seed)
    retids <- sample(rentrez_search$ids, n_retmax, replace = FALSE)
  } else {
    retids <- rentrez_search$ids
  }
  
  #-----------------------------------------------------------------#
  # Retrieve using entrez_fetch() function (rettype = "fasta")
  # Save as FASTA file
  #-----------------------------------------------------------------#
  entrez_fasta <- entrez_fetch(db = "nucleotide", id = retids,
                               rettype="fasta", parsed = FALSE)
  write.table(entrez_fasta, file = sprintf("%s/entrez_fasta.fa", output_dir),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  #-----------------------------------------------------------------#
  # Compile the sequence file using seqkit in Shell
  #-----------------------------------------------------------------#
  shell_command1 <- sprintf("seqkit seq -w 0 %s/entrez_fasta.fa > %s/download.fa",
                            output_dir, output_dir)
  system(shell_command1)
  # Delete a temporal file
  system(sprintf("rm %s/entrez_fasta.fa", output_dir))
  
  #-----------------------------------------------------------------#
  # Generate degenerate primer list
  #-----------------------------------------------------------------#
  # Check degenerated bases
  deg_in_F <- !all(strsplit(F_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))
  deg_in_R <- !all(strsplit(R_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))
  
  if((deg_in_F | deg_in_R) & (n_mismatch > 0)) {
    # Expand degenerate primers
    if(deg_in_F) expanded_F_primer <- expand_degenerate_primer(F_primer) else expanded_F_primer <- F_primer
    if(deg_in_R) expanded_R_primer <- expand_degenerate_primer(R_primer) else expanded_R_primer <- R_primer
    # Generate primer list
    expanded_primer_list <- expand.grid(expanded_F_primer, expanded_R_primer)
    rownames(expanded_primer_list) <- sprintf("p%s", 1:nrow(expanded_primer_list))
    write.table(expanded_primer_list,
                file = sprintf("%s/expanded_primer_list.tsv", output_dir),
                col.names = FALSE, row.names = TRUE, sep = "\t", quote = FALSE)
    
    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -p %s/expanded_primer_list.tsv -m %s %s/download.fa | seqkit rmdup -w 0 > %s/amplified.fa",
                              output_dir,
                              n_mismatch,
                              output_dir, output_dir)
    system(shell_command2)
  } else {
    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -F %s -R %s -w 0 -m %s %s/download.fa > %s/amplified.fa",
                              F_primer, R_primer,
                              n_mismatch,
                              output_dir, output_dir)
    system(shell_command2)
  }

  #-----------------------------------------------------------------#
  # Save query and stats
  #-----------------------------------------------------------------#
  # Collect parameters
  parameter_list <- matrix(c(paste("NCBI search query:", search_query),
                             paste("Forward primer:", F_primer),
                             paste("Reverse primer:", R_primer),
                             paste("N of retrieved sequences:", length(retids)),
                             paste("Maximum N of retrieved sequences:", n_retmax),
                             paste("N of retrieved IDs:", length(rentrez_search$ids)),
                             paste("Maximum N of retrieved IDs:", n_retidmax),
                             paste("Maximum N of mismatches:", n_mismatch),
                             paste("Random sampling:", random_sampling),
                             paste("Random sampling seed:", random_sampling_seed),
                             paste("Output directory:", output_dir)),
                           ncol = 1)
  if (save_parameter) {
    write.table(parameter_list, sprintf("%s/parameter_list.txt", output_dir),
                quote = FALSE, col.names = FALSE, row.names = FALSE) 
  }
  
  if (save_stat) {
    system(sprintf("seqkit stats -T %s/*.fa > %s/stat.tsv", output_dir, output_dir))
    system(sprintf("cat %s/stat.tsv", output_dir))
  } else {
    system(sprintf("seqkit stats -T %s/*.fa", output_dir))
  }
  
  # Generate output message
  time_elapsed <- (proc.time() - time_start)[3]
  message(paste("\nNCBI search query:", search_query))
  message(paste("Primer set:", F_primer, "-", R_primer, "\n"))
  message(paste(round(time_elapsed, 2), "sec elapsed for downloading and extracting sequences\n"))
}



#-----------------------------------------------------------------#
# doamp_custom(): Extract (possible) amplicon sequences from
#                 user-specified FASTA file
#-----------------------------------------------------------------#
doamp_custom <- function (target_fasta,
                          F_primer,
                          R_primer,
                          n_mismatch = 0,
                          output_dir = "rDoAMP_Out",
                          save_parameter = TRUE,
                          save_stat = TRUE,
                          overwrite_output_dir = FALSE) {
  #-----------------------------------------------------------------#
  # Creat output directory
  #-----------------------------------------------------------------#
  time_start <- proc.time() # Measure elapsed time
  # Check and create output directory
  if(overwrite_output_dir) {
    dir.create(output_dir, showWarnings = FALSE)
    if(length(list.files(sprintf("%s/", output_dir))) > 0) system(sprintf("rm %s/*", output_dir)) # Clean files in output_dir
  } else {
    if(dir.exists(output_dir)) {
      stop("Output directory already exists")
    } else {
      dir.create(output_dir)
    }
  }
  # Copy target fasta to output_dir
  system(sprintf("cp %s %s/custom_db0.fa", target_fasta, output_dir))
  
  #-----------------------------------------------------------------#
  # Compile the sequence file using seqkit in Shell
  #-----------------------------------------------------------------#
  shell_command1 <- sprintf("seqkit seq -w 0 %s/custom_db0.fa > %s/custom_db.fa",
                            output_dir, output_dir)
  system(shell_command1)
  # Delete a temporal file
  system(sprintf("rm %s/custom_db0.fa", output_dir))
  
  #-----------------------------------------------------------------#
  # Generate degenerate primer list
  #-----------------------------------------------------------------#
  # Check degenerated bases
  deg_in_F <- !all(strsplit(F_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))
  deg_in_R <- !all(strsplit(R_primer, split = NULL)[[1]] %in% c("A", "T", "G", "C"))
  
  if((deg_in_F | deg_in_R) & (n_mismatch > 0)){
    # Expand degenerate primers
    if(deg_in_F) expanded_F_primer <- expand_degenerate_primer(F_primer) else expanded_F_primer <- F_primer
    if(deg_in_R) expanded_R_primer <- expand_degenerate_primer(R_primer) else expanded_R_primer <- R_primer
    # Generate primer list
    expanded_primer_list <- expand.grid(expanded_F_primer, expanded_R_primer)
    rownames(expanded_primer_list) <- sprintf("p%s", 1:nrow(expanded_primer_list))
    write.table(expanded_primer_list,
                file = sprintf("%s/expanded_primer_list.tsv", output_dir),
                col.names = FALSE, row.names = TRUE, sep = "\t", quote = FALSE)
    
    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -p %s/expanded_primer_list.tsv -m %s %s/custom_db.fa | seqkit rmdup -w 0 > %s/amplified.fa",
                              output_dir,
                              n_mismatch,
                              output_dir, output_dir)
    system(shell_command2)
  } else {
    #-----------------------------------------------------------------#
    # Extract the barcoding region using seqkit amplicon
    #-----------------------------------------------------------------#
    shell_command2 <- sprintf("seqkit amplicon -F %s -R %s -w 0 -m %s %s/custom_db.fa > %s/amplified.fa",
                              F_primer, R_primer,
                              n_mismatch,
                              output_dir, output_dir)
    system(shell_command2)
  }
  
  #-----------------------------------------------------------------#
  # Save query and stats
  #-----------------------------------------------------------------#
  # Collect parameters
  parameter_list <- matrix(c(paste("Target FASTA file:", target_fasta),
                             paste("Forward primer:", F_primer),
                             paste("Reverse primer:", R_primer),
                             paste("N of maximum mismatch:", n_mismatch),
                             paste("Output directory:", output_dir)),
                           ncol = 1)
  if (save_parameter) {
    write.table(parameter_list, sprintf("%s/parameter_list.txt", output_dir),
                quote = FALSE, col.names = FALSE, row.names = FALSE) 
  }
  
  if (save_stat) {
    system(sprintf("seqkit stats -T %s/*.fa > %s/stat.tsv", output_dir, output_dir))
    system(sprintf("cat %s/stat.tsv", output_dir))
  } else {
    system(sprintf("seqkit stats -T %s/*.fa", output_dir))
  }
  
  # Generate output message
  time_elapsed <- (proc.time() - time_start)[3]
  message(paste("\nTarget FASTA file:", target_fasta))
  message(paste("Primer set:", F_primer, "-", R_primer, "\n"))
  message(paste(round(time_elapsed, 2), "sec elapsed for extracting amplicons"))
}



#-----------------------------------------------------------------#
# Popular primer sets
#-----------------------------------------------------------------#
if (FALSE) {
  # Primer sequences
  MiFish_U <- list(forward = "GTCGGTAAAACTCGTGCCAGC",
                   reverse = "CATAGTGGGGTATCTAATCCCAGTTTG")
  Prok16S_515F_806R <- list(forward = "GTGYCAGCMGCCGCGGTAA",
                            reverse = "GGACTACNVGGGTWTCTAAT")
  Prok16S_515F_926R <- list(forward = "GTGYCAGCMGCCGCGGTAA",
                            reverse = "CCGYCAATTYMTTTRAGTTT")
  FungiITS_ITS1_F_KYO1_ITS2_KYO2 <- list(forward = "CTHGGTCATTTAGAGGAASTAA",
                                         reverse = "TTYRCTRCGTTCTTCATC")
  Euk18S_1391F_EukBr <- list(forward = "GTACACACCGCCCGTC",
                             reverse = "TGATCCTTCTGCAGGTTCACCTAC")
  Euk18S_1389F_1510R <- list(forward = "TTGTACACACCGCCC",
                             reverse = "CCTTCYGCAGGTTCACCTAC")
  AnimalCOI_mlCOIintF_HCO2198 <- list(forward = "GGWACWGGWTGAACWGTWTAYCCYCC",
                                      reverse = "TAAACTTCAGGGTGACCAAAAAATCA")
  Euk_nuSSU1334_50_nuSSU1648_30 <- list(forward = "CGATAACGAACGAGACCT",
                                        reverse = "ANCCATTCAATCGGTANT")
  
  # Combine as a list
  popular_primer_set <- list(MiFish_U,
                             Prok16S_515F_806R,
                             Prok16S_515F_926R,
                             FungiITS_ITS1_F_KYO1_ITS2_KYO2,
                             Euk18S_1391F_EukBr,
                             Euk18S_1389F_1510R,
                             AnimalCOI_mlCOIintF_HCO2198,
                             Euk_nuSSU1334_50_nuSSU1648_30)
  names(popular_primer_set) <- c("MiFish_U",
                                 "Prok16S_515F_806R",
                                 "Prok16S_515F_926R",
                                 "FungiITS_ITS1_F_KYO1_ITS2_KYO2",
                                 "Euk18S_1391F_EukBr",
                                 "Euk18S_1389F_1510R",
                                 "AnimalCOI_mlCOIintF_HCO2198",
                                 "Euk_nuSSU1334_50_nuSSU1648_30")
  
  # Save primer list
  saveRDS(popular_primer_set, "popular_primer_set.obj")
}

