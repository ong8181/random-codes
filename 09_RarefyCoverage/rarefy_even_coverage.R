####
#### rarefy_even_coverage.R
#### script to rarefy samples based on a user-specified coverage
####

# Required packages
require(phyloseq)
require(vegan)
require(tidyverse)
require(purrrlyr)


# Define function
rarefy_even_coverage <-  function(ps_obj,
                                  coverage = 0.97,
                                  rarefy_step = 10,
                                  remove_not_rarefied = FALSE,
                                  show_plot = FALSE,
                                  minimum_reads = 10,
                                  slope_range = 0.003
                                  ){
  ## Main function
  if(min(rowSums(otu_table(ps_obj))) > minimum_reads){
    # Convert ps_obj to OTU table
    com_mat <- as.matrix(otu_table(ps_obj))

    ## Get coverage-based rarefied otu_table
    rare_slopelist <- com_mat %>% array_tree(1) %>%
      map(function(x) rareslope(x, round(seq(1, sum(x), by = rarefy_step))))
    
    # Specify the coverage (slope = 0.03 ==> coverage = 97%)
    slope <- 1 - coverage

    # Define function to obtain a value that is close to the coverage 
    # (Here error = 0.003 is specified, but you may change the value)
    find_closest_id <- function(x) max(which(abs(x - slope) < slope_range))
    coverage_id0 <- rare_slopelist %>% map(find_closest_id) %>% unlist

    # Select samples with finite values
    coverage_id <- which(is.finite(coverage_id0))
    
    # Collect values
    coverage_slope <- coverage_reads <- rep(NA, length(coverage_id0))
    names(coverage_slope) <- names(coverage_reads) <- names(coverage_id0)
    for (i in 1:length(coverage_id0)) {
      if (is.finite(coverage_id0)[i]) {
        coverage_slope[i] <- rare_slopelist[[i]][coverage_id0[i]]
        coverage_reads[i] <- rare_slopelist[[i]][coverage_id0[i]] %>%
          names %>% str_split(pattern = "N") %>%
          unlist %>% .[2] %>% as.numeric
      }
    }
    
    # Inspect raw v.s. rarefied reads
    if (show_plot) {
      plot(rowSums(com_mat[coverage_id,]),
           coverage_reads[coverage_id],
           xlab = "Total sequence reads",
           ylab = "Rarefied sequence reads", las = 1,
           main = sprintf("Coverage: %s", coverage))
    }
    
    # Get rarefied counts
    rarefied_count_list <-
      com_mat[coverage_id,] %>% array_tree(1) %>%
      list(x = ., y = coverage_reads[coverage_id] %>% array_tree) %>%
      pmap(function(x,y) rrarefy(x, y))
    # Combined as data.frame
    rarefied_count <- as.data.frame(do.call(rbind, rarefied_count_list))
    #rowSums(rarefied_count)
    rownames(rarefied_count) <- names(coverage_reads[coverage_id])

    #dim(covraredat); rownames(covraredat)
    # Replace original data
    com_mat2 <- com_mat
    com_mat2[rownames(rarefied_count),] <- as.matrix(rarefied_count)
    # Import to the original phyloseq object
    ps_rare <- ps_obj
    otu_table(ps_rare) <- com_mat2
    # Add rarefied/not-rarefied information
    sample_data(ps_rare) <- sample_data(ps_rare) %>%
      data.frame %>%
      mutate(rarefied = as.logical(is.finite(coverage_id0)),
             sum_reads_original = sample_sums(ps_obj),
             sum_reads_rarefied = sample_sums(ps_rare))

    # Output message
    if (any(is.infinite(coverage_id0))) {
      message1 <- sprintf("Following samples were NOT rarefied: %s",
                          names(coverage_id0)[is.infinite(coverage_id0)] %>%
                            paste(collapse=" "))
    } else {
      message1 <- "All samples were successfully rarefied!"
    }
    # Show message
    message(message1)
    
    # Add rarefied OR not_rarefied information to the sample data
    if (remove_not_rarefied) {
      # Remove not-rarefied samples
      ps_rare <- ps_rare %>% prune_samples(sample_names(.)[coverage_id], .)
      # Return rarefied phyloseq object
      message("Rarefied samples were removed from output as you specified")
      return(ps_rare)
    } else {
      # Return rarefied phyloseq object
      message("Rarefied/not-rarefied samples were kept in the phyloseq object.\nPlease check \'rarefied\' column of the sample_data()")
      return(ps_rare)
    }
    
  }else{
    error_message <- sprintf("Please remove samples with sequence reads < %s. They will not be rarefied.", minimum_reads)
    stop(error_message)
  }
}


