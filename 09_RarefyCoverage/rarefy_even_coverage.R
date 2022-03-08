####
#### rarefy_even_coverage.R
#### Functions to rarefy samples based on a user-specified coverage
####

# Required packages
require(tidyverse)
require(phyloseq)
require(vegan)

#ps_obj = ps_sample
#coverage = 0.97
#rarefy_step = 10
#remove_not_rarefied = FALSE
#show_plot = FALSE
#minimum_reads = 10
#slope_range = 0.003

# ---------------------------------------------- #
# Coverage-based rarefaction
# ---------------------------------------------- #
rarefy_even_coverage <-  function(ps_obj,
                                  coverage = 0.97,
                                  rarefy_step = 10,
                                  remove_not_rarefied = FALSE,
                                  show_plot = FALSE,
                                  minimum_reads = 10,
                                  slope_range = 0.003,
                                  ran_seed = 1234
                                  ){
  ## Check phyloseq object structure
  if (!dim(otu_table(ps_obj))[2] == dim(tax_table(ps_obj))[1]) stop("\'Taxa\' must be in \'column\' in the phyloseq object!")
  
  ## Main function
  if(min(rowSums(otu_table(ps_obj))) > minimum_reads){
    # Set random seed
    set.seed(ran_seed)
    # Convert ps_obj to OTU table
    com_mat <- as.matrix(otu_table(ps_obj))

    # Get coverage-based rarefied otu_table
    rare_slopelist <- com_mat %>% array_tree(1) %>%
      map(function(x) rareslope(x, round(seq(1, sum(x), by = rarefy_step))))
    # Record rarefied reads
    rare_readslist <- com_mat %>% array_tree(1) %>%
      map(function(x) round(seq(1, sum(x), by = rarefy_step)))
    
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
        coverage_reads[i] <- rare_readslist[[i]][coverage_id0[i]]
        #coverage_reads[i] <- rare_slopelist[[i]][coverage_id0[i]] %>%
        #  names %>% str_split(pattern = "N") %>%
        #  unlist %>% .[2] %>% as.numeric
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
    rrlist <- com_mat[coverage_id,] %>% array_tree(1) %>%
      list(x = ., y = coverage_reads[coverage_id] %>% array_tree)
    rarefied_count_list <- rrlist %>% pmap(function(x,y) rrarefy(x, y))
    ## Repeat three rarefactions to mitigate random sampling effects
    ## This increases computation time
    #rarefied_count_list1 <- rrlist %>% pmap(function(x,y) rrarefy(x, y))
    #rarefied_count_list2 <- rrlist %>% pmap(function(x,y) rrarefy(x, y))
    #rarefied_count_list3 <- rrlist %>% pmap(function(x,y) rrarefy(x, y))
    #rarefied_count_list <- rarefied_count_list1 # temproal object
    #for (j in 1:length(rrlist[[1]])) {
    #  rarefied_count_list_tmp <- round(colMeans(rbind(rarefied_count_list1[[j]], 
    #                                                    rarefied_count_list2[[j]], 
    #                                                    rarefied_count_list3[[j]]))-0.3)
    #  rarefied_count_list[[j]] <- rarefied_count_list_tmp
    #}
    
    # Combined as data.frame
    rarefied_count <- as.data.frame(do.call(rbind, rarefied_count_list))
    #rowSums(rarefied_count)
    rownames(rarefied_count) <- names(coverage_reads[coverage_id])
    rarefied_sp <- rowSums(rarefied_count > 0)

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
             n_taxa_original = rowSums(otu_table(ps_obj) > 0),
             sum_reads_rarefied = sample_sums(ps_rare),
             n_taxa_rarefied = rowSums(otu_table(ps_rare) > 0),
             rarefied_slope = coverage_slope)

    # Output message
    if (any(is.infinite(coverage_id0))) {
      message1 <- sprintf("Following samples were NOT rarefied: %s\n",
                          names(coverage_id0)[is.infinite(coverage_id0)] %>%
                            paste(collapse=" "))
    } else {
      message1 <- "All samples were successfully rarefied!\n"
    }
    # Show message
    message(message1)
    
    # Add rarefied OR not_rarefied information to the sample data
    if (remove_not_rarefied) {
      # Remove not-rarefied samples
      ps_rare <- ps_rare %>% prune_samples(sample_names(.)[coverage_id], .)
      # Return rarefied phyloseq object
      message("Rarefied samples were removed from output as you specified.")
      return(ps_rare)
    } else {
      # Return rarefied phyloseq object
      message2 <- "Rarefied/not-rarefied samples were kept in the phyloseq object."
      message3 <- "Sequence reads of the not-rarefied samples were not changed."
      message4 <- "Please check \'rarefied\' column of sample_data()."
      message(message2); message(message3); message(message4)
      return(ps_rare)
    }
    
  }else{
    error_message <- sprintf("Please remove samples with sequence reads < %s. They will not be rarefied.", minimum_reads)
    stop(error_message)
  }
}



# ---------------------------------------------- #
# Visualize rarefaction curve
# ---------------------------------------------- #
#ps_obj <- ps_sample
#ps_obj_rare <- ps_rare
#rarefy_step = 10
plot_rarefy <- function (ps_obj,
                         ps_obj_rare,
                         rarefy_step = 10,
                         plot_rarefied_point = TRUE,
                         ran_seed = 1234) {
  # Convert ps_obj to OTU table
  com_mat <- as.matrix(otu_table(ps_obj))
  # Get coverage-based rarefied otu_table
  set.seed(ran_seed)
  rare_xy <- com_mat %>% array_tree(1) %>%
    map(function(x) rarefy(x, round(seq(1, sum(x), by = rarefy_step)))) %>%
    map(function(i) data.frame(x = attributes(i)$Subsample, y = as.numeric(i)))
  # Label list
  names(rare_xy) <- rownames(sample_data(ps_obj))
  # Convert to the data frame
  rare_df <- map_dfr(rare_xy, function(x) return(data.frame(x)), .id = "sample")
  
  # Compile rarefied data
  rare_df2 <- data.frame(sample_data(ps_obj_rare)) %>%
    rename(x = sum_reads_rarefied, y = n_taxa_rarefied) %>% mutate(sample = rownames(.))
  rare_df2_slope <- rare_df2$rarefied_slope
  rare_df2_intercept <- rare_df2$y - rare_df2$rarefied_slope * rare_df2$x
  
  # Visualize with ggplot
  g_rare <- ggplot(rare_df, aes(x = x, y = y, color = sample)) +
    geom_line(size = 0.5) +
    xlab("Sequence reads") +
    ylab("The number of species") +
    NULL
  
  if (plot_rarefied_point) {
    # Add rarefied information
    g_rare <- g_rare +
      geom_point(data = rare_df2, aes(x = x, y = y, color = sample),
               size = 3, shape = 18) +
      geom_abline(slope = rare_df2_slope,
                  intercept = rare_df2_intercept,
                  linetype = 3, size = 0.5, alpha = 0.8) +
      #xlim(0,10000) +
      NULL
  }
  
  # Return ggplot object
  return(g_rare)
}




