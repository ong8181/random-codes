####
#### F02. Helper functions
####

# Taxa name summarize function
taxa_name_summarize <- function(ps_object,
                                taxa_rank,
                                top_taxa_n = 10,
                                taxa_rank_names = colnames(tax_table(ps_object))){
  tax_df <- as.data.frame(tax_table(ps_object))
  if(is.null(tax_df$rep_tax)) tax_df$rep_tax <- "Undetermined"

  # Search Others and Undetermined taxa
  rep_tax_cond1 <- tax_df[,taxa_rank] == "" & !is.na(tax_df[,taxa_rank])
  tax_col1 <- which(colnames(tax_df) == taxa_rank)
  tax_col2 <- which(colnames(tax_df) == "species")
  rep_tax_cond2 <- apply(tax_df[,(tax_col1+1):tax_col2] == "", 1, sum) == (tax_col2 - tax_col1)
  
  # Replace taxa names
  tax_df[!rep_tax_cond1, "rep_tax"] <- as.character(tax_df[!rep_tax_cond1, taxa_rank])
  tax_df[rep_tax_cond1 & !rep_tax_cond2, "rep_tax"] <- "Others"
  
  # Re-import phyloseq object with revised tax_table
  ps_object2 <- phyloseq(otu_table(ps_object), sample_data(ps_object), tax_table(as.matrix(tax_df)))
  
  # Repalce low abundance taxa name with Others
  taxa_abundance_rank <- aggregate(taxa_sums(ps_object2), by = list(tax_table(ps_object2)[,"rep_tax"]), sum)
  taxa_abundance_rank <- taxa_abundance_rank[order(taxa_abundance_rank$x, decreasing = T),]
  taxa_top <- taxa_abundance_rank[1:top_taxa_n,]
  
  low_tax <- is.na(match(tax_table(ps_object2)[,"rep_tax"], as.character(taxa_top[,1])))
  tax_table(ps_object2)[low_tax,"rep_tax"] <- "Others"
  
  return(ps_object2)
}


# Inspect library size
inspect_libsize <- function (phyloseq_obj, neg_id = "sample_nc") {
  df1 <- as.data.frame(sample_data(phyloseq_obj)) # Put sample_data into a ggplot-friendly data.frame
  df1$LibrarySize <- sample_sums(phyloseq_obj)
  df1 <- df1[order(df1$LibrarySize),]; df1$Index <- seq(nrow(df1))
  
  # Return ggplot object
  return(df1 %>%
           ggplot(aes_string(x = "Index", y = "LibrarySize", color = as.name(neg_id))) +
           geom_point())
}


# Visualize decontam prevalence results
inspect_prevalence <- function (phyloseq_obj,
                                  decontam_result,
                                  neg_id = "is_neg_decontam") {
  # Extract pos/neg information
  ps_taxsum_neg <- phyloseq_obj %>% transform_sample_counts(function(x) 1*(x>0)) %>%
    prune_samples(sample_data(.) %>% pull(as.name(neg_id)), .) %>% taxa_sums
  ps_taxsum_pos <- phyloseq_obj %>% transform_sample_counts(function(x) 1*(x>0)) %>%
    prune_samples(!sample_data(.) %>% pull(as.name(neg_id)), .) %>% taxa_sums
  # Make data.frame of prevalence in positive and negative samples
  df1_pa <- data.frame(pa_pos = ps_taxsum_pos,
                       pa_neg = ps_taxsum_neg,
                       contaminant=decontam_result$contaminant)
  
  # Return ggplot object
  return(ggplot(data=df1_pa, aes(x = pa_neg, y = pa_pos, color = contaminant)) + geom_point() +
           xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)"))
}


