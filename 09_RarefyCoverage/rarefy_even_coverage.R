####
#### rarefy_even_coverage.R
#### script to rarefy samples based on the coverage
#### Reference: https://www.fifthdimension.jp/wiki.cgi?page=%BC%AB%CD%B3%BD%B8%B2%F12016%A1%A7%A5%E1%A5%BF%A5%D0%A1%BC%A5%B3%A1%BC%A5%C7%A5%A3%A5%F3%A5%B0%A1%A6%B4%C4%B6%ADDNA%A5%D0%A1%BC%A5%B3%A1%BC%A5%C7%A5%A3%A5%F3%A5%B0%B2%F2%C0%CF%A4%CE%B5%BB%CB%A1&file=Rscript%5FKadowaki%2Ezip&action=ATTACH
####

# Set working directory
#setwd("09_RarefyCoverage/")

# Define function
rarefy_even_coverage <-  function(ps_obj, coverage = 0.97){
  #ps_obj <- ps_sample
  #coverage <- 0.97
  
  if(min(rowSums(otu_table(ps_obj))) > 10){
    # Convert ps_obj to OTU table
    comm1 <- as.matrix(otu_table(ps_obj))
    
    ## Get coverage-based rarefied otu_table
    rare_slopelist <- list()
    for(i in 1:nrow(comm1)){
      rare_slopelist[[i]] <- comm1[i,] %>%
        vegan::rareslope(round(seq(1, sum(.), by = 10)))
    }
    
    # Specify the coverage (slope = 0.03 ==> coverage 97%)
    slope <- 1 - coverage
    # Define function to obtain a value that is close to the coverage 
    # (Here error = 0.003 is specified, but you may change the value)
    # (The results of rarefaction may probably not depend too much on the value)
    cvr_rare <- lapply(rare_slopelist,
                       function(x) {
                         max(which(abs(x - slope) < 0.003))
                       }
    ) %>% unlist()
    #print(cvrrare)
    # Select samples with finite values
    cvrno <- which(is.finite(cvr_rare))
    # Inspect raw v.s. rarefied reads
    if (F) {
      plot(rowSums(comm1[cvrno,]),
           cvr_rare[cvrno],
           xlab = "Total No of reads",
           ylab = "Rarefied No of reads")
    }
    
    # Get rarefied counts
    covraredat <- data.frame()
    for(i in 1:length(cvrno)){
      covraredat <- rbind(covraredat, vegan::rrarefy(comm1[cvrno[i],], cvr_rare[cvrno[i]]))
    }
    #dim(covraredat); rownames(covraredat)
    # Replace original data
    comm2 <- comm1
    comm2[rownames(covraredat),] <- as.matrix(covraredat)
    #  i <- 30; plot(as.numeric(comm2[i,]), as.numeric(comm1[i,]))
    
    otu_table(ps_obj) <- comm2
    print(is.finite(cvr_rare))
    ps_obj <- subset_samples(ps_obj, is.finite(cvr_rare))
    # plot_richness(phyloseq, measures = "Observed")
    return(ps_obj)
  }else{
    print("Error: please remove samples with no. reads < 10")
  }
}



