####
#### Generate all possible combinations of degenerate primers
####

#------------------------------------------------------#
# expand_degenerate_primer(character)
#------------------------------------------------------#
expand_degenerate_primer <- function(input_primer = "ATGCN") {
  all_degenerate_bases <- c("R", "M", "W", "S", "Y", "K", "H", "B", "D", "V", "N")
  base2_degenerate_bases <- c("R", "M", "W", "S", "Y", "K")
  base3_degenerate_bases <- c("H", "B", "D", "V")
  R = c("A",      "G"     )
  M = c("A",           "C")
  W = c("A", "T"          )
  S = c(          "G", "C")
  Y = c(     "T",      "C")
  K = c(     "T", "G"     )
  H = c("A", "T",      "C")
  B = c(     "T", "G", "C")
  D = c("A", "T", "G"     )
  V = c("A",      "G", "C")
  N = c("A", "T", "G", "C")
  degenerate_base_list <- list(R, M, W, S, Y, K, H, B, D, V, N)
  names(degenerate_base_list) <- all_degenerate_bases
  
  # Split primer
  split_primer <- strsplit(input_primer, split = NULL)[[1]]
  
  # Identify the number of degenerate bases
  not_ATGC_base <- !(split_primer %in% c("A", "T", "G", "C"))
  base_all_list <- degenerate_base_list[split_primer[not_ATGC_base]]
  expand_combination <- expand.grid(base_all_list)
  expand_combination <- apply(expand_combination, 2, as.character)
  
  # Prepare output objects
  primer_list <- matrix(rep(split_primer, nrow(expand_combination)),
                        ncol = length(split_primer), byrow = T)
  primer_list[,not_ATGC_base] <- expand_combination
  primer_list_onerow <- apply(primer_list, 1, function(x) paste(x, collapse = ""))
  primer_list_onerow <- matrix(primer_list_onerow, ncol = 1)
  
  # Return expand primers
  return(primer_list_onerow)
  # Output expanded primer list
  #write.table(primer_list_onerow, "primer_list.txt",
  #            col.names = FALSE, row.names = FALSE, quote = FALSE)  
}



