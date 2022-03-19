# Convenient functions for coverage-based rarefaction and its visualization
This repository includes convenient functions to perform coverage-based rarefaction. `phyloseq` object can be easily rarefied based on a user-specified coverage by the following command. Functions immplemented in `iNEXT` package is used. `rarefy_even_coverage` returns identical results with `phyloseq_coverage_raref` function in `metagMisc`, but a phyloseq object with `taxa_are_rows = TRUE` or `taxa_are_rows = FALSE` is accepted. In addition, `rarefy_even_coverage` returns a rarefaction curve for visualization.


```{r}
# Calculate rarefied matrix and rarefaction curve simultaneously
ps_rare_raw <- rarefy_even_coverage(ps_sample,
                                    coverage = 0.97,
                                    include_iNEXT_results = TRUE)
# Extract phyloseq object
ps_rare <- ps_rare_raw[[1]]                      
```

In addition, results of the coverage-based rarefaction can be checked by visualizing the rarefaction curves.

```{r}
plot_rarefy(ps_rare_raw)
```

<img src="img/rarefy_plot.png" width="800px">

For more detail, please run `demo_rarefy.R`.