# Convenient functions for coverage-based rarefaction and its visualization
This repository includes convenient functions to perform coverage-based rarefaction. `phyloseq` object can be easily rarefied based on a user-specified coverage by the following command.

Currently, `vegan` version and `iNEXT` version are implemented, but `vegan` version will be deleted in near future.

```{r}
# vegan::rarefy version
ps_rare <- rarefy_even_coverage(ps_sample, coverage = 0.97, rarefy_step = 10)

# iNEXT version
ps_rare2 <- rarefy_coverage_inext(ps_sample, coverage = 0.97)
```

In addition, results of the coverage-based rarefaction can be checked by visualizing the rarefaction curves.

```{r}
# vegan::rarefy version
plot_rarefy(ps_sample, ps_rare)

# iNEXT version
ps_rare_inext <- rarefy_coverage_inext(ps_sample, coverage = 0.97,
                                       include_iNEXT_results = TRUE,
                                       knots = 200, nboot = 1)
plot_rarefy_inext(ps_rare_inext)
```

<img src="img/rarefy_plot.png" width="800px">

For more detail, please run `demo_rarefy.R`.