# Manual implementation of the regularized S-map
The regularized S-map has been implemented according to Cenci et al. (2019) Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13150)

Currently, the regularization relies on "glmnet" package of R.

## Test results
Regularized S-map functions are in "functions" folder, and they are tested in 01_RegularizedSmap.R. Accordign to a test using a subset of the data in Ushio et al. (2018) (see "data" folder), ridge (or lasso) S-map shows a better performance than the standard S-map, which is consistent with Cenci et al. (2019) (run 01_RegularizedSmap.R or see "figs" folder). Current implementations of the regularized and standard S-map are much slower than the standard S-map in rEDM package of R.

## Usage
Regularized S-map can be performed using extended_lnlp() function. The default setting of extended_lnlp() function is as follows (similar to block_lnlp() function in rEDM package):

extended_lnlp(block_time, lib = c(1, NROW(block_time)), pred = lib, tp = 1, target_column = 1, lib_column = 1:NCOL(block_time), num_neighbors = NCOL(block_time) + 1, theta = 0, method = "s-map", regularized = FALSE, lambda = NULL, alpha = 0, glmnet_parallel = FALSE, random_seed = NULL, no_parallel = glmnet_parallel, save_smap_coefficients = FALSE)

Here are descriptions of some important parameters:
- block_time: Dataframe or matrix.
- method: Specify method to make predictions (currently only "s-map" is valid).
- theta: Nonlinearity parameter for S-map.
- regularized: Specify whether regularized or standard S-map is used.
- lambda: Tuning parameter for regularization. If lambda = 0, there will be no regularization.
- alpha: Tuning parameter for elastic net. If alpha = 0, ridge regression will be performed. If alpha = 1, lasso regression will be performed.
- random_seed: Set this seed to reproduce results of the parallel computation mode.
- no_parallel: If FALSE, then the computation will be parallel.
