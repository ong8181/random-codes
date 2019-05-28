####
#### Test regularized S-map performance
####
#### 2019.5.15 Ushio
#### R 3.5.2
####

# Set random seeds (for reproduction)
ran_seed <- 123
set.seed(ran_seed)

# Load packages
library(rEDM); packageVersion("rEDM") # 0.7.4, 2019.5.15
library(glmnet); packageVersion("glmnet") # 2.0.16, 2019.5.15
library(doParallel); packageVersion("doParallel")

# Load pacakges for visualization
library(ggplot2); packageVersion("ggplot2") # 3.1.0, 2019.5.15
library(cowplot); packageVersion("cowplot") # 0.9.4, 2019.5.15

# Load defined functions
source("functions/Extended_SSR_v1.R")
source("functions/Extended_Smap_v1.R")

# Load workspace
#d <- read.table("data/data_2sp.txt") # Model time series
# A subet of fish time series used in Ushio et al. (2018)
d <- read.csv("data/data_four_fish_species.csv")


##--------------- Test implemented functions ---------------##
#---------- Preparations ----------#
# Normalize and check model time series
y1 <- as.numeric(scale(d[,1]))
plot(y1, type = "l", col = 2)

# Estimate best embeding dimension
simp_res1 <- simplex(y1, E = 1:20, silent = T)
(Ey1 <- simp_res1[which.min(simp_res1$rmse), "E"])
block_y1 <- make_block(y1, max_lag = Ey1)[,2:(Ey1+1)]

# Set parameter ranges
theta_test <- c(0, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2)
lambda_test <- c(0, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2)


#---------- rEDM implementation of the S-map ----------#
y1_redm <- block_lnlp(block_y1, theta = theta_test, method = "s-map", silent = T) # rEDM implementation


#---------- Manual implementation of the regularized S-map ----------#
# Register cores
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

# Prepare result objects
y1_smap_stats <- y1_ridge_stats <- y1_lasso_stats <-
  data.frame(theta = NA, lambda = NA, N = NA, rho = NA, mae = NA, rmse = NA)

# Check parameter dependence (take time, depending on your environment)
# If "lambda" value is not provided, 10-fold cross validation will be automatically applied to determine the best lambda.
for(i in 1:length(theta_test))
{
  for(j in 1:length(lambda_test))
  {
    y1_smap <- extended_lnlp(block_y1, theta = theta_test[i], lambda = lambda_test[j],
                             method = "s-map", regularized = F, random_seed = ran_seed, no_parallel = FALSE)
    y1_ridge <- extended_lnlp(block_y1, theta = theta_test[i], lambda = lambda_test[j],
                              method = "s-map", regularized = T, alpha = 0, random_seed = ran_seed, no_parallel = FALSE)
    y1_lasso <- extended_lnlp(block_y1, theta = theta_test[i], lambda = lambda_test[j],
                              method = "s-map", regularized = T, alpha = 1, random_seed = ran_seed, no_parallel = FALSE)
    
    # Summarize results
    y1_smap_stats[((i-1)*8+j),] <- data.frame(theta = theta_test[i], lambda = lambda_test[j], y1_smap$stats)
    y1_ridge_stats[((i-1)*8+j),] <- data.frame(theta = theta_test[i], lambda = lambda_test[j], y1_ridge$stats)
    y1_lasso_stats[((i-1)*8+j),] <- data.frame(theta = theta_test[i], lambda = lambda_test[j], y1_lasso$stats)
  }
}

parallel::stopCluster(cl)

#---------- Compare results ----------#
# Check rEDM results
# Resutls of rEDM and manual S-map are identical
y1_smap_redm <-cbind(y1_redm[,c("theta")],
                     data.frame(lambda = NA),
                     y1_redm[,c("num_pred", "rho", "mae", "rmse")],
                     data.frame(method = "rEDM"))
colnames(y1_smap_redm)[1] <- "theta"
y1_smap_stats[seq(1,57,8),]; y1_smap_redm

# Combine results
y1_smap_stats$method <- "smap_manual"
y1_ridge_stats$method <- "ridge"
y1_lasso_stats$method <- "lasso"
stats_all <- rbind(y1_smap_stats[seq(1,57,8),], y1_ridge_stats, y1_lasso_stats)
stats_all$lambda <- factor(stats_all$lambda)
             
# ggplot visualization      
g1 <- ggplot(stats_all, aes(x = theta, y = rmse, color = lambda, facet = method, group = lambda))
g1 <- g1 + geom_point() + geom_line() + facet_grid(.~method)
g1 <- g1 + xlab(expression(theta)) + ylab("RMSE")
g1 <- g1 +  scale_color_discrete(name = expression(lambda))

pdf(file = "fig/RegularizedSmapTest.pdf", width = 10, height = 5)
g1; dev.off()

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("sessionInfo/RegularizedSmap_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
