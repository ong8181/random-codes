####
#### Implementation of regularized S-map
#### Originally introduced by Cenci et al. (2019) Methods in Ecology and Evolution
#### see https://doi.org/10.1111/2041-210X.13150
####
#### 2019.05.15 Masayuki Ushio
####

#---------- History ----------#
# 2019.05.15: Fisrt implementation of Regularized S-map
#             Rely on the use of "glmnet" package of R
# 2019.05.16: "weights" option in "glmnet" and "cv.glmnet" functions revised
# 2019.05.23: CV criterion, "lambda.min", has been changed to "lambda.1se"
# 2019.05.28: Parallel version implemented
# 2019.08.28: Bug fixed
#-----------------------------#

extended_smap <- function(vectors,
                          target,
                          lib_indices,
                          pred_indices,
                          theta,
                          regularized = FALSE,
                          lambda = NULL,
                          alpha = 0, # default is the ridge regression.
                                     # If alpha = 1, it is the lasso regression
                          glmnet_parallel = FALSE,
                          random_seed = NULL,
                          no_parallel = glmnet_parallel,
                          save_smap_coefficients = FALSE)
{
  #require(glmnet)

  # set E here
  E <- NCOL(vectors)

  # setup output
  pred_vals <- rep.int(NaN, times = length(target))
  smap_coefficient_vals <- matrix(NaN, ncol = NCOL(vectors) + 1, nrow = length(target))

  # exclude libs that contains NaN in target
  lib_nona <- !apply(cbind(vectors, target), 1, function(x) any(is.na(x)))
  lib_indices <- lib_nona & lib_indices

  # Make new pred_indices (exclude indices that contains NaN)
  pred_indices <- !apply(vectors, 1, function(x) any(is.na(x))) & pred_indices

  # Add NaN if vectors contains NaN #<--- should be deleted?
  pred_vals[!pred_indices] <- NaN
  smap_coefficient_vals[!pred_indices,] <- NaN

  # Make random seed array
  if (!is.null(random_seed))
  {
    set.seed(random_seed)
    p_seeds <- as.integer(runif(length(pred_indices)) * 32768)
    next_seed <- as.integer(runif(1) * 32768)
  }

  #-------------------- Main loop to make predictions --------------------#
    make_one_predict <- function(p)
    {
      if (!is.null(random_seed)) set.seed(p_seeds[p])

      lib_indices_p <- lib_indices
      lib_indices_p[p] <- FALSE
      libs <- which(lib_indices_p)

      # compute distances of temporal information
      q <- matrix(rep(vectors[p,], length(libs)), nrow = length(libs), byrow = T)
      distances <- sqrt(rowSums((vectors[libs,] - q) ^ 2))

      # compute temporal weights
      d_bar <- mean(distances, na.rm = TRUE)
      c_ws <- exp(-theta * distances / d_bar)

      if (regularized)
      {
        # do regularized S-map
        # Currently rely on "glmnet" package of R
        if (E == 1)
        {
          A <- cbind(vectors[libs,], 1)
          B <- cbind(target[libs])
          if(is.null(lambda)){
            # make prediction
            fit <- cv.glmnet(A, B, weights = c_ws, type.measure = "mae", alpha = alpha, family = "gaussian", nfolds = 10,
                             parallel = glmnet_parallel, intercept = TRUE)
            pred_val <- predict(fit, s = fit$lambda.1se, newx = matrix(c(vectors[p,], 1), nrow = 1))
            smap_coefficient_val <- matrix(t(coef(fit, s = fit$lambda.1se)), nrow = 1)[c(2, 1)]
          }else{
            # make prediction
            fit <- glmnet(A, B, weights = c_ws, alpha = alpha, family = "gaussian", lambda = lambda, intercept = TRUE)
            pred_val <- predict(fit, s = fit$lambda, newx = matrix(c(vectors[p,], 1), nrow = 1))
            smap_coefficient_val <- matrix(t(coef(fit, s = fit$lambda)), nrow = 1)[c(2, 1)]
          }
        } else {
          A <- cbind(vectors[libs,])
          B <- cbind(target[libs])
          if(is.null(lambda)){
            # make prediction
            fit <- cv.glmnet(A, B, weights = c_ws, type.measure = "mae", alpha = alpha, family = "gaussian", nfolds = 10,
                             parallel = glmnet_parallel, intercept = TRUE)
            pred_val <- predict(fit, s = fit$lambda.1se, newx = matrix(c(vectors[p,]), nrow = 1))
            smap_coefficient_val <- matrix(t(coef(fit, s = fit$lambda.1se)), nrow = 1)[c(2:(NCOL(vectors) + 1), 1)]
          }else{
            # make prediction
            fit <- glmnet(A, B, weights = c_ws, alpha = alpha, family = "gaussian", lambda = lambda, intercept = TRUE)
            pred_val <- predict(fit, s = fit$lambda, newx = matrix(c(vectors[p,]), nrow = 1))
            smap_coefficient_val <- matrix(t(coef(fit, s = fit$lambda)), nrow = 1)[c(2:(NCOL(vectors) + 1), 1)]
          }
        }
      } else {
        # do singular-value decomposition
        A <- cbind(vectors[libs,], 1) * c_ws
        A_svd <- svd(A)

        # remove singular values that are too small
        s <- A_svd$d
        s_inv <- diag(ifelse(s >= max(s) * 1e-5, 1 / s, 0))

        # perform back-substitute to solve        
        map <- A_svd$v %*% s_inv %*% t(A_svd$u) %*% (c_ws * target[libs])

        # make prediction
        pred_val <- sum(map * c(vectors[p,], 1))
        smap_coefficient_val <- t(map)
      }

      list(p = p, pred_val = pred_val, smap_coefficient_val = smap_coefficient_val)
    }
    
    # Perform parallel computing
    if (!no_parallel && foreach::getDoParWorkers() > 1)
    {
      `%dopar%` <- foreach::`%dopar%`
      all_predicted <- foreach::foreach(p = which(pred_indices), .packages = "glmnet") %dopar%
                                  make_one_predict(p)
      for (predicted in all_predicted)
      {
        pred_vals[predicted$p] <- predicted$pred_val
        smap_coefficient_vals[predicted$p,] <- predicted$smap_coefficient_val
      }
    } else {
      for (p in which(pred_indices))
      {
        predicted <- make_one_predict(p)
        pred_vals[p] <- predicted$pred_val
        smap_coefficient_vals[p,] <- predicted$smap_coefficient_val
      }
    }
  #-------------------- Main loop to finished --------------------#

  if (!is.null(random_seed)) set.seed(next_seed)

  # Add colnames for smap_coefficient_vals
  colnames(smap_coefficient_vals) <- c(paste0("c_", 1:NCOL(vectors)), "c_0")

  # return output & stats
  if (save_smap_coefficients)
  {
    return(list(pred = pred_vals,
                stats = compute_stats_SSR(target[pred_indices], pred_vals[pred_indices]),
                smap_coefficients = smap_coefficient_vals))
  } else {
    return(list(pred = pred_vals,
                stats = compute_stats_SSR(target[pred_indices], pred_vals[pred_indices])))
  }
}

