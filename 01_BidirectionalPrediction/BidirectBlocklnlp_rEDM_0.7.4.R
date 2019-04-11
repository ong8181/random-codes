####
#### Function: Bidirectional forecasting
#### 2018.12.23 Ushio
####

#---------- History ----------#
# v1, 2018.12.23, implemented in R3.4.4
# v2, 2019.3.26, "E" in "xxxx_backword" has been revised as "E*tau" in R3.5.2
# v3, 2019.3.28, "compute_stats" has been revised as "compute_stats2" to avoid conflict
# v4, 2019.4.2, "lib = c(1, NROW(time_series))" in simplex_forward and simplex_backward is revised to "lib = lib"
#     2019.4.2, "pred = lib" in simplex_forward and simplex_backward is revised to "pred = pred"
# v5, 2019.4.3, Re-write "complete_case_only == FALSE". Bugs fixed
#     2019.4.3, Output style modified ("tp", "theta", "embedding" added)

#---------- Usage ----------#
## CAUTION: Currently I assume the case only when tp = 1
## The function may work for different tp and tau, but need your own validation!

require(rEDM)

# Define function
bidirect_block_lnlp <- function(time_series,
                                complete_case_only = FALSE, # if TRUE, averaged predictions will be calculated only for complete pairs
                                lib = c(1, NROW(time_series)),
                                pred = lib,
                                tau = 1,
                                tp_forward = 1,
                                tp_backward = -NCOL(time_series)*tau - tp_forward + 1,
                                theta = 8,
                                method = "s-map",
                                num_neighbors = 0,
                                stats_only = TRUE,
                                silent = TRUE,
                                exclusion_radius = NULL,
                                epsilon = NULL){
  # Perform forward and backward simplex projection
  pred_forward <- block_lnlp(time_series, tp = tp_forward,
                             lib = lib, pred = pred,
                             method = method, num_neighbors = num_neighbors,
                             stats_only = F, silent = silent, theta = theta)
  pred_backward <- block_lnlp(time_series, tp = tp_backward,
                              lib = lib, pred = pred,
                              method = method, num_neighbors = num_neighbors,
                              stats_only = F, silent = silent, theta = theta)
  
  # Averaging forward and backward predictions
  if(complete_case_only){
    # Match time indices of the model output
    time_shared <- as.numeric(na.omit(intersect(pred_forward$model_output[[1]]$time,
                                                pred_backward$model_output[[1]]$time)))
    # Extract model output
    time_forward <- pred_forward$model_output[[1]][match(time_shared, pred_forward$model_output[[1]]$time),]
    time_backward <- pred_backward$model_output[[1]][match(time_shared, pred_backward$model_output[[1]]$time),]
    
    averaged_prediction <- (time_forward$pred + time_backward$pred)/2
    
    pred_bidirect <- data.frame(time = time_forward$time,
                                obs = time_forward$obs,
                                pred = averaged_prediction)
  }else{
    # Match time indices of the model output
    pred_forward_nona <- pred_forward$model_output[[1]][!is.na(pred_forward$model_output[[1]]$time),]
    pred_backward_nona <- pred_backward$model_output[[1]][!is.na(pred_backward$model_output[[1]]$time),]
    time_covered <- sort(unique(c(pred_forward_nona$time, pred_backward_nona$time)))
    
    # Identify time index of simplex_forward and simplex_backward
    time_forward_id <- match(pred_forward_nona$time, time_covered)
    time_backward_id <- match(pred_backward_nona$time, time_covered)
    
    # Create simplex_bidirect object
    pred_bidirect0 <- data.frame(time = time_covered)
    pred_bidirect0$forward_pred <- pred_bidirect0$forward_obs <- pred_bidirect0$forward_time <- NaN
    pred_bidirect0$backward_pred <- pred_bidirect0$backward_obs <- pred_bidirect0$backward_time <- NaN
    pred_bidirect0[time_forward_id,]$forward_time <- pred_forward_nona$time
    pred_bidirect0[time_forward_id,]$forward_obs <- pred_forward_nona$obs
    pred_bidirect0[time_forward_id,]$forward_pred <- pred_forward_nona$pred
    pred_bidirect0[time_backward_id,]$backward_time <- pred_backward_nona$time
    pred_bidirect0[time_backward_id,]$backward_obs <- pred_backward_nona$obs
    pred_bidirect0[time_backward_id,]$backward_pred <- pred_backward_nona$pred
    
    # Sanity check
    if(F){
      sum(pred_bidirect0$forward_time - pred_bidirect0$backward_time, na.rm = T) == 0
      sum(pred_bidirect0$time - pred_bidirect0$backward_time, na.rm = T) == 0
      sum(pred_bidirect0$forward_obs - pred_bidirect0$backward_obs, na.rm = T) == 0
    }
    
    # Calculate average prediction
    obs_combined <- apply(cbind(pred_bidirect0$forward_obs, pred_bidirect0$backward_obs), 1, function(x) mean(x, na.rm = T))
    averaged_prediction <- apply(cbind(pred_bidirect0$forward_pred, pred_bidirect0$backward_pred), 1, function(x) mean(x, na.rm = T))
    
    # Summary
    pred_bidirect <- data.frame(time = pred_bidirect0$time,
                                   obs = obs_combined,
                                   pred = averaged_prediction)
    }
  
  # Collecting parameters
  params_forward <- pred_forward[1:16]
  params_backward <- pred_backward[1:16]
  params_all <- cbind(data.frame(method = c("forward_prediction","backward_prediction")),
                     rbind(params_forward, params_backward))
  
  # Calculate prediction accuracy
  stats_forward <- pred_forward[,c("num_pred", "embedding", "theta", "tp", "rho", "mae", "rmse")]
  stats_backward <- pred_backward[,c("num_pred", "embedding", "theta", "tp", "rho", "mae", "rmse")]
  stats_bidirection <- compute_stats2(pred_bidirect$obs, pred_bidirect$pred)
  stats_bidirection$embedding <- pred_forward$embedding
  stats_bidirection$theta <- pred_forward$theta
  stats_bidirection$tp <- sprintf("%s, %s", pred_forward$tp, pred_backward$tp)
  stats_bidirection <- stats_bidirection[,c("num_pred", "embedding", "theta", "tp", "rho", "mae", "rmse")]
  stats_all <- cbind(data.frame(method = c("bidirect_prediction","forward_prediction","backward_prediction")),
                     rbind(stats_bidirection, stats_forward, stats_backward))      

  # Return output
  if(stats_only){
    pred_bidirect_all <- data.frame(stats_all)
  }else{
    pred_bidirect_all <- list(params = params_all, model_output = pred_bidirect, stats = stats_all)
  }
  return(pred_bidirect_all)
}

compute_stats2 <- function(obs, pred){
  # computes performance metrics for how well predictions match observations
  # obs = vector of observations
  # pred = vector of prediction
  
  num_pred <- sum(is.finite(obs) & is.finite(pred))
  rho <- cor(obs, pred, use = "pairwise.complete.obs")
  mae <- mean(abs(obs-pred), na.rm = TRUE)
  rmse <- sqrt(mean((obs-pred)^2, na.rm = TRUE))
  return(data.frame(num_pred = num_pred, rho = rho, mae = mae, rmse = rmse))
}

# Calculate best theta using bidirectional block_lnlp
bestT_bidirect_block_lnlp <- function(time_series,
                                      theta_range = c(0, 1e-04, 3e-04,
                                                      0.001, 0.003, 0.01, 0.03,
                                                      0.1, 0.3, 0.5, 0.75,
                                                      1, 1.5, 2, 3, 4, 6, 8),
                                      lib = c(1, NROW(time_series)),
                                      tp_forward = 1,
                                      tp_backward =  -NCOL(time_series)*tau - tp_forward + 1,
                                      tau = 1,
                                      complete_case_only = F,
                                      criteria = "rmse",
                                      show_fig = F,
                                      save_stats = F){
  pred_res <- bidirect_block_lnlp(time_series, lib = lib, tp_forward = tp_forward, tp_backward = tp_backward, stats_only = T, complete_case_only = complete_case_only, tau = tau, theta = theta_range[1])[1,]
  for(i in 2:length(theta_range)) pred_res[i,] <- bidirect_block_lnlp(time_series, lib = lib, tp_forward = tp_forward, tp_backward = tp_backward, stats_only = T, complete_case_only = complete_case_only, tau = tau, theta = theta_range[i])[1,]
  best_T <- theta_range[which.min(pred_res[,criteria])]
  
  if(show_fig) plot(theta_range, pred_res[,criteria], type = "b", ylab = criteria, xlab = expression(theta)) 
  if(save_stats){
    all_res <- list(E = best_T, stats = pred_res)
  }else{
    all_res <- best_T
  }
  return(all_res)
}