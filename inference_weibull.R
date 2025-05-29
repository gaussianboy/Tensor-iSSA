# Parallel kernels specification (Using 50 cores)
library(parallel)
cl <- makeCluster(50)

# Specifying working directory
setwd("PROVIDE YOUR WORKING DIRECTORY")
load("gam_simulation_weibull.RData")

# Exporting objects to the nodes
clusterExport(cl = cl, varlist = c("my_data", "landscape"), envir = .GlobalEnv)

# Iterative function
iteration_list <- as.list(1:100)
gam_objects <- parLapply(cl, iteration_list, function(i) {
  library(dplyr)
  library(tidyverse)
  library(lubridate)
  library(raster)
  library(rgdal)
  library(sf)
  library(rgeos)
  library(sp)
  library(deldir)
  library(maptools)
  library(spatstat)
  library(circular)
  library(mgcv)
  library(circular)
  library(copula)
  library(cylcop)
  library(amt)
  library(distr)

  ### Animal tracks data preparation ###
  my_df <- my_data[[i]]
  my_df$t <- as.POSIXlt(my_df$t, format = "%Y-%m-%d %H:%M:%OS")
  my_df <- my_df %>% dplyr::filter(!(is.na(t)))
  my_df <- my_df %>% arrange(t)

  my_df <- make_track(my_df, .x = x, .y = y, .t = t)

  my_df <- my_df %>% track_resample(rate = minutes(6), tolerance = minutes(2))

  df1 <- as.data.frame(my_df)

  # Creating step lengths
  xyset <- cbind(df1$x, df1$y)
  step_lengths <- sqrt(rowSums((xyset[-1, ] - xyset[-nrow(xyset), ])^2))
  step_lengths <- ifelse(step_lengths == 0, 1e-6, step_lengths)

  # Calculating maximum step length
  max_sl <- max(step_lengths)

  #### NON-UNIFORM INTEGRATION POINTS ####
  # Sampling random steps from an empirical movement kernel
  issa_df_MC <- steps_by_burst(my_df)
  issa_df_MC <- issa_df_MC %>% random_steps(n_control = 300)

  issa_df_MC <- issa_df_MC %>%
    extract_covariates(landscape)

  issa_df_MC$count <- as.numeric(issa_df_MC$case_)

  issa_df_MC$my_times <- 1

  issa_model_classic <- gam(
    cbind(my_times, step_id_) ~ -1 +
      x1 +
      x2 +
      cen +
      sl_ +
      log(sl_) +
      cos(ta_),
    data = issa_df_MC,
    family = cox.ph,
    weights = count
  )

  #### UNIFORMLY SAMPLED INTEGRATION POINTS ####
  issa_df <- steps_by_burst(my_df)

  set.seed(999)

  issa_df <- issa_df %>% random_steps(
    n_control = 300,
    rand_sl = sqrt(random_numbers(make_unif_distr(min = 0, max = max_sl^2), n = 1e5)),
    angle = 0,
    rand_ta = random_numbers(make_unif_distr(), n = 1e+05)
  )

  issa_df <- issa_df %>%
    extract_covariates(landscape)

  issa_df$count <- as.numeric(issa_df$case_)

  # 1) iSSA model with uniform integration points
  issa_df$my_times <- 1

  issa <- gam(
    cbind(my_times, step_id_) ~ -1 +
      x1 +
      x2 +
      cen +
      sl_ +
      log(sl_) +
      cos(ta_),
    data = issa_df,
    family = cox.ph,
    weights = count
  )

  # 2) iSSA model without x2
  issa_x2 <- gam(
    cbind(my_times, step_id_) ~ -1 +
      x1 +
      cen +
      sl_ +
      log(sl_) +
      cos(ta_),
    data = issa_df,
    family = cox.ph,
    weights = count
  )

  # 3) iSSA model without x2 but accounting for missing spatial variation
  issa_x2_spatial <- gam(
    cbind(my_times, step_id_) ~ -1 +
      x1 +
      cen +
      sl_ +
      log(sl_) +
      cos(ta_) +
      s(x2_, y2_, k = 100),
    data = issa_df,
    family = cox.ph,
    weights = count
  )

  # 4) GAM-iSSA model using x1, x2 and cen in the linear predictor
  gam_issa <- gam(
    cbind(my_times, step_id_) ~ -1 +
      x1 +
      x2 +
      cen +
      te(sl_, ta_, bs = c("bs", "cc")),
    data = issa_df,
    family = cox.ph,
    weights = count
  )

  # 5) GAM-iSSA without x2
  gam_issa_x2 <- gam(
    cbind(my_times, step_id_) ~ -1 +
      x1 +
      cen +
      te(sl_, ta_, bs = c("bs", "cc")),
    data = issa_df,
    family = cox.ph,
    weights = count
  )

  # 6) GAM-iSSA without x2 but accounting for it
  gam_issa_x2_spatial <- gam(
    cbind(my_times, step_id_) ~ -1 +
      x1 +
      cen +
      te(sl_, ta_, bs = c("bs", "cc")) +
      s(x2_, y2_, k = 100),
    data = issa_df,
    family = cox.ph,
    weights = count
  )

  ### MSE calculation based on the initial Equation of Forester et al (2009) ###
  times <- as.list(as.numeric(issa_df$step_id_))

  ## Sampling different integration points for test set
  issa_df_test <- steps_by_burst(my_df)

  set.seed(666)

  issa_df_test <- issa_df_test %>% random_steps(
    n_control = 300,
    rand_sl = sqrt(random_numbers(make_unif_distr(min = 0, max = max_sl^2), n = 1e5)),
    angle = 0,
    rand_ta = random_numbers(make_unif_distr(), n = 1e+05)
  )


  issa_df_test <- issa_df_test %>%
    extract_covariates(landscape)

  issa_df_test$count <- as.numeric(issa_df_test$case_)
  issa_df_test$my_times <- 1


  MSE_fun <- function(time, model) {
    quadrature_points <- issa_df_test %>% filter(step_id_ == time)
    ta_sl <- cbind(quadrature_points$ta_, quadrature_points$sl_)


    ssf_values <- exp(1.5 * quadrature_points$x1 - quadrature_points$x2 -
      0.02 * quadrature_points$cen) * dweibull(ta_sl[, 2], shape = 4, scale = 4) *
      dvonmises(circular(ta_sl[, 1]), mu = circular(0), kappa = 1)

    real_lik <- ssf_values / sum(ssf_values)

    pred_ssf_values <- exp(predict(model, newdata = quadrature_points))
    pred_lik <- pred_ssf_values / sum(pred_ssf_values)

    SE <- mean((real_lik - pred_lik)^2)

    return(SE)
  }

  SE_model_iSSA_classic <- mean(sapply(times, FUN = MSE_fun, model = issa_model_classic), na.rm = TRUE)
  SE_model_iSSA <- mean(sapply(times, FUN = MSE_fun, model = issa), na.rm = TRUE)
  SE_model_iSSA_x2 <- mean(sapply(times, FUN = MSE_fun, model = issa_x2), na.rm = TRUE)
  SE_model_iSSA_x2_spatial <- mean(sapply(times, FUN = MSE_fun, model = issa_x2_spatial), na.rm = TRUE)
  SE_model_gam_iSSA <- mean(sapply(times, FUN = MSE_fun, model = gam_issa), na.rm = TRUE)
  SE_model_gam_iSSA_x2 <- mean(sapply(times, FUN = MSE_fun, model = gam_issa_x2), na.rm = TRUE)
  SE_model_gam_iSSA_x2_spatial <- mean(sapply(times, FUN = MSE_fun, model = gam_issa_x2_spatial), na.rm = TRUE)


  return(list(
    issa_model_classic, issa, issa_x2, issa_x2_spatial,
    gam_issa, gam_issa_x2, gam_issa_x2_spatial,
    SE_model_iSSA_classic, SE_model_iSSA, SE_model_iSSA_x2, SE_model_iSSA_x2_spatial,
    SE_model_gam_iSSA, SE_model_gam_iSSA_x2, SE_model_gam_iSSA_x2_spatial
  ))
})

save.image(file = "inference_weibull_mises.RData")

parallel::stopCluster(cl)
