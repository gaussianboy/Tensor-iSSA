### Function from MoveHMM package (in order to calculate the turning angles later) ###
turnAngle <- function(x, y, z, LLangle) {
  if (all(x == y) | all(y == z)) {
    return(NA)
  }

  if (LLangle) {
    angle <- (bearing(x, y) - bearing(y, z)) / 180 * pi
  } else {
    v <- c(y[1] - x[1], y[2] - x[2])
    w <- matrix(c(z[, 1] - y[1], z[, 2] - y[2]), ncol = 2, byrow = FALSE)
    angle <- atan2(w[, 2], w[, 1]) - atan2(v[2], v[1])
  }

  for (i in 1:length(angle)) {
    while (angle[i] <= -pi) {
      angle[i] <- angle[i] + 2 * pi
    }
    while (angle[i] >= pi) {
      angle[i] <- angle[i] - 2 * pi
    }
  }

  return(angle)
}

### Function used to generate sequential movement decisions
simdata_generic <- function(parameter, case, center, Ts, seed, vzero, landscape) {
  loc <- matrix(NA, ncol = 2, nrow = Ts) # location of the animal
  colnames(loc) <- c("x", "y")
  loc <- as.data.frame(loc)
  # Note: need to convert to data frame so that plotting the points in raster works
  steps <- numeric(Ts) # save selected step length
  cells <- ncell(landscape) # number of cells in the raster
  # Setting seed
  set.seed(seed)
  # Initial 2 locations are random
  ext_hr <- 20
  loc[1, ] <- c(
    sample((center[1] - ext_hr):(center[1] + ext_hr), 1),
    sample((center[2] - ext_hr):(center[2] + ext_hr), 1)
  )

  loc[2, ] <- c(
    sample((center[1] - ext_hr):(center[1] + ext_hr), 1),
    sample((center[2] - ext_hr):(center[2] + ext_hr), 1)
  )
  #
  xy <- xyFromCell(landscape, 1:length(landscape$x1))

  #### Different movement kernel scenarios ####

  ## Gamma(bimodal-von Mises) ##
  if (case == "bimodal") {
    dist_mix <- UnivarMixingDistribution(distr::Gammad(shape = 10, scale = 1), distr::Gammad(shape = 50, scale = 1),
      mixCoeff = c(1 / 4, 3 / 4)
    )

    my_fun_sl <- distr::d(dist_mix)

    for (ts in 3:Ts) {
      # faster version using just vectors
      sl1 <- distanceFromPoints(landscape, loc[ts - 1, ]) # distance to next point
      # (set value of current location to a small number because otherwise gamma is NA (at zero))

      ta <- turnAngle(as.numeric(loc[ts - 2, ]), as.numeric(loc[ts - 1, ]), xy, LLangle = FALSE)

      timepoint <- ts - 2

      sl1[cellFromXY(landscape, loc[ts - 1, ])] <- vzero
      ta[cellFromXY(landscape, loc[ts - 1, ])] <- vzero


      # calculate kernel from this (dgamma)
      kern1 <- my_fun_sl(values(sl1)) *
        dvonmises(circular(ta), mu = circular(0), kappa = 1) /
        (values(sl1))

      # calculate hsf (habitat selection function)
      hsf <- exp(as.matrix(values(landscape)) %*% parameter$beta)
      # calculate step probability
      step.prob <- (kern1 * hsf) / sum(kern1 * hsf)
      # sample a new step
      temp <- sample(1:cells, 1, prob = step.prob)
      loc[ts, ] <- xyFromCell(landscape, temp)
      steps[ts - 1] <- sl1[cellFromXY(landscape, loc[ts, ])]
    }
    return(cbind(loc))
  } else if (case == "gamma_mises") {
    for (ts in 3:Ts) {
      # faster version using just vectors
      sl1 <- distanceFromPoints(landscape, loc[ts - 1, ]) # distance to next point

      ta <- turnAngle(as.numeric(loc[ts - 2, ]), as.numeric(loc[ts - 1, ]), xy, LLangle = FALSE)

      timepoint <- ts - 2

      sl1[cellFromXY(landscape, loc[ts - 1, ])] <- vzero
      ta[cellFromXY(landscape, loc[ts - 1, ])] <- vzero

      # calculate kernel from this (dgamma)
      kern1 <- dgamma(values(sl1), shape = 10, rate = 1) *
        dvonmises(circular(ta), mu = circular(0), kappa = 1) /
        (values(sl1))

      # calculate hsf (habitat selection function)
      hsf <- exp(as.matrix(values(landscape)) %*% parameter$beta)
      # calculate step probability
      step.prob <- (kern1 * hsf) / sum(kern1 * hsf)
      # sample a new step
      temp <- sample(1:cells, 1, prob = step.prob)
      loc[ts, ] <- xyFromCell(landscape, temp)
      steps[ts - 1] <- sl1[cellFromXY(landscape, loc[ts, ])]
    }

    return(cbind(loc))
  } else if (case == "weibull") {
    for (ts in 3:Ts) {
      # faster version using just vectors
      sl1 <- distanceFromPoints(landscape, loc[ts - 1, ]) # distance to next point
      # (set value of current location to a small number because otherwise gamma is NA (at zero))

      ta <- turnAngle(as.numeric(loc[ts - 2, ]), as.numeric(loc[ts - 1, ]), xy, LLangle = FALSE)

      timepoint <- ts - 2

      sl1[cellFromXY(landscape, loc[ts - 1, ])] <- vzero
      ta[cellFromXY(landscape, loc[ts - 1, ])] <- vzero


      # calculate kernel from this (weibull)
      kern1 <- dweibull(values(sl1), shape = 4, scale = 4) *
        dvonmises(circular(ta), mu = circular(0), kappa = 1) /
        (values(sl1))

      # calculate hsf (habitat selection function)
      hsf <- exp(as.matrix(values(landscape)) %*% parameter$beta)
      # calculate step probability
      step.prob <- (kern1 * hsf) / sum(kern1 * hsf)
      # sample a new step
      temp <- sample(1:cells, 1, prob = step.prob)
      loc[ts, ] <- xyFromCell(landscape, temp)
      steps[ts - 1] <- sl1[cellFromXY(landscape, loc[ts, ])]
    }
    return(cbind(loc))
  } else if (case == "copula") {
    # copula
    cop <- cyl_quadsec()
    marginal_1 <- list(name = "wrappedcauchy", coef = list(location = 0, scale = 0.3))
    marginal_2 <- list(name = "weibull", coef = list(shape = 3, scale = 7))

    set.seed(666)
    sample <- rjoint(10000, cop, marginal_1, marginal_2)

    marginal_1 <- fit_angle(theta = sample[, 1], parametric = FALSE)
    marginal_2 <- fit_steplength(x = sample[, 2], parametric = "lnorm")

    for (ts in 3:Ts) {
      # faster version using just vectors
      sl1 <- distanceFromPoints(landscape, loc[ts - 1, ]) # distance to next point
      # (set value of current location to a small number because otherwise gamma is NA (at zero))

      ta <- turnAngle(as.numeric(loc[ts - 2, ]), as.numeric(loc[ts - 1, ]), xy, LLangle = FALSE)

      timepoint <- ts - 2

      sl1[cellFromXY(landscape, loc[ts - 1, ])] <- vzero
      ta[cellFromXY(landscape, loc[ts - 1, ])] <- vzero


      # calculate kernel from this
      ta_sl <- cbind(ta, values(sl1))
      kern1 <- djoint(x = ta_sl, copula = cop, marginal_1 = marginal_1, marginal_2 = marginal_2)

      # calculate hsf (habitat selection function)
      hsf <- exp(as.matrix(values(landscape)) %*% parameter$beta)
      # calculate step probability
      step.prob <- (kern1 * hsf) / sum(kern1 * hsf)
      # sample a new step
      temp <- sample(1:cells, 1, prob = step.prob)
      loc[ts, ] <- xyFromCell(landscape, temp)
      steps[ts - 1] <- sl1[cellFromXY(landscape, loc[ts, ])]
    }
    return(cbind(loc))
  } else if (case == "copula_2") {
    # copula
    cop <- cyl_quadsec()
    marginal_1 <- list(name = "wrappedcauchy", coef = list(location = 0, scale = 0.5))
    marginal_2 <- list(name = "weibull", coef = list(shape = 2, scale = 9))
    
    set.seed(666)
    sample <- rjoint(10000, cop, marginal_1, marginal_2)
    
    marginal_1 <- fit_angle(theta = sample[, 1], parametric = FALSE)
    marginal_2 <- fit_steplength(x = sample[, 2], parametric = "lnorm")
    
    for (ts in 3:Ts) {
      # faster version using just vectors
      sl1 <- distanceFromPoints(landscape, loc[ts - 1, ]) # distance to next point
      # (set value of current location to a small number because otherwise gamma is NA (at zero))
      
      ta <- turnAngle(as.numeric(loc[ts - 2, ]), as.numeric(loc[ts - 1, ]), xy, LLangle = FALSE)
      
      timepoint <- ts - 2
      
      sl1[cellFromXY(landscape, loc[ts - 1, ])] <- vzero
      ta[cellFromXY(landscape, loc[ts - 1, ])] <- vzero
      
      
      # calculate kernel from this
      ta_sl <- cbind(ta, values(sl1))
      kern1 <- djoint(x = ta_sl, copula = cop, marginal_1 = marginal_1, marginal_2 = marginal_2)
      
      # calculate hsf (habitat selection function)
      hsf <- exp(as.matrix(values(landscape)) %*% parameter$beta)
      # calculate step probability
      step.prob <- (kern1 * hsf) / sum(kern1 * hsf)
      # sample a new step
      temp <- sample(1:cells, 1, prob = step.prob)
      loc[ts, ] <- xyFromCell(landscape, temp)
      steps[ts - 1] <- sl1[cellFromXY(landscape, loc[ts, ])]
    }
    return(cbind(loc))
  } 
}
