library("copula")

generate_synthetic_data <- function(copula, n, inverse_marginals = c(NA), 
  rounding_precision = c(NA), omission_probability = c(0), x_default = c(0), 
  x_lower = c(NA), x_upper = c(NA)) {
  
  # Initial processing of the parameter of the data generating algorithm
  d <- data.frame(cop=dim(copula))
  d$im <- length(inverse_marginals)
  d$rp <- length(rounding_precision)
  d$op <- length(omission_probability)
  d$xl <- length(x_lower)
  d$xu <- length(x_upper)
  d$xd <- length(x_default)
  
  # Validate the input of the inverse marginals.
  if (! d$im %in% c(1, d$cop)) {
    warning(paste0(
      "Count of inverse marginals ", d$im, " is not 1 or ", d$cop, ". ",
      "Missing functions are set to uniform marginals. ", 
      "Additional functions provided will be ignored."
    ));
    inverse_marginals <- 
      c(inverse_marginals, rep(NA, max(d$cop - d$im,0)))[1:d$cop]
  } else {
    inverse_marginals <- rep(inverse_marginals, d$cop / d$im)
  }

  # Validate the input of the rounding precisions.
  if (! d$rp %in% c(1, d$cop)) {
    warning(paste0(
      "Count of rounding precisions ", d$rp, " is not 1 or ", d$cop, ". ",
      "Missing parameter will be set to maximum precision. ", 
      "Additional parameter will be ignored."
    ));
    rounding_precision <- 
      c(rounding_precision, rep(NA, max(d$cop - d$rp,0)))[1:d$cop]
  } else {
    rounding_precision <- rep(rounding_precision, d$cop / d$rp)
  }
  
  # Validate the input of the omission probabilities.
  if (! d$op %in% c(1, d$cop)) {
    warning(paste0(
      "Count of omission probabilities ", d$op, " is not 1 or ", d$cop, ". ",
      "Missing parameter will be set to zero. ", 
      "Additional parameter will be ignored."
    ));
    omission_probability <- 
      c(omission_probability, rep(0, max(d$cop - d$op,0)))[1:d$cop]
  } else {
    omission_probability <- rep(omission_probability, d$cop / d$op)
  }
  
  # Validate the input of the lower boundary vector.
  if (! d$xl %in% c(1, d$cop)) {
    warning(paste0(
      "Dimension of the lower boundary ", d$xl, " is not 1 or ", d$cop, ". ",
      "Missing parameter will be set to NA / -infinity. ", 
      "Additional parameter will be ignored."
    ));
    x_lower <- c(x_lower, rep(NA, max(d$cop - d$xl,0)))[1:d$cop]
  } else {
    x_lower <- rep(x_lower, d$cop / d$xl)
  }
  
  # Validate the input of the upper boundary vector.
  if (! d$xu %in% c(1, d$cop)) {
    warning(paste0(
      "Dimension of the upper boundary ", d$xu, " is not 1 or ", d$cop, ". ",
      "Missing parameter will be set to NA / +infinity. ",  
      "Additional parameter will be ignored."
    ));
    x_upper <- c(x_upper, rep(NA, max(d$cop - d$xu,0)))[1:d$cop]
  } else {
    x_upper <- rep(x_upper, d$cop / d$xu)
  }
  
  # Validate the input of the default value boundary vector.
  if (! d$xd %in% c(1, d$cop)) {
    warning(paste0(
      "Dimension of the default vector ", d$xd, " is not 1 or ", d$cop, ". ",
      "Missing parameter will be set to zero. ",
      "Additional parameter will be ignored."
    ));
    x_default <- c(x_default, rep(0, max(d$cop - d$xd,0)))[1:d$cop]
  } else {
    x_default <- rep(x_default, d$cop / d$xd)
  }
  
  # 1. Generate a sample u ~ copula of size n and initialize the two output
  # variables x (modified) and y (unmodified).
  u <- copula::rCopula(n, copula)
  x <- y <- u
  
  # 2. Quantile Transformation (ITS) w. r. t. the inverse marginals.
  for (i in 1:d$cop) { if (is.function(inverse_marginals[[i]])) {
    x[, i] <- inverse_marginals[[i]](u[, i])
  }}
  
  # 3. Rounding of data.
  for (i in 1:d$cop) { if (!is.na(rounding_precision[i])) {
    x[, i] <- round(x[, i] / rounding_precision[i]) * rounding_precision[i]
  }}
  
  # 4. Random omission of entries.
  for (i in 1:d$cop) { if (!is.na(omission_probability[i])) {
    if (0 < omission_probability[i] && omission_probability[i] < 1) {
      b <- rbinom(n,1,omission_probability[i])
      x[, i] <- (1 - b) * x[, i] + b * x_default[i]
    }
  }}
  
  # 5. Limiting of values w. r. t. the boundaries.
  for (i in 1:d$cop) { if (!is.na(x_lower[i])) {
    x[, i][x[, i] <= x_lower[i]] <- x_lower[i]
  }}
  for (i in 1:d$cop) { if (!is.na(x_upper[i])) {
    x[, i][x[, i] >= x_upper[i]] <- x_upper[i]
  }}

  # Output of the generated data.
  out <- list()
  out$modified <- x
  out$unmodified <- y
  
  return(out)
}

