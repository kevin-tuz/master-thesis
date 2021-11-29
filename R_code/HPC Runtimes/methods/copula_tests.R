# --------------------------------------------------------------------------- #
# Bivariate radial symmetry test - Genest et Nešlehová (2014)
# Implementation based on the MATLAB code by J.F. Quessy

radSymTest_genestEtNeslehova <- function(X, ties.method="max", b=1, N=1000) {
  # Calculate the pseudo observations U w.r.t. to the sample X. Measure the 
  # sample size n and construct the size-adapted binwidth bn.
  U <- copula::pobs(X, ties.method = ties.method);
  n <- nrow(X);
  bn <- b / sqrt(n);
  
  # Calculation of the test statistic Rn.
  E <- matrix(nrow = n, ncol = n);
  for (i.1 in 1:n) { for (i.2 in 1:n) {
    E[i.1, i.2] <- prod(U[i.1, ] <= U[i.2, ]) - prod(1 - U[i.1, ] <= U[i.2, ]);
  }}
  Rn <- sum(E %*% t(E)) / (n^2);
  
  # Computation of the multiplier bootstrap replicates (deterministic part).
  dCdu1 <- rep(0, n);
  dCdu2 <- rep(0, n);
  
  for (i in 1:n) { 
    C.upr.1 <- mean((U[, 1] <= (U[i, 1] + bn)) * (U[, 2] <= U[i, 2]));
    C.lwr.1 <- mean((U[, 1] <= (U[i, 1] - bn)) * (U[, 2] <= U[i, 2]));
    C.upr.2 <- mean((U[, 1] <= U[i, 1]) * (U[, 2] <= (U[i, 2] + bn)));
    C.lwr.2 <- mean((U[, 1] <= U[i, 1]) * (U[, 2] <= (U[i, 2] - bn)));
    
    dCdu1[i] <- (C.upr.1 - C.lwr.1) / (2 * bn);
    dCdu2[i] <- (C.upr.2 - C.lwr.2) / (2 * bn);
  }
  
  E <- matrix(nrow = n, ncol = n);
  for (i.1 in 1:n) { for (i.2 in 1:n) {
    Pn <- prod(U[i.1, ] <= U[i.2, ]) - prod((1 - U[i.1, ]) <= U[i.2, ]);
    Pn.1 <- (U[i.1, 1] <= U[i.2, 1]) - (1 - U[i.1, 1] <= U[i.2, 1]);
    Pn.2 <- (U[i.1, 2] <= U[i.2, 2]) - (1 - U[i.1, 2] <= U[i.2, 2]);
    dCdu.Pn.1 <- dCdu1[i.2] * Pn.1;
    dCdu.Pn.2 <- dCdu2[i.2] * Pn.2;
    E[i.1, i.2] <- Pn - dCdu.Pn.1 - dCdu.Pn.2;
  }}
  E2 <- E %*% t(E);
  
  # Computation of the multiplier bootstrap replicates (random part).
  Rn_bootstrap = rep(0, N);
  for (k in 1:N) {
    xi <- rexp(n);
    Xi <- xi / mean(xi) - 1;
    Rn_bootstrap[k] <- t(Xi) %*% E2 %*% Xi / (n^2);
  }
  
  # Output of the results
  structure(class = "htest", list(
    statistic = c(statistic = Rn),
    p.value = (sum(Rn_bootstrap > Rn) + 0.5) / (N + 1)
  ));
}

# --------------------------------------------------------------------------- #
# Bivariate exchangeability test - Genest, Nešlehová et Quessy (2012)

exchTest_genestNeslehovaQuessy <- function(X, ties.method="max", b=1, N=1000) {
  # Calculate the pseudo observations U w.r.t. to the sample X. Measure the 
  # sample size and construct the size-adapted binwidth bn.
  U <- copula::pobs(X, ties.method = ties.method);
  n <- nrow(X);
  bn <- b / sqrt(n);
  
  # Calculation of the statistic An.
  E <- matrix(nrow = n, ncol = n);
  for (i.1 in 1:n) { for (i.2 in 1:n) {
    E[i.1, i.2] <- prod(U[i.1,] <= U[i.2, 1:2]) - prod(U[i.1,] <= U[i.2, 2:1]);
  }}
  An <- sum(E %*% t(E)) / (n^2);
  
  # Computation of the multiplier bootstrap replicates (deterministic part).
  dCdu1 <- rep(0, n);
  dCdu2 <- rep(0, n);
  
  for (i in 1:n) { 
    C.upr.1 <- mean((U[, 1] <= (U[i, 1] + bn)) * (U[, 2] <= U[i, 2]));
    C.lwr.1 <- mean((U[, 1] <= (U[i, 1] - bn)) * (U[, 2] <= U[i, 2]));
    C.upr.2 <- mean((U[, 1] <= U[i, 1]) * (U[, 2] <= (U[i, 2] + bn)));
    C.lwr.2 <- mean((U[, 1] <= U[i, 1]) * (U[, 2] <= (U[i, 2] - bn)));
    
    dCdu1[i] <- (C.upr.1 - C.lwr.1) / (2 * bn);
    dCdu2[i] <- (C.upr.2 - C.lwr.2) / (2 * bn);
  }
  
  E <- matrix(nrow = n, ncol = n);
  for (i.1 in 1:n) { for (i.2 in 1:n) {
    Pn <- prod(U[i.1, ] <= U[i.2, 1:2]) - prod(U[i.1, ] <= U[i.2, 2:1]);
    Pn.1 <- (U[i.1, 1] <= U[i.2, 1]) - (U[i.1, 2] <= U[i.2, 1]);
    Pn.2 <- (U[i.1, 2] <= U[i.2, 2]) - (U[i.1, 1] <= U[i.2, 2]);
    dCdu.Pn.1 <- dCdu1[i.2] * Pn.1;
    dCdu.Pn.2 <- dCdu2[i.2] * Pn.2;
    E[i.1, i.2] <- Pn - dCdu.Pn.1 - dCdu.Pn.2;
  }}
  E2 <- E %*% t(E);
  
  # Computation of the multiplier bootstrap replicates (random part).
  An_bootstrap = rep(0, N);
  for (k in 1:N) {
    xi <- rexp(n);
    Xi <- xi / mean(xi) - 1;
    An_bootstrap[k] <- t(Xi) %*% E2 %*% Xi / (n^2);
  }
  
  # Output of the results
  structure(class = "htest", list(
    statistic = c(statistic = An),
    p.value = (sum(An_bootstrap > An) + 0.5) / (N + 1)
  ));
}

# --------------------------------------------------------------------------- #
# Bivariate interval censored MPLE - Implementation based on Li et al. (2020),
# under modified notations.

censoredMPLE <- function(
  copula, method, X, start = NULL, lower = NULL, upper = NULL, 
  optim.control = list(maxit = 1000), estimate.variance = NA, 
  hideWarnings = TRUE, bound.eps = .Machine$double.eps^0.5
) {
  # Model fitting via maximum pseudo likelihood estimation based on a censored 
  # log pseudo likelihood function CensorLogLikelihood.
  
  U <- copula::pobs(X);
  
  # Validation of inputs.
  stopifnot(is.numeric(d <- ncol(U)), d >= 2);
  if (copula@dimension != d) { 
    stop("The copula and sample dimensions are not equal.");
  }
  if (is.null(start)) { 
    start <- copula:::fitCopStart(copula, U);
  }
  if (any(is.na(start))) { 
    stop("The parameter 'start' contains NA values.");
  }

  if (length(copula@parameters) != length(start)) {
    stop(paste0(
      "The vector 'start' (", length(start),") is not feasible. A vector of ",
      "length ", length(copula@parameters), " for the copula parameter is ",
      "required."
    ));
  }
  
  # Setup of the optimization control.
  control <- c(optim.control, fnscale = -1);
  control <- control[!vapply(control, is.null, NA)];
  if (!is.null(optim.control[[1]])) { control <- c(control, optim.control); }
  
  meth.has.bounds <- method %in% c("Brent", "L-BFGS-B");
  lower <- ifelse (meth.has.bounds, copula@param.lowbnd + bound.eps, -Inf );
  upper <- ifelse (meth.has.bounds, copula@param.upbnd - bound.eps, Inf );
  (if (hideWarnings) { suppressWarnings } else { identity })(
    fit <- optim(
      start, CensorLogLikelihood, lower = lower, upper = upper,
      method = method, copula = copula, X = X, control = control
    )
  )
  return(fit);
}

CensorLogLikelihood <- function(param, X, copula) {
  copula@parameters <- param;
  duplicatesIndex <- function(X) {
    return <- duplicated(X) | duplicated(X, fromLast=TRUE);
  }
  maxRanks <- apply(X, 2, rank, ties.method = "max") / (nrow(X) + 1);
  minRanks <- apply(X, 2, rank, ties.method = "min") / (nrow(X) + 1);
  censored <- apply(maxRanks, 2, duplicatesIndex);
  n <- nrow(X)
  
  # t = tied => censorship true
  # f = censorship false
  tt <- ff <- tf <- ft <- NULL;
  TT <- FF <- TF <- FT <- 1;
  
  for (i in 1:n) {
    c <- censored[i, ];
    minR <- minRanks[i, ];
    maxR <- maxRanks[i, ];
    
    # Check the type of censoring intervals required.
    if (!c[1] & !c[2]) { ff <- rbind(ff, c(minR, maxR)); } 
    else if (c[1] & c[2]) { tt <- rbind(tt, c(minR, maxR)); } 
    else if (!c[1] & c[2]) { ft <- rbind(ft, c(minR, maxR)); }
    else if (c[1] & !c[2]) { tf <- rbind(tf, c(minR, maxR)); }
  }
  
  # Likelihood terms for the four cases stated above.
  if (!is.null(tt)) {
    TT <- copula::pCopula(tt[, 3:4],copula) + 
      copula::pCopula(tt[, 1:2],copula) -
      copula::pCopula(cbind(tt[, 1], tt[, 4]), copula) - 
      copula::pCopula(cbind(tt[, 3], tt[, 2]), copula); 
  }
  if (!is.null(ff)) { 
    FF <- copula::dCopula(ff[, 1:2], copula);
  } 
  if (!is.null(tf)) { 
    TF <- copula:::dCdu(copula, matrix(tf[, 3:4], ncol=2))[, 2] - 
      copula:::dCdu(copula, matrix(tf[, 1:2], ncol=2))[, 2]; 
  }
  if (!is.null(ft)) {
    FT <- copula:::dCdu(copula, matrix(ft[, 3:4], ncol=2))[, 1] -
      copula:::dCdu(copula, matrix(ft[, 1:2], ncol=2))[, 1]; 
  }
  
  # Return the censored log pseudo likelihood
  return((sum(log(TT)) + sum(log(TF)) + sum(log(FT)) + sum(log(FF))))
}

# --------------------------------------------------------------------------- #
# Bivariate interval censored goodness of fit test.

gofTest_censored_normalCopula <- function(X, ties.method="max", N=250, start=NULL) {
  copula <- copula::normalCopula();
  U <- copula::pobs(X);
  n <- nrow(X);
  
  F.1.inv <- function(u) { return(quantile(U[, 1], u, type = 1, names = FALSE)); }
  F.2.inv <- function(u) { return(quantile(U[, 2], u, type = 1, names = FALSE)); }
  
  # First fit for the test statistic calcula
  theta_n <- copula@parameters;
  try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X, start)$par }, silent = TRUE);
  
  # Stop, if the fit fails.
  if (any(is.na(theta_n))) { 
    return(structure(class="htest", list(statistic=NA, p.value=NA))); 
    stop();
  }
  copula.theta_n <- copula::normalCopula(dim=2, theta_n);
  
  # Calculate the main test statistic
  eval.fitCopula <- copula::pCopula(U, copula.theta_n);
  eval.empCopula <- copula::C.n(U,U);
  
  Gn <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE); 
  
  # Estimate the p value with bootstrap replicates:
  Gn_bootstrap = rep(0, N);
  for (k in 1:N) {
    U_bootstrap <- copula::rCopula(n, copula.theta_n);
    X_bootstrap <- U_bootstrap;
    X_bootstrap[, 1] <- F.1.inv(U_bootstrap[, 1]);
    X_bootstrap[, 2] <- F.2.inv(U_bootstrap[, 2]);
    
    theta_n <- copula@parameters;
    try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X_bootstrap, start)$par }, silent = TRUE);
    
    if (any(is.na(theta_n))) { 
      Gn_bootstrap[k] <- NA;
    } else {
      copula.theta_n <- copula::normalCopula(dim=2, theta_n);
      eval.fitCopula <- copula::pCopula(U_bootstrap, copula.theta_n);
      eval.empCopula <- copula::C.n(U_bootstrap, U_bootstrap);
      Gn_bootstrap[k] <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE);
    }
  }
  
  # Output of the results
  structure(class = "htest", list(
    statistic = c(statistic = Gn),
    p.value = (sum(Gn_bootstrap > Gn, na.rm = TRUE) + 0.5) / (N + 1)
  ));
}

gofTest_censored_tCopula <- function(X, ties.method="max", N=250, start=NULL) {
  copula <- copula::tCopula();
  U <- copula::pobs(X);
  n <- nrow(X);
  
  F.1.inv <- function(u) { return(quantile(U[, 1], u, type = 1, names = FALSE)); }
  F.2.inv <- function(u) { return(quantile(U[, 2], u, type = 1, names = FALSE)); }
  
  # First fit for the test statistic calcula
  theta_n <- copula@parameters;
  try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X, start)$par }, silent = TRUE);
  
  # Stop, if the fit fails.
  if (any(is.na(theta_n))) { 
    return(structure(class="htest", list(statistic=NA, p.value=NA))); 
    stop();
  }
  
  copula.theta_n <- copula::tCopula(dim=2, theta_n[1], df = 4);
  
  # Calculate the main test statistic
  eval.fitCopula <- copula::pCopula(U, copula.theta_n);
  eval.empCopula <- copula::C.n(U,U);
  
  Gn <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE); 
  
  # Estimate the p value with bootstrap replicates:
  Gn_bootstrap = rep(0, N);
  for (k in 1:N) {
    U_bootstrap <- copula::rCopula(n, copula.theta_n);
    X_bootstrap <- U_bootstrap;
    X_bootstrap[, 1] <- F.1.inv(U_bootstrap[, 1]);
    X_bootstrap[, 2] <- F.2.inv(U_bootstrap[, 2]);
    
    theta_n <- copula@parameters;
    try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X_bootstrap, start)$par }, silent = TRUE);
    
    if (any(is.na(theta_n))) { 
      Gn_bootstrap[k] <- NA;
    } else {
      copula.theta_n <- copula::tCopula(dim=2, theta_n[1], df = 4);
      eval.fitCopula <- copula::pCopula(U_bootstrap, copula.theta_n);
      eval.empCopula <- copula::C.n(U_bootstrap, U_bootstrap);
      Gn_bootstrap[k] <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE);
    }
  }
  
  # Output of the results
  structure(class = "htest", list(
    statistic = c(statistic = Gn),
    p.value = (sum(Gn_bootstrap > Gn, na.rm = TRUE) + 0.5) / (N + 1)
  ));
}

gofTest_censored_gumbelCopula <- function(X, ties.method="max", N=250, start=NULL) {
  copula <- copula::gumbelCopula();
  U <- copula::pobs(X);
  n <- nrow(X);
  
  F.1.inv <- function(u) { return(quantile(U[, 1], u, type = 1, names = FALSE)); }
  F.2.inv <- function(u) { return(quantile(U[, 2], u, type = 1, names = FALSE)); }
  
  # First fit for the test statistic calcula
  theta_n <- copula@parameters;
  try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X, start)$par }, silent = TRUE);
  
  # Stop, if the fit fails.
  if (any(is.na(theta_n))) { 
    return(structure(class="htest", list(statistic=NA, p.value=NA))); 
    stop();
  }
  copula.theta_n <- copula::gumbelCopula(dim=2, theta_n);
  
  # Calculate the main test statistic
  eval.fitCopula <- copula::pCopula(U, copula.theta_n);
  eval.empCopula <- copula::C.n(U,U);
  
  Gn <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE); 
  
  # Estimate the p value with bootstrap replicates:
  Gn_bootstrap = rep(0, N);
  for (k in 1:N) {
    U_bootstrap <- copula::rCopula(n, copula.theta_n);
    X_bootstrap <- U_bootstrap;
    X_bootstrap[, 1] <- F.1.inv(U_bootstrap[, 1]);
    X_bootstrap[, 2] <- F.2.inv(U_bootstrap[, 2]);
    
    theta_n <- copula@parameters;
    try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X_bootstrap, start)$par }, silent = TRUE);
    
    if (is.na(theta_n)) { 
      Gn_bootstrap[k] <- NA;
    } else {
      copula.theta_n <- copula::gumbelCopula(dim=2, theta_n);
      eval.fitCopula <- copula::pCopula(U_bootstrap, copula.theta_n);
      eval.empCopula <- copula::C.n(U_bootstrap, U_bootstrap);
      Gn_bootstrap[k] <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE);
    }
  }
  
  # Output of the results
  structure(class = "htest", list(
    statistic = c(statistic = Gn),
    p.value = (sum(Gn_bootstrap > Gn, na.rm = TRUE) + 0.5) / (N + 1)
  ));
}

gofTest_censored_claytonCopula <- function(X, ties.method="max", N=250, start=NULL) {
  copula <- copula::claytonCopula();
  U <- copula::pobs(X);
  n <- nrow(X);
  
  F.1.inv <- function(u) { return(quantile(U[, 1], u, type = 1, names = FALSE)); }
  F.2.inv <- function(u) { return(quantile(U[, 2], u, type = 1, names = FALSE)); }
  
  # First fit for the test statistic calcula
  theta_n <- copula@parameters;
  try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X, start)$par }, silent = TRUE);
  
  # Stop, if the fit fails.
  if (any(is.na(theta_n))) { 
    return(structure(class="htest", list(statistic=NA, p.value=NA))); 
    stop();
  }
  
  copula.theta_n <- copula::claytonCopula(dim=2, theta_n);
  
  # Calculate the main test statistic
  eval.fitCopula <- copula::pCopula(U, copula.theta_n);
  eval.empCopula <- copula::C.n(U,U);
  
  Gn <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE); 
  
  # Estimate the p value with bootstrap replicates:
  Gn_bootstrap = rep(0, N);
  for (k in 1:N) {
    U_bootstrap <- copula::rCopula(n, copula.theta_n);
    X_bootstrap <- U_bootstrap;
    X_bootstrap[, 1] <- F.1.inv(U_bootstrap[, 1]);
    X_bootstrap[, 2] <- F.2.inv(U_bootstrap[, 2]);
    
    theta_n <- copula@parameters;
    try({ theta_n <- censoredMPLE(copula=copula, method="BFGS", X_bootstrap, start)$par }, silent = TRUE);
    
    if (any(is.na(theta_n))) { 
      Gn_bootstrap[k] <- NA;
    } else {
      copula.theta_n <- copula::claytonCopula(dim=2, theta_n);
      eval.fitCopula <- copula::pCopula(U_bootstrap, copula.theta_n);
      eval.empCopula <- copula::C.n(U_bootstrap, U_bootstrap);
      Gn_bootstrap[k] <- sum((eval.fitCopula - eval.empCopula) ** 2, na.rm = TRUE);
    }
  }
  
  # Output of the results
  structure(class = "htest", list(
    statistic = c(statistic = Gn),
    p.value = (sum(Gn_bootstrap > Gn, na.rm = TRUE) + 0.5) / (N + 1)
  ));
}