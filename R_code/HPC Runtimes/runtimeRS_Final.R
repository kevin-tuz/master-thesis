# !/usr/bin/env Rscript
# Implementation: Kevin Tuz [kevin.tuz@hhu.de / ketuz100]

  library(optparse)

# Runtime options
  option_list = list(
    make_option(c("--copula"), type = "character", default = "normal"),
    make_option(c("--tau"),    type = "double",    default = 0.75),
    make_option(c("--iter"),   type = "integer",   default = 500),
    make_option(c("--sample"), type = "double",    default = 250),
    make_option(c("--prec"),   type = "double",    default = NA),
    make_option(c("--omit"),   type = "double",    default = NA),
    make_option(c("--map"),    type = "character", default = NA),
    make_option(c("--path"),   type = "character", default = "./outputs/")
  )
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
# Report the parameter
  print(paste0("Starting the simulation for radial symmetrie at ", Sys.time(), " with following parameters:"))
  print(paste0(" - Copula ", opt$copula, ", with Kendall's tau ", opt$tau, "."))
  print(paste0(" - ", opt$iter, " iterations; sample size: ", opt$sample, "."))
  if (!is.na(opt$prec)) { print(paste0(" - Rounding with 10E-", opt$round, " precision.")) }
  if (!is.na(opt$omit)) { print(paste0(" - Omission with ", opt$omit*100, " % probability.")) }
  if (!is.na(opt$map)) { print(paste0(" - Empirical structure mapping: '", opt$map, "'.")) }

# Loading of libraries and files
  library(copula)
  source("methods/copula_tests.R")
  source("datagenerating_process/datagenerator.R")
  source("methods/sample_reductions.R")

# Loading and setup of empirical functions
  empInvCDF <- function(p, x) { return(quantile(x, p, names=F, na.rm=T)) }
  
  x.BC1 <- data.frame(read.csv("empirical_data_for_transformation/BC1.csv"))$rt
  x.BC2 <- data.frame(read.csv("empirical_data_for_transformation/BC2.csv"))$rt
  x.BC3 <- data.frame(read.csv("empirical_data_for_transformation/BC3.csv"))$rt
  x.AC1 <- data.frame(read.csv("empirical_data_for_transformation/AC1.csv"))$rt
  x.AC2 <- data.frame(read.csv("empirical_data_for_transformation/AC2.csv"))$rt
  x.AC3 <- data.frame(read.csv("empirical_data_for_transformation/AC3.csv"))$rt
  x.WC1 <- data.frame(read.csv("empirical_data_for_transformation/WC1.csv"))$rt
  x.WC2 <- data.frame(read.csv("empirical_data_for_transformation/WC2.csv"))$rt
  x.WC3 <- data.frame(read.csv("empirical_data_for_transformation/WC3.csv"))$rt
  
  q.BC1 <- function(p) {return(empInvCDF(p, x.BC1)) }
  q.BC2 <- function(p) {return(empInvCDF(p, x.BC2)) }
  q.BC3 <- function(p) {return(empInvCDF(p, x.BC3)) }
  q.AC1 <- function(p) {return(empInvCDF(p, x.AC1)) }
  q.AC2 <- function(p) {return(empInvCDF(p, x.AC2)) }
  q.AC3 <- function(p) {return(empInvCDF(p, x.AC3)) }
  q.WC1 <- function(p) {return(empInvCDF(p, x.WC1)) }
  q.WC2 <- function(p) {return(empInvCDF(p, x.WC2)) }
  q.WC3 <- function(p) {return(empInvCDF(p, x.WC3)) }
  
  BC <- c(q.BC1, q.BC2, q.BC3)
  AC <- c(q.AC1, q.AC2, q.AC3)
  WC <- c(q.WC1, q.WC2, q.WC3)

# Translation of the input flags and model calibration
  if (opt$copula == "clayton") {
    copula.uncalibrated <- copula::claytonCopula(dim = 2)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::claytonCopula(dim = 2, param = theta)
  }
  if (opt$copula == "gumbel") {
    copula.uncalibrated <- copula::gumbelCopula(dim = 2)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::gumbelCopula(dim = 2, param = theta)
  }
  if (opt$copula == "normal") {
    copula.uncalibrated <- copula::normalCopula(dim = 2)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::normalCopula(dim = 2, param = theta)
  }
  if (opt$copula == "t") {
    copula.uncalibrated <- copula::tCopula(dim = 2, df = 4, df.fixed = T)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::tCopula(dim = 2, df = 4, df.fixed = T, param = theta)
  }
  if (opt$copula == "khoclayton") {
    copula.uncalibrated <- copula::claytonCopula(dim = 2)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::khoudrajiCopula(
        copula1 = copula::indepCopula(), 
        copula2 = copula::claytonCopula(dim = 2, param = theta),
        shapes = c(0.20, 0.95)
    )
  }
  if (opt$copula == "khogumbel") {
    copula.uncalibrated <- copula::gumbelCopula(dim = 2)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::khoudrajiCopula(
      copula1 = copula::indepCopula(), 
      copula2 = copula::gumbelCopula(dim = 2, param = theta),
      shapes = c(0.20, 0.95)
    )
  }
  if (opt$copula == "khonormal") {
    copula.uncalibrated <- copula::normalCopula(dim = 2)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::khoudrajiCopula(
      copula1 = copula::indepCopula(), 
      copula2 = copula::normalCopula(dim = 2, param = theta),
      shapes = c(0.20, 0.95)
    )
  }
  if (opt$copula == "khot") {
    copula.uncalibrated <- copula::tCopula(dim = 2, df = 4, df.fixed = T)
    theta <- copula::iTau(copula.uncalibrated, opt$tau)
    copula <- copula::khoudrajiCopula(
      copula1 = copula::indepCopula(), 
      copula2 = copula::tCopula(dim = 2, df = 4, df.fixed = T, param = theta),
      shapes = c(0.20, 0.95)
    )
  }

# Abbreviated variable names for the runtime parameter
  N <- opt$iter         # Iteration count N
  n <- opt$sample       # Sample size n
  rp <- 10^(-opt$prec)  # Precision
  rp.jitter <- if (!is.na(rp)) { rp } else { 10^(-6) }
  op <- opt$omit        # Omission probability
  em <- list()          # Empirical mapping
  
  if (is.na(opt$map)) {
    em$x <- em$y <- c(NA, NA, NA)
  } else {
    if (opt$map %in% c("WxW", "WxA", "WxB")) { em$x <- WC }
    if (opt$map %in% c("AxW", "AxA", "AxB")) { em$x <- AC }
    if (opt$map %in% c("BxW", "BxA", "BxB")) { em$x <- BC }
    if (opt$map %in% c("WxW", "AxW", "BxW")) { em$y <- WC }  
    if (opt$map %in% c("WxA", "AxA", "BxA")) { em$y <- AC }
    if (opt$map %in% c("WxB", "AxB", "BxB")) { em$y <- BC }
  }
  
  n.greedy.10 <- 0.10 * n
  n.greedy.25 <- 0.25 * n

# Preparation of output structure
  results <- data.frame(
    index = integer(),
    setup.copula = character(),
    setup.tau = numeric(),
    setup.sample = integer(),
    setup.prec = numeric(),
    setup.omit = numeric(),
    setup.map = character(),
    
    data.unm.tied_full = numeric(),
    data.unm.tied_part = numeric(),
    data.mod.tied_full = numeric(),
    data.mod.tied_part = numeric(),
    
    m1.Reference.unm.stat = numeric(),
    m1.Reference.mod.stat = numeric(),
    m1.Reference.unm.pval = numeric(),
    m1.Reference.mod.pval = numeric(),
    m1.Reference.benchmark = numeric(),
    
    m2.VReduct.unm.stat = numeric(),
    m2.VReduct.mod.stat = numeric(),
    m2.VReduct.unm.pval = numeric(),
    m2.VReduct.mod.pval = numeric(),
    m2.VReduct.benchmark = numeric(),
  
    m3.GReduct25.unm.stat = numeric(),
    m3.GReduct25.mod.stat = numeric(),
    m3.GReduct25.unm.pval = numeric(),
    m3.GReduct25.mod.pval = numeric(),
    m3.GReduct25.benchmark = numeric(),
    
    m4.RndRanks.unm.stat = numeric(),
    m4.RndRanks.mod.stat = numeric(),
    m4.RndRanks.unm.pval = numeric(),
    m4.RndRanks.mod.pval = numeric(),
    m4.RndRanks.benchmark = numeric(),
    
    m5.PiJitter.unm.stat = numeric(),
    m5.PiJitter.mod.stat = numeric(),
    m5.PiJitter.unm.pval = numeric(),
    m5.PiJitter.mod.pval = numeric(),
    m5.PiJitter.benchmark = numeric(),
    
    m6.TauJitter.unm.stat = numeric(),
    m6.TauJitter.mod.stat = numeric(),
    m6.TauJitter.unm.pval = numeric(),
    m6.TauJitter.mod.pval = numeric(),
    m6.TauJitter.benchmark = numeric(),
    
    m7.Kojadinovic.unm.stat = numeric(),
    m7.Kojadinovic.mod.stat = numeric(),
    m7.Kojadinovic.unm.pval = numeric(),
    m7.Kojadinovic.mod.pval = numeric(), 
    m7.Kojadinovic.benchmark = numeric(),
    
    m8.VReductRndRanks.unm.stat = numeric(),
    m8.VReductRndRanks.mod.stat = numeric(),
    m8.VReductRndRanks.unm.pval = numeric(),
    m8.VReductRndRanks.mod.pval = numeric(), 
    m8.VReductRndRanks.benchmark = numeric(),
    
    m9.GReduct25RndRanks.unm.stat = numeric(),
    m9.GReduct25RndRanks.mod.stat = numeric(),
    m9.GReduct25RndRanks.unm.pval = numeric(),
    m9.GReduct25RndRanks.mod.pval = numeric(),
    m9.GReduct25RndRanks.benchmark = numeric()
    
  )

# Simulation runtime
ssts <- Sys.time()
deltas <- c()
for (i in 1:N) {
  
  inverse_marginals <- c(NA, NA)
  im.x <- sample(1:3,1)
  im.y <- sample(1:3,1)
  inverse_marginals <- c(em$x[[im.x]], em$y[[im.y]])
  
  data <- generate_synthetic_data(copula, n, inverse_marginals, rp, op)
  x.unm <- copula::pobs(data$unmodified)
  x.mod <- copula::pobs(data$modified)
  
  data.unm.tied_full <- sum(apply(x.unm, 1, function(y) {
    sum((x.unm[,1] == y[1]) & (x.unm[,2] == y[2])) - 1
  }) > 0)
  data.mod.tied_full <- sum(apply(x.mod, 1, function(y) {
    sum((x.mod[,1] == y[1]) & (x.mod[,2] == y[2])) - 1
  }) > 0)
  
  data.unm.tied_part <- sum(apply(x.unm, 1, function(y) {
    sum((x.unm[,1] == y[1]) | (x.unm[,2] == y[2])) - 1
  }) > 0)
  data.mod.tied_part <- sum(apply(x.mod, 1, function(y) {
    sum((x.mod[,1] == y[1]) | (x.mod[,2] == y[2])) - 1
  }) > 0)
  
  # Running simulation (Method 1)
    tm <- Sys.time()
    m1.Reference.unm <- list(statistic = NA, p.value = NA)
    m1.Reference.mod <- list(statistic = NA, p.value = NA)
    try({m1.Reference.unm <- radSymTest_genestEtNeslehova(x.unm)})
    try({m1.Reference.mod <- radSymTest_genestEtNeslehova(x.mod)})
    m1.Reference.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

    # Running simulation (Method 2)
    tm <- Sys.time()
    m2.VReduct.unm <- list(statistic = NA, p.value = NA)
    m2.VReduct.mod <- list(statistic = NA, p.value = NA)
    x.unm.VReduct <- fulltie_samplereduction(x.unm)
    x.mod.VReduct <- fulltie_samplereduction(x.mod)
    try({m2.VReduct.unm <- radSymTest_genestEtNeslehova(x.unm.VReduct)})
    try({m2.VReduct.mod <- radSymTest_genestEtNeslehova(x.mod.VReduct)})
    m2.VReduct.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

    # Running simulation (Method 3) 
    tm <- Sys.time()
    m3.GReduct25.unm <- list(statistic = NA, p.value = NA)
    m3.GReduct25.mod <- list(statistic = NA, p.value = NA)
    try({
      x.unm.GReduct.25 <- greedy_samplereduction_BIV(x.unm, n.greedy.25)
      m3.GReduct25.unm <- radSymTest_genestEtNeslehova(x.unm.GReduct.25)
    })
    try({
      x.mod.GReduct.25 <- greedy_samplereduction_BIV(x.mod, n.greedy.25)
      m3.GReduct25.mod <- radSymTest_genestEtNeslehova(x.mod.GReduct.25)
    })
    m3.GReduct25.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
    
    # Running simulation (Method 4)
    tm <- Sys.time()
    m4.RndRanks.unm <- list(statistic = NA, p.value = NA)
    m4.RndRanks.mod <- list(statistic = NA, p.value = NA)
    try({m4.RndRanks.unm <- radSymTest_genestEtNeslehova(x.unm, ties.method = "random")})
    try({m4.RndRanks.mod <- radSymTest_genestEtNeslehova(x.mod, ties.method = "random")})
    m4.RndRanks.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

    # Running simulation (Method 5)
    tm <- Sys.time()
    m5.PiJitter.unm <- list(statistic = NA, p.value = NA)
    m5.PiJitter.mod <- list(statistic = NA, p.value = NA)
    noise.PiJitter <- copula::rCopula(n, copula::indepCopula(2)) * rp.jitter
    try({m5.PiJitter.unm <- radSymTest_genestEtNeslehova(x.unm + noise.PiJitter)})
    try({m5.PiJitter.mod <- radSymTest_genestEtNeslehova(x.mod + noise.PiJitter)})
    m5.PiJitter.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
    
    # Running simulation (Method 6)
    tm <- Sys.time()
    copList <- c(copula::normalCopula(0,2), copula::upfhCopula(2))
    empKTau.unm <- empKTau.mod <- 0
    try({empKTau.unm <- min(max(copula::corKendall(x.unm)[1,2],0),1)})
    try({empKTau.mod <- min(max(copula::corKendall(x.mod)[1,2],0),1)})

    try({
      mixCopula.unm <- copula::mixCopula(copList, c(empKTau.unm, 1 - empKTau.unm))
      mixCopula.mod <- copula::mixCopula(copList, c(empKTau.mod, 1 - empKTau.mod))
      noise.TauJitter.unm <- copula::rCopula(n, mixCopula.unm) * rp.jitter
      noise.TauJitter.mod <- copula::rCopula(n, mixCopula.mod) * rp.jitter
    })
    m6.TauJitter.unm <- list(statistic = NA, p.value = NA)
    m6.TauJitter.mod <- list(statistic = NA, p.value = NA)
    try({m6.TauJitter.unm <- radSymTest_genestEtNeslehova(x.unm + noise.TauJitter.unm)})
    try({m6.TauJitter.mod <- radSymTest_genestEtNeslehova(x.mod + noise.TauJitter.mod)})
    m6.TauJitter.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
    
    # Running simulation (Method 7)
    tm <- Sys.time()  
    m7.Kojadinovic.unm <- list(statistic = NA, p.value = NA)
    m7.Kojadinovic.mod <- list(statistic = NA, p.value = NA)
    try({m7.Kojadinovic.unm <- copula::radSymTest(x.unm, ties = TRUE)})
    try({m7.Kojadinovic.mod <- copula::radSymTest(x.mod, ties = TRUE)})
    m7.Kojadinovic.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
    
    # Running simulation (Method 8)
    tm <- Sys.time()  
    m8.VReductRndRanks.unm <- list(statistic = NA, p.value = NA)
    m8.VReductRndRanks.mod <- list(statistic = NA, p.value = NA)
    try({
      x.unm.VReduct <- fulltie_samplereduction(x.unm)
      m8.VReductRndRanks.unm <- radSymTest_genestEtNeslehova(
        x.unm.VReduct, ties.method = "random")
    })
    try({
      x.mod.VReduct <- fulltie_samplereduction(x.mod)
      m8.VReductRndRanks.mod <- radSymTest_genestEtNeslehova(
        x.mod.VReduct, ties.method = "random")
    })
    m8.VReductRndRanks.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

    # Running simulation (Method 9)
    tm <- Sys.time()
    m9.GReduct25RndRanks.unm <- list(statistic = NA, p.value = NA)
    m9.GReduct25RndRanks.mod <- list(statistic = NA, p.value = NA)
    try({
      x.unm.GReduct.25 <- greedy_samplereduction_BIV(x.unm, n.greedy.25)
      m9.GReduct25RndRanks.unm <- radSymTest_genestEtNeslehova(
        x.unm.GReduct.25, ties.method = "random")
    })
    try({
    x.mod.GReduct.25 <- greedy_samplereduction_BIV(x.mod, n.greedy.25)
    m9.GReduct25RndRanks.mod <- radSymTest_genestEtNeslehova(
      x.mod.GReduct.25, ties.method = "random")
    })
    m9.GReduct25RndRanks.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
    
  iresult <- list(
      i, opt$copula, opt$tau, opt$sample, rp, op, opt$map,
      data.unm.tied_full, data.unm.tied_part, 
      data.mod.tied_full, data.mod.tied_part,
      
      m1.Reference.unm$statistic,
      m1.Reference.mod$statistic,
      m1.Reference.unm$p.value,
      m1.Reference.mod$p.value,
      m1.Reference.benchmark,
      
      m2.VReduct.unm$statistic,
      m2.VReduct.mod$statistic,
      m2.VReduct.unm$p.value,
      m2.VReduct.mod$p.value,
      m2.VReduct.benchmark,
      
      m3.GReduct25.unm$statistic,
      m3.GReduct25.mod$statistic,
      m3.GReduct25.unm$p.value,
      m3.GReduct25.mod$p.value,
      m3.GReduct25.benchmark,
      
      m4.RndRanks.unm$statistic,
      m4.RndRanks.mod$statistic,
      m4.RndRanks.unm$p.value,
      m4.RndRanks.mod$p.value,
      m4.RndRanks.benchmark,
      
      m5.PiJitter.unm$statistic,
      m5.PiJitter.mod$statistic,
      m5.PiJitter.unm$p.value,
      m5.PiJitter.mod$p.value,
      m5.PiJitter.benchmark,
      
      m6.TauJitter.unm$statistic,
      m6.TauJitter.mod$statistic,
      m6.TauJitter.unm$p.value,
      m6.TauJitter.mod$p.value,
      m6.TauJitter.benchmark,
      
      m7.Kojadinovic.unm$statistic,
      m7.Kojadinovic.mod$statistic,
      m7.Kojadinovic.unm$p.value,
      m7.Kojadinovic.mod$p.value, 
      m7.Kojadinovic.benchmark,
      
      m8.VReductRndRanks.unm$statistic,
      m8.VReductRndRanks.mod$statistic,
      m8.VReductRndRanks.unm$p.value,
      m8.VReductRndRanks.mod$p.value,
      m8.VReductRndRanks.benchmark,
      
      m9.GReduct25RndRanks.unm$statistic,
      m9.GReduct25RndRanks.mod$statistic,
      m9.GReduct25RndRanks.unm$p.value,
      m9.GReduct25RndRanks.mod$p.value,
      m9.GReduct25RndRanks.benchmark
    )
    
    results[nrow(results) + 1, ] <- iresult

    if (i %% ceiling(N/10) == 0) {
      delta <- difftime(Sys.time(), ssts, units = "secs")[[1]]
      deltas <- c(deltas, delta)

      spi <- sum(deltas) / i
      remaining <- spi * (N - i)
      ssts <- Sys.time()

      print(paste0(i, " iterations in ", round(delta, 2), " s. ETA: ", (ssts + remaining)))
      
      write.csv(x <- results, file <- paste0(
          opt$path2, "temp_", round(i / N * 100), "perc_RS_", opt$copula, "_", opt$tau, "_N", N,
          "_n", n,"_rp", rp,"_op", op, "_em", opt$map, ".csv")
      )
    }
}

# Simulation finished. --------
  print("Saving the results ...")
  write.csv(
    x <- results,
    file <- paste0(
      opt$path, "RS_", opt$copula, "_", opt$tau, "_N", N,
      "_n", n,"_rp", rp,"_op", op, "_em", opt$map, ".csv")
  )

  print(paste0(
    "Stopping the simulation at ", Sys.time(), "."
  ))