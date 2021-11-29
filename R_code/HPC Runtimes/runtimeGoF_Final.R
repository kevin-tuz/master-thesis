# !/usr/bin/env Rscript
# Implementation: Kevin Tuz [kevin.tuz@hhu.de / ketuz100]

library(optparse)

# Runtime options
option_list <- list(
  make_option(c("--copula"), type = "character", default = "gumbel"),
  make_option(c("--tau"),    type = "double",    default = 0.75),
  make_option(c("--iter"),   type = "integer",   default = 250),
  make_option(c("--sample"), type = "double",    default = 250),
  make_option(c("--prec"),   type = "double",    default = NA),
  make_option(c("--omit"),   type = "double",    default = NA),
  make_option(c("--map"),    type = "character", default = NA),
  make_option(c("--path"),   type = "character", default = "./outputs/"),
  make_option(c("--path2"),   type = "character", default = "./temp/")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
  
# Report the parameter
print(paste0("Starting the sim for fit & GoF at ", Sys.time(), ":"))
print(paste0("- Runtime code: ", opt$copula, "(", opt$tau, "), N=", opt$iter, ", n=", opt$sample, "."))
if (!is.na(opt$prec)) { print(paste0("- Rounding with 10E-", opt$prec, " precision.")) }
if (!is.na(opt$omit)) { print(paste0("- Omission with ", opt$omit*100, " % probability.")) }
if (!is.na(opt$map)) { print(paste0("- Empirical structure mapping: '", opt$map, "'.")) }

# Loading of libraries and files
library(copula)
source("methods/copula_tests.R")
source("datagenerating_process/datagenerator.R")
source("methods/sample_reductions.R")

# Loading and setup of empirical functions
empInvCDF <- function(p, x) { return(quantile(x, p, names = FALSE, na.rm = TRUE)) }
  
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


# Model calibration -----
if (opt$copula == "clayton") {
  copula.uncalibrated <- copula::claytonCopula(dim = 2)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::claytonCopula(dim = 2, param = theta)
  initial.theta <- 2
}
if (opt$copula == "gumbel") {
  copula.uncalibrated <- copula::gumbelCopula(dim = 2)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::gumbelCopula(dim = 2, param = theta)
  initial.theta <- 2
}
if (opt$copula == "normal") {
  copula.uncalibrated <- copula::normalCopula(dim = 2)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::normalCopula(dim = 2, param = theta)
  initial.theta <- 0.7071068
}
if (opt$copula == "t") {
  copula.uncalibrated <- copula::tCopula(dim = 2, df = 4, df.fixed = T)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::tCopula(dim = 2, df = 4, df.fixed = T, param = theta)
  initial.theta <- c(0.7071068, 4)
}
if (opt$copula == "khoclayton") {
  copula.uncalibrated <- copula::claytonCopula(dim = 2)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::khoudrajiCopula(
    copula1 = copula::indepCopula(), 
    copula2 = copula::claytonCopula(dim = 2, param = theta),
    shapes = c(0.20, 0.95)
  )
  initial.theta <- 2
}
if (opt$copula == "khogumbel") {
  copula.uncalibrated <- copula::gumbelCopula(dim = 2)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::khoudrajiCopula(
    copula1 = copula::indepCopula(), 
    copula1 = copula::gumbelCopula(dim = 2, param = theta),
    shapes = c(0.20, 0.95)
  )
  initial.theta <- 2
}
if (opt$copula == "khonormal") {
  copula.uncalibrated <- copula::normalCopula(dim = 2)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::khoudrajiCopula(
    copula1 = copula::indepCopula(), 
    copula2 = copula::normalCopula(dim = 2, param = theta),
    shapes = c(0.20, 0.95)
  )
  initial.theta <- 0.7071068
}
if (opt$copula == "khot") {
  copula.uncalibrated <- copula::tCopula(dim = 2, df = 4, df.fixed = T)
  theta <- copula::iTau(copula.uncalibrated, opt$tau)
  copula <- copula::khoudrajiCopula(
    copula1 = copula::indepCopula(), 
    copula2 = copula::tCopula(dim = 2, df = 4, df.fixed = T, param = theta),
    shapes = c(0.20, 0.95)
  )
  initial.theta <- c(0.7071068,4)
}

N <- opt$iter         # Iteration count N
N.gof <- 250
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

n.greedy.25 <- 0.10 * n
n.greedy.25 <- ceiling(0.25 * n)

# Preparation of output structure
results <- data.frame(
  index = integer(),
  setup.copula = character(),
  setup.tau = numeric(),
  setup.sample = integer(),
  setup.prec = numeric(),
  setup.omit = numeric(),
  setup.map = character(),

  setup.truetheta = numeric(),

  data.unm.tied_full = numeric(),
  data.unm.tied_part = numeric(),
  data.mod.tied_full = numeric(),
  data.mod.tied_part = numeric(),

  # Method 1 - Reference - Fit & GoF
  m1.Reference.estimate = numeric(),
  m1.Reference.stderror = numeric(),
  m1.Reference.CI95.lwr = numeric(),
  m1.Reference.CI95.upr = numeric(),
  m1.Reference.fit.benchmark = numeric(),
  m1.Reference.normal.stat = numeric(),
  m1.Reference.normal.pval = numeric(),
  m1.Reference.t.stat = numeric(),
  m1.Reference.t.pval = numeric(),
  m1.Reference.clayton.stat = numeric(),
  m1.Reference.clayton.pval = numeric(),
  m1.Reference.gumbel.stat = numeric(),
  m1.Reference.gumbel.pval = numeric(),
  m1.Reference.gof.benchmark = numeric(),

  # Method 2 - V-Reduction - Fit & GoF
  m2.VReduct.estimate = numeric(),
  m2.VReduct.stderror = numeric(),
  m2.VReduct.CI95.lwr = numeric(),
  m2.VReduct.CI95.upr = numeric(),
  m2.VReduct.fit.benchmark = numeric(),
  m2.VReduct.normal.stat = numeric(),
  m2.VReduct.normal.pval = numeric(),
  m2.VReduct.t.stat = numeric(),
  m2.VReduct.t.pval = numeric(),
  m2.VReduct.clayton.stat = numeric(),
  m2.VReduct.clayton.pval = numeric(),
  m2.VReduct.gumbel.stat = numeric(),
  m2.VReduct.gumbel.pval = numeric(),
  m2.VReduct.gof.benchmark = numeric(),

  # Method 3 - G-Reduction - Fit & GoF
  m3.GReduct.estimate = numeric(),
  m3.GReduct.stderror = numeric(),
  m3.GReduct.CI95.lwr = numeric(),
  m3.GReduct.CI95.upr = numeric(),
  m3.GReduct.fit.benchmark = numeric(),
  m3.GReduct.normal.stat = numeric(),
  m3.GReduct.normal.pval = numeric(),
  m3.GReduct.t.stat = numeric(),
  m3.GReduct.t.pval = numeric(),
  m3.GReduct.clayton.stat = numeric(),
  m3.GReduct.clayton.pval = numeric(),
  m3.GReduct.gumbel.stat = numeric(),
  m3.GReduct.gumbel.pval = numeric(),
  m3.GReduct.gof.benchmark = numeric(),

  # Method 4 - Random Ranks - Fit & GoF
  m4.RndRanks.estimate = numeric(),
  m4.RndRanks.stderror = numeric(),
  m4.RndRanks.CI95.lwr = numeric(),
  m4.RndRanks.CI95.upr = numeric(),
  m4.RndRanks.fit.benchmark = numeric(),
  m4.RndRanks.normal.stat = numeric(),
  m4.RndRanks.normal.pval = numeric(),
  m4.RndRanks.t.stat = numeric(),
  m4.RndRanks.t.pval = numeric(),
  m4.RndRanks.clayton.stat = numeric(),
  m4.RndRanks.clayton.pval = numeric(),
  m4.RndRanks.gumbel.stat = numeric(),
  m4.RndRanks.gumbel.pval = numeric(),
  m4.RndRanks.gof.benchmark = numeric(),

  # Method 5 - Pi-Jitter - Fit & GoF
  m5.PiJitter.estimate = numeric(),
  m5.PiJitter.stderror = numeric(),
  m5.PiJitter.CI95.lwr = numeric(),
  m5.PiJitter.CI95.upr = numeric(),
  m5.PiJitter.fit.benchmark = numeric(),
  m5.PiJitter.normal.stat = numeric(),
  m5.PiJitter.normal.pval = numeric(),
  m5.PiJitter.t.stat = numeric(),
  m5.PiJitter.t.pval = numeric(),
  m5.PiJitter.clayton.stat = numeric(),
  m5.PiJitter.clayton.pval = numeric(),
  m5.PiJitter.gumbel.stat = numeric(),
  m5.PiJitter.gumbel.pval = numeric(),
  m5.PiJitter.gof.benchmark = numeric(),

  # Method 6 - Tau-Jitter - Fit & GoF
  m6.TauJitter.estimate = numeric(),
  m6.TauJitter.stderror = numeric(),
  m6.TauJitter.CI95.lwr = numeric(),
  m6.TauJitter.CI95.upr = numeric(),
  m6.TauJitter.fit.benchmark = numeric(),
  m6.TauJitter.normal.stat = numeric(),
  m6.TauJitter.normal.pval = numeric(),
  m6.TauJitter.t.stat = numeric(),
  m6.TauJitter.t.pval = numeric(),
  m6.TauJitter.clayton.stat = numeric(),
  m6.TauJitter.clayton.pval = numeric(),
  m6.TauJitter.gumbel.stat = numeric(),
  m6.TauJitter.gumbel.pval = numeric(),
  m6.TauJitter.gof.benchmark = numeric(),

  # Method 7 - Censoring - Fit & GoF
  m7.Censoring.estimate = numeric(),
  m7.Censoring.stderror = numeric(),
  m7.Censoring.CI95.lwr = numeric(),
  m7.Censoring.CI95.upr = numeric(),
  m7.Censoring.fit.benchmark = numeric(),
  m7.Censoring.normal.stat = numeric(),
  m7.Censoring.normal.pval = numeric(),
  m7.Censoring.t.stat = numeric(),
  m7.Censoring.t.pval = numeric(),
  m7.Censoring.clayton.stat = numeric(),
  m7.Censoring.clayton.pval = numeric(),
  m7.Censoring.gumbel.stat = numeric(),
  m7.Censoring.gumbel.pval = numeric(),
  m7.Censoring.gof.benchmark = numeric(),

  # Method 8 - Kojadinovic - Fit & GoF
  m8.Kojadinovic.estimate = numeric(),
  m8.Kojadinovic.stderror = numeric(),
  m8.Kojadinovic.CI95.lwr = numeric(),
  m8.Kojadinovic.CI95.upr = numeric(),
  m8.Kojadinovic.fit.benchmark = numeric(),
  m8.Kojadinovic.normal.stat = numeric(),
  m8.Kojadinovic.normal.pval = numeric(),
  m8.Kojadinovic.t.stat = numeric(),
  m8.Kojadinovic.t.pval = numeric(),
  m8.Kojadinovic.clayton.stat = numeric(),
  m8.Kojadinovic.clayton.pval = numeric(),
  m8.Kojadinovic.gumbel.stat = numeric(),
  m8.Kojadinovic.gumbel.pval = numeric(),
  m8.Kojadinovic.gof.benchmark = numeric(),

  # Method 9 - V-Reduct-Random-Ranks - Fit & GoF
  m9.VReductRndRanks.estimate = numeric(),
  m9.VReductRndRanks.stderror = numeric(),
  m9.VReductRndRanks.CI95.lwr = numeric(),
  m9.VReductRndRanks.CI95.upr = numeric(),
  m9.VReductRndRanks.fit.benchmark = numeric(),
  m9.VReductRndRanks.normal.stat = numeric(),
  m9.VReductRndRanks.normal.pval = numeric(),
  m9.VReductRndRanks.t.stat = numeric(),
  m9.VReductRndRanks.t.pval = numeric(),
  m9.VReductRndRanks.clayton.stat = numeric(),
  m9.VReductRndRanks.clayton.pval = numeric(),
  m9.VReductRndRanks.gumbel.stat = numeric(),
  m9.VReductRndRanks.gumbel.pval = numeric(),
  m9.VReductRndRanks.gof.benchmark = numeric(),

  # Method 10 - G-Reduct-Random-Ranks - Fit & GoF
  m10.GReductRndRanks.estimate = numeric(),
  m10.GReductRndRanks.stderror = numeric(),
  m10.GReductRndRanks.CI95.lwr = numeric(),
  m10.GReductRndRanks.CI95.upr = numeric(),
  m10.GReductRndRanks.fit.benchmark = numeric(),
  m10.GReductRndRanks.normal.stat = numeric(),
  m10.GReductRndRanks.normal.pval = numeric(),
  m10.GReductRndRanks.t.stat = numeric(),
  m10.GReductRndRanks.t.pval = numeric(),
  m10.GReductRndRanks.clayton.stat = numeric(),
  m10.GReductRndRanks.clayton.pval = numeric(),
  m10.GReductRndRanks.gumbel.stat = numeric(),
  m10.GReductRndRanks.gumbel.pval = numeric(),
  m10.GReductRndRanks.gof.benchmark = numeric(),

  # Method 11 - V-Reduct-Censoring - Fit & GoF
  m11.VReductCensoring.estimate = numeric(),
  m11.VReductCensoring.stderror = numeric(),
  m11.VReductCensoring.CI95.lwr = numeric(),
  m11.VReductCensoring.CI95.upr = numeric(),
  m11.VReductCensoring.fit.benchmark = numeric(),
  m11.VReductCensoring.normal.stat = numeric(),
  m11.VReductCensoring.normal.pval = numeric(),
  m11.VReductCensoring.t.stat = numeric(),
  m11.VReductCensoring.t.pval = numeric(),
  m11.VReductCensoring.clayton.stat = numeric(),
  m11.VReductCensoring.clayton.pval = numeric(),
  m11.VReductCensoring.gumbel.stat = numeric(),
  m11.VReductCensoring.gumbel.pval = numeric(),
  m11.VReductCensoring.gof.benchmark = numeric(),

  # Method 12 - G-Reduct-Censoring - Fit & GoF
  m12.GReductCensoring.estimate = numeric(),
  m12.GReductCensoring.stderror = numeric(),
  m12.GReductCensoring.CI95.lwr = numeric(),
  m12.GReductCensoring.CI95.upr = numeric(),
  m12.GReductCensoring.fit.benchmark = numeric(),
  m12.GReductCensoring.normal.stat = numeric(),
  m12.GReductCensoring.normal.pval = numeric(),
  m12.GReductCensoring.t.stat = numeric(),
  m12.GReductCensoring.t.pval = numeric(),
  m12.GReductCensoring.clayton.stat = numeric(),
  m12.GReductCensoring.clayton.pval = numeric(),
  m12.GReductCensoring.gumbel.stat = numeric(),
  m12.GReductCensoring.gumbel.pval = numeric(),
  m12.GReductCensoring.gof.benchmark = numeric()
)

# Running simulation -------
ssts <- Sys.time()
deltas <- c()
for (i in 1:N) {
  # Running simulation (Generate)
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

  # Method 1 - Reference
  tm <- Sys.time()
  m1.Reference.estimate <- NA
  m1.Reference.stderror <- NA
  m1.Reference.CI95.lwr <- NA
  m1.Reference.CI95.upr <- NA
  
  try({ fit.mod <- copula:::fitCopula(copula, x.mod, method = "mpl") }, silent = TRUE)
  try({ m1.Reference.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m1.Reference.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m1.Reference.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m1.Reference.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m1.Reference.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  tm <- Sys.time()
  m1.Reference.normal <- list(statistic = NA, p.value = NA)
  m1.Reference.t <- list(statistic = NA, p.value = NA)
  m1.Reference.clayton <- list(statistic = NA, p.value = NA)
  m1.Reference.gumbel <- list(statistic = NA, p.value = NA)
  try({ m1.Reference.normal <- copula:::gofCopula(copula::normalCopula(), x.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m1.Reference.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m1.Reference.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m1.Reference.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  m1.Reference.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  # Method 2: V-Reduction
  tm <- Sys.time()
  m2.VReduct.estimate <- NA
  m2.VReduct.stderror <- NA
  m2.VReduct.CI95.lwr <- NA
  m2.VReduct.CI95.upr <- NA
  x.mod.VReduct <- x.mod
  try({ x.mod.VReduct <- fulltie_samplereduction(x.mod) }, silent = TRUE)
  try({ fit.mod <- copula:::fitCopula(copula, x.mod.VReduct, method = "mpl") }, silent = TRUE)
  try({ m2.VReduct.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m2.VReduct.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m2.VReduct.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m2.VReduct.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m2.VReduct.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  tm <- Sys.time()
  m2.VReduct.normal <- list(statistic = NA, p.value = NA)
  m2.VReduct.t <- list(statistic = NA, p.value = NA)
  m2.VReduct.clayton <- list(statistic = NA, p.value = NA)
  m2.VReduct.gumbel <- list(statistic = NA, p.value = NA)
  try({ m2.VReduct.normal <- copula:::gofCopula(copula::normalCopula(), x.mod.VReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m2.VReduct.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod.VReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m2.VReduct.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod.VReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m2.VReduct.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod.VReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  m2.VReduct.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  # Method 3: G-Reduction
  tm <- Sys.time()
  m3.GReduct.estimate <- NA
  m3.GReduct.stderror <- NA
  m3.GReduct.CI95.lwr <- NA
  m3.GReduct.CI95.upr <- NA
  x.mod.GReduct <- x.mod
  try({ x.mod.GReduct <- greedy_samplereduction_BIV(x.mod, n.greedy.25) }, silent = TRUE)
  try({ fit.mod <- copula:::fitCopula(copula, x.mod.GReduct, method = "mpl") }, silent = TRUE)
  try({ m3.GReduct.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m3.GReduct.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m3.GReduct.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m3.GReduct.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m3.GReduct.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  tm <- Sys.time()
  m3.GReduct.normal <- list(statistic = NA, p.value = NA)
  m3.GReduct.t <- list(statistic = NA, p.value = NA)
  m3.GReduct.clayton <- list(statistic = NA, p.value = NA)
  m3.GReduct.gumbel <- list(statistic = NA, p.value = NA)
  try({ m3.GReduct.normal <- copula:::gofCopula(copula::normalCopula(), x.mod.GReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({m3.GReduct.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod.GReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m3.GReduct.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod.GReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m3.GReduct.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod.GReduct, N = N.gof, ties = FALSE) }, silent = TRUE)
  m3.GReduct.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  # Method 4: Random-Ranks
  tm <- Sys.time()
  m4.RndRanks.estimate <- NA
  m4.RndRanks.stderror <- NA
  m4.RndRanks.CI95.lwr <- NA
  m4.RndRanks.CI95.upr <- NA
  try({ fit.mod <- copula:::fitCopula(copula, copula:::pobs(x.mod, ties.method = "random"), method = "mpl") }, silent = TRUE)
  try({ m4.RndRanks.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m4.RndRanks.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m4.RndRanks.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m4.RndRanks.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m4.RndRanks.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  tm <- Sys.time()
  m4.RndRanks.normal <- list(statistic = NA, p.value = NA)
  m4.RndRanks.t <- list(statistic = NA, p.value = NA)
  m4.RndRanks.clayton <- list(statistic = NA, p.value = NA)
  m4.RndRanks.gumbel <- list(statistic = NA, p.value = NA)
  try({ m4.RndRanks.normal <- copula:::gofCopula(copula::normalCopula(), x.mod, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m4.RndRanks.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m4.RndRanks.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m4.RndRanks.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  m4.RndRanks.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  # Method 5: Pi-Jitter
  tm <- Sys.time()
  m5.PiJitter.estimate <- NA
  m5.PiJitter.stderror <- NA
  m5.PiJitter.CI95.lwr <- NA
  m5.PiJitter.CI95.upr <- NA
  noise.PiJitter <- copula::rCopula(n, copula::indepCopula(2)) * rp.jitter
  try({ fit.mod <- copula:::fitCopula(copula, x.mod + noise.PiJitter, method = "mpl") }, silent = TRUE)
  try({ m5.PiJitter.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m5.PiJitter.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m5.PiJitter.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m5.PiJitter.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m5.PiJitter.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  tm <- Sys.time()
  noise.PiJitter <- copula::rCopula(n, copula::indepCopula(2)) * rp.jitter
  m5.PiJitter.normal <- list(statistic = NA, p.value = NA)
  m5.PiJitter.t <- list(statistic = NA, p.value = NA)
  m5.PiJitter.clayton <- list(statistic = NA, p.value = NA)
  m5.PiJitter.gumbel <- list(statistic = NA, p.value = NA)
  try({ m5.PiJitter.normal <- copula:::gofCopula(copula::normalCopula(), x.mod + noise.PiJitter, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m5.PiJitter.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod + noise.PiJitter, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m5.PiJitter.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod + noise.PiJitter, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m5.PiJitter.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod + noise.PiJitter, N = N.gof, ties = FALSE) }, silent = TRUE)
  m5.PiJitter.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  # Method 6: Tau-Jitter
  tm <- Sys.time()
  m6.TauJitter.estimate <- NA
  m6.TauJitter.stderror <- NA
  m6.TauJitter.CI95.lwr <- NA
  m6.TauJitter.CI95.upr <- NA
  copList <- c(copula::normalCopula(0,2), copula::upfhCopula(2))
  empKTau.mod <- 0
  try({ empKTau.mod <- min(max(copula::corKendall(x.mod)[1,2],0),1) }, silent = TRUE)
  try({ mixCopula.mod <- copula::mixCopula(copList, c(empKTau.mod, 1 - empKTau.mod)) }, silent = TRUE)
  try({ noise.TauJitter.mod <- copula::rCopula(n, mixCopula.mod) * rp.jitter }, silent = TRUE)
  try({ fit.mod <- copula:::fitCopula(copula, x.mod + noise.TauJitter.mod, method = "mpl") }, silent = TRUE)
  try({ m6.TauJitter.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m6.TauJitter.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m6.TauJitter.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m6.TauJitter.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m6.TauJitter.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  tm <- Sys.time()
  copList <- c(copula::normalCopula(0,2), copula::upfhCopula(2))
  empKTau.mod <- 0
  try({ empKTau.mod <- min(max(copula::corKendall(x.mod)[1,2],0),1) }, silent = TRUE)
  try({ mixCopula.mod <- copula::mixCopula(copList, c(empKTau.mod, 1 - empKTau.mod)) }, silent = TRUE)
  try({ noise.TauJitter.mod <- copula::rCopula(n, mixCopula.mod) * rp.jitter }, silent = TRUE)
  m6.TauJitter.normal <- list(statistic = NA, p.value = NA)
  m6.TauJitter.t <- list(statistic = NA, p.value = NA)
  m6.TauJitter.clayton <- list(statistic = NA, p.value = NA)
  m6.TauJitter.gumbel <- list(statistic = NA, p.value = NA)
  try({ m6.TauJitter.normal <- copula:::gofCopula(copula::normalCopula(), x.mod + noise.TauJitter.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m6.TauJitter.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod + noise.TauJitter.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m6.TauJitter.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod + noise.TauJitter.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  try({ m6.TauJitter.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod + noise.TauJitter.mod, N = N.gof, ties = FALSE) }, silent = TRUE)
  m6.TauJitter.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  # Method 7: Censoring
  tm <- Sys.time()
  m7.Censoring.estimate <- NA
  m7.Censoring.stderror <- NA
  m7.Censoring.CI95.lwr <- NA
  m7.Censoring.CI95.upr <- NA
  try({ fit.mod <- censoredMPLE(copula = copula, method="BFGS", x.mod, start=initial.theta) }, silent = TRUE)
  try({ m7.Censoring.estimate <- fit.mod$par }, silent = TRUE)
  try({ m7.Censoring.stderror <- NA }, silent = TRUE)
  try({ m7.Censoring.CI95.lwr <- NA }, silent = TRUE)
  try({ m7.Censoring.CI95.upr <- NA }, silent = TRUE)
  m7.Censoring.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  tm <- Sys.time()
  m7.Censoring.normal <- list(statistic = NA, p.value = NA)
  m7.Censoring.t <- list(statistic = NA, p.value = NA)
  m7.Censoring.clayton <- list(statistic = NA, p.value = NA)
  m7.Censoring.gumbel <- list(statistic = NA, p.value = NA)
  try({ m7.Censoring.normal <- gofTest_censored_normalCopula(x.mod, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  try({ m7.Censoring.t <- gofTest_censored_tCopula(x.mod, N = N.gof, start = c(initial.theta[1], 4)) }, silent = TRUE)
  try({ m7.Censoring.clayton <- gofTest_censored_claytonCopula(x.mod, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  try({ m7.Censoring.gumbel <- gofTest_censored_gumbelCopula(x.mod, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  m7.Censoring.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  # Method 8: Kojadinovic
  tm <- Sys.time()
  m8.Kojadinovic.estimate <- NA
  m8.Kojadinovic.stderror <- NA
  m8.Kojadinovic.CI95.lwr <- NA
  m8.Kojadinovic.CI95.upr <- NA
  try({ fit.mod <- copula:::fitCopula(copula, x.mod, method = "mpl") }, silent = TRUE)
  try({ m8.Kojadinovic.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m8.Kojadinovic.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m8.Kojadinovic.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m8.Kojadinovic.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m8.Kojadinovic.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  tm <- Sys.time()
  m8.Kojadinovic.normal <- list(statistic = NA, p.value = NA)
  m8.Kojadinovic.t <- list(statistic = NA, p.value = NA)
  m8.Kojadinovic.clayton <- list(statistic = NA, p.value = NA)
  m8.Kojadinovic.gumbel <- list(statistic = NA, p.value = NA)
  try({ m8.Kojadinovic.normal <- copula:::gofCopula(copula::normalCopula(), x.mod, N = N.gof, ties = TRUE) }, silent = TRUE)
  try({ m8.Kojadinovic.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod, N = N.gof, ties = TRUE) }, silent = TRUE)
  try({ m8.Kojadinovic.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod, N = N.gof, ties = TRUE) }, silent = TRUE)
  try({ m8.Kojadinovic.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod, N = N.gof, ties = TRUE) }, silent = TRUE)
  m8.Kojadinovic.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  # Method 9: V-Reduct + Rnd Ranks
  tm <- Sys.time()
  m9.VReductRndRanks.estimate <- NA
  m9.VReductRndRanks.stderror <- NA
  m9.VReductRndRanks.CI95.lwr <- NA
  m9.VReductRndRanks.CI95.upr <- NA
  x.mod.VReduct <- x.mod
  try({ x.mod.VReduct <- fulltie_samplereduction(x.mod) }, silent = TRUE)
  try({ fit.mod <- copula:::fitCopula(copula, copula:::pobs(x.mod.VReduct, ties.method = "random"), method = "mpl") }, silent = TRUE)
  try({ m9.VReductRndRanks.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m9.VReductRndRanks.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m9.VReductRndRanks.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m9.VReductRndRanks.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m9.VReductRndRanks.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  tm <- Sys.time()
  x.mod.VReduct <- x.mod
  try({ x.mod.VReduct <- fulltie_samplereduction(x.mod) }, silent = TRUE)
  m9.VReductRndRanks.normal <- list(statistic = NA, p.value = NA)
  m9.VReductRndRanks.t <- list(statistic = NA, p.value = NA)
  m9.VReductRndRanks.clayton <- list(statistic = NA, p.value = NA)
  m9.VReductRndRanks.gumbel <- list(statistic = NA, p.value = NA)
  try({ m9.VReductRndRanks.normal <- copula:::gofCopula(copula::normalCopula(), x.mod.VReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m9.VReductRndRanks.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod.VReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m9.VReductRndRanks.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod.VReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m9.VReductRndRanks.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod.VReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  m9.VReductRndRanks.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

# Method 10: G-Reduct + Rnd Ranks
  tm <- Sys.time()
  m10.GReductRndRanks.estimate <- NA
  m10.GReductRndRanks.stderror <- NA
  m10.GReductRndRanks.CI95.lwr <- NA
  m10.GReductRndRanks.CI95.upr <- NA

  x.mod.GReduct <- x.mod
  try({ x.mod.GReduct <- greedy_samplereduction_BIV(x.mod, n.greedy.25) }, silent = TRUE)
  try({ fit.mod <- copula:::fitCopula(copula, copula:::pobs(x.mod.GReduct, ties.method = "random"), method = "mpl") }, silent = TRUE)
  try({ m10.GReductRndRanks.estimate <- fit.mod@estimate }, silent = TRUE)
  try({ m10.GReductRndRanks.stderror <- sqrt(fit.mod@var.est) }, silent = TRUE)
  try({ m10.GReductRndRanks.CI95.lwr <- confint(fit.mod)[1] }, silent = TRUE)
  try({ m10.GReductRndRanks.CI95.upr <- confint(fit.mod)[2] }, silent = TRUE)
  m10.GReductRndRanks.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  tm <- Sys.time()
  x.mod.GReduct <- x.mod
  try({ x.mod.GReduct <- greedy_samplereduction_BIV(x.mod, n.greedy.25) }, silent = TRUE)
  m10.GReductRndRanks.normal <- list(statistic = NA, p.value = NA)
  m10.GReductRndRanks.t <- list(statistic = NA, p.value = NA)
  m10.GReductRndRanks.clayton <- list(statistic = NA, p.value = NA)
  m10.GReductRndRanks.gumbel <- list(statistic = NA, p.value = NA)
  try({ m10.GReductRndRanks.normal <- copula:::gofCopula(copula::normalCopula(), x.mod.GReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m10.GReductRndRanks.t <- copula:::gofCopula(copula::tCopula(df=4, df.fixed=T), x.mod.GReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m10.GReductRndRanks.clayton <- copula:::gofCopula(copula::claytonCopula(), x.mod.GReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  try({ m10.GReductRndRanks.gumbel <- copula:::gofCopula(copula::gumbelCopula(), x.mod.GReduct, N = N.gof, ties = FALSE, ties.method = "random") }, silent = TRUE)
  m10.GReductRndRanks.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  # Method 11: V-Reduct + Censoring
  tm <- Sys.time()
  m11.VReductCensoring.estimate <- NA
  m11.VReductCensoring.stderror <- NA
  m11.VReductCensoring.CI95.lwr <- NA
  m11.VReductCensoring.CI95.upr <- NA
  x.mod.VReduct <- x.mod
  try({ x.mod.VReduct <- fulltie_samplereduction(x.mod) }, silent = TRUE)
  try({ fit.mod <- censoredMPLE(copula = copula, method="BFGS", x.mod.VReduct, start=initial.theta) }, silent = TRUE)
  try({ m11.VReductCensoring.estimate <- fit.mod$par }, silent = TRUE)
  try({ m11.VReductCensoring.stderror <- NA }, silent = TRUE)
  try({ m11.VReductCensoring.CI95.lwr <- NA }, silent = TRUE)
  try({ m11.VReductCensoring.CI95.upr <- NA }, silent = TRUE)
  m11.VReductCensoring.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  tm <- Sys.time()
  x.mod.VReduct <- x.mod
  try({ x.mod.VReduct <- fulltie_samplereduction(x.mod) }, silent = TRUE)
  m11.VReductCensoring.normal <- list(statistic = NA, p.value = NA)
  m11.VReductCensoring.t <- list(statistic = NA, p.value = NA)
  m11.VReductCensoring.clayton <- list(statistic = NA, p.value = NA)
  m11.VReductCensoring.gumbel <- list(statistic = NA, p.value = NA)
  try({ m11.VReductCensoring.normal <- gofTest_censored_normalCopula(x.mod.VReduct, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  try({ m11.VReductCensoring.t <- gofTest_censored_tCopula(x.mod.VReduct, N = N.gof, start = c(initial.theta[1], 4)) }, silent = TRUE)
  try({ m11.VReductCensoring.clayton <- gofTest_censored_claytonCopula(x.mod.VReduct, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  try({ m11.VReductCensoring.gumbel <- gofTest_censored_gumbelCopula(x.mod.VReduct, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  m11.VReductCensoring.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  # Method 12: G-Reduct + Censoring
  tm <- Sys.time()
  m12.GReductCensoring.estimate <- NA
  m12.GReductCensoring.stderror <- NA
  m12.GReductCensoring.CI95.lwr <- NA
  m12.GReductCensoring.CI95.upr <- NA
  x.mod.GReduct <- x.mod
  try({ x.mod.GReduct <- greedy_samplereduction_BIV(x.mod, n.greedy.25) }, silent = TRUE)
  try({ fit.mod <- censoredMPLE(copula=copula, method="BFGS", x.mod.GReduct, start=initial.theta) }, silent = TRUE)
  try({ m12.GReductCensoring.estimate <- fit.mod$par }, silent = TRUE)
  try({ m12.GReductCensoring.stderror <- NA }, silent = TRUE)
  try({ m12.GReductCensoring.CI95.lwr <- NA }, silent = TRUE)
  try({ m12.GReductCensoring.CI95.upr <- NA }, silent = TRUE)
  m12.GReductCensoring.fit.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]
  
  tm <- Sys.time()
  x.mod.GReduct <- x.mod
  try({ x.mod.GReduct <- greedy_samplereduction_BIV(x.mod, n.greedy.25) }, silent = TRUE)
  m12.GReductCensoring.normal <- list(statistic = NA, p.value = NA)
  m12.GReductCensoring.t <- list(statistic = NA, p.value = NA)
  m12.GReductCensoring.clayton <- list(statistic = NA, p.value = NA)
  m12.GReductCensoring.gumbel <- list(statistic = NA, p.value = NA)
  try({ m12.GReductCensoring.normal <- gofTest_censored_normalCopula(x.mod.GReduct, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  try({ m12.GReductCensoring.t <- gofTest_censored_tCopula(x.mod.GReduct, N = N.gof, start = c(initial.theta[1],4)) }, silent = TRUE)
  try({ m12.GReductCensoring.clayton <- gofTest_censored_claytonCopula(x.mod.GReduct, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  try({ m12.GReductCensoring.gumbel <- gofTest_censored_gumbelCopula(x.mod.GReduct, N = N.gof, start = initial.theta[1]) }, silent = TRUE)
  m12.GReductCensoring.gof.benchmark <- difftime(Sys.time(), tm, units = "secs")[[1]]

  iresult <- list(
    i, opt$copula, opt$tau, opt$sample, rp, op, opt$map,
    initial.theta[1],
    data.unm.tied_full, data.unm.tied_part, 
    data.mod.tied_full, data.mod.tied_part,
  
    m1.Reference.estimate,
    m1.Reference.stderror,
    m1.Reference.CI95.lwr,
    m1.Reference.CI95.upr,
    m1.Reference.fit.benchmark,
    m1.Reference.normal$statistic,
    m1.Reference.normal$p.value,
    m1.Reference.t$statistic,
    m1.Reference.t$p.value,
    m1.Reference.clayton$statistic,
    m1.Reference.clayton$p.value,
    m1.Reference.gumbel$statistic,
    m1.Reference.gumbel$p.value,
    m1.Reference.gof.benchmark,
    
    m2.VReduct.estimate,
    m2.VReduct.stderror,
    m2.VReduct.CI95.lwr,
    m2.VReduct.CI95.upr,
    m2.VReduct.fit.benchmark,
    m2.VReduct.normal$statistic,
    m2.VReduct.normal$p.value,
    m2.VReduct.t$statistic,
    m2.VReduct.t$p.value,
    m2.VReduct.clayton$statistic,
    m2.VReduct.clayton$p.value,
    m2.VReduct.gumbel$statistic,
    m2.VReduct.gumbel$p.value,
    m2.VReduct.gof.benchmark,

    m3.GReduct.estimate,
    m3.GReduct.stderror,
    m3.GReduct.CI95.lwr,
    m3.GReduct.CI95.upr,
    m3.GReduct.fit.benchmark,
    m3.GReduct.normal$statistic,
    m3.GReduct.normal$p.value,
    m3.GReduct.t$statistic,
    m3.GReduct.t$p.value,
    m3.GReduct.clayton$statistic,
    m3.GReduct.clayton$p.value,
    m3.GReduct.gumbel$statistic,
    m3.GReduct.gumbel$p.value,
    m3.GReduct.gof.benchmark,

    m4.RndRanks.estimate,
    m4.RndRanks.stderror,
    m4.RndRanks.CI95.lwr,
    m4.RndRanks.CI95.upr,
    m4.RndRanks.fit.benchmark,
    m4.RndRanks.normal$statistic,
    m4.RndRanks.normal$p.value,
    m4.RndRanks.t$statistic,
    m4.RndRanks.t$p.value,
    m4.RndRanks.clayton$statistic,
    m4.RndRanks.clayton$p.value,
    m4.RndRanks.gumbel$statistic,
    m4.RndRanks.gumbel$p.value,
    m4.RndRanks.gof.benchmark,   

    m5.PiJitter.estimate,
    m5.PiJitter.stderror,
    m5.PiJitter.CI95.lwr,
    m5.PiJitter.CI95.upr,
    m5.PiJitter.fit.benchmark,
    m5.PiJitter.normal$statistic,
    m5.PiJitter.normal$p.value,
    m5.PiJitter.t$statistic,
    m5.PiJitter.t$p.value,
    m5.PiJitter.clayton$statistic,
    m5.PiJitter.clayton$p.value,
    m5.PiJitter.gumbel$statistic,
    m5.PiJitter.gumbel$p.value,
    m5.PiJitter.gof.benchmark,   
    
    m6.TauJitter.estimate,
    m6.TauJitter.stderror,
    m6.TauJitter.CI95.lwr,
    m6.TauJitter.CI95.upr,
    m6.TauJitter.fit.benchmark,
    m6.TauJitter.normal$statistic,
    m6.TauJitter.normal$p.value,
    m6.TauJitter.t$statistic,
    m6.TauJitter.t$p.value,
    m6.TauJitter.clayton$statistic,
    m6.TauJitter.clayton$p.value,
    m6.TauJitter.gumbel$statistic,
    m6.TauJitter.gumbel$p.value,
    m6.TauJitter.gof.benchmark,   
    
    m7.Censoring.estimate,
    m7.Censoring.stderror,
    m7.Censoring.CI95.lwr,
    m7.Censoring.CI95.upr,
    m7.Censoring.fit.benchmark,
    m7.Censoring.normal$statistic,
    m7.Censoring.normal$p.value,
    m7.Censoring.t$statistic,
    m7.Censoring.t$p.value,
    m7.Censoring.clayton$statistic,
    m7.Censoring.clayton$p.value,
    m7.Censoring.gumbel$statistic,
    m7.Censoring.gumbel$p.value,
    m7.Censoring.gof.benchmark,   
    
    m8.Kojadinovic.estimate,
    m8.Kojadinovic.stderror,
    m8.Kojadinovic.CI95.lwr,
    m8.Kojadinovic.CI95.upr,
    m8.Kojadinovic.fit.benchmark,
    m8.Kojadinovic.normal$statistic,
    m8.Kojadinovic.normal$p.value,
    m8.Kojadinovic.t$statistic,
    m8.Kojadinovic.t$p.value,
    m8.Kojadinovic.clayton$statistic,
    m8.Kojadinovic.clayton$p.value,
    m8.Kojadinovic.gumbel$statistic,
    m8.Kojadinovic.gumbel$p.value,
    m8.Kojadinovic.gof.benchmark,   
    
    m9.VReductRndRanks.estimate,
    m9.VReductRndRanks.stderror,
    m9.VReductRndRanks.CI95.lwr,
    m9.VReductRndRanks.CI95.upr,
    m9.VReductRndRanks.fit.benchmark,  
    m9.VReductRndRanks.normal$statistic,
    m9.VReductRndRanks.normal$p.value,
    m9.VReductRndRanks.t$statistic,
    m9.VReductRndRanks.t$p.value,
    m9.VReductRndRanks.clayton$statistic,
    m9.VReductRndRanks.clayton$p.value,
    m9.VReductRndRanks.gumbel$statistic,
    m9.VReductRndRanks.gumbel$p.value,
    m9.VReductRndRanks.gof.benchmark,   
    
    m10.GReductRndRanks.estimate,
    m10.GReductRndRanks.stderror,
    m10.GReductRndRanks.CI95.lwr,
    m10.GReductRndRanks.CI95.upr,
    m10.GReductRndRanks.fit.benchmark,
    m10.GReductRndRanks.normal$statistic,
    m10.GReductRndRanks.normal$p.value,
    m10.GReductRndRanks.t$statistic,
    m10.GReductRndRanks.t$p.value,
    m10.GReductRndRanks.clayton$statistic,
    m10.GReductRndRanks.clayton$p.value,
    m10.GReductRndRanks.gumbel$statistic,
    m10.GReductRndRanks.gumbel$p.value,
    m10.GReductRndRanks.gof.benchmark,   
    
    m11.VReductCensoring.estimate,
    m11.VReductCensoring.stderror,
    m11.VReductCensoring.CI95.lwr,
    m11.VReductCensoring.CI95.upr,
    m11.VReductCensoring.fit.benchmark,
    m11.VReductCensoring.normal$statistic,
    m11.VReductCensoring.normal$p.value,
    m11.VReductCensoring.t$statistic,
    m11.VReductCensoring.t$p.value,
    m11.VReductCensoring.clayton$statistic,
    m11.VReductCensoring.clayton$p.value,
    m11.VReductCensoring.gumbel$statistic,
    m11.VReductCensoring.gumbel$p.value,
    m11.VReductCensoring.gof.benchmark,

    m12.GReductCensoring.estimate,
    m12.GReductCensoring.stderror,
    m12.GReductCensoring.CI95.lwr,
    m12.GReductCensoring.CI95.upr,
    m12.GReductCensoring.fit.benchmark,
    m12.GReductCensoring.normal$statistic,
    m12.GReductCensoring.normal$p.value,
    m12.GReductCensoring.t$statistic,
    m12.GReductCensoring.t$p.value,
    m12.GReductCensoring.clayton$statistic,
    m12.GReductCensoring.clayton$p.value,
    m12.GReductCensoring.gumbel$statistic,
    m12.GReductCensoring.gumbel$p.value,
    m12.GReductCensoring.gof.benchmark
  )

  try({ try(results[nrow(results) + 1, ] <- iresult) })

  if (i %% ceiling(N/10) == 0) {
    delta <- difftime(Sys.time(), ssts, units = "secs")[[1]]
    deltas <- c(deltas, delta)

    spi <- sum(deltas) / i
    remaining <- spi * (N - i)
    ssts <- Sys.time()

    print(paste0(i, " iterations in ", round(delta, 2), " s. ETA: ", (ssts + remaining)))
      
    write.csv(x <- results, file <- paste0(
        opt$path2, "temp_", round(i / N * 100), "perc_GoF_", opt$copula, "_", opt$tau, "_N", N,
        "_n", n,"_rp", rp,"_op", op, "_em", opt$map, ".csv")
    )
  }
}

print("Saving the results ...")
write.csv(
  x <- results,
  file <- paste0(
      opt$path, "GoF_", opt$copula, "_", opt$tau, "_N", N,
      "_n", n,"_rp", rp,"_op", op, "_em", opt$map, ".csv")
)

print(paste0(
  "Stopping the simulation at ", Sys.time(), "."
))