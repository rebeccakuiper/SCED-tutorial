################################################################################
# SCED tutorial - GORICA evaluation of trajectory hypotheses for participant P7
# Multivariate regression (two outcomes: PCL-5 and CLL) fitted with lm()
################################################################################


## 0. Packages ----------------------------------------------------------------
# Install once (uncomment if needed):
# install.packages("restriktor")
library(restriktor)  # goric(), benchmark()
packageVersion("restriktor")  # report this version in the paper

# The data file is read relative to the working directory: 
# put this script in a folder with a subfolder 'data' containing 'ENET07.dat' (check with getwd()).


## 1. Data preparation ---------------------------------------------------------
data_ENET_P7 <- read.table("data/ENET07.dat", header = TRUE)  # first line = variable names
head(data_ENET_P7)

# Missing values are coded -999
data_ENET_P7[data_ENET_P7 == -999] <- NA

# 'event' is a 0/1 indicator (0 = no critical incident, 1 = one or more).
# Kept numeric: for a 0/1 variable this gives the same model as factor().

# Time: measured in days since the first intervention session.
# (a) Optional: Express it in weeks. In days, sqTime runs up to ~16,000 and its
#     coefficients are ~1e-4, which makes the covariance matrix nearly singular.
#     Rescaling changes no weights, only the units of the coefficients.
# (b) Centre it at one-third of the intervention phase. 
#     Btw Centring does NOT determine where the maximum/minimum of the parabola lies; 
#     it only makes beta_1 the slope at that time point (beta_2 is unaffected).
#
# (a):
#data_ENET_P7$Time   <- data_ENET_P7$Time / 7
# (b):
data_ENET_P7$Time   <- data_ENET_P7$Time - quantile(data_ENET_P7$Time, 1/3, na.rm = TRUE)
data_ENET_P7$sqTime <- data_ENET_P7$Time^2


## 2. Fit the unconstrained model (Step 1) --------------------------------------
# Multivariate multiple regression: 
# two outcomes, and for each the same three predictors.
lm_fit_P7 <- lm(cbind(PCL5, CLL) ~ Time + sqTime + event, data = data_ENET_P7)
#nobs(lm_fit_P7)  # number of observations used (rows with NA are dropped)
#
# lm_fit_P7
# summary(lm_fit_P7)
# Do not inspect the estimates before the hypotheses are specified.

# When creating this file, goric() cannot yet take a multivariate lm object, 
# so pass the estimates and their covariance matrix. 
# Label them by name, not by position:
est <- coef(lm_fit_P7)  # 4 x 2 matrix: rows = predictors, columns = outcomes
est <- setNames(c(est), c(paste0("beta_PCL5_", 0:3), paste0("beta_CLL_", 0:3)))
# 0 = intercept, 1 = Time (linear), 2 = sqTime (quadratic), 3 = event
VCOV <- vcov(lm_fit_P7)
#stopifnot(identical(length(est), nrow(VCOV)))
#dimnames(VCOV) <- list(names(est), names(est))


## 3. Hypotheses (Step 2) -------------------------------------------------------
# Restrictions can be separated by ',', ';', and/or '&'.
# Here, ',' is used to separate the restrictions within one outcome,
#       '&' is used to separate the restrictions of the two outcomes.

H_NoTrend <- "beta_PCL5_1 = 0, beta_PCL5_2 = 0 &
              beta_CLL_1  = 0, beta_CLL_2  = 0"
# no linear and no quadratic effect of Time on either outcome

H_Linear <- "beta_PCL5_1 < 0, beta_PCL5_2 = 0 &
             beta_CLL_1  > 0, beta_CLL_2  = 0"
# PCL-5 decreases linearly, CLL increases linearly

H_Parabolic <- "beta_PCL5_2 < 0 &
                 beta_CLL_1 > 0, beta_CLL_2 = 0"
# PCL-5: inverted-U (negative quadratic effect); beta_PCL5_1 is unrestricted.
#        NB: this also covers a parabola whose peak lies outside the observed period,
#        i.e. a PCL-5 trend that only rises or only falls, 
#        and (on the boundary beta_PCL5_2 = 0) a purely linear rise. 
#        See H_RiseFall for a stricter version.
# CLL: linear increase.

# Optional, stricter version of 'first deterioration, then improvement':
# the PCL-5 slope at time t is beta_PCL5_1 + 2*t*beta_PCL5_2, 
# so require it to be positive at the start and 
# negative at the end of the intervention phase.
t_start <- min(data_ENET_P7$Time, na.rm = TRUE)
t_end   <- max(data_ENET_P7$Time, na.rm = TRUE)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0,
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0 &
                       beta_CLL_1 > 0, beta_CLL_2 = 0", 2 * t_start, 2 * t_end)

# The unconstrained hypothesis is added by default as fail-safe.


## 4. GORICA evaluation (Step 3) -------------------------------------------------
# The penalty is computed by simulation, so set a seed for reproducibility;
# re-run with another seed to check that the penalties are stable
# (if not, increase the number of simulation iterations; see ?goric).
set.seed(123)
output_P7 <- goric(est, VCOV = VCOV,
                   hypotheses = list(H_NoTrend    = H_NoTrend,
                                     H_Linear     = H_Linear,
                                     H_Parabolic = H_Parabolic))
output_P7
# summary(output_P7)   # more detail
# output_P7$ratio.gw   # matrix with ratios of GORICA weights

# Always check the log-likelihood (loglik) values and/or the loglik weights
# of the preferred hypothesis (see 'Guidelines output GORIC.html' and
# 'Decision tree - Interpretation.pdf' on https://github.com/rebeccakuiper/Tutorials):
# - The unconstrained hypothesis always has the highest loglik.
# - If the preferred hypothesis has (about) the same loglik (weight) as the
#   unconstrained one, its restrictions are in line with the data.
# - If its loglik (weight) is (clearly) lower, it is preferred (partly) because
#   it is less complex. Since the unconstrained (fail-safe) hypothesis is
#   always included in the set, the issue is then that we may not have enough
#   power to detect that the restrictions do not hold, that is, to prefer the
#   fail-safe hypothesis.
#   Inspect this further with the order-restricted estimates (or-mle's): they
#   show which restrictions are active, i.e., which estimates are on the
#   boundary.
#
# 1. loglik values and loglik weights:
output_P7$result[, c("model", "loglik", "loglik.weights")]
# 2. or-mle's per hypothesis (last row: the unconstrained estimates):
slopes <- c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")
signif(as.matrix(coef(output_P7))[, slopes], 3)
#
# Here (P7): the unconstrained hypothesis is preferred, so all three
# informative hypotheses are weak. Even the best-fitting one (H_Linear) has a
# clearly lower loglik (-5.534 vs. 23.568). None of the hypothesized trajectories
# describes the data well; comparing the or-mle's with the unconstrained
# estimates (last row) shows which restrictions conflict with the data.
# Another trajectory (not among the informative hypotheses) may better
# describe the data.

# Sensitivity check with the stricter rise-then-fall hypothesis:
# set.seed(123)
# output_P7_rf <- goric(est, VCOV = VCOV,
#                       hypotheses = list(H_NoTrend    = H_NoTrend,
#                                         H_Linear     = H_Linear,
#                                         H_Parabolic  = H_Parabolic,
#                                         H_RiseFall   = H_RiseFall))
# output_P7_rf
#
# Because H_RiseFall is probably an alternative specification of H_Parabolic,
# use the set with only H_RiseFall (not both):
set.seed(123)
output_P7_rf <- goric(est, VCOV = VCOV,
                      hypotheses = list(H_NoTrend    = H_NoTrend,
                                        H_Linear     = H_Linear,
                                        H_RiseFall   = H_RiseFall))
output_P7_rf


## 5. Extra: benchmarks -----------------------------------------------------------
# How (un)usual is the observed support for the preferred hypothesis under
# hypothetical populations? Parameters are set by name, not by position.
slopes   <- c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")
est_no   <- replace(est, slopes, 0)                         # no trend
est_lin  <- replace(est, c("beta_PCL5_2", "beta_CLL_2"), 0) # linear trend only
my_pop_est <- rbind(NoTrend  = est_no,
                    Linear   = est_lin,
                    Observed = est)

set.seed(123)
benchmarks_P7 <- benchmark(output_P7, pop_est = my_pop_est,
                           iter = 100)  # 100 to save time; use >= 1000 for reporting
benchmarks_P7
plot(benchmarks_P7, log_scale = TRUE)


## 6.
sessionInfo()  # record R and package versions
