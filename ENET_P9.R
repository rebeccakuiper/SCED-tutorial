################################################################################
# SCED tutorial - GORICA evaluation of trajectory hypotheses for participant P9
# Multivariate regression (two outcomes: PCL-5 and CLL) fitted with lm()
################################################################################


## 0. Packages ----------------------------------------------------------------
# Install once (uncomment if needed):
# install.packages("restriktor")
library(restriktor)  # goric(), benchmark()
packageVersion("restriktor")  # report this version in the paper

# The data file is read relative to the working directory: 
# put this script in a folder with a subfolder 'data' containing 'ENET09.dat' (check with getwd()).


## 1. Data preparation ---------------------------------------------------------
data_ENET_P9 <- read.table("data/ENET09.dat", header = TRUE)  # first line = variable names
head(data_ENET_P9)

# Missing values are coded -999
data_ENET_P9[data_ENET_P9 == -999] <- NA

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
#data_ENET_P9$Time   <- data_ENET_P9$Time / 7
# (b):
data_ENET_P9$Time   <- data_ENET_P9$Time - quantile(data_ENET_P9$Time, 1/3)
data_ENET_P9$sqTime <- data_ENET_P9$Time^2


## 2. Fit the unconstrained model (Step 1) --------------------------------------
# Multivariate multiple regression: 
# two outcomes, and for each the same three predictors.
lm_fit_P9 <- lm(cbind(PCL5, CLL) ~ Time + sqTime + event, data = data_ENET_P9)
#nobs(lm_fit_P9)  # number of observations used (rows with NA are dropped)
#
# lm_fit_P9
# summary(lm_fit_P9)
# Do not inspect the estimates before the hypotheses are specified.

# When creating this file, goric() cannot yet take a multivariate lm object, 
# so pass the estimates and their covariance matrix. 
# Label them by name, not by position:
est <- coef(lm_fit_P9)  # 4 x 2 matrix: rows = predictors, columns = outcomes
est <- setNames(c(est), c(paste0("beta_PCL5_", 0:3), paste0("beta_CLL_", 0:3)))
# 0 = intercept, 1 = Time (linear), 2 = sqTime (quadratic), 3 = event
VCOV <- vcov(lm_fit_P9)
#stopifnot(identical(length(est), nrow(VCOV)))
#dimnames(VCOV) <- list(names(est), names(est))


## 3. Hypotheses (Step 2) -------------------------------------------------------
# Restrictions are separated by ';', ',', and/or '&'.
# Here, '&' is used to separate the restrictions per outcome,
#       ';' is used to separate the restrictions within one outcome.

H_NoTrend <- "beta_PCL5_1 = 0; beta_PCL5_2 = 0 &
              beta_CLL_1  = 0; beta_CLL_2  = 0"
# no linear and no quadratic effect of Time on either outcome

H_Linear <- "beta_PCL5_1 < 0; beta_PCL5_2 = 0 &
             beta_CLL_1  > 0; beta_CLL_2  = 0"
# PCL-5 decreases linearly, CLL increases linearly

H_Parabolic <- "beta_PCL5_2 < 0 &
                 beta_CLL_1 > 0; beta_CLL_2 = 0"
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
t_start <- min(data_ENET_P9$Time)
t_end   <- max(data_ENET_P9$Time)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0;
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0;
                       beta_CLL_1 > 0; beta_CLL_2 = 0", 2 * t_start, 2 * t_end)

# The unconstrained hypothesis is added by default as fail-safe.


## 4. GORICA evaluation (Step 3) -------------------------------------------------
# The penalty is computed by simulation, so set a seed for reproducibility;
# re-run with another seed to check that the penalties are stable
# (if not, increase the number of simulation iterations; see ?goric).
set.seed(123)
output_P9 <- goric(est, VCOV = VCOV,
                   hypotheses = list(H_NoTrend    = H_NoTrend,
                                     H_Linear     = H_Linear,
                                     H_Parabolic = H_Parabolic))
output_P9
# summary(output_P9)   # more detail
# output_P9$ratio.gw   # matrix with ratios of GORICA weights

# Always check what the preferred hypothesis looks like: 
# the order-restricted estimates show whether the restrictions are active 
# (estimate on the boundary).
round(coef(output_P9)[, c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")], 4)
# E.g., if H_Parabolic is preferred but has beta_PCL5_2 = 0 (boundary),
# it does not describe deterioration followed by improvement.

# Sensitivity check with the stricter rise-then-fall hypothesis:
# set.seed(123)
# output_P9_rf <- goric(est, VCOV = VCOV,
#                       hypotheses = list(H_NoTrend    = H_NoTrend,
#                                         H_Linear     = H_Linear,
#                                         H_Parabolic  = H_Parabolic,
#                                         H_RiseFall   = H_RiseFall))
# output_P9_rf
#
# Because H_RiseFall is probably an alternative specification of H_Parabolic,
# use the set with only H_RiseFall (not both):
set.seed(123)
output_P9_rf <- goric(est, VCOV = VCOV,
                      hypotheses = list(H_NoTrend    = H_NoTrend,
                                        H_Linear     = H_Linear,
                                        H_RiseFall   = H_RiseFall))
output_P9_rf


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
benchmarks_P9 <- benchmark(output_P9, pop_est = my_pop_est,
                           iter = 100)  # 100 to save time; use >= 1000 for reporting
benchmarks_P9
plot(benchmarks_P9, log_scale = TRUE)


## 6.
sessionInfo()  # record R and package versions
