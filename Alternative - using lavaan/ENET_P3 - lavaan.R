################################################################################
# SCED tutorial - GORICA evaluation for participant P3, model fitted with lavaan
# (alternative to the lm() script; same model, same hypotheses)
################################################################################


## 0. Packages ----------------------------------------------------------------
# install.packages(c("restriktor", "lavaan"))
library(restriktor)
library(lavaan)


## 1. Data preparation (identical to the lm() script) ------------------------------
data_ENET_P3 <- read.table("data/ENET03.dat", header = TRUE)
data_ENET_P3[data_ENET_P3 == -999] <- NA

# Keep 'event' numeric (0/1). 
# lavaan converts a factor to 1/2 codes, 
# which shifts the intercepts relative to lm().
#
# Time in weeks (numerical stability), centered at 1/3:
# (a) Express it in weeks. In days, sqTime runs up to ~16,000 and its
#     coefficients are ~1e-4, which makes the covariance matrix nearly singular.
#     Rescaling changes no weights, only the units of the coefficients.
data_ENET_P3$Time   <- data_ENET_P3$Time / 7
# (b) Centre it at one-third of the intervention phase. 
#     Btw Centring does NOT determine where the maximum/minimum of the parabola lies; 
#     it only makes beta_1 the slope at that time point (beta_2 is unaffected).
data_ENET_P3$Time   <- data_ENET_P3$Time - quantile(data_ENET_P3$Time, 1/3, na.rm = TRUE)
data_ENET_P3$sqTime <- data_ENET_P3$Time^2


## 2. Fit the unconstrained model (Step 1) --------------------------------------
model <- '
  PCL5 ~ beta_PCL5_1*Time + beta_PCL5_2*sqTime + beta_PCL5_3*event
  CLL  ~ beta_CLL_1*Time  + beta_CLL_2*sqTime  + beta_CLL_3*event
'
# 1 = Time (linear), 2 = sqTime (quadratic), 3 = event.
#
# Residuals of PCL5 and CLL are correlated by default (as in the lm() model).
fit_P3 <- sem(model, data = data_ENET_P3, meanstructure = TRUE)


## 3. Hypotheses (Step 2) -------------------------------------------------------
# (',' separates restrictions within one outcome, '&' those of the two outcomes)
H_NoTrend    <- "beta_PCL5_1 = 0, beta_PCL5_2 = 0 & beta_CLL_1 = 0, beta_CLL_2 = 0"
H_Linear     <- "beta_PCL5_1 < 0, beta_PCL5_2 = 0 & beta_CLL_1 > 0, beta_CLL_2 = 0"
H_Parabolic  <- "beta_PCL5_2 < 0 & beta_CLL_1 > 0, beta_CLL_2 = 0"
# Optional:
t_start <- min(data_ENET_P3$Time, na.rm = TRUE)
t_end   <- max(data_ENET_P3$Time, na.rm = TRUE)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0,
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0 &
                       beta_CLL_1 > 0, beta_CLL_2 = 0", 2 * t_start, 2 * t_end)


## 4. GORICA evaluation (Step 3) -------------------------------------------------
# Only call goric(fit_P3, ...) with restriktor version 0.7 or higher.
# Workaround: extract the labelled estimates and their covariance matrix.
betas <- c("beta_PCL5_1", "beta_PCL5_2", "beta_PCL5_3",
           "beta_CLL_1",  "beta_CLL_2",  "beta_CLL_3")
est  <- coef(fit_P3)[betas]
VCOV <- vcov(fit_P3)[betas, betas]
# Note: lavaan (ML) divides by N instead of N - k, 
# so its standard errors are smaller than lm()'s (here variances x 21/25) 
# and the weights differ somewhat from the lm() results reported in the paper.

set.seed(123)
output_P3 <- goric(est, VCOV = VCOV,
                   hypotheses = list(H_NoTrend    = H_NoTrend,
                                     H_Linear     = H_Linear,
                                     H_Parabolic  = H_Parabolic))
output_P3

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
output_P3$result[, c("model", "loglik", "loglik.weights")]
# 2. or-mle's per hypothesis (last row: the unconstrained estimates):
slopes <- c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")
signif(as.matrix(coef(output_P3))[, slopes], 3)
#
# Here (P3): H_NoTrend is preferred, but its loglik (weight) is clearly lower
# than that of the unconstrained hypothesis (loglik 6.050 vs. 8.460; loglik weights 0.061 vs. 0.682).
# H_NoTrend only contains equality restrictions (all trend parameters 0),
# so its or-mle's are 0 by definition. It is preferred because it is the
# least complex hypothesis, not because it fits (about) as well as the
# unconstrained one: the data show some change over time, but too little
# (given the few observations) to support a more complex hypothesis.
# So we may not have enough power to detect a (small) trend, that is, to
# prefer the fail-safe unconstrained hypothesis.

#Optional / Alternative:
set.seed(123)
output_P3_rf <- goric(est, VCOV = VCOV,
                      hypotheses = list(H_NoTrend    = H_NoTrend,
                                        H_Linear     = H_Linear,
                                        H_RiseFall   = H_RiseFall))
output_P3_rf


## 5. Extra: benchmarks (parameters set by name) ----------------------------------
est_no  <- replace(est, c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2"), 0)
est_lin <- replace(est, c("beta_PCL5_2", "beta_CLL_2"), 0)
my_pop_est <- rbind(NoTrend = est_no, Linear = est_lin, Observed = est)

set.seed(123)
benchmarks_P3 <- benchmark(output_P3, pop_est = my_pop_est, iter = 100)
benchmarks_P3
