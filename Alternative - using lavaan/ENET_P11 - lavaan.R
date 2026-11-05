################################################################################
# SCED tutorial - GORICA evaluation for participant P11, model fitted with lavaan
# (alternative to the lm() script; same model, same hypotheses)
################################################################################


## 0. Packages ----------------------------------------------------------------
# install.packages(c("restriktor", "lavaan"))
library(restriktor)
library(lavaan)


## 1. Data preparation (identical to the lm() script) ------------------------------
data_ENET_P11 <- read.table("data/ENET11.dat", header = TRUE)
data_ENET_P11[data_ENET_P11 == -999] <- NA

# Keep 'event' numeric (0/1). 
# lavaan converts a factor to 1/2 codes, 
# which shifts the intercepts relative to lm().
#
# Time in weeks (numerical stability), centered at 1/3:
# (a) Express it in weeks. In days, sqTime runs up to ~16,000 and its
#     coefficients are ~1e-4, which makes the covariance matrix nearly singular.
#     Rescaling changes no weights, only the units of the coefficients.
data_ENET_P11$Time   <- data_ENET_P11$Time / 7
# (b) Centre it at one-third of the intervention phase. 
#     Btw Centring does NOT determine where the maximum/minimum of the parabola lies; 
#     it only makes beta_1 the slope at that time point (beta_2 is unaffected).
data_ENET_P11$Time   <- data_ENET_P11$Time - quantile(data_ENET_P11$Time, 1/3, na.rm = TRUE)
data_ENET_P11$sqTime <- data_ENET_P11$Time^2


## 2. Fit the unconstrained model (Step 1) --------------------------------------
model <- '
  PCL5 ~ beta_PCL5_1*Time + beta_PCL5_2*sqTime + beta_PCL5_3*event
  CLL  ~ beta_CLL_1*Time  + beta_CLL_2*sqTime  + beta_CLL_3*event
'
# 1 = Time (linear), 2 = sqTime (quadratic), 3 = event.
#
# Residuals of PCL5 and CLL are correlated by default (as in the lm() model).
fit_P11 <- sem(model, data = data_ENET_P11, meanstructure = TRUE)


## 3. Hypotheses (Step 2) -------------------------------------------------------
# (',' separates restrictions within one outcome, '&' those of the two outcomes)
H_NoTrend    <- "beta_PCL5_1 = 0, beta_PCL5_2 = 0 & beta_CLL_1 = 0, beta_CLL_2 = 0"
H_Linear     <- "beta_PCL5_1 < 0, beta_PCL5_2 = 0 & beta_CLL_1 > 0, beta_CLL_2 = 0"
H_Parabolic  <- "beta_PCL5_2 < 0 & beta_CLL_1 > 0, beta_CLL_2 = 0"
# Optional:
t_start <- min(data_ENET_P11$Time, na.rm = TRUE)
t_end   <- max(data_ENET_P11$Time, na.rm = TRUE)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0,
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0 &
                       beta_CLL_1 > 0, beta_CLL_2 = 0", 2 * t_start, 2 * t_end)


## 4. GORICA evaluation (Step 3) -------------------------------------------------
# Only call goric(fit_P11, ...) with restriktor version 0.7 or higher.
# Workaround: extract the labelled estimates and their covariance matrix.
betas <- c("beta_PCL5_1", "beta_PCL5_2", "beta_PCL5_3",
           "beta_CLL_1",  "beta_CLL_2",  "beta_CLL_3")
est  <- coef(fit_P11)[betas]
VCOV <- vcov(fit_P11)[betas, betas]
# Note: lavaan (ML) divides by N instead of N - k, 
# so its standard errors are smaller than lm()'s (here variances x 21/25) 
# and the weights differ somewhat from the lm() results reported in the paper.

set.seed(123)
output_P11 <- goric(est, VCOV = VCOV,
                   hypotheses = list(H_NoTrend    = H_NoTrend,
                                     H_Linear     = H_Linear,
                                     H_Parabolic  = H_Parabolic))
output_P11

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
output_P11$result[, c("model", "loglik", "loglik.weights")]
# 2. or-mle's per hypothesis (last row: the unconstrained estimates):
slopes <- c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")
signif(as.matrix(coef(output_P11))[, slopes], 3)
#
# Here (P11): H_Linear is preferred. Its loglik (weight) is lower than that of
# the unconstrained hypothesis (loglik 11.275 vs. 12.176; loglik weights 0.219 vs. 0.540).
# The or-mle's show that its inequality restrictions are not active
# (beta_PCL5_1 = -0.465, beta_CLL_1 = 0.0363; not on the boundary), so the expected
# directions are in line with the data. The loss in fit comes from the
# equality restriction(s) (no quadratic effects for PCL-5 and CLL), which set those
# estimates to zero.
# Since this loss is small, H_Linear describes the data well.

#Optional / Alternative:
set.seed(123)
output_P11_rf <- goric(est, VCOV = VCOV,
                      hypotheses = list(H_NoTrend    = H_NoTrend,
                                        H_Linear     = H_Linear,
                                        H_RiseFall   = H_RiseFall))
output_P11_rf


## 5. Extra: benchmarks (parameters set by name) ----------------------------------
est_no  <- replace(est, c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2"), 0)
est_lin <- replace(est, c("beta_PCL5_2", "beta_CLL_2"), 0)
my_pop_est <- rbind(NoTrend = est_no, Linear = est_lin, Observed = est)

set.seed(123)
benchmarks_P11 <- benchmark(output_P11, pop_est = my_pop_est, iter = 100)
benchmarks_P11
