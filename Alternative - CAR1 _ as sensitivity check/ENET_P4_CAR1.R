################################################################################
# SCED tutorial - GORICA evaluation of trajectory hypotheses for participant P4
# Sensitivity analysis: same multivariate regression as in ENET_P4.R, but with
# continuous-time AR(1) (autocorrelated) residuals, estimated by REML.
#
# Only the residual covariance changes; the trend parameters (and thus the
# hypotheses) have the same meaning as in the lm() analysis.
################################################################################


## 0. Packages and helper function ----------------------------------------------
library(restriktor)            # goric()
source("CAR1_functions.R")     # fit_car1_biv(); same folder as this script


## 1. Data preparation (identical to ENET_P4.R) -----------------------------------
data_ENET_P4 <- read.table("data/ENET04.dat", header = TRUE)
data_ENET_P4[data_ENET_P4 == -999] <- NA
#
data_ENET_P4$Time   <- data_ENET_P4$Time - quantile(data_ENET_P4$Time, 1/3, na.rm = TRUE)
data_ENET_P4$sqTime <- data_ENET_P4$Time^2


## 2. Fit the unconstrained model with CAR(1) residuals (Step 1) ------------------
# Residual covariance: Cov(e_js, e_kt) = Sigma[j,k] * phi^|t - s|, with the lag
# in weeks (time_unit = 7 days). Time itself is still in days, as in ENET_P4.R.
car1_fit_P4 <- fit_car1_biv(~ Time + sqTime + event, data = data_ENET_P4,
                            outcomes = c("PCL5", "CLL"), time = "Time",
                            time_unit = 7)
car1_fit_P4$nobs        # observations used (complete cases)
car1_fit_P4$phi         # autocorrelation at a lag of one week
car1_fit_P4$rho         # contemporaneous correlation between PCL5 and CLL residuals
car1_fit_P4$LR_AR       # LR statistic: CAR(1) versus independent residuals
#                         (phi = 0 lies on the boundary: conservative test)

# Estimates and covariance matrix, labelled as in ENET_P4.R
est  <- car1_fit_P4$coef
names(est) <- c(paste0("beta_PCL5_", 0:3), paste0("beta_CLL_", 0:3))
VCOV <- car1_fit_P4$vcov
dimnames(VCOV) <- list(names(est), names(est))


## 3. Hypotheses (Step 2; identical to ENET_P4.R) ---------------------------------
H_NoTrend <- "beta_PCL5_1 = 0, beta_PCL5_2 = 0 &
              beta_CLL_1  = 0, beta_CLL_2  = 0"
H_Linear <- "beta_PCL5_1 < 0, beta_PCL5_2 = 0 &
             beta_CLL_1  > 0, beta_CLL_2  = 0"
H_Parabolic <- "beta_PCL5_2 < 0 &
                beta_CLL_1 > 0, beta_CLL_2 = 0"
#
# Start and end of the intervention phase (all sessions, as in ENET_P4.R)
t_start <- min(data_ENET_P4$Time, na.rm = TRUE)
t_end   <- max(data_ENET_P4$Time, na.rm = TRUE)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0,
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0 &
                       beta_CLL_1 > 0, beta_CLL_2 = 0", 2 * t_start, 2 * t_end)


## 4. GORICA evaluation (Step 3) ---------------------------------------------------
set.seed(123)
output_P4_CAR1 <- goric(est, VCOV = VCOV,
                        hypotheses = list(H_NoTrend   = H_NoTrend,
                                          H_Linear    = H_Linear,
                                          H_Parabolic = H_Parabolic))
output_P4_CAR1

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
output_P4_CAR1$result[, c("model", "loglik", "loglik.weights")]
# 2. or-mle's per hypothesis (last row: the unconstrained estimates):
slopes <- c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")
signif(as.matrix(coef(output_P4_CAR1))[, slopes], 3)
#
# Here (P4): H_Parabolic is preferred. Its loglik (weight) is clearly lower than that of
# the unconstrained hypothesis (loglik 17.626 vs. 19.469; loglik weights 0.135 vs. 0.854).
# The or-mle's show that its inequality restrictions are not active
# (beta_PCL5_2 = -0.00215, beta_CLL_1 = 0.00596; not on the boundary), so the expected
# directions are in line with the data. The loss in fit comes from the
# equality restriction(s) (no quadratic effect for CLL), which set those
# estimates to zero.
# Since this loss is considerable, H_Parabolic is preferred mainly because it is
# less complex: the expected directions are supported, but the hypothesis
# may not fully capture the trajectory (e.g., a quadratic effect).

# Alternative set with the stricter rise-then-fall hypothesis
set.seed(123)
output_P4_CAR1_rf <- goric(est, VCOV = VCOV,
                           hypotheses = list(H_NoTrend  = H_NoTrend,
                                             H_Linear   = H_Linear,
                                             H_RiseFall = H_RiseFall))
output_P4_CAR1_rf
