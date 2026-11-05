################################################################################
# SCED tutorial - GORICA evaluation of trajectory hypotheses for participant P1
# Sensitivity analysis: same multivariate regression as in ENET_P1.R, but with
# continuous-time AR(1) (autocorrelated) residuals, estimated by REML.
#
# Only the residual covariance changes; the trend parameters (and thus the
# hypotheses) have the same meaning as in the lm() analysis.
################################################################################


## 0. Packages and helper function ----------------------------------------------
library(restriktor)            # goric()
source("CAR1_functions.R")     # fit_car1_biv(); same folder as this script


## 1. Data preparation (identical to ENET_P1.R) -----------------------------------
data_ENET_P1 <- read.table("data/ENET01.dat", header = TRUE)
data_ENET_P1[data_ENET_P1 == -999] <- NA
#
data_ENET_P1$Time   <- data_ENET_P1$Time - quantile(data_ENET_P1$Time, 1/3, na.rm = TRUE)
data_ENET_P1$sqTime <- data_ENET_P1$Time^2


## 2. Fit the unconstrained model with CAR(1) residuals (Step 1) ------------------
# Residual covariance: Cov(e_js, e_kt) = Sigma[j,k] * phi^|t - s|, with the lag
# in weeks (time_unit = 7 days). Time itself is still in days, as in ENET_P1.R.
car1_fit_P1 <- fit_car1_biv(~ Time + sqTime + event, data = data_ENET_P1,
                            outcomes = c("PCL5", "CLL"), time = "Time",
                            time_unit = 7)
car1_fit_P1$nobs        # observations used (complete cases)
car1_fit_P1$phi         # autocorrelation at a lag of one week
car1_fit_P1$rho         # contemporaneous correlation between PCL5 and CLL residuals
car1_fit_P1$LR_AR       # LR statistic: CAR(1) versus independent residuals
#                         (phi = 0 lies on the boundary: conservative test)

# Estimates and covariance matrix, labelled as in ENET_P1.R
est  <- car1_fit_P1$coef
names(est) <- c(paste0("beta_PCL5_", 0:3), paste0("beta_CLL_", 0:3))
VCOV <- car1_fit_P1$vcov
dimnames(VCOV) <- list(names(est), names(est))


## 3. Hypotheses (Step 2; identical to ENET_P1.R) ---------------------------------
H_NoTrend <- "beta_PCL5_1 = 0, beta_PCL5_2 = 0 &
              beta_CLL_1  = 0, beta_CLL_2  = 0"
H_Linear <- "beta_PCL5_1 < 0, beta_PCL5_2 = 0 &
             beta_CLL_1  > 0, beta_CLL_2  = 0"
H_Parabolic <- "beta_PCL5_2 < 0 &
                beta_CLL_1 > 0, beta_CLL_2 = 0"
#
# Start and end of the intervention phase (all sessions, as in ENET_P1.R)
t_start <- min(data_ENET_P1$Time, na.rm = TRUE)
t_end   <- max(data_ENET_P1$Time, na.rm = TRUE)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0,
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0 &
                       beta_CLL_1 > 0, beta_CLL_2 = 0", 2 * t_start, 2 * t_end)


## 4. GORICA evaluation (Step 3) ---------------------------------------------------
set.seed(123)
output_P1_CAR1 <- goric(est, VCOV = VCOV,
                        hypotheses = list(H_NoTrend   = H_NoTrend,
                                          H_Linear    = H_Linear,
                                          H_Parabolic = H_Parabolic))
output_P1_CAR1

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
output_P1_CAR1$result[, c("model", "loglik", "loglik.weights")]
# 2. or-mle's per hypothesis (last row: the unconstrained estimates):
slopes <- c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")
signif(as.matrix(coef(output_P1_CAR1))[, slopes], 3)
#
# Here (P1): H_Parabolic is preferred, but its loglik (weight) is lower than
# that of the unconstrained hypothesis (loglik 19.743 vs. 20.094; loglik weights 0.385 vs. 0.547).
# So H_Parabolic may be preferred mainly because it is less complex, and we
# may not have enough power to detect that its restrictions do not hold (that
# is, to prefer the fail-safe unconstrained hypothesis). The or-mle's show
# that the restriction beta_PCL5_2 < 0 is active: under H_Parabolic,
# beta_PCL5_2 = 0 (on the boundary; the unconstrained estimate, 0.000489, is positive, i.e., a U-shape),
# with beta_PCL5_1 = 0.0612 > 0. Hence, the best-fitting trajectory within
# H_Parabolic is a linear INCREASE in PCL-5, not a deterioration followed by
# improvement: H_Parabolic is preferred because it is the only hypothesis
# that allows PCL-5 to increase, not because the data show an inverted U.
# The stricter H_RiseFall (below) evaluates a rise followed by a fall.

# Alternative set with the stricter rise-then-fall hypothesis
set.seed(123)
output_P1_CAR1_rf <- goric(est, VCOV = VCOV,
                           hypotheses = list(H_NoTrend  = H_NoTrend,
                                             H_Linear   = H_Linear,
                                             H_RiseFall = H_RiseFall))
output_P1_CAR1_rf
