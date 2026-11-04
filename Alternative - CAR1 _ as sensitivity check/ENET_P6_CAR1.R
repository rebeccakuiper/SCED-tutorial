################################################################################
# SCED tutorial - GORICA evaluation of trajectory hypotheses for participant P6
# Sensitivity analysis: same multivariate regression as in ENET_P6.R, but with
# continuous-time AR(1) (autocorrelated) residuals, estimated by REML.
#
# Only the residual covariance changes; the trend parameters (and thus the
# hypotheses) have the same meaning as in the lm() analysis.
################################################################################


## 0. Packages and helper function ----------------------------------------------
library(restriktor)            # goric()
source("ENET_CAR1_sensitivity/CAR1_functions.R")     # fit_car1_biv(); same folder as this script


## 1. Data preparation (identical to ENET_P6.R) -----------------------------------
data_ENET_P6 <- read.table("data/ENET06.dat", header = TRUE)
data_ENET_P6[data_ENET_P6 == -999] <- NA
#
data_ENET_P6$Time   <- data_ENET_P6$Time - quantile(data_ENET_P6$Time, 1/3)
data_ENET_P6$sqTime <- data_ENET_P6$Time^2


## 2. Fit the unconstrained model with CAR(1) residuals (Step 1) ------------------
# Residual covariance: Cov(e_js, e_kt) = Sigma[j,k] * phi^|t - s|, with the lag
# in weeks (time_unit = 7 days). Time itself is still in days, as in ENET_P6.R.
car1_fit_P6 <- fit_car1_biv(~ Time + sqTime + event, data = data_ENET_P6,
                            outcomes = c("PCL5", "CLL"), time = "Time",
                            time_unit = 7)
car1_fit_P6$nobs        # observations used (complete cases)
car1_fit_P6$phi         # autocorrelation at a lag of one week
car1_fit_P6$rho         # contemporaneous correlation between PCL5 and CLL residuals
car1_fit_P6$LR_AR       # LR statistic: CAR(1) versus independent residuals
#                         (phi = 0 lies on the boundary: conservative test)

# Estimates and covariance matrix, labelled as in ENET_P6.R
est  <- car1_fit_P6$coef
names(est) <- c(paste0("beta_PCL5_", 0:3), paste0("beta_CLL_", 0:3))
VCOV <- car1_fit_P6$vcov
dimnames(VCOV) <- list(names(est), names(est))


## 3. Hypotheses (Step 2; identical to ENET_P6.R) ---------------------------------
H_NoTrend <- "beta_PCL5_1 = 0; beta_PCL5_2 = 0 &
              beta_CLL_1  = 0; beta_CLL_2  = 0"
H_Linear <- "beta_PCL5_1 < 0; beta_PCL5_2 = 0 &
             beta_CLL_1  > 0; beta_CLL_2  = 0"
H_Parabolic <- "beta_PCL5_2 < 0 &
                beta_CLL_1 > 0; beta_CLL_2 = 0"
#
# Start and end of the intervention phase (all sessions, as in ENET_P6.R)
t_start <- min(data_ENET_P6$Time)
t_end   <- max(data_ENET_P6$Time)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0;
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0;
                       beta_CLL_1 > 0; beta_CLL_2 = 0", 2 * t_start, 2 * t_end)


## 4. GORICA evaluation (Step 3) ---------------------------------------------------
set.seed(123)
output_P6_CAR1 <- goric(est, VCOV = VCOV,
                        hypotheses = list(H_NoTrend   = H_NoTrend,
                                          H_Linear    = H_Linear,
                                          H_Parabolic = H_Parabolic))
output_P6_CAR1
round(coef(output_P6_CAR1)[, c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")], 4)

# Alternative set with the stricter rise-then-fall hypothesis
set.seed(123)
output_P6_CAR1_rf <- goric(est, VCOV = VCOV,
                           hypotheses = list(H_NoTrend  = H_NoTrend,
                                             H_Linear   = H_Linear,
                                             H_RiseFall = H_RiseFall))
output_P6_CAR1_rf
