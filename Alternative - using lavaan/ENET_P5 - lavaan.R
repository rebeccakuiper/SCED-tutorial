################################################################################
# SCED tutorial - GORICA evaluation for participant P5, model fitted with lavaan
# (alternative to the lm() script; same model, same hypotheses)
################################################################################


## 0. Packages ----------------------------------------------------------------
# install.packages(c("restriktor", "lavaan"))
library(restriktor)
library(lavaan)


## 1. Data preparation (identical to the lm() script) ------------------------------
data_ENET_P5 <- read.table("data/ENET05.dat", header = TRUE)
data_ENET_P5[data_ENET_P5 == -999] <- NA

# Keep 'event' numeric (0/1). 
# lavaan converts a factor to 1/2 codes, 
# which shifts the intercepts relative to lm().
#
# Time in weeks (numerical stability), centered at 1/3:
# (a) Express it in weeks. In days, sqTime runs up to ~16,000 and its
#     coefficients are ~1e-4, which makes the covariance matrix nearly singular.
#     Rescaling changes no weights, only the units of the coefficients.
data_ENET_P5$Time   <- data_ENET_P5$Time / 7
# (b) Centre it at one-third of the intervention phase. 
#     Btw Centring does NOT determine where the maximum/minimum of the parabola lies; 
#     it only makes beta_1 the slope at that time point (beta_2 is unaffected).
data_ENET_P5$Time   <- data_ENET_P5$Time - quantile(data_ENET_P5$Time, 1/3)
data_ENET_P5$sqTime <- data_ENET_P5$Time^2


## 2. Fit the unconstrained model (Step 1) --------------------------------------
model <- '
  PCL5 ~ beta_PCL5_1*Time + beta_PCL5_2*sqTime + beta_PCL5_3*event
  CLL  ~ beta_CLL_1*Time  + beta_CLL_2*sqTime  + beta_CLL_3*event
'
# 1 = Time (linear), 2 = sqTime (quadratic), 3 = event.
#
# Residuals of PCL5 and CLL are correlated by default (as in the lm() model).
fit_P5 <- sem(model, data = data_ENET_P5, meanstructure = TRUE)


## 3. Hypotheses (Step 2) -------------------------------------------------------
H_NoTrend    <- "beta_PCL5_1 = 0; beta_PCL5_2 = 0; beta_CLL_1 = 0; beta_CLL_2 = 0"
H_Linear     <- "beta_PCL5_1 < 0; beta_PCL5_2 = 0; beta_CLL_1 > 0; beta_CLL_2 = 0"
H_Parabolic  <- "beta_PCL5_2 < 0; beta_CLL_1 > 0; beta_CLL_2 = 0"
# Optional:
t_start <- min(data_ENET_P5$Time)
t_end   <- max(data_ENET_P5$Time)
H_RiseFall <- sprintf("beta_PCL5_1 + %.3f*beta_PCL5_2 > 0;
                       beta_PCL5_1 + %.3f*beta_PCL5_2 < 0;
                       beta_CLL_1 > 0; beta_CLL_2 = 0", 2 * t_start, 2 * t_end)


## 4. GORICA evaluation (Step 3) -------------------------------------------------
# Only call goric(fit_P5, ...) with restriktor version 0.7 or higher.
# Workaround: extract the labelled estimates and their covariance matrix.
betas <- c("beta_PCL5_1", "beta_PCL5_2", "beta_PCL5_3",
           "beta_CLL_1",  "beta_CLL_2",  "beta_CLL_3")
est  <- coef(fit_P5)[betas]
VCOV <- vcov(fit_P5)[betas, betas]
# Note: lavaan (ML) divides by N instead of N - k, 
# so its standard errors are smaller than lm()'s (here variances x 21/25) 
# and the weights differ somewhat from the lm() results reported in the paper.

set.seed(123)
output_P5 <- goric(est, VCOV = VCOV,
                   hypotheses = list(H_NoTrend    = H_NoTrend,
                                     H_Linear     = H_Linear,
                                     H_Parabolic  = H_Parabolic))
output_P5
#round(coef(output_P5)[, c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2")], 4)

#Optional / Alternative:
set.seed(123)
output_P5_rf <- goric(est, VCOV = VCOV,
                      hypotheses = list(H_NoTrend    = H_NoTrend,
                                        H_Linear     = H_Linear,
                                        H_RiseFall   = H_RiseFall))
output_P5_rf


## 5. Extra: benchmarks (parameters set by name) ----------------------------------
est_no  <- replace(est, c("beta_PCL5_1", "beta_PCL5_2", "beta_CLL_1", "beta_CLL_2"), 0)
est_lin <- replace(est, c("beta_PCL5_2", "beta_CLL_2"), 0)
my_pop_est <- rbind(NoTrend = est_no, Linear = est_lin, Observed = est)

set.seed(123)
benchmarks_P5 <- benchmark(output_P5, pop_est = my_pop_est, iter = 100)
benchmarks_P5
