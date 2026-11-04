################################################################################
# Sensitivity analysis: GORICA weights with independent residuals (lm(); as in
# the paper) versus continuous-time AR(1) residuals (ENET_Pk_CAR1.R)
#
# Uses output_Pk (from ENET_Pk.R) and output_Pk_CAR1 / car1_fit_Pk (from
# ENET_Pk_CAR1.R); scripts whose results are not in the workspace are run.
# Works with restriktor 0.6-x (CRAN) and 0.7+ (columns selected by name).
################################################################################

P_names <- c("P1", "P3", "P4", "P5", "P6", "P7", "P9", "P10", "P11")

for (p in P_names) {
  if (!exists(paste0("output_", p)))           source(sprintf("ENET_%s.R", p))
  if (!exists(paste0("output_", p, "_CAR1")))  source(sprintf("ENET_%s_CAR1.R", p))
}

weights_of <- function(out) {
  res <- as.data.frame(out$result)
  setNames(res[["gorica.weights"]], res$model)
}


## Per participant: autocorrelation and preferred hypothesis -----------------------
overview <- do.call(rbind, lapply(P_names, function(p) {
  w_lm  <- weights_of(get(paste0("output_", p)))
  w_ar  <- weights_of(get(paste0("output_", p, "_CAR1")))
  fit   <- get(paste0("car1_fit_", p))
  data.frame(Participant    = p,
             n              = fit$nobs,
             phi_1week      = round(fit$phi, 3),
             LR_AR          = round(fit$LR_AR, 2),
             Preferred_lm   = names(which.max(w_lm)),
             Weight_lm      = round(max(w_lm), 3),
             Preferred_CAR1 = names(which.max(w_ar)),
             Weight_CAR1    = round(max(w_ar), 3))
}))
overview
# phi_1week: residual autocorrelation at a lag of one week;
# LR_AR    : likelihood-ratio statistic CAR(1) vs independent residuals
#            (boundary: compare with a 50:50 mixture of chi2(0) and chi2(1),
#             i.e. critical value 2.71 at alpha = .05).


## All weights side by side -------------------------------------------------------
all_weights <- do.call(rbind, lapply(P_names, function(p) {
  w_lm <- weights_of(get(paste0("output_", p)))
  w_ar <- weights_of(get(paste0("output_", p, "_CAR1")))
  data.frame(Participant = p, Hypothesis = names(w_lm),
             Weight_lm = round(w_lm, 3), Weight_CAR1 = round(w_ar[names(w_lm)], 3),
             row.names = NULL)
}))
all_weights

write.table(all_weights, "Results_lm_vs_CAR1.xls", row.names = FALSE, sep = "\t", dec = ",")
