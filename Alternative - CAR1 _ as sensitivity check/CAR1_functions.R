################################################################################
# Bivariate regression with continuous-time AR(1) residuals, estimated by REML
#
#   y_jt = x_t' beta_j + e_jt,   j = 1, 2 (e.g., PCL5 and CLL)
#
# Residual covariance (separable):
#   Cov(e_js, e_kt) = Sigma[j, k] * phi^|t - s|
#   - Sigma : 2 x 2 contemporaneous residual covariance matrix
#             (variances of both outcomes and their correlation)
#   - phi   : autocorrelation at a lag of one time unit (here: one week);
#             continuous time, so unequally spaced measurements and gaps due to
#             missing data are handled correctly (phi^lag, lag in weeks).
#
# The mean structure (and thus the trend parameters) is the same as in
# lm(cbind(y1, y2) ~ ...), so the same informative hypotheses can be evaluated.
# For phi -> 0 the model reduces to the multivariate lm() model: then the
# estimates and their covariance matrix are identical to those of lm().
################################################################################

fit_car1_biv <- function(formula, data, outcomes = c("PCL5", "CLL"),
                         time = "Time", time_unit = 7,
                         phi_bounds = c(1e-4, 0.99)) {

  # Complete cases on outcomes, predictors and time (as lm() does)
  vars <- unique(c(outcomes, all.vars(formula), time))
  d <- data[stats::complete.cases(data[, vars]), , drop = FALSE]

  X  <- stats::model.matrix(formula, d)              # n x p, same for both outcomes
  Y  <- as.matrix(d[, outcomes])                     # n x 2
  n  <- nrow(X); p <- ncol(X); q <- length(outcomes)
  tt <- d[[time]] / time_unit                        # time in units for phi
  L  <- abs(outer(tt, tt, "-"))                      # lags

  # Stacked form: y = vec(Y) (outcome 1 first), Xs = I_q (x) X
  y  <- c(Y)
  Xs <- kronecker(diag(q), X)

  # Build Sigma (x) R(phi) from unconstrained parameters
  make_V <- function(par) {
    s   <- exp(par[1:2]); rho <- tanh(par[3]); phi <- plogis(par[4])
    Sig <- matrix(c(s[1]^2, rho * s[1] * s[2], rho * s[1] * s[2], s[2]^2), 2)
    R   <- phi^L
    list(V = kronecker(Sig, R), Sig = Sig, phi = phi)
  }

  # Restricted (REML) log-likelihood, with beta profiled out
  reml <- function(par) {
    m <- make_V(par)
    cV <- tryCatch(chol(m$V), error = function(e) NULL)
    if (is.null(cV)) return(-1e10)
    Vi_X <- backsolve(cV, forwardsolve(t(cV), Xs))
    Vi_y <- backsolve(cV, forwardsolve(t(cV), y))
    XtViX <- crossprod(Xs, Vi_X)
    cA <- tryCatch(chol(XtViX), error = function(e) NULL)
    if (is.null(cA)) return(-1e10)
    b  <- backsolve(cA, forwardsolve(t(cA), crossprod(Xs, Vi_y)))
    r  <- y - Xs %*% b
    rVir <- sum(r * backsolve(cV, forwardsolve(t(cV), r)))
    N <- length(y); k <- ncol(Xs)
    -0.5 * ((N - k) * log(2 * pi) + 2 * sum(log(diag(cV))) +
              2 * sum(log(diag(cA))) + rVir)
  }

  # Starting values from the multivariate OLS fit
  ols <- stats::lm.fit(X, Y)
  S0  <- crossprod(ols$residuals) / (n - p)
  start <- c(log(sqrt(diag(S0))), atanh(stats::cov2cor(S0)[1, 2]), qlogis(0.3))
  lo <- c(-Inf, -Inf, -Inf, qlogis(phi_bounds[1]))
  up <- c( Inf,  Inf,  Inf, qlogis(phi_bounds[2]))

  opt <- stats::optim(start, function(p) -reml(p), method = "L-BFGS-B",
                      lower = lo, upper = up)
  # Model without autocorrelation (phi at lower bound) for comparison
  opt0 <- stats::optim(start[1:3], function(p) -reml(c(p, lo[4])), method = "BFGS")

  m  <- make_V(opt$par)
  cV <- chol(m$V)
  Vi_X  <- backsolve(cV, forwardsolve(t(cV), Xs))
  XtViX <- crossprod(Xs, Vi_X)
  vcov_b <- solve(XtViX)
  b <- c(vcov_b %*% crossprod(Xs, backsolve(cV, forwardsolve(t(cV), y))))

  nm <- c(outer(colnames(X), outcomes, function(a, o) paste0(o, ":", a)))
  names(b) <- nm; dimnames(vcov_b) <- list(nm, nm)

  list(coef = b, vcov = vcov_b, phi = m$phi, Sigma = m$Sig,
       rho = stats::cov2cor(m$Sig)[1, 2],
       logLik_REML = -opt$value, logLik_REML_noAR = -opt0$value,
       # REML likelihoods are comparable here (same fixed effects);
       # phi = 0 lies on the boundary, so the LR test is conservative.
       LR_AR = 2 * (opt0$value - opt$value),
       phi_at_bound = abs(opt$par[4] - up[4]) < 1e-6 || abs(opt$par[4] - lo[4]) < 1e-6,
       nobs = n, convergence = opt$convergence)
}
