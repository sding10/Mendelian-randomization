overall_start <- Sys.time()
cat("Script started at:", format(overall_start), "\n")


## Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(gt)
library(sandwich)
library(gridExtra)
library(grid)
library(readxl)
library(purrr)
library(tibble)

### Five methods

#### Data preparation

expit <- function(x) {
  1 / (1 + base::exp(-x))
}

logit <- function(p) {
  base::log(p / (1 - p))
}

solve_psd <- function(A, tol = 1e-10) {
  eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < tol) {
    A <- A + diag(tol, nrow(A))
  }
  solve(A)
}

trim_prob <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}

as_binary <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  as.numeric(x)
}

wald_glm <- function(fit, term, alpha = 0.05) {
  coef_tab <- summary(fit)$coefficients
  
  if (!(term %in% rownames(coef_tab))) {
    return(list(
      est = NA_real_,
      se = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      z = NA_real_,
      p = NA_real_
    ))
  }
  
  est <- unname(coef_tab[term, "Estimate"])
  se  <- unname(coef_tab[term, "Std. Error"])
  z   <- unname(coef_tab[term, "z value"])
  p   <- unname(coef_tab[term, "Pr(>|z|)"])
  zcrit <- stats::qnorm(1 - alpha / 2)
  
  list(
    est = est,
    se = se,
    lcl = est - zcrit * se,
    ucl = est + zcrit * se,
    z = z,
    p = p
  )
}

build_mr_data <- function(dat, z, x, y) {
  vars <- unique(c(z, x, y))
  dat <- dat[, vars, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  
  if (nrow(dat) == 0) {
    stop("No complete cases remain after filtering.")
  }
  
  for (v in c(x, y, z)) {
    dat[[v]] <- as_binary(dat[[v]])
  }
  
  if (!all(dat[[x]] %in% c(0, 1))) {
    stop(sprintf("%s must be coded 0/1.", x))
  }
  
  if (!all(dat[[y]] %in% c(0, 1))) {
    stop(sprintf("%s must be coded 0/1.", y))
  }
  
  for (v in z) {
    if (!all(dat[[v]] %in% c(0, 1, 2))) {
      stop(sprintf("IV %s must be coded 0/1/2.", v))
    }
  }
  
  dat
}


#### Instrument strength model

fit_first_stage <- function(dat, z, x) {
  fml <- stats::as.formula(
    paste(x, "~", paste(z, collapse = " + "))
  )
  stats::glm(fml, data = dat, family = stats::binomial())
}


#### Wald

est_wald <- function(dat, z, x, y, alpha = 0.05, eps = 1e-8) {
  if (length(z) != 1) {
    stop("est_wald() supports exactly one IV.")
  }
  
  z1 <- z[1]
  
  out_fail <- c(
    method = "Wald",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0
  )
  
  fit_x <- tryCatch(
    stats::glm(
      stats::as.formula(paste(x, "~", z1)),
      data = dat,
      family = stats::binomial()
    ),
    error = function(e) NULL
  )
  
  fit_y <- tryCatch(
    stats::glm(
      stats::as.formula(paste(y, "~", z1)),
      data = dat,
      family = stats::binomial()
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_x) || is.null(fit_y)) return(out_fail)
  
  coef_x <- summary(fit_x)$coefficients
  coef_y <- summary(fit_y)$coefficients
  
  if (!(z1 %in% rownames(coef_x)) || !(z1 %in% rownames(coef_y))) {
    return(out_fail)
  }
  
  a_hat <- unname(coef_x[z1, "Estimate"])
  se_a  <- unname(coef_x[z1, "Std. Error"])
  g_hat <- unname(coef_y[z1, "Estimate"])
  se_g  <- unname(coef_y[z1, "Std. Error"])
  
  if (!is.finite(a_hat) || !is.finite(g_hat) ||
      !is.finite(se_a) || !is.finite(se_g) ||
      abs(a_hat) < eps) {
    return(out_fail)
  }
  
  est <- g_hat / a_hat
  
  scores_x <- tryCatch(estfun(fit_x), error = function(e) NULL)
  scores_y <- tryCatch(estfun(fit_y), error = function(e) NULL)
  
  if (is.null(scores_x) || is.null(scores_y)) return(out_fail)
  if (nrow(scores_x) != nrow(scores_y)) return(out_fail)
  
  V_x <- tryCatch(sandwich(fit_x), error = function(e) NULL)
  V_y <- tryCatch(sandwich(fit_y), error = function(e) NULL)
  
  if (is.null(V_x) || is.null(V_y)) return(out_fail)
  
  if (!(z1 %in% rownames(V_x)) || !(z1 %in% rownames(V_y))) {
    return(out_fail)
  }
  
  var_a <- unname(V_x[z1, z1])
  var_g <- unname(V_y[z1, z1])
  
  cross_cov <- tryCatch(
    stats::vcov(fit_x) %*% base::crossprod(scores_x, scores_y) %*% stats::vcov(fit_y),
    error = function(e) NULL
  )
  
  if (is.null(cross_cov) ||
      !(z1 %in% rownames(cross_cov)) ||
      !(z1 %in% colnames(cross_cov))) {
    return(out_fail)
  }
  
  cov_ag <- unname(cross_cov[z1, z1])
  
  var <- var_g / (a_hat^2) +
    (g_hat^2 * var_a) / (a_hat^4) -
    (2 * g_hat * cov_ag) / (a_hat^3)
  
  if (!is.finite(var) || var <= 0) {
    return(out_fail)
  }
  
  se <- sqrt(var)
  
  if (!is.finite(est) || !is.finite(se) || se <= eps) {
    return(out_fail)
  }
  
  z_stat <- est / se
  p_val  <- 2 * (1 - stats::pnorm(abs(z_stat)))
  zcrit  <- stats::qnorm(1 - alpha / 2)
  lcl    <- est - zcrit * se
  ucl    <- est + zcrit * se
  
  c(
    method = "Wald",
    est = est,
    se = se,
    lcl = lcl,
    ucl = ucl,
    z = z_stat,
    p = p_val,
    conv = 1
  )
}


#### 2SPS
est_2sps <- function(dat, z, x, y, eps = 1e-8) {
  out_fail <- c(
    method = "2SPS",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0
  )
  
  fml_1 <- stats::as.formula(
    paste(x, "~", paste(z, collapse = " + "))
  )
  
  fit_1 <- tryCatch(
    stats::glm(fml_1, data = dat, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(fit_1)) return(out_fail)
  
  dat$x_hat <- tryCatch(
    stats::predict(fit_1, type = "response"),
    error = function(e) rep(NA_real_, nrow(dat))
  )
  
  if (any(!is.finite(dat$x_hat))) return(out_fail)
  
  fit_2 <- tryCatch(
    stats::glm(
      stats::as.formula(paste(y, "~ x_hat")),
      data = dat,
      family = stats::binomial()
    ),
    error = function(e) NULL
  )
  if (is.null(fit_2)) return(out_fail)
  
  out <- wald_glm(fit_2, "x_hat")
  
  est <- out$est
  se  <- out$se
  
  if (!is.finite(est) || !is.finite(se) || se <= eps) {
    return(out_fail)
  }
  
  c(method = "2SPS", out, conv = 1)
}


#### 2SRI

est_2sri <- function(dat, z, x, y, eps = 1e-8) {
  out_fail <- c(
    method = "2SRI",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0
  )
  
  fml_1 <- stats::as.formula(
    paste(x, "~", paste(z, collapse = " + "))
  )
  
  fit_1 <- tryCatch(
    stats::glm(fml_1, data = dat, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(fit_1)) return(out_fail)
  
  dat$x_hat <- tryCatch(
    stats::predict(fit_1, type = "response"),
    error = function(e) rep(NA_real_, nrow(dat))
  )
  if (any(!is.finite(dat$x_hat))) return(out_fail)
  
  dat$r_1 <- dat[[x]] - dat$x_hat
  
  fit_2 <- tryCatch(
    stats::glm(
      stats::as.formula(paste(y, "~", x, "+ r_1")),
      data = dat,
      family = stats::binomial()
    ),
    error = function(e) NULL
  )
  if (is.null(fit_2)) return(out_fail)
  
  out <- wald_glm(fit_2, x)
  
  est <- out$est
  se  <- out$se
  
  if (!is.finite(est) || !is.finite(se) || se <= eps) {
    return(out_fail)
  }
  
  c(method = "2SRI", out, conv = 1)
}


#### GMM (Corrected)

est_gmm <- function(dat, z, x, y,
                    start = c(0, 0),
                    alpha = 0.05,
                    maxit = 1000,
                    eps = 1e-8) {
  y_vec <- dat[[y]]
  X <- cbind(1, dat[[x]])
  Z <- cbind(1, as.matrix(dat[, z, drop = FALSE]))
  n <- nrow(dat)
  
  if (length(y_vec) != n || nrow(X) != n || nrow(Z) != n) {
    stop("Input dimensions do not match.")
  }
  
  p <- ncol(X)
  
  if (length(start) != p) {
    stop("Length of 'start' must equal the number of regression parameters.")
  }
  
  out_fail <- c(
    method = "GMM",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0,
    obj = NA_real_
  )
  
  ZZ <- crossprod(Z) / n
  W <- tryCatch(solve_psd(ZZ), error = function(e) NULL)
  
  if (is.null(W) || any(!is.finite(W))) {
    return(out_fail)
  }
  
  g_bar <- function(b) {
    eta <- drop(X %*% b)
    mu  <- expit(eta)
    r   <- y_vec - mu
    drop(crossprod(Z, r) / n)
  }
  
  obj_gmm <- function(b) {
    g <- g_bar(b)
    
    if (any(!is.finite(g))) {
      return(1e12)
    }
    
    val <- as.numeric(t(g) %*% W %*% g)
    if (!is.finite(val)) 1e12 else val
  }
  
  opt <- tryCatch(
    stats::optim(
      par = start,
      fn = obj_gmm,
      method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-10)
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || is.null(opt$par) || any(!is.finite(opt$par))) {
    return(out_fail)
  }
  
  b_hat <- opt$par
  
  if (opt$convergence != 0) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  eta_hat <- drop(X %*% b_hat)
  mu_hat  <- expit(eta_hat)
  r_hat   <- y_vec - mu_hat
  
  if (any(!is.finite(mu_hat)) || any(mu_hat <= 0) || any(mu_hat >= 1)) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  G <- Z * as.numeric(r_hat)
  S_hat <- crossprod(G) / n
  
  w_hat <- as.numeric(mu_hat * (1 - mu_hat))
  D_hat <- -crossprod(Z, X * w_hat) / n
  
  B <- t(D_hat) %*% W %*% D_hat
  B_inv <- tryCatch(solve_psd(B), error = function(e) NULL)
  
  if (is.null(B_inv) || any(!is.finite(B_inv))) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  M <- t(D_hat) %*% W %*% S_hat %*% W %*% D_hat
  
  if (any(!is.finite(M))) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  vcov_b <- (B_inv %*% M %*% B_inv) / n
  
  if (any(!is.finite(vcov_b))) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  var_b <- diag(vcov_b)
  var_b[var_b < 0] <- NA_real_
  se_b <- sqrt(var_b)
  
  est <- b_hat[2]
  se  <- se_b[2]
  
  if (!is.finite(est) ||
      !is.finite(se) || is.na(se) || se <= eps) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  z_stat <- est / se
  p_val  <- 2 * (1 - stats::pnorm(abs(z_stat)))
  zcrit  <- stats::qnorm(1 - alpha / 2)
  lcl    <- est - zcrit * se
  ucl    <- est + zcrit * se
  
  c(
    method = "GMM",
    est = est,
    se = se,
    lcl = lcl,
    ucl = ucl,
    z = z_stat,
    p = p_val,
    conv = 1,
    obj = opt$value
  )
}


#### IV-MVB (Corrected)
est_ivmvb <- function(dat, z, x, y,
                      start = NULL,
                      alpha = 0.05,
                      maxit = 1000,
                      eps = 1e-8) {
  x_obs <- dat[[x]]
  y_obs <- dat[[y]]
  Z <- as.matrix(dat[, z, drop = FALSE])
  k <- ncol(Z)
  
  rho_obs <- suppressWarnings(stats::cor(x_obs, y_obs, use = "complete.obs"))
  
  if (!is.finite(rho_obs)) {
    rho_obs <- NA_real_
  }
  
  if (is.finite(rho_obs)) {
    rho_obs <- max(min(rho_obs, 1 - eps), -1 + eps)
  }
  
  if (is.null(start)) {
    start <- rep(0, k + 3)
  }
  
  out_fail <- c(
    method = "IV-MVB",
    est = NA_real_,
    se = NA_real_,
    lcl = NA_real_,
    ucl = NA_real_,
    z = NA_real_,
    p = NA_real_,
    conv = 0,
    obj = NA_real_,
    rho_obs = rho_obs
  )
  
  if (!is.finite(rho_obs)) {
    return(out_fail)
  }
  
  nll <- function(par) {
    a0 <- par[1]
    a  <- par[2:(k + 1)]
    b0 <- par[k + 2]
    b1 <- par[k + 3]
    
    mu_x <- as.numeric(a0 + Z %*% a)
    pi_x <- expit(mu_x)
    
    mu_y <- b0 + b1 * pi_x
    pi_y <- expit(mu_y)
    
    pi_x <- trim_prob(pi_x, eps)
    pi_y <- trim_prob(pi_y, eps)
    
    denom <- sqrt(pi_x * (1 - pi_x) * pi_y * (1 - pi_y))
    theta <- rho_obs / denom
    
    p11 <- pi_x * pi_y * (1 + theta * (1 - pi_x) * (1 - pi_y))
    p10 <- pi_x * (1 - pi_y) * (1 - theta * (1 - pi_x) * pi_y)
    p01 <- (1 - pi_x) * pi_y * (1 - theta * pi_x * (1 - pi_y))
    p00 <- (1 - pi_x) * (1 - pi_y) * (1 + theta * pi_x * pi_y)
    
    if (any(!is.finite(p11)) || any(!is.finite(p10)) ||
        any(!is.finite(p01)) || any(!is.finite(p00))) {
      return(1e12)
    }
    
    if (any(p11 <= 0) || any(p10 <= 0) || any(p01 <= 0) || any(p00 <= 0)) {
      return(1e12)
    }
    
    p_sum <- p11 + p10 + p01 + p00
    if (any(abs(p_sum - 1) > 1e-6)) {
      return(1e12)
    }
    
    ll <- x_obs * y_obs * log(p11) +
      x_obs * (1 - y_obs) * log(p10) +
      (1 - x_obs) * y_obs * log(p01) +
      (1 - x_obs) * (1 - y_obs) * log(p00)
    
    -sum(ll)
  }
  
  opt <- tryCatch(
    stats::optim(
      par = start,
      fn = nll,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = maxit, reltol = 1e-10)
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || is.null(opt$par)) {
    return(out_fail)
  }
  
  conv <- as.integer(opt$convergence == 0)
  b_hat <- opt$par
  
  if (conv == 0 || any(!is.finite(b_hat))) {
    out <- out_fail
    out["obj"] <- if (!is.null(opt$value)) opt$value else NA_real_
    return(out)
  }
  
  se_b <- rep(NA_real_, length(b_hat))
  if (!is.null(opt$hessian)) {
    V <- tryCatch(solve_psd(opt$hessian), error = function(e) NULL)
    if (!is.null(V)) {
      var_b <- diag(V)
      var_b[var_b < 0] <- NA_real_
      se_b <- sqrt(var_b)
    }
  }
  
  est <- b_hat[k + 3]
  se  <- se_b[k + 3]
  
  if (!is.finite(est) || !is.finite(se) || is.na(se) || se <= eps) {
    out <- out_fail
    out["obj"] <- opt$value
    return(out)
  }
  
  z_stat <- est / se
  p_val  <- 2 * (1 - stats::pnorm(abs(z_stat)))
  zcrit  <- stats::qnorm(1 - alpha / 2)
  lcl    <- est - zcrit * se
  ucl    <- est + zcrit * se
  
  c(
    method = "IV-MVB",
    est = est,
    se = se,
    lcl = lcl,
    ucl = ucl,
    z = z_stat,
    p = p_val,
    conv = 1,
    obj = opt$value,
    rho_obs = rho_obs
  )
}

### Wrapper

MRbinary <- function(dat,
                     z,
                     x,
                     y,
                     alpha = 0.05,
                     maxit = 1000) {
  
  dat <- build_mr_data(dat, z = z, x = x, y = y)
  
  safe_est <- function(expr, method) {
    tryCatch(
      expr,
      error = function(e) c(
        method = method,
        est = NA_real_,
        se = NA_real_,
        lcl = NA_real_,
        ucl = NA_real_,
        z = NA_real_,
        p = NA_real_,
        conv = 0
      )
    )
  }
  
  res <- list()
  
  res$stage1 <- tryCatch(
    fit_first_stage(dat, z, x),
    error = function(e) NULL
  )
  
  res$Wald     <- safe_est(est_wald(dat, z, x, y, alpha = alpha), "Wald")
  res$`2SPS`   <- safe_est(est_2sps(dat, z, x, y), "2SPS")
  res$`2SRI`   <- safe_est(est_2sri(dat, z, x, y), "2SRI")
  res$GMM      <- safe_est(est_gmm(dat, z, x, y, alpha = alpha, maxit = maxit), "GMM")
  res$`IV-MVB` <- safe_est(est_ivmvb(dat, z, x, y, alpha = alpha, maxit = maxit), "IV-MVB")
  
  tab <- do.call(
    rbind,
    lapply(res[c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB")], function(o) {
      data.frame(
        method = as.character(o["method"]),
        est = as.numeric(o["est"]),
        se = as.numeric(o["se"]),
        lcl = as.numeric(o["lcl"]),
        ucl = as.numeric(o["ucl"]),
        z = as.numeric(o["z"]),
        p = as.numeric(o["p"]),
        conv = if ("conv" %in% names(o)) as.numeric(unname(o["conv"])) else 1,
        stringsAsFactors = FALSE
      )
    })
  )
  
  zcrit <- stats::qnorm(1 - alpha / 2)
  ok <- tab$conv == 1 & is.finite(tab$est) & is.finite(tab$se) & tab$se > 0
  
  tab$lcl[ok] <- tab$est[ok] - zcrit * tab$se[ok]
  tab$ucl[ok] <- tab$est[ok] + zcrit * tab$se[ok]
  tab$z[ok]   <- tab$est[ok] / tab$se[ok]
  tab$p[ok]   <- 2 * (1 - stats::pnorm(abs(tab$z[ok])))
  
  tab$lcl[!ok] <- NA_real_
  tab$ucl[!ok] <- NA_real_
  tab$z[!ok]   <- NA_real_
  tab$p[!ok]   <- NA_real_
  
  list(
    dat = dat,
    stage1 = res$stage1,
    summary = tab,
    raw = res
  )
}


### Simulation

#### Simulation once
sim_once <- function(n,
                     a1,
                     b1,
                     c_x = 1.5,
                     c_y = 1.5,
                     p_z = 0.2,
                     sigma_u = 0.5,
                     a0 = 0,
                     b0 = 0,
                     maxit = 1000) {
  
  Z <- rbinom(n, size = 1, prob = p_z)
  U <- rnorm(n, mean = 0, sd = sigma_u)
  
  p_x <- expit(a0 + a1 * Z + c_x * U)
  p_x <- trim_prob(p_x)
  X <- rbinom(n, size = 1, prob = p_x)
  
  p_y <- expit(b0 + b1 * X + c_y * U)
  p_y <- trim_prob(p_y)
  Y <- rbinom(n, size = 1, prob = p_y)
  
  dat <- data.frame(
    Z = Z,
    X = X,
    Y = Y,
    U = U
  )
  
  fit <- tryCatch(
    MRbinary(
      dat = dat[, c("Z", "X", "Y")],
      z = "Z",
      x = "X",
      y = "Y",
      maxit = maxit
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || is.null(fit$summary)) {
    return(data.frame(
      method = c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB"),
      est = NA_real_,
      se = NA_real_,
      lcl = NA_real_,
      ucl = NA_real_,
      z = NA_real_,
      p = NA_real_,
      conv = 0,
      b1 = b1,
      stringsAsFactors = FALSE
    ))
  }
  
  tab <- fit$summary
  
  if ("est" %in% names(tab)) {
    tab$est[!is.finite(tab$est)] <- NA_real_
  }
  
  if ("se" %in% names(tab)) {
    tab$se[!is.finite(tab$se)] <- NA_real_
    tab$se[tab$se <= 0] <- NA_real_
  }
  
  if ("lcl" %in% names(tab)) {
    tab$lcl[!is.finite(tab$lcl)] <- NA_real_
  }
  
  if ("ucl" %in% names(tab)) {
    tab$ucl[!is.finite(tab$ucl)] <- NA_real_
  }
  
  if ("z" %in% names(tab)) {
    tab$z[!is.finite(tab$z)] <- NA_real_
  }
  
  if ("p" %in% names(tab)) {
    tab$p[!is.finite(tab$p)] <- NA_real_
  }
  
  needed <- c("method", "est", "se", "lcl", "ucl", "z", "p", "conv")
  for (v in needed) {
    if (!v %in% names(tab)) {
      tab[[v]] <- NA
    }
  }
  
  tab$b1 <- b1
  
  tab
}


#### Simulation setup
run_sim <- function(n_sim,
                    n,
                    a1,
                    b1,
                    c_x = 1.5,
                    c_y = 1.5,
                    p_z = 0.2,
                    sigma_u = 0.5,
                    scen = NA_character_,
                    methods = c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB"),
                    seed = NULL,
                    maxit = 1000) {
  
  if (!is.null(seed)) set.seed(seed)
  
  res_list <- vector("list", n_sim)
  run_start <- Sys.time()
  
  cat("\n--------------------------------------------------\n")
  cat("Starting run_sim()",
      "| scen =", scen,
      "| n =", n,
      "| a1 =", a1,
      "| b1 =", b1,
      "| n_sim =", n_sim, "\n")
  cat("Start time:", format(run_start), "\n")
  
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  on.exit(close(pb), add = TRUE)
  
  
  make_fail <- function(i) {
    data.frame(
      sim    = i,
      scen   = scen,
      n      = n,
      a1     = a1,
      b1     = b1,
      method = methods,
      est    = NA_real_,
      se     = NA_real_,
      lcl    = NA_real_,
      ucl    = NA_real_,
      z      = NA_real_,
      p      = NA_real_,
      conv   = 0,
      stringsAsFactors = FALSE
    )
  }
  
  for (i in seq_len(n_sim)) {
    
    out <- tryCatch(
      sim_once(
        n       = n,
        a1      = a1,
        b1      = b1,
        c_x     = c_x,
        c_y     = c_y,
        p_z     = p_z,
        sigma_u = sigma_u,
        maxit   = maxit
      ),
      error = function(e) NULL
    )
    
    if (is.null(out)) {
      res_list[[i]] <- make_fail(i)
      next
    }
    
    tab <- out
    
    if (!"method" %in% names(tab) && "Method" %in% names(tab)) {
      names(tab)[names(tab) == "Method"] <- "method"
    }
    if (!"est" %in% names(tab) && "Estimate" %in% names(tab)) {
      names(tab)[names(tab) == "Estimate"] <- "est"
    }
    if (!"se" %in% names(tab) && "SE" %in% names(tab)) {
      names(tab)[names(tab) == "SE"] <- "se"
    }
    if (!"lcl" %in% names(tab) && "LCL" %in% names(tab)) {
      names(tab)[names(tab) == "LCL"] <- "lcl"
    }
    if (!"ucl" %in% names(tab) && "UCL" %in% names(tab)) {
      names(tab)[names(tab) == "UCL"] <- "ucl"
    }
    if (!"z" %in% names(tab)) {
      tab$z <- NA_real_
    }
    if (!"p" %in% names(tab)) {
      tab$p <- NA_real_
    }
    if (!"conv" %in% names(tab)) {
      tab$conv <- NA_real_
    }
    
    tab$sim  <- i
    tab$scen <- scen
    tab$n    <- n
    tab$a1   <- a1
    tab$b1   <- b1
    
    miss <- setdiff(methods, tab$method)
    if (length(miss) > 0) {
      tab_miss <- data.frame(
        sim    = i,
        scen   = scen,
        n      = n,
        a1     = a1,
        b1     = b1,
        method = miss,
        est    = NA_real_,
        se     = NA_real_,
        lcl    = NA_real_,
        ucl    = NA_real_,
        z      = NA_real_,
        p      = NA_real_,
        conv   = 0,
        stringsAsFactors = FALSE
      )
      tab <- rbind(tab, tab_miss)
    }
    
    keep <- c("sim", "scen", "n", "a1", "b1",
              "method", "est", "se", "lcl", "ucl", "z", "p", "conv")
    
    for (v in keep) {
      if (!v %in% names(tab)) tab[[v]] <- NA
    }
    
    tab <- tab[, keep, drop = FALSE]
    tab <- tab[match(methods, tab$method), , drop = FALSE]
    
    res_list[[i]] <- tab
  }
  
  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  
  res
}


### Alpha grid
a1_grid <- c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)


### Scenario A
Scenario_A <- function(n_sim = 8000) {
  res_list <- vector("list", length(a1_grid))
  
  for (i in seq_along(a1_grid)) {
    a1 <- a1_grid[i]
    
    res_list[[i]] <- run_sim(
      n_sim = n_sim,
      n = 1000,
      a1 = a1,
      b1 = 0,
      c_x = 1.5,
      c_y = 1.5,
      p_z = 0.2,
      sigma_u = 0.5,
      scen = "A",
      seed = 1000 + i
    )
  }
  
  dplyr::bind_rows(res_list)
}


### Scenario B
Scenario_B <- function(n_sim = 8000) {
  res_list <- vector("list", length(a1_grid))
  
  for (i in seq_along(a1_grid)) {
    a1 <- a1_grid[i]
    
    res_list[[i]] <- run_sim(
      n_sim = n_sim,
      n = 1000,
      a1 = a1,
      b1 = 1,
      c_x = 1.5,
      c_y = 1.5,
      p_z = 0.2,
      sigma_u = 0.5,
      scen = "B",
      seed = 2000 + i
    )
  }
  
  dplyr::bind_rows(res_list)
}


### Combine A+B
sim_res_ab <- dplyr::bind_rows(
  Scenario_A(n_sim = 8000),
  Scenario_B(n_sim = 8000)
)


### Sample size grid
n_grid <- c(50, 200, 500, 1000, 1500, 2000, 2500)


### Scenario C

Scenario_C <- function(n_sim = 8000) {
  res_list <- vector("list", length(n_grid))
  
  for (i in seq_along(n_grid)) {
    n_i <- n_grid[i]
    
    res_list[[i]] <- run_sim(
      n_sim = n_sim,
      n = n_i,
      a1 = 0.7,
      b1 = 0,
      c_x = 1.5,
      c_y = 1.5,
      p_z = 0.2,
      sigma_u = 0.5,
      scen = "C",
      seed = 3000 + i
    )
  }
  
  dplyr::bind_rows(res_list)
}


### Scenario D

Scenario_D <- function(n_sim = 8000) {
  res_list <- vector("list", length(n_grid))
  
  for (i in seq_along(n_grid)) {
    n_i <- n_grid[i]
    
    res_list[[i]] <- run_sim(
      n_sim = n_sim,
      n = n_i,
      a1 = 0.7,
      b1 = 1,
      c_x = 1.5,
      c_y = 1.5,
      p_z = 0.2,
      sigma_u = 0.5,
      scen = "D",
      seed = 4000 + i
    )
  }
  
  dplyr::bind_rows(res_list)
}


### Combine C+D
sim_res_cd <- dplyr::bind_rows(
  Scenario_C(n_sim = 8000),
  Scenario_D(n_sim = 8000)
)


#### Combine A+B+C+D

sim_res_all <- dplyr::bind_rows(sim_res_ab, sim_res_cd) %>%
  dplyr::mutate(
    scen = factor(scen, levels = c("A", "B", "C", "D")),
    method = factor(method, levels = c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB"))
  )



output_dir  <- "results"
figures_dir <- file.path(output_dir, "figures")
tables_dir  <- file.path(output_dir, "tables")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)


#### Figure 2
plot_dat <- sim_res_all %>%
  dplyr::filter(conv == 1, !is.na(est), !is.na(b1)) %>%
  dplyr::mutate(
    bias   = est - b1,
    scen   = as.character(scen),
    method = factor(method, levels = c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB"))
  )

method_levels <- c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB")

fill_cols <- c(
  "Wald"   = "black",
  "2SPS"   = "grey55",
  "2SRI"   = "grey80",
  "GMM"    = "white",
  "IV-MVB" = "cyan"
)

border_col    <- "grey45"
zero_line_col <- "darkseagreen3"
median_col    <- "red"


make_bias_panel <- function(dat, scen_id, xvar, x_levels,
                            xlab_text, panel_letter) {
  
  par(
    mar = c(4.8, 5.0, 2.0, 1.0),
    mgp = c(2.5, 0.7, 0),
    tcl = -0.3,
    bg  = "grey92"
  )
  
  dat_sub <- dat %>%
    dplyr::filter(scen == scen_id) %>%
    dplyr::mutate(
      xgrp = if (xvar == "a1") as.character(a1) else as.character(n),
      xgrp = factor(xgrp, levels = as.character(x_levels)),
      method = factor(method, levels = method_levels)
    )
  
  n_x     <- length(x_levels)
  centers <- seq_len(n_x) * 3
  offsets <- c(-0.8, -0.4, 0, 0.4, 0.8)
  
  box_list <- list()
  at_vec   <- c()
  col_vec  <- c()
  medians  <- c()
  
  k <- 1
  for (i in seq_along(x_levels)) {
    for (j in seq_along(method_levels)) {
      
      vals <- dat_sub %>%
        dplyr::filter(
          xgrp == as.character(x_levels[i]),
          method == method_levels[j]
        ) %>%
        dplyr::pull(bias)
      
      box_list[[k]] <- vals
      at_vec[k]     <- centers[i] + offsets[j]
      col_vec[k]    <- fill_cols[method_levels[j]]
      medians[k]    <- if (length(vals) > 0) median(vals, na.rm = TRUE) else NA_real_
      k <- k + 1
    }
  }
  
  boxplot(
    box_list,
    at      = at_vec,
    col     = col_vec,
    border  = border_col,
    outline = FALSE,
    xaxt    = "n",
    yaxt    = "n",
    xlab    = "",
    ylab    = "",
    ylim    = c(-10, 10),
    xlim    = c(min(centers) - 1.5, max(centers) + 1.5),
    pars    = list(boxwex = 0.22)
  )
  
  abline(h = 0, col = zero_line_col, lwd = 1)
  
  for (i in seq_along(at_vec)) {
    if (is.finite(medians[i])) {
      segments(
        at_vec[i] - 0.10, medians[i],
        at_vec[i] + 0.10, medians[i],
        col = median_col, lwd = 1.1
      )
    }
  }
  
  axis(1, at = centers, labels = x_levels, cex.axis = 0.85)
  axis(2, at = c(-10, -5, 0, 5, 10), las = 1, cex.axis = 0.85)
  
  mtext(xlab_text, side = 1, line = 2.8, font = 2, cex = 0.88)
  mtext("Estimated Bias in Causal Parameter", side = 2, line = 3.5, font = 2, cex = 0.82)
  mtext(panel_letter, side = 3, line = 0.3, adj = -0.12, cex = 1.05)
  
  box()
}


png(
  filename = file.path(figures_dir, "fig_bias_all.png"),
  width  = 1000,
  height = 780,
  res    = 120,
  bg     = "grey92"
)

layout(
  matrix(c(1, 2,
           3, 4,
           5, 5), nrow = 3, byrow = TRUE),
  heights = c(1, 1, 0.16)
)

make_bias_panel(plot_dat, "A", "a1",
                c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 1),
                "Instrument Strength", "(a)")

make_bias_panel(plot_dat, "B", "a1",
                c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 1),
                "Instrument Strength", "(b)")

make_bias_panel(plot_dat, "C", "n",
                c(50, 200, 500, 1000, 1500, 2000, 2500),
                "Sample Size", "(c)")

make_bias_panel(plot_dat, "D", "n",
                c(50, 200, 500, 1000, 1500, 2000, 2500),
                "Sample Size", "(d)")

par(mar = c(0, 0, 0, 0), bg = "grey92")
plot.new()

legend(
  "center",
  legend = method_levels,
  horiz  = TRUE,
  bty    = "n",
  pt.bg  = fill_cols[method_levels],
  pch    = 22,
  col    = border_col
)

dev.off()


#### Figure 3

plot_cov_dat <- sim_res_all %>%
  filter(conv == 1, !is.na(lcl), !is.na(ucl), !is.na(b1)) %>%
  mutate(
    scen   = as.character(scen),
    method = factor(as.character(method), levels = method_levels),
    cover  = as.numeric(lcl <= b1 & ucl >= b1)
  ) %>%
  group_by(scen, a1, n, method) %>%
  summarise(
    cov = mean(cover, na.rm = TRUE),
    .groups = "drop"
  )


legend_labels <- c("Wald", "2SPS", "2SRI", "GMM", "MVB")

line_cols <- c(
  "Wald"   = "black",
  "2SPS"   = "grey50",
  "2SRI"   = "grey65",
  "GMM"    = "grey80",
  "IV-MVB" = "cyan3"
)

pch_vals <- c(
  "Wald"   = 16,
  "2SPS"   = 8,
  "2SRI"   = 18,
  "GMM"    = 15,
  "IV-MVB" = 17
)

lty_vals <- c(
  "Wald"   = 2,
  "2SPS"   = 2,
  "2SRI"   = 2,
  "GMM"    = 2,
  "IV-MVB" = 1
)


make_cov_panel <- function(dat, scen_id, xvar, xlab_text, panel_letter) {
  
  par(
    mar = c(4.8, 5.0, 2.0, 1.0),
    mgp = c(2.5, 0.7, 0),
    tcl = -0.3,
    bg  = "grey92"
  )
  
  dat_sub <- dat %>%
    filter(scen == scen_id) %>%
    mutate(x = if (xvar == "a1") a1 else n)
  
  if (xvar == "a1") {
    xlim_use  <- c(0, 1.02)
    xtick_at  <- c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
    xtick_lab <- c("0.01", "0.1", "0.2", "0.3", "0.5", "0.7", "1.0")
  } else {
    xlim_use  <- c(0, 2600)
    xtick_at  <- c(50, 200, 500, 1000, 1500, 2000, 2500)
    xtick_lab <- c("50", "200", "500", "1000", "1500", "2000", "2500")
  }
  
  plot(
    NA, NA,
    xlim = xlim_use,
    ylim = c(-0.02, 1.05),
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    bty  = "l"
  )
  
  axis(1, at = xtick_at, labels = xtick_lab, cex.axis = 0.85)
  axis(
    2,
    at = seq(0, 1, by = 0.2),
    labels = formatC(seq(0, 1, by = 0.2), format = "f", digits = 1),
    las = 1,
    cex.axis = 0.85
  )
  
  mtext(xlab_text, side = 1, line = 2.8, font = 2, cex = 0.88)
  mtext("Coverage Probability", side = 2, line = 3.5, font = 2, cex = 0.82)
  mtext(panel_letter, side = 3, line = 0.3, adj = -0.12, cex = 1.05)
  
  for (m in method_levels) {
    dd <- dat_sub %>%
      filter(method == m) %>%
      arrange(x)
    
    if (nrow(dd) == 0) next
    
    lines(
      dd$x, dd$cov,
      col = line_cols[m],
      lty = lty_vals[m],
      lwd = 1.4
    )
    
    points(
      dd$x, dd$cov,
      col = line_cols[m],
      pch = pch_vals[m],
      cex = 0.85
    )
  }
  
  box(bty = "l")
}


png(
  filename = file.path(figures_dir, "figure_3_coverage_oldstyle.png"),
  width  = 1000,
  height = 780,
  res    = 120,
  bg     = "grey92"
)

layout(
  matrix(c(1, 2,
           3, 4,
           5, 5), nrow = 3, byrow = TRUE),
  heights = c(1, 1, 0.16)
)

make_cov_panel(
  dat          = plot_cov_dat,
  scen_id      = "A",
  xvar         = "a1",
  xlab_text    = "Instrument Strength",
  panel_letter = "(a)"
)

make_cov_panel(
  dat          = plot_cov_dat,
  scen_id      = "B",
  xvar         = "a1",
  xlab_text    = "Instrument Strength",
  panel_letter = "(b)"
)

make_cov_panel(
  dat          = plot_cov_dat,
  scen_id      = "C",
  xvar         = "n",
  xlab_text    = "Sample Size",
  panel_letter = "(c)"
)

make_cov_panel(
  dat          = plot_cov_dat,
  scen_id      = "D",
  xvar         = "n",
  xlab_text    = "Sample Size",
  panel_letter = "(d)"
)


par(mar = c(0, 0, 0, 0), bg = "grey92")
plot.new()

legend(
  "center",
  legend    = legend_labels,
  horiz     = TRUE,
  bty       = "n",
  lty       = lty_vals[method_levels],
  lwd       = 1.4,
  col       = line_cols[method_levels],
  pch       = pch_vals[method_levels],
  pt.cex    = 0.85,
  cex       = 0.95,
  x.intersp = 0.8,
  y.intersp = 0.8
)

dev.off()

# Supplementary Table 2
scen_levels   <- c("A", "B", "C", "D")

tab_bias_dat <- sim_res_all %>%
  filter(conv == 1, !is.na(est), !is.na(b1)) %>%
  mutate(
    scen = factor(as.character(scen), levels = scen_levels),
    method = factor(as.character(method), levels = method_levels),
    bias = est - b1
  ) %>%
  group_by(scen, a1, n, method) %>%
  summarise(
    n_conv = dplyr::n(),
    med_bias = median(bias, na.rm = TRUE),
    q1_bias = quantile(bias, probs = 0.25, na.rm = TRUE),
    q3_bias = quantile(bias, probs = 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    set_lab = case_when(
      scen %in% c("A", "B") ~ paste0("a1 = ", a1),
      scen %in% c("C", "D") ~ paste0("n = ", n)
    ),
    set_val = case_when(
      scen %in% c("A", "B") ~ a1,
      scen %in% c("C", "D") ~ n
    ),
    bias_iqr = sprintf("%.2f (%.2f, %.2f)", med_bias, q1_bias, q3_bias)
  ) %>%
  select(
    Scenario = scen,
    Setting = set_lab,
    sort_value = set_val,
    Method = method,
    `Number Converged` = n_conv,
    `Median Bias (Q1, Q3)` = bias_iqr
  ) %>%
  arrange(Scenario, sort_value, Method)

tab_bias_save <- tab_bias_dat %>%
  select(-sort_value)

write.csv(
  tab_bias_save,
  file.path(output_dir, "supp_table_s2_bias.csv"),
  row.names = FALSE
)

tab_bias_grob <- gridExtra::tableGrob(
  tab_bias_save,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 10,
    padding = unit(c(4, 4), "pt")
  )
)

ggsave(
  filename = file.path(tables_dir, "supp_table_s2_bias.pdf"),
  plot = tab_bias_grob,
  width = 11,
  height = 14,
  units = "in"
)

tab_cov_dat <- sim_res_all %>%
  filter(conv == 1, !is.na(lcl), !is.na(ucl), !is.na(b1)) %>%
  mutate(
    scen = factor(as.character(scen), levels = scen_levels),
    method = factor(as.character(method), levels = method_levels),
    cover = (lcl <= b1 & ucl >= b1)
  ) %>%
  group_by(scen, a1, n, method) %>%
  summarise(
    n_conv = dplyr::n(),
    cov = mean(cover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    set_lab = case_when(
      scen %in% c("A", "B") ~ paste0("a1 = ", a1),
      scen %in% c("C", "D") ~ paste0("n = ", n)
    ),
    set_val = case_when(
      scen %in% c("A", "B") ~ a1,
      scen %in% c("C", "D") ~ n
    ),
    cov_pct = sprintf("%.1f%%", 100 * cov)
  ) %>%
  select(
    Scenario = scen,
    Setting = set_lab,
    sort_value = set_val,
    Method = method,
    `Number Converged` = n_conv,
    `Coverage Probability` = cov_pct
  ) %>%
  arrange(Scenario, sort_value, Method)

tab_cov_save <- tab_cov_dat %>%
  select(-sort_value)

write.csv(
  tab_cov_save,
  file.path(output_dir, "supp_table_s3_coverage.csv"),
  row.names = FALSE
)

tab_cov_grob <- gridExtra::tableGrob(
  tab_cov_save,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 10,
    padding = unit(c(4, 4), "pt")
  )
)

ggsave(
  filename = file.path(tables_dir, "supp_table_s3_coverage.pdf"),
  plot = tab_cov_grob,
  width = 11,
  height = 14,
  units = "in"
)


###########CLEAR DATASET####################
run_clear_analysis <- FALSE
if (run_clear_analysis) {
  path_clear_main <- file.path("data", "clear_data_final.csv")
  path_clear_iv   <- file.path("data", "clear_IVSNPs_final.csv")
  
  clear_main <- readr::read_csv(path_clear_main, show_col_types = FALSE)
  clear_iv_info <- readr::read_csv(path_clear_iv, show_col_types = FALSE)


y_var <- "rastatus_n"
x_var <- "eversmoker_n"
z_vars <- clear_iv_info$SNPvars

clear_mr_data <- clear_main %>%
  dplyr::select(dplyr::all_of(c(y_var, x_var, z_vars)))

clear_mr_base <- clear_mr_data %>%
  dplyr::filter(
    !is.na(.data[[y_var]]),
    !is.na(.data[[x_var]])
  )


iv_missing_summary <- data.frame(
  SNP = z_vars,
  n_missing = sapply(clear_mr_base[z_vars], function(v) sum(is.na(v))),
  n_nonmissing = sapply(clear_mr_base[z_vars], function(v) sum(!is.na(v)))
) %>%
  dplyr::arrange(dplyr::desc(n_missing))

iv_missing_summary
summary(iv_missing_summary$n_nonmissing)

run_single_iv <- function(dat, snp, x_var, y_var) {
  
  dat_sub <- dat %>%
    dplyr::select(dplyr::all_of(c(y_var, x_var, snp))) %>%
    na.omit()
  
  n_obs <- nrow(dat_sub)
  
  fit <- MRbinary(
    dat = dat_sub,
    z = snp,
    x = x_var,
    y = y_var
  )
  
  res <- fit$summary
  res$SNP <- snp
  res$n_obs <- n_obs
  
  res
}

res_single_iv <- map_dfr(
  z_vars,
  ~run_single_iv(clear_mr_base, .x, x_var, y_var)
)

dim(res_single_iv)
table(res_single_iv$method)
summary(res_single_iv$n_obs)


#### Figure 4

snp_order <- c(
  "rs4822102",
  "rs9627272",
  "rs17194666",
  "rs8093502",
  "rs67777803",
  "rs79211555",
  "rs16969968",
  "rs4310804",
  "rs56902655",
  "rs62007780",
  "rs1811739",
  "rs7141236",
  "rs933006",
  "rs11048674",
  "rs3782464",
  "rs7956611",
  "rs117837188",
  "rs2959084",
  "rs61886926",
  "rs2796793",
  "rs17481720",
  "rs10807839",
  "rs2056475",
  "rs2214107",
  "rs6968380",
  "rs10044618",
  "rs13109980",
  "rs55900829",
  "rs4580516",
  "rs904592",
  "rs9288999",
  "rs1869073",
  "rs1901477",
  "rs62180324",
  "rs6738833",
  "rs7600835",
  "rs11162019",
  "rs12142143"
)

setequal(snp_order, clear_iv_info$SNP)

figure4_data <- res_single_iv %>%
  dplyr::rename(SNP_var = SNP) %>%
  dplyr::left_join(
    clear_iv_info %>%
      dplyr::select(SNPvars, SNP),
    by = c("SNP_var" = "SNPvars")
  ) %>%
  dplyr::mutate(
    SNP_label = factor(SNP, levels = snp_order),
    method = factor(
      method,
      levels = c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB")
    )
  ) %>%
  dplyr::filter(
    !is.na(SNP_label),
    is.finite(est),
    is.finite(lcl),
    is.finite(ucl)
  ) %>%
  dplyr::mutate(
    est_plot = pmax(pmin(est, 10), -10),
    lcl_plot = pmax(pmin(lcl, 10), -10),
    ucl_plot = pmax(pmin(ucl, 10), -10)
  )

method_colors <- c(
  "IV-MVB" = "#19B5FE",
  "GMM"    = "#D9D9D9",
  "2SRI"   = "#BDBDBD",
  "2SPS"   = "#969696",
  "Wald"   = "black"
)

method_shapes <- c(
  "IV-MVB" = 17,
  "GMM"    = 15,
  "2SRI"   = 16,
  "2SPS"   = 18,
  "Wald"   = 16
)

figure4_plot <- ggplot(
  figure4_data,
  aes(x = est_plot, y = method, color = method, shape = method)
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    color = "black"
  ) +
  geom_errorbarh(
    aes(xmin = lcl_plot, xmax = ucl_plot),
    height = 0.14,
    linewidth = 0.4
  ) +
  geom_point(
    size = 1.7,
    stroke = 0.2
  ) +
  facet_wrap(
    ~ SNP_label,
    ncol = 4,
    scales = "fixed"
  ) +
  scale_x_continuous(
    limits = c(-10, 10),
    breaks = c(-10, -6, -2, 2, 6, 10)
  ) +
  scale_y_discrete(
    limits = c("IV-MVB", "GMM", "2SRI", "2SPS", "Wald")
  ) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(
    x = "Causal Parameter Estimate (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    strip.text = element_text(
      face = "bold",
      size = 8
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.title.x = element_text(face = "bold", size = 11),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

figure4_plot

ggsave(
  filename = file.path(figures_dir, "figure4_clear_single_iv.png"),
  plot = figure4_plot,
  width = 11,
  height = 14,
  dpi = 300
)

ggsave(
  filename = file.path(figures_dir, "figure4_clear_single_iv.pdf"),
  plot = figure4_plot,
  width = 11,
  height = 14
)


# #### Supplementary Table 4
snp_order <- c(
  "rs4822102",
  "rs9627272",
  "rs17194666",
  "rs8093502",
  "rs67777803",
  "rs79211555",
  "rs16969968",
  "rs4310804",
  "rs56902655",
  "rs62007780",
  "rs1811739",
  "rs7141236",
  "rs933006",
  "rs11048674",
  "rs3782464",
  "rs7956611",
  "rs117837188",
  "rs2959084",
  "rs61886926",
  "rs2796793",
  "rs17481720",
  "rs10807839",
  "rs2056475",
  "rs2214107",
  "rs6968380",
  "rs10044618",
  "rs13109980",
  "rs55900829",
  "rs4580516",
  "rs904592",
  "rs9288999",
  "rs1869073",
  "rs1901477",
  "rs62180324",
  "rs6738833",
  "rs7600835",
  "rs11162019",
  "rs12142143"
)

table_s4 <- map_dfr(z_vars, function(snp) {
  
  dat_sub <- clear_mr_base %>%
    dplyr::select(all_of(c(snp, x_var, y_var))) %>%
    na.omit()
  
  n_obs <- nrow(dat_sub)
  z <- dat_sub[[snp]]
  x <- dat_sub[[x_var]]
  
  allele_freq <- mean(z, na.rm = TRUE) / 2
  observed_maf <- min(allele_freq, 1 - allele_freq)
  
  fit_x <- glm(
    formula = as.formula(paste(x_var, "~", snp)),
    data = dat_sub,
    family = binomial()
  )
  
  fit_x_coef <- summary(fit_x)$coefficients[2, ]
  
  a1_est <- unname(fit_x_coef["Estimate"])
  a1_se  <- unname(fit_x_coef["Std. Error"])
  a1_p   <- unname(fit_x_coef["Pr(>|z|)"])
  
  tbl <- table(z, x)
  chisq_p <- suppressWarnings(chisq.test(tbl)$p.value)
  
  x_success <- tapply(x, z, sum)
  x_total   <- tapply(x, z, length)
  scores    <- as.numeric(names(x_success))
  
  ca_p <- prop.trend.test(
    x = as.numeric(x_success),
    n = as.numeric(x_total),
    score = scores
  )$p.value
  
  tibble(
    SNP_var = snp,
    N_obs = n_obs,
    Observed_MAF = observed_maf,
    IV_strength = a1_est,
    SE_a1 = a1_se,
    p_value = a1_p,
    chisq_p_value = chisq_p,
    cochran_armitage_p_value = ca_p
  )
})

table_s4_final <- table_s4 %>%
  left_join(
    clear_iv_info %>%
      dplyr::select(SNPvars, SNP, risk_allele, nonrisk_allele, MAF),
    by = c("SNP_var" = "SNPvars")
  ) %>%
  rename(
    GWAS_MAF = MAF,
    `Risk Allele` = risk_allele,
    `Non-Risk Allele` = nonrisk_allele,
    `Observed MAF` = Observed_MAF,
    `N Obs` = N_obs,
    `IV Strength (α̂1)` = IV_strength,
    `SE(α̂1)` = SE_a1,
    `p-value` = p_value,
    `χ² p-value` = chisq_p_value,
    `Cochran-Armitage p-value` = cochran_armitage_p_value
  ) %>%
  mutate(
    SNP = factor(SNP, levels = snp_order)
  ) %>%
  arrange(SNP) %>%
  mutate(
    SNP = as.character(SNP),
    GWAS_MAF = round(GWAS_MAF, 3),
    `Observed MAF` = round(`Observed MAF`, 3),
    `IV Strength (α̂1)` = round(`IV Strength (α̂1)`, 3),
    `SE(α̂1)` = round(`SE(α̂1)`, 3),
    `p-value` = ifelse(
      `p-value` < 1e-4, "<0.0001",
      formatC(`p-value`, digits = 4, format = "f")
    ),
    `χ² p-value` = ifelse(
      `χ² p-value` < 1e-4, "<0.0001",
      formatC(`χ² p-value`, digits = 4, format = "f")
    ),
    `Cochran-Armitage p-value` = ifelse(
      `Cochran-Armitage p-value` < 1e-4, "<0.0001",
      formatC(`Cochran-Armitage p-value`, digits = 4, format = "f")
    )
  ) %>%
  dplyr::select(
    SNP,
    `Risk Allele`,
    `Non-Risk Allele`,
    GWAS_MAF,
    `Observed MAF`,
    `N Obs`,
    `IV Strength (α̂1)`,
    `SE(α̂1)`,
    `p-value`,
    `χ² p-value`,
    `Cochran-Armitage p-value`
  )

table_s4_final

write.csv(
  table_s4_final,
  file.path(tables_dir, "table_s4.csv"),
  row.names = FALSE
)

table_s4_gt <- table_s4_final %>%
  gt() %>%
  tab_header(
    title = "Supplementary Table 4",
    subtitle = "Relevant information and statistics for the 38 SNPs used as instruments in the Mendelian randomization analyses"
  )

gtsave(
  data = table_s4_gt,
  filename = file.path(tables_dir, "table_s4.html")
)



# Supplementary Table 5

method_levels <- c("Wald", "2SPS", "2SRI", "GMM", "IV-MVB")

format_p <- function(x) {
  ifelse(
    is.na(x),
    "N/A",
    ifelse(x < 1e-4, "<0.0001", formatC(x, digits = 4, format = "f"))
  )
}

format_est_ci_partial <- function(est, lcl, ucl) {
  dplyr::case_when(
    is.na(est) ~ "N/A",
    !is.na(est) & !is.na(lcl) & !is.na(ucl) ~ sprintf("%.2f (%.2f, %.2f)", est, lcl, ucl),
    TRUE ~ sprintf("%.2f", est)
  )
}


table_s5_weighted <- res_single_iv %>%
  filter(
    conv == 1,
    is.finite(est),
    is.finite(n_obs),
    n_obs > 0
  ) %>%
  group_by(method) %>%
  summarise(
    k = n(),
    sum_w = sum(n_obs),
    sum_w2 = sum(n_obs^2),
    n_total = sum(n_obs),
    est_weighted = weighted.mean(est, w = n_obs),
    
    var_between = {
      mu <- weighted.mean(est, w = n_obs)
      denom <- sum(n_obs) - sum(n_obs^2) / sum(n_obs)
      
      if (n() > 1 && denom > 0) {
        sum(n_obs * (est - mu)^2) / denom
      } else {
        NA_real_
      }
    },
    
    n_eff = ifelse(sum_w2 > 0, (sum_w^2) / sum_w2, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    se_weighted = ifelse(
      !is.na(var_between) & !is.na(n_eff) & n_eff > 0,
      sqrt(var_between / n_eff),
      NA_real_
    ),
    lcl_weighted = ifelse(
      !is.na(se_weighted),
      est_weighted - 1.96 * se_weighted,
      NA_real_
    ),
    ucl_weighted = ifelse(
      !is.na(se_weighted),
      est_weighted + 1.96 * se_weighted,
      NA_real_
    ),
    z_weighted = ifelse(
      !is.na(se_weighted) & se_weighted > 0,
      est_weighted / se_weighted,
      NA_real_
    ),
    p_weighted = ifelse(
      !is.na(z_weighted),
      2 * (1 - pnorm(abs(z_weighted))),
      NA_real_
    )
  ) %>%
  dplyr::select(
    method,
    k,
    n_total,
    est_weighted,
    se_weighted,
    lcl_weighted,
    ucl_weighted,
    p_weighted
  )

table_s5_weighted


fit_multi_38 <- MRbinary(
  dat = clear_mr_base,
  z = z_vars,
  x = x_var,
  y = y_var
)

table_s5_multi <- fit_multi_38$summary %>%
  transmute(
    method,
    est_multi = est,
    se_multi = se,
    lcl_multi = lcl,
    ucl_multi = ucl,
    p_multi = p,
    conv_multi = conv
  )

table_s5_multi

table_s5_combined <- dplyr::full_join(
  table_s5_weighted,
  table_s5_multi,
  by = "method"
) %>%
  dplyr::mutate(
    method = factor(trimws(as.character(method)), levels = method_levels)
  ) %>%
  dplyr::arrange(method) %>%
  dplyr::mutate(
    method = as.character(method)
  )

table_s5_combined


table_s5_final <- table_s5_combined %>%
  dplyr::transmute(
    `MR Method` = method,
    `Avg. Causal Parameter Estimate (95% CI)` =
      format_est_ci_partial(est_weighted, lcl_weighted, ucl_weighted),
    `p-value (Weighted)` = format_p(p_weighted),
    `Causal Parameter Estimate (95% CI)` =
      format_est_ci_partial(est_multi, lcl_multi, ucl_multi),
    `p-value (Additive)` = format_p(p_multi)
  )

table_s5_final


csv_out  <- file.path(tables_dir, "table_s5.csv")
html_out <- file.path(tables_dir, "table_s5.html")
pdf_out  <- file.path(tables_dir, "table_s5.pdf")

write.csv(
  table_s5_final,
  csv_out,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


table_s5_gt <- table_s5_final %>%
  gt() %>%
  tab_header(
    title = "Supplementary Table 5",
    subtitle = paste(
      "Results of existing MR analysis methods and the proposed method",
      "using data from the CLEAR study to estimate the causal effect of smoking on RA risk."
    )
  ) %>%
  cols_label(
    `MR Method` = "MR Method",
    `Avg. Causal Parameter Estimate (95% CI)` =
      md("1) Weighted Average of 38 Single IV Models<br>Avg. Causal Parameter Estimate (95% CI)"),
    `p-value (Weighted)` =
      md("1) Weighted Average of 38 Single IV Models<br>p-value"),
    `Causal Parameter Estimate (95% CI)` =
      md("2) Single Additive Model<br>Causal Parameter Estimate (95% CI)"),
    `p-value (Additive)` =
      md("2) Single Additive Model<br>p-value")
  )

gtsave(data = table_s5_gt, filename = html_out)
gtsave(data = table_s5_gt, filename = pdf_out)


#### Supplementary Table 6

snp_24_paper <- c(
  "rs10044618", "rs10807839", "rs11048674", "rs117837188",
  "rs12142143", "rs13109980", "rs16969968", "rs17481720",
  "rs1811739", "rs1869073", "rs6738833", "rs7956611",
  "rs4580516", "rs6968380", "rs933006", "rs9288999",
  "rs3782464", "rs4310804", "rs55900829", "rs56902655",
  "rs7600835", "rs67777803", "rs7141236", "rs79211555"
)

z_vars_24_legacy <- c(
  "rs10044618_C",
  "rs10807839_G",
  "rs11048674_A",
  "rs117837188_A",
  "rs12142143_C",
  "rs13109980_A",
  "rs16969968_A",
  "rs17481720_T",
  "rs1811739_A",
  "rs1869073_T",
  "rs6738833_T",
  "rs7956611_C",
  "rs4580516_A",
  "rs6968380_G",
  "rs933006_G",
  "rs9288999_A",
  "rs3782464_C",
  "rs4310804_G",
  "rs55900829_A",
  "rs56902655_G",
  "rs7600835_A",
  "rs67777803_T",
  "rs7141236_T",
  "rs79211555_T"
)

length(z_vars_24_legacy)
z_vars_24_legacy

missing_legacy_cols <- setdiff(z_vars_24_legacy, names(clear_main))
if (length(missing_legacy_cols) > 0) {
  stop(
    paste(
      "These legacy Supplementary Table 6 SNP columns are missing from clear_main:",
      paste(missing_legacy_cols, collapse = ", ")
    )
  )
}

clear_mr_base_24_legacy <- clear_main %>%
  dplyr::select(dplyr::all_of(c(y_var, x_var, z_vars_24_legacy))) %>%
  dplyr::filter(
    !is.na(.data[[y_var]]),
    !is.na(.data[[x_var]])
  )

fit_multi_24_legacy <- MRbinary(
  dat = clear_mr_base_24_legacy,
  z = z_vars_24_legacy,
  x = x_var,
  y = y_var
)

table_s6_raw <- fit_multi_24_legacy$summary
table_s6_raw

table_s6_final <- table_s6_raw %>%
  dplyr::filter(method %in% c("2SPS", "2SRI", "GMM", "IV-MVB")) %>%
  dplyr::mutate(
    method = trimws(as.character(method)),
    est = ifelse(conv == 0, NA, est),
    lcl = ifelse(conv == 0, NA, lcl),
    ucl = ifelse(conv == 0, NA, ucl),
    p   = ifelse(conv == 0, NA, p)
  ) %>%
  dplyr::mutate(
    method = factor(method, levels = c("2SPS", "2SRI", "GMM", "IV-MVB"))
  ) %>%
  dplyr::arrange(method) %>%
  dplyr::transmute(
    `MR Method` = as.character(method),
    `Causal Parameter Estimate (95% CI)` = format_est_ci_partial(est, lcl, ucl),
    `p-value` = format_p(p)
  )

table_s6_final

write.csv(
  table_s6_final,
  file.path(tables_dir, "table_s6.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

table_s6_gt <- table_s6_final %>%
  gt() %>%
  tab_header(
    title = "Supplementary Table 6",
    subtitle = paste(
      "Multiple IV results of existing MR analysis methods and the proposed method",
      "using data from the CLEAR study to estimate the causal effect of smoking on RA risk."
    )
  ) %>%
  cols_label(
    `MR Method` = "MR Method",
    `Causal Parameter Estimate (95% CI)` = "Causal Parameter Estimate (95% CI)",
    `p-value` = "p-value"
  ) %>%
  tab_source_note(
    source_note = paste(
      "† SNPs include:",
      paste(snp_24_paper, collapse = ", ")
    )
  ) %>%
  tab_source_note(
    source_note = "‡ If the GMM algorithm did not converge in this multiple-IV scenario, results are displayed as N/A."
  )

gtsave(
  data = table_s6_gt,
  filename = file.path(tables_dir, "table_s6.html")
)

gtsave(
  data = table_s6_gt,
  filename = file.path(tables_dir, "table_s6.pdf")
)

}
overall_end <- Sys.time()
cat("Script finished at:", format(overall_end), "\n")
cat("Total elapsed time:",
    round(as.numeric(difftime(overall_end, overall_start, units = "mins")), 2),
    "minutes\n")

