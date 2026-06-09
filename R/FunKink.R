# Functional kink model:
# Y_i(t) = beta_1(t) X_i + beta_2(t) (X_i - gamma)_+
#          + Z_i' beta_3(t) + e_i(t)
#
# Intercept convention:
# - For simulation code matching the original scripts, pass a column of ones in Z
#   and keep intercept = FALSE.
# - For real data, either include that column in Z yourself or set
#   intercept = TRUE to let the package add it internally.

prep_kink_data <- function(Y, x, Z = NULL, t = NULL,
                           scale_xz = TRUE,
                           center_y = FALSE,
                           standardize_y = FALSE,
                           intercept = FALSE,
                           remove_na = TRUE) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  m <- ncol(Y)
  x <- as.numeric(x)
  Z <- .as_z(Z, n)

  if (length(x) != n || nrow(Z) != n) {
    stop("`Y`, `x`, and `Z` must have the same number of rows.")
  }
  if (is.null(t)) t <- seq(0, 1, length.out = m)
  if (length(t) != m) stop("`t` must have length ncol(Y).")

  keep <- rep(TRUE, n)
  if (remove_na) {
    keep <- stats::complete.cases(data.frame(x = x, Z, Y, check.names = FALSE))
    Y <- Y[keep, , drop = FALSE]
    x <- x[keep]
    Z <- Z[keep, , drop = FALSE]
  }

  Y_raw <- Y
  x_raw <- x
  Z_raw <- Z

  if (scale_xz) {
    x_center <- mean(x_raw)
    x_scale <- stats::sd(x_raw)
    if (!is.finite(x_scale) || x_scale <= 0) x_scale <- 1
    x <- (x_raw - x_center) / x_scale

    if (ncol(Z_raw) > 0) {
      Z_center <- colMeans(Z_raw)
      Z_scale <- apply(Z_raw, 2, stats::sd)
      Z_scale[!is.finite(Z_scale) | Z_scale <= 0] <- 1
      Z <- sweep(Z_raw, 2, Z_center, "-")
      Z <- sweep(Z, 2, Z_scale, "/")
    } else {
      Z_center <- numeric(0)
      Z_scale <- numeric(0)
    }
  } else {
    x_center <- 0
    x_scale <- 1
    Z_center <- rep(0, ncol(Z_raw))
    Z_scale <- rep(1, ncol(Z_raw))
  }

  if (standardize_y) {
    Y_center <- colMeans(Y)
    Y_scale <- apply(Y, 2, stats::sd)
    Y_scale[!is.finite(Y_scale) | Y_scale <= 0] <- 1
    Y <- sweep(Y, 2, Y_center, "-")
    Y <- sweep(Y, 2, Y_scale, "/")
  } else if (center_y) {
    Y_center <- colMeans(Y)
    Y_scale <- rep(1, ncol(Y))
    Y <- sweep(Y, 2, Y_center, "-")
  } else {
    Y_center <- rep(0, ncol(Y))
    Y_scale <- rep(1, ncol(Y))
  }

  out <- list(
    Y = Y,
    x = x,
    Z = Z,
    t = as.numeric(t),
    intercept = isTRUE(intercept),
    keep = keep,
    Y_raw = Y_raw,
    x_raw = x_raw,
    Z_raw = Z_raw,
    scale_xz = isTRUE(scale_xz),
    center_y = isTRUE(center_y),
    standardize_y = isTRUE(standardize_y),
    Y_center = Y_center,
    Y_scale = Y_scale,
    x_center = x_center,
    x_scale = x_scale,
    Z_center = Z_center,
    Z_scale = Z_scale
  )
  class(out) <- "kink_data"
  out
}

raw_gamma <- function(gamma, data) {
  if (!isTRUE(data$scale_xz)) return(gamma)
  gamma * data$x_scale + data$x_center
}

fit_kink <- function(data,
                     gamma_grid = NULL,
                     grid_len = 300,
                     nbasis = 10,
                     rho = 1e-3,
                     degree = 3,
                     diff_order = 2) {
  Y <- as.matrix(data$Y)
  x <- as.numeric(data$x)
  Z <- .as_z(data$Z, length(x))
  m <- ncol(Y)

  if (is.null(gamma_grid)) {
    gamma_grid <- seq(
      stats::quantile(x, 0.05, names = FALSE),
      stats::quantile(x, 0.95, names = FALSE),
      length.out = grid_len
    )
  }

  cache <- .kink_cache(
    x = x,
    Z = Z,
    m = m,
    gamma_grid = gamma_grid,
    nbasis = nbasis,
    rho = rho,
    degree = degree,
    diff_order = diff_order,
    intercept = isTRUE(data$intercept)
  )

  rss <- numeric(length(gamma_grid))
  best <- 1L
  best_fit <- NULL

  for (j in seq_along(gamma_grid)) {
    fit <- .fit_gamma(Y, cache, j)
    rss[j] <- fit$rss
    if (is.null(best_fit) || fit$rss < best_fit$rss) {
      best <- j
      best_fit <- fit
    }
  }

  out <- list(
    model = "right_kink",
    gamma = gamma_grid[best],
    beta = best_fit$beta,
    coef_names = cache$coef_names,
    rss = best_fit$rss,
    resid = best_fit$resid,
    fitted = best_fit$fitted,
    theta = best_fit$theta,
    rss_path = rss,
    gamma_grid = gamma_grid,
    best_idx = best,
    cache = cache,
    settings = list(nbasis = nbasis, rho = rho, degree = degree,
                    diff_order = diff_order)
  )
  class(out) <- "kink_fit"
  out
}

test_kink <- function(data,
                      B = 500,
                      nbasis = 12,
                      grid_len = 100,
                      seed = 123) {
  set.seed(seed)

  Y <- as.matrix(data$Y)
  x <- as.numeric(data$x)
  Z <- .as_z(data$Z, length(x))
  n <- length(x)
  m <- ncol(Y)
  intercept <- isTRUE(data$intercept)

  Bmat <- .basis_matrix(m, nbasis, degree = 3)
  X0 <- if (intercept) cbind("(Intercept)" = 1, X = x, Z) else cbind(X = x, Z)
  p0 <- ncol(X0)
  pp <- p0 * nbasis

  Dbig <- matrix(0, n * m, pp)
  for (v in seq_len(p0)) {
    idx <- ((v - 1) * nbasis + 1):(v * nbasis)
    Dbig[, idx] <- X0[, v] %x% Bmat
  }

  Kn <- crossprod(Dbig) / n
  Kn_inv <- solve(Kn + diag(1e-6, pp))

  Yvec <- as.vector(t(Y))
  delta <- Kn_inv %*% (crossprod(Dbig, Yvec) / n)
  resid <- matrix(Yvec - Dbig %*% delta, nrow = n, ncol = m, byrow = TRUE)

  grid <- seq(
    stats::quantile(x, 0.10, names = FALSE),
    stats::quantile(x, 0.90, names = FALSE),
    length.out = grid_len
  )

  psi1 <- array(0, dim = c(n, m, length(grid)))
  psis <- array(0, dim = c(n, m, length(grid)))
  Vhat <- matrix(0, length(grid), m)

  psi2 <- matrix(0, n, pp)
  for (i in seq_len(n)) {
    Di <- Dbig[((i - 1) * m + 1):(i * m), , drop = FALSE]
    psi2[i, ] <- crossprod(Di, resid[i, ])
  }

  for (gidx in seq_along(grid)) {
    g <- grid[gidx]
    hinge <- pmax(x - g, 0)

    for (j in seq_len(m)) {
      Dt <- Dbig[seq(j, n * m, by = m), , drop = FALSE]
      H <- -crossprod(Dt, hinge) / n
      p1 <- resid[, j] * hinge
      ps <- p1 + psi2 %*% (Kn_inv %*% H)

      psi1[, j, gidx] <- p1
      psis[, j, gidx] <- ps
      Vhat[gidx, j] <- mean(ps^2)
    }
  }

  Vinv <- 1 / pmax(n * m * Vhat, 1e-8)
  Tgrid <- sapply(seq_along(grid), function(gidx) {
    sum(colSums(psi1[, , gidx, drop = FALSE][,,1])^2 * Vinv[gidx, ])
  })
  Tobs <- max(Tgrid)

  boot_stat <- function(w) {
    vals <- sapply(seq_along(grid), function(gidx) {
      sum(colSums(w * psis[, , gidx, drop = FALSE][,,1])^2 * Vinv[gidx, ])
    })
    max(vals)
  }

  Tboot <- replicate(B, boot_stat(stats::rnorm(n)))

  out <- list(
    model = "right_kink_test",
    p_value = mean(Tboot >= Tobs),
    T_obs = Tobs,
    T_boot = Tboot,
    gamma_grid = grid,
    settings = list(B = B, nbasis = nbasis, grid_len = grid_len, seed = seed)
  )
  class(out) <- "kink_test"
  out
}

predict_kink <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted)
  x <- as.numeric(newdata$x)
  Z <- .as_z(newdata$Z, length(x))
  hinge <- pmax(x - object$gamma, 0)
  W <- if (isTRUE(object$cache$intercept)) {
    cbind("(Intercept)" = 1, X = x, hinge = hinge, Z)
  } else {
    cbind(X = x, hinge = hinge, Z)
  }
  W %*% t(object$beta)
}

coef_kink <- function(object, ...) {
  beta <- object$beta
  colnames(beta) <- object$coef_names
  beta
}

print.kink_fit <- function(x, ...) {
  cat("Functional kink fit\n")
  cat("  gamma:", format(x$gamma), "\n")
  cat("  rss:", format(x$rss), "\n")
  invisible(x)
}

print.kink_test <- function(x, ...) {
  cat("Functional kink test\n")
  cat("  T_obs:", format(x$T_obs), "\n")
  cat("  p_value:", format(x$p_value), "\n")
  invisible(x)
}

.as_z <- function(Z, n) {
  if (is.null(Z)) matrix(numeric(0), nrow = n, ncol = 0) else as.matrix(Z)
}

.basis_matrix <- function(m, nbasis, degree = 3) {
  splines::bs(seq(0, 1, length.out = m),
              df = nbasis,
              degree = degree,
              intercept = TRUE)
}

.penalty <- function(K, p, diff_order = 2) {
  D <- diff(diag(K), differences = diff_order)
  kronecker(diag(p), crossprod(D))
}

.kink_cache <- function(x, Z, m, gamma_grid,
                        nbasis, rho, degree, diff_order, intercept) {
  n <- length(x)
  B <- .basis_matrix(m, nbasis, degree)
  K <- ncol(B)
  hinge_name <- "hinge"
  coef_names <- if (isTRUE(intercept)) {
    c("(Intercept)", "X", hinge_name, colnames_or_default(Z))
  } else {
    c("X", hinge_name, colnames_or_default(Z))
  }
  p <- length(coef_names)
  BtB <- crossprod(B)
  P <- .penalty(K, p, diff_order)

  W_list <- vector("list", length(gamma_grid))
  Ainv_list <- vector("list", length(gamma_grid))

  for (j in seq_along(gamma_grid)) {
    g <- gamma_grid[j]
    hinge <- pmax(x - g, 0)
    W <- if (isTRUE(intercept)) {
      cbind("(Intercept)" = 1, X = x, hinge = hinge, Z)
    } else {
      cbind(X = x, hinge = hinge, Z)
    }
    XtX <- kronecker(crossprod(W), BtB)
    A <- XtX + n * m * rho * P + 1e-8 * diag(nrow(XtX))
    W_list[[j]] <- W
    Ainv_list[[j]] <- chol2inv(chol(A))
  }

  list(
    x = x,
    Z = Z,
    m = m,
    n = n,
    B = B,
    K = K,
    p = p,
    intercept = intercept,
    coef_names = coef_names,
    gamma_grid = gamma_grid,
    W_list = W_list,
    Ainv_list = Ainv_list
  )
}

.fit_gamma <- function(Y, cache, idx) {
  W <- cache$W_list[[idx]]
  G <- Y %*% cache$B
  rhs <- as.vector(t(crossprod(W, G)))
  theta <- cache$Ainv_list[[idx]] %*% rhs
  coef <- matrix(theta, nrow = cache$K, ncol = cache$p)
  beta <- cache$B %*% coef
  colnames(beta) <- cache$coef_names
  fitted <- W %*% t(beta)
  resid <- Y - fitted

  list(
    gamma = cache$gamma_grid[idx],
    beta = beta,
    rss = sum(resid^2),
    resid = resid,
    fitted = fitted,
    theta = theta,
    idx = idx
  )
}

colnames_or_default <- function(Z) {
  if (ncol(Z) == 0) return(character(0))
  nm <- colnames(Z)
  if (is.null(nm)) nm <- paste0("Z", seq_len(ncol(Z)))
  nm
}
