test_that("right-side kink fit and test run", {
  set.seed(1)
  n <- 30
  m <- 8
  tt <- seq(0, 1, length.out = m)
  x <- rnorm(n)
  z <- rnorm(n)
  Z <- cbind(1, z)

  beta1 <- cos(2 * pi * tt)
  beta2 <- sin(2 * pi * tt)
  beta30 <- sqrt(tt + 0.01)
  beta31 <- tt

  Y <- x %*% t(beta1) +
    pmax(x, 0) %*% t(beta2) +
    Z[, 1] %*% t(beta30) +
    Z[, 2] %*% t(beta31) +
    matrix(rnorm(n * m, sd = 0.1), n, m)

  dat <- prep_kink_data(Y, x, Z = Z, t = tt, scale_xz = FALSE)
  fit <- fit_kink(dat, grid_len = 15, nbasis = 5)
  expect_s3_class(fit, "kink_fit")
  expect_true(is.finite(fit$gamma))
  expect_equal(nrow(predict_kink(fit)), n)
  expect_true("hinge" %in% colnames(coef_kink(fit)))

  tst <- test_kink(dat, B = 5, nbasis = 5, grid_len = 10, seed = 1)
  expect_s3_class(tst, "kink_test")
  expect_true(is.finite(tst$p_value))
})
