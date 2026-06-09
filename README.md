# FunKink

`FunKink` is a minimal R package for functional kink model analysis in
function-on-scalar regression. The model is designed for structural changes that
remain continuous while allowing a kink in covariate effects:

```tex
Y_i(t) = \beta_1(t)X_i + \beta_2(t)(X_i-\gamma)_+
       + \mathbf Z_i^\top \boldsymbol\beta_3(t) + e_i(t).
```

It contains only estimation and testing code. Plotting, table formatting, and
real-data reporting utilities are intentionally left out.

## Model convention

The kink term is right-side:

```r
pmax(x - gamma, 0)
```

The intercept is optional and user-controlled.

- Simulation convention: pass `Z = cbind(1, z)` and use `intercept = FALSE`.
- Real-data convention: either include a column of ones in `Z`, or set
  `intercept = TRUE` in `prep_kink_data()`.

## Usage

```r
library(FunKink)

set.seed(1)
n <- 100
m <- 15
t <- seq(0, 1, length.out = m)
x <- rnorm(n)
z <- rnorm(n)
Z <- cbind(1, z)

beta1 <- 4 * cos(2 * pi * t)
beta2 <- 4 * cos(2 * pi * t)
beta30 <- sqrt(2 * t)
beta31 <- sin(2 * pi * t)
Y <- x %*% t(beta1) +
  pmax(x - 0, 0) %*% t(beta2) +
  Z[, 1] %*% t(beta30) +
  Z[, 2] %*% t(beta31) +
  matrix(rnorm(n * m), n, m)

dat <- prep_kink_data(Y, x, Z = Z, t = t, intercept = FALSE)

fit <- fit_kink(dat, grid_len = 100)
fit$gamma
raw_gamma(fit$gamma, dat)

test <- test_kink(dat, B = 199)
test$p_value
```

## Install locally

From the directory containing the package folder:

```bash
R CMD build FunKink
R CMD INSTALL FunKink_0.1.0.tar.gz
```
