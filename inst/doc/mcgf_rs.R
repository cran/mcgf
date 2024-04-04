## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----regime-------------------------------------------------------------------
library(mcgf)
K <- 2
N <- 5000
lag <- 5

set.seed(123)
tran_mat <- matrix(rnorm(K^2, mean = 0.06 / (K - 1), sd = 0.01), nrow = K)
diag(tran_mat) <- rnorm(K, mean = 0.94, sd = 0.1)
tran_mat <- sweep(abs(tran_mat), 1, rowSums(tran_mat), `/`)
tran_mat

regime <- rep(NA, N)
regime[1] <- 1

for (n in 2:(N)) {
    regime[n] <- sample(1:K, 1, prob = tran_mat[regime[n - 1], ])
}
table(regime)

## ----locations----------------------------------------------------------------
set.seed(123)
x <- stats::rnorm(25, -110)
y <- stats::rnorm(25, 50)
locations_all <- cbind(lon = x, lat = y)
locations <- locations_all[1:20, ]
locations_new <- locations_all[-c(1:20), ]
locations_all

## ----data---------------------------------------------------------------------
# simulate RS MCGF
par_base <- list(
    par_s = list(nugget = 0, c = 0.005, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.2)
)
par_lagr1 <- list(v1 = 100, v2 = 100, k = 2)
par_lagr2 <- list(v1 = 50, v2 = 100, k = 2)

h <- find_dists_new(locations, locations_new)

set.seed(123)
data_all <- mcgf_rs_sim(
    N = N, label = regime,
    base_ls = list("sep"),
    lagrangian_ls = list("lagr_tri"),
    par_base_ls = list(par_base),
    par_lagr_ls = list(par_lagr1, par_lagr2),
    lambda_ls = list(0.3, 0.3),
    lag_ls = list(lag, lag),
    dists_ls = list(h, h)
)
data_all <- data_all[-c(1:(lag + 1)), ]
rownames(data_all) <- 1:nrow(data_all)
head(data_all)

## ----hold---------------------------------------------------------------------
data_old <- data_all[, 1:21]
data_new <- data_all[, -c(1:21)]

## ----mcgf---------------------------------------------------------------------
data_mcgf <- mcgf_rs(data_old[, -1],
    locations = locations, longlat = TRUE,
    label = data_old[, 1]
)
data_mcgf <- add_acfs(data_mcgf, lag_max = lag)
data_mcgf <- add_ccfs(data_mcgf, lag_max = lag)

# If multiple cores are available
# data_mcgf <- add_ccfs(data_mcgf, lag_max = lag, ncores = 8)

## ----acfs---------------------------------------------------------------------
acfs(data_mcgf)

## ----ccfs---------------------------------------------------------------------
# Regime 1
ccfs(data_mcgf)$ccfs_rs$`Regime 1`[1:6, 1:6, 2]

# Regime 2
ccfs(data_mcgf)$ccfs_rs$`Regime 2`[1:6, 1:6, 2]

## ----spatial, message=F-------------------------------------------------------
fit_spatial <- fit_base(
    data_mcgf,
    model_ls = "spatial",
    lag_ls = lag,
    par_init_ls = list(c(c = 0.000001, gamma = 0.5)),
    par_fixed_ls = list(c(nugget = 0)),
    rs = FALSE
)
fit_spatial[[1]]$fit

## ----temporal-----------------------------------------------------------------
fit_temporal <- fit_base(
    data_mcgf,
    model = "temporal",
    lag_ls = lag,
    par_init_ls = list(c(a = 1, alpha = 0.5)),
    rs = FALSE
)
fit_temporal[[1]]$fit

## ----sep----------------------------------------------------------------------
data_mcgf <- add_base(data_mcgf,
    fit_s = fit_spatial,
    fit_t = fit_temporal,
    sep = T
)

## -----------------------------------------------------------------------------
fit_sep <- fit_base(
    data_mcgf,
    model_ls = "sep",
    lag_ls = lag,
    par_init_ls = list(c(c = 0.000001, gamma = 0.5, a = 1, alpha = 0.5)),
    par_fixed_ls = list(c(nugget = 0)),
    control = list(list(iter.max = 10000, eval.max = 10000)),
    rs = FALSE
)
fit_sep[[1]]$fit

## ----base---------------------------------------------------------------------
data_mcgf <- add_base(data_mcgf, fit_base = fit_sep, old = TRUE)
model(data_mcgf, model = "base", old = TRUE)

## ----westerly-----------------------------------------------------------------
fit_lagr_rs <- fit_lagr(data_mcgf,
    model_ls = list("lagr_tri"),
    par_init_ls = list(list(lambda = 0.1, v1 = 100, v2 = 100, k = 2))
)
lapply(fit_lagr_rs[1:2], function(x) x$fit)

## ----stat---------------------------------------------------------------------
data_mcgf <- add_lagr(data_mcgf, fit_lagr = fit_lagr_rs)
model(data_mcgf, old = TRUE)

## ----krige--------------------------------------------------------------------
krige_base_new <- krige_new(
    x = data_mcgf,
    locations_new = locations_new,
    model = "base",
    interval = TRUE
)
krige_stat_new <- krige_new(
    x = data_mcgf,
    locations_new = locations_new,
    model = "all",
    interval = TRUE
)

## ----rmse---------------------------------------------------------------------
# RMSE
rmse_base <-
    sqrt(colMeans((data_new - krige_base_new$fit[, -c(1:20)])^2, na.rm = T))
rmse_stat <-
    sqrt(colMeans((data_new - krige_stat_new$fit[, -c(1:20)])^2, na.rm = T))
rmse <- c(
    "Base" = mean(rmse_base),
    "STAT" = mean(rmse_stat)
)
rmse

## ----MAE----------------------------------------------------------------------
mae_base <- colMeans(abs(data_new - krige_base_new$fit[, -c(1:20)]), na.rm = T)
mae_stat <- colMeans(abs(data_new - krige_stat_new$fit[, -c(1:20)]), na.rm = T)
mae <- c(
    "Base" = mean(mae_base),
    "STAT" = mean(mae_stat)
)
mae

## ----popi---------------------------------------------------------------------
# POPI
popi_base <- colMeans(
    data_new < krige_base_new$lower[, -c(1:20)] | data_new > krige_base_new$upper[, -c(1:20)],
    na.rm = T
)
popi_stat <- colMeans(
    data_new < krige_stat_new$lower[, -c(1:20)] | data_new > krige_stat_new$upper[, -c(1:20)],
    na.rm = T
)
popi <- c(
    "Base" = mean(popi_base),
    "STAT" = mean(popi_stat)
)
popi

