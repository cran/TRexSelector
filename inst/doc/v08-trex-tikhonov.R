## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----data---------------------------------------------------------------------
set.seed(2718)
n <- 120
p <- 60
true_support <- 21:30
beta <- numeric(p)
beta[true_support] <- 1
X <- matrix(rnorm(n * p), n, p)
y <- as.numeric(X %*% beta + rnorm(n))

## ----K------------------------------------------------------------------------
gamma <- diff(diag(p))                # first differences, (p - 1) x p
if (requireNamespace("Matrix", quietly = TRUE)) {
  K <- trex_gamma_to_k(gamma)
  class(K)
} else {
  K <- crossprod(gamma)
}
K[1:4, 1:4]

## ----run----------------------------------------------------------------------
tik <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, seed = 4, verbose = FALSE,
                                tik_control = trex_tikhonov_control(tikhonov_K = K, lambda_2 = 1),
                                control = trex_control(K = 10))
tik$select()
tik$selected_indices
trex_fdp(tik$selected_indices, true_support)
trex_tpp(tik$selected_indices, true_support)
tik$lambda2_used

## ----pilot--------------------------------------------------------------------
pilot <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, seed = 4, verbose = FALSE,
                                  tik_control = trex_tikhonov_control(
                                    tikhonov_K = K,
                                    pilot_grid = c(0, 1, 10), pilot_K = 4, pilot_T = 2),
                                  control = trex_control(K = 10))
pilot$select()
pilot$lambda2_used
pilot$lambda2_pilot$grid
pilot$lambda2_pilot$counts

## ----power-rule---------------------------------------------------------------
profile <- trex_tikhonov_lambda2_power_rule(X, y, tFDR = 0.1, seed = 4,
                                            tik_control = trex_tikhonov_control(tikhonov_K = K),
                                            control = trex_control(K = 10), grid = c(0, 1, 10))
profile$lambda2
profile$counts

## ----anchor-------------------------------------------------------------------
en <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, seed = 4, verbose = FALSE,
                               tik_control = trex_tikhonov_control(
                                 tikhonov_K = diag(p), lambda_2 = 1,
                                 cluster_dummies = FALSE, fold_dummy_coupling = FALSE),
                               control = trex_control(K = 10))
en$select()
en$selected_indices

## ----centering----------------------------------------------------------------
groups <- rep(1:6, each = 10)
gc <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, seed = 4, verbose = FALSE,
                               tik_control = trex_tikhonov_control(
                                 penalty_structure = "GROUP_CENTERING", prior_groups = groups,
                                 lambda_2 = 1),
                               control = trex_control(K = 10))
gc$select()
gc$selected_indices

## ----tiknet-------------------------------------------------------------------
pw <- TRexTikhonovSelector$new(X, y, tFDR = 0.1, seed = 4, verbose = FALSE,
                               tik_control = trex_tikhonov_control(tikhonov_K = K, lambda_2 = 1),
                               control = trex_control(solver = "TIKNET", K = 10))
pw$select()
pw$selected_indices

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

