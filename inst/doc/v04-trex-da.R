## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----data---------------------------------------------------------------------
set.seed(321)
n <- 100
p <- 100
rho <- 0.7
X <- matrix(0, n, p)
X[, 1] <- rnorm(n)
for (j in 2:p) X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
true_support <- c(10, 30, 50, 70, 90)
beta <- numeric(p)
beta[true_support] <- 2
y <- as.numeric(X %*% beta + rnorm(n))

## ----compare------------------------------------------------------------------
plain <- TRexSelector$new(X, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                          control = trex_control(K = 10))
plain$select()
da <- TRexDASelector$new(X, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                         da_control = trex_da_control(da_method = "BT"),
                         control = trex_control(K = 10))
da$select()
data.frame(selector = c("T-Rex", "T-Rex+DA"),
           selected = c(length(plain$selected_indices), length(da$selected_indices)),
           FDP = c(trex_fdp(plain$selected_indices, true_support),
                   trex_fdp(da$selected_indices, true_support)),
           TPP = c(trex_tpp(plain$selected_indices, true_support),
                   trex_tpp(da$selected_indices, true_support)))

## ----fields-------------------------------------------------------------------
da$da_method_used
da$rho_thresh
round(da$rho_grid, 3)

## ----ar1----------------------------------------------------------------------
ar1 <- TRexDASelector$new(X, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                          da_control = trex_da_control(da_method = "AR1"),
                          control = trex_control(K = 10))
ar1$select()
ar1$selected_indices
round(ar1$cor_coef, 3)

## ----prior--------------------------------------------------------------------
blocks <- rep(1:10, each = 10)
pg <- TRexDASelector$new(X, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                         da_control = trex_da_control(prior_groups = list(blocks)),
                         control = trex_control(K = 10))
pg$select()
pg$selected_indices

## ----greedy-------------------------------------------------------------------
greedy <- TRexDASelector$new(X, y, tFDR = 0.1, seed = 7, verbose = FALSE,
                             da_control = trex_da_control(da_method = "BT"),
                             control = trex_control(solver = "TOMP", K = 10,
                                                    exch_tie_alpha = 0.25))
greedy$select()
greedy$selected_indices

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

