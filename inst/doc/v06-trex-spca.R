## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----data---------------------------------------------------------------------
set.seed(77)
n <- 200
p <- 60
Z <- matrix(rnorm(n * 2), n, 2)
V <- matrix(0, p, 2)
V[1:10, 1] <- 1 / sqrt(10)
V[31:40, 2] <- 1 / sqrt(10)
X <- 6 * Z[, 1] %*% t(V[, 1]) + 3 * Z[, 2] %*% t(V[, 2]) + matrix(rnorm(n * p), n, p)

## ----run----------------------------------------------------------------------
spca <- TRexSPCASelector$new(X, M = 2, tFDR = 0.1, seed = 9, verbose = FALSE,
                             control = trex_spca_control(K = 10))
spca$select()
spca$active_sets

## ----fields-------------------------------------------------------------------
dim(spca$Z)
dim(spca$V)
round(spca$adjusted_ev, 3)
round(spca$cumulative_ev, 3)

## ----loadings, fig.cap = "The two sparse loading vectors. The true supports are variables 1 to 10 and 31 to 40."----
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (m in 1:2) {
  plot(spca$V[, m], type = "h", lwd = 2, col = "#336C68", xlab = "variable",
       ylab = sprintf("loading %d", m))
}
par(op)

## ----thresholded--------------------------------------------------------------
thr <- TRexSPCASelector$new(X, M = 2, tFDR = 0.1, seed = 9, verbose = FALSE,
                            control = trex_spca_control(mode = "Thresholded", K = 10))
thr$select()
thr$active_sets
round(thr$cumulative_ev, 3)

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

