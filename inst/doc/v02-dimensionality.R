## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
options(width = 80)
library(TRexSelector)
# Vignettes build on at most two threads (the CRAN rule)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----regimes-map, echo = FALSE, fig.cap = "Figure 1: The three regimes in the plane of p and n.", out.width = "85%"----
knitr::include_graphics("figures/dimensionality_map.png")

## ----regimes------------------------------------------------------------------
set.seed(7)
n <- 100
num_act <- 5
trials <- 10
true_support <- seq_len(num_act)
for (p in c(50, 200, 1000)) {
  fdr <- 0
  tpr <- 0
  for (trial in seq_len(trials)) {
    X <- matrix(rnorm(n * p), n, p)
    beta <- c(rep(1, num_act), rep(0, p - num_act))
    y <- as.vector(X %*% beta + rnorm(n))
    ctrl <- trex_control(K = 10)
    sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 6 + trial, verbose = FALSE, control = ctrl)
    sel$select()
    fdr <- fdr + trex_fdp(sel$selected_indices, true_support) / trials
    tpr <- tpr + trex_tpp(sel$selected_indices, true_support) / trials
  }
  cat(sprintf("p = %4d  n/p = %5.2f  FDR = %.2f  TPR = %.2f\n", p, n / p, fdr, tpr))
}

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

