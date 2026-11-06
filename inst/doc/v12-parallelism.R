## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----threads------------------------------------------------------------------
trex_get_num_threads()
trex_get_max_threads()

## ----data---------------------------------------------------------------------
set.seed(12)
n <- 100
p <- 150
beta <- c(rep(1, 5), rep(0, p - 5))
X <- matrix(rnorm(n * p), n, p)
y <- as.numeric(X %*% beta + rnorm(n))

## ----parallel-----------------------------------------------------------------
serial <- TRexSelector$new(X, y, tFDR = 0.1, seed = 42, verbose = FALSE,
                           control = trex_control(K = 10))
serial$select()
parallel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 42, verbose = FALSE,
                             control = trex_control(K = 10, parallel_rnd_experiments = TRUE,
                                                    max_outer_threads = 2))
parallel$select()
budget <- TRexSelector$new(X, y, tFDR = 0.1, seed = 42, verbose = FALSE,
                           control = trex_control(K = 10, max_inner_threads = 2))
budget$select()
identical(serial$selected_indices, parallel$selected_indices)
identical(serial$selected_indices, budget$selected_indices)

## ----mc, eval = FALSE---------------------------------------------------------
# trials <- parallel::mclapply(seq_len(40), function(mc) {
#   trex_set_num_threads(1)
#   set.seed(mc)
#   X <- matrix(rnorm(n * p), n, p)
#   y <- as.numeric(X %*% beta + rnorm(n))
#   sel <- TRexSelector$new(X, y, tFDR = 0.1, verbose = FALSE)   # seed = -1 per trial
#   sel$select()
#   trex_fdp(sel$selected_indices, which(beta != 0))
# }, mc.cores = 4)
# mean(unlist(trials))

## ----cap, eval = FALSE--------------------------------------------------------
# old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))
# # ... the work ...
# trex_set_num_threads(old_threads)

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

