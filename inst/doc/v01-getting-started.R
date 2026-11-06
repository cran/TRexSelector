## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4)
options(width = 80)
library(TRexSelector)
# Vignettes build on at most two threads (the CRAN rule)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----eval = FALSE-------------------------------------------------------------
# install.packages(c("Rcpp", "RcppEigen", "BH", "Rcereal", "R6"))

## ----eval = FALSE-------------------------------------------------------------
# install.packages("remotes")
# remotes::install_github("trexselector/TRexSelector", subdir = "R/TRexSelector")

## ----first-run----------------------------------------------------------------
data(Gauss_data)
X <- Gauss_data$X
y <- Gauss_data$y
true_support <- which(Gauss_data$support)
true_support

## ----first-select-------------------------------------------------------------
sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 1234, verbose = FALSE)
sel$select()
sel$selected_indices

## ----metrics------------------------------------------------------------------
trex_fdp(sel$selected_indices, true_support)
trex_tpp(sel$selected_indices, true_support)

## ----fields-------------------------------------------------------------------
sel$num_selected   # number of selected variables
sel$T_stop         # stopping time: dummies included when the T-loop stopped
sel$L              # calibrated dummy multiplier, the selector used L * p dummies
sel$v_thresh       # voting level the selection was read at
sel$K              # random experiments

## ----matrices-----------------------------------------------------------------
dim(sel$fdp_hat_mat)
head(sel$voting_grid)
round(sel$fdp_hat_mat[sel$T_stop, ], 3)

## ----solver-choice------------------------------------------------------------
ctrl <- trex_control(solver = "TLASSO", K = 10)
sel2 <- TRexSelector$new(X, y, tFDR = 0.1, seed = 1234, verbose = FALSE, control = ctrl)
sel2$select()
sel2$selected_indices

## ----eval = FALSE-------------------------------------------------------------
# citation("TRexSelector")

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

