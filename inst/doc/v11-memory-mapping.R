## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----agnostic, echo = FALSE, out.width = "90%"--------------------------------
knitr::include_graphics("figures/memory_agnostic.png")

## ----create-------------------------------------------------------------------
set.seed(5)
n <- 100
p <- 200
X <- matrix(rnorm(n * p), n, p)
path <- tempfile(fileext = ".bin")
mm <- mem_matrix_from(X, path)
class(mm)
dim(mm)
file.size(path) == 8 * n * p

## ----ops----------------------------------------------------------------------
round(mm[1:2, 1:3], 3)
identical(mm[1:2, 1:3], X[1:2, 1:3])
mm[1, 1] <- 42
mm[1, 1]
mm[1, 1] <- X[1, 1]
max(abs(as.matrix(mm) - X))
ro <- mem_matrix(path, rows = n, cols = p, mode = "readonly")
tryCatch(ro[1, 1] <- 0, error = function(e) cat("error:", conditionMessage(e), "\n"))

## ----stream-------------------------------------------------------------------
path2 <- tempfile(fileext = ".bin")
con <- file(path2, open = "wb")
for (j in seq_len(p)) {
  set.seed(1000 + j)
  writeBin(rnorm(n), con)
}
close(con)
big <- mem_matrix(path2, rows = n, cols = p, mode = "readwrite")
set.seed(1000 + 3)
max(abs(as.numeric(big[, 3]) - rnorm(n)))

## ----select-------------------------------------------------------------------
beta <- c(rep(1.5, 4), rep(0, p - 4))
y <- as.numeric(X %*% beta + rnorm(n))
sel_mem <- TRexSelector$new(X, y, tFDR = 0.1, seed = 3, verbose = FALSE,
                            control = trex_control(K = 10))
sel_mem$select()
sel_map <- TRexSelector$new(mm, y, tFDR = 0.1, seed = 3, verbose = FALSE,
                            control = trex_control(K = 10))
sel_map$select()
identical(sel_mem$selected_indices, sel_map$selected_indices)
max(abs(as.matrix(mm) - X))      # the file holds the original values again

## ----solver-------------------------------------------------------------------
D <- matrix(rnorm(n * p), n, p)
s <- TLARS_Solver$new(mm, D, y)
s$execute_step(T_stop = 2)
s$get_active_predictor_indices()
s$restore()
z <- MLZScoreScaler$new()
invisible(z$fit_transform_inplace(mm))
round(mean(mm[, 1]), 12)
invisible(z$inverse_transform_inplace(mm))
en <- MLElasticNet$new(low_memory = TRUE)
en$fit(mm, y, alpha = 1, n_lambda = 20)
dim(en$get_coef())

## ----cleanup, include = FALSE-------------------------------------------------
rm(mm, ro, big, s)
invisible(gc())
unlink(c(path, path2))
trex_set_num_threads(old_threads)

