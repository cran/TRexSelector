## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----double-------------------------------------------------------------------
Xi <- matrix(sample(0:2, 60 * 20, replace = TRUE), 60, 20)
storage.mode(Xi)
storage.mode(Xi) <- "double"     # one explicit conversion, then zero-copy from here on
storage.mode(Xi)

## ----selector-----------------------------------------------------------------
set.seed(1)
n <- 80
p <- 40
X <- matrix(rnorm(n * p), n, p)
y <- as.numeric(X[, 1] * 2 + rnorm(n))
X_before <- X + 0
sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 2, verbose = FALSE, control = trex_control(K = 10))
round(sqrt(colSums(X^2))[1:3], 6)     # the working design: unit columns
sel$select()
max(abs(X - X_before))                 # the original values are back

## ----writers------------------------------------------------------------------
D <- matrix(rnorm(n * p), n, p)
w1 <- TLARS_Solver$new(X, D, y)
tryCatch(TLASSO_Solver$new(X, D, y), error = function(e) cat("error:", conditionMessage(e), "\n"))
w1$restore()

## ----readers------------------------------------------------------------------
sx <- MLLpNormScaler$new(norm_type = 2, center = TRUE)
invisible(sx$fit_transform_inplace(X))
sd <- MLLpNormScaler$new(norm_type = 2, center = TRUE)
invisible(sd$fit_transform_inplace(D))
y_c <- y - mean(y)
opts <- tsolver_options(normalize = FALSE, intercept = FALSE, precentered = TRUE)
r1 <- TLARS_Solver$new(X, D, y_c, options = opts)
r2 <- TLASSO_Solver$new(X, D, y_c, options = opts)
r1$execute_step(T_stop = 2)
r2$execute_step(T_stop = 2)
invisible(sx$inverse_transform_inplace(X))
invisible(sd$inverse_transform_inplace(D))
max(abs(X - X_before))

## ----low-memory---------------------------------------------------------------
cv_a <- MLRidgeCV$new()
cv_a$fit(X, y, n_folds = 5, n_lambda = 30)
cv_b <- MLRidgeCV$new(low_memory = TRUE)
cv_b$fit(X, y, n_folds = 5, n_lambda = 30)
all.equal(cv_a$get_cv_errors(), cv_b$get_cv_errors())

## ----bytes--------------------------------------------------------------------
r1$get_memory_bytes()

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

