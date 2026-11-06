## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----data---------------------------------------------------------------------
set.seed(11)
n <- 100
p <- 100
group_size <- 5
groups <- rep(seq_len(p / group_size), each = group_size)
latent <- matrix(rnorm(n * max(groups)), n, max(groups))
X <- 0.9 * latent[, groups] + sqrt(1 - 0.9^2) * matrix(rnorm(n * p), n, p)
true_support <- which(groups %in% 1:2)
beta <- numeric(p)
beta[true_support] <- 1
y <- as.numeric(X %*% beta + rnorm(n))

## ----en-----------------------------------------------------------------------
gvs <- TRexGVSSelector$new(X, y, tFDR = 0.1, seed = 5, verbose = FALSE,
                           gvs_control = trex_gvs_control(gvs_type = "EN", lambda_2 = 1),
                           control = trex_control(K = 10))
gvs$select()
sort(gvs$selected_indices)
trex_fdp(gvs$selected_indices, true_support)
trex_tpp(gvs$selected_indices, true_support)

## ----fields-------------------------------------------------------------------
gvs$gvs_type
gvs$max_clusters
gvs$hc_method_used
gvs$lambda2_used
table(gvs$groups)[1:5]

## ----cv-----------------------------------------------------------------------
cv <- TRexGVSSelector$new(X, y, tFDR = 0.1, seed = 5, verbose = FALSE,
                          gvs_control = trex_gvs_control(gvs_type = "EN", cv_n_lambda = 100),
                          control = trex_control(K = 10))
cv$select()
cv$lambda2_used
str(cv$lambda2_cv_diagnostics)

## ----ien----------------------------------------------------------------------
ien <- TRexGVSSelector$new(X, y, tFDR = 0.1, seed = 5, verbose = FALSE,
                           gvs_control = trex_gvs_control(gvs_type = "IEN", groups = groups,
                                                          lambda_2 = 1),
                           control = trex_control(solver = "TIENET", K = 10))
ien$select()
sort(ien$selected_indices)
trex_fdp(ien$selected_indices, true_support)
trex_tpp(ien$selected_indices, true_support)

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

