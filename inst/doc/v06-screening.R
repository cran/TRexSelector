## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----data---------------------------------------------------------------------
set.seed(2024)
n <- 100
p <- 300
true_support <- c(3, 50, 120, 200, 280)
beta <- numeric(p)
beta[true_support] <- 2
X <- matrix(rnorm(n * p), n, p)
y <- as.numeric(X %*% beta + rnorm(n))

## ----screen-------------------------------------------------------------------
scr <- TRexScreeningSelector$new(X, y, seed = 3, verbose = FALSE,
                                 control = trex_control(K = 20))
scr$select()
scr$selected_indices
scr$estimated_FDR

## ----mats---------------------------------------------------------------------
dim(scr$beta_mat)
dim(scr$dummy_betas)

## ----boot---------------------------------------------------------------------
boot <- TRexScreeningSelector$new(X, y, seed = 3, verbose = FALSE,
                                  screen_control = trex_screen_control(use_bootstrap_CI = TRUE,
                                                                       R_boot = 200),
                                  control = trex_control(K = 20))
boot$select()
boot$used_bootstrap
boot$confidence_level
boot$estimated_FDR
boot$selected_indices

## ----da-----------------------------------------------------------------------
ar <- TRexScreeningSelector$new(X, y, seed = 3, verbose = FALSE,
                                screen_control = trex_screen_control(trex_method = "TREX_DA_AR1"),
                                control = trex_control(K = 20))
ar$select()
round(ar$estimated_correlation, 3)

## ----biobank------------------------------------------------------------------
Y <- cbind(y, as.numeric(X %*% (0.5 * beta) + rnorm(n)), rnorm(n))
colnames(Y) <- c("strong", "weak", "null")
bio <- TRexBiobankScreeningSelector$new(X, Y, tFDR = 0.1, seed = 3, verbose = FALSE,
                                        biobank_control = trex_biobank_control(0.05, 0.15),
                                        control = trex_control(K = 20))
res <- bio$select()
length(res)

## ----records------------------------------------------------------------------
summary_tab <- do.call(rbind, lapply(res, function(r)
  data.frame(phenotype = r$phenotype_index, selected = length(r$selected_indices),
             estimated_FDR = round(r$estimated_FDR, 3), method = r$method_used,
             fallback = r$used_fallback_trex)))
summary_tab
res[[1]]$selected_indices

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

