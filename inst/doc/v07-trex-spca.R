## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
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

## ----loadings-----------------------------------------------------------------
library(ggplot2)
loadings <- data.frame(variable = rep(seq_len(p), 2), loading = c(spca$V[, 1], spca$V[, 2]),
                       component = rep(c("loading 1", "loading 2"), each = p))
on_support <- (loadings$component == "loading 1" & loadings$variable %in% 1:10) |
  (loadings$component == "loading 2" & loadings$variable %in% 31:40)
loadings$kind <- ifelse(on_support, "true support", "other")
ggplot(loadings, aes(variable, loading, colour = kind)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_segment(aes(xend = variable, yend = 0), linewidth = 0.8) +
  facet_wrap(~ component) +
  scale_colour_manual(values = c("true support" = "#336C68", "other" = "#D1495B")) +
  labs(x = "Variable", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "white", colour = "grey60"))

## ----thresholded--------------------------------------------------------------
thr <- TRexSPCASelector$new(X, M = 2, tFDR = 0.1, seed = 9, verbose = FALSE,
                            control = trex_spca_control(mode = "Thresholded", K = 10))
thr$select()
thr$active_sets
round(thr$cumulative_ev, 3)

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

