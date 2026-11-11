## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----framework, echo = FALSE, fig.cap = "Figure 1: The T-Rex framework.", out.width = "85%"----
knitr::include_graphics("figures/T-Rex_framework.png")

## ----enlarged, echo = FALSE, fig.cap = "Figure 2: The enlarged design (predictors and dummies).", out.width = "65%"----
knitr::include_graphics("figures/predictor_matrix_with_dummies.png")

## ----data---------------------------------------------------------------------
set.seed(123)
n <- 100
p <- 150
beta <- c(rep(1, 5), rep(0, p - 5))
true_support <- which(beta != 0)
X <- matrix(rnorm(n * p), n, p)
y <- as.numeric(X %*% beta + rnorm(n))

## ----run----------------------------------------------------------------------
sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 42, verbose = FALSE)
sel$select()
sel$selected_indices

## ----K------------------------------------------------------------------------
sel$voting_grid

## ----fields-------------------------------------------------------------------
sel$T_stop
sel$L
sel$v_thresh

## ----fdp-hat------------------------------------------------------------------
fdp <- sel$fdp_hat_mat
dim(fdp)
round(fdp[sel$T_stop, ], 3)
sel$r_mat

## ----phi, fig.cap = "Relative occurrences after the last stopping time. The true support stands out."----
library(ggplot2)
occurrence <- data.frame(variable = seq_len(p), phi = sel$phi_mat[sel$T_stop, ],
                         kind = ifelse(seq_len(p) %in% true_support, "true support", "other"))
threshold <- data.frame(v = sel$v_thresh, rule = "voting threshold")
ggplot(occurrence, aes(variable, phi, colour = kind)) +
  geom_segment(aes(xend = variable, yend = 0), linewidth = 0.8) +
  geom_hline(data = threshold, aes(yintercept = v, linetype = rule), colour = "#D1495B") +
  scale_colour_manual(values = c("true support" = "#336C68", "other" = "grey60")) +
  scale_linetype_manual(values = c("voting threshold" = "dashed")) +
  labs(x = "Variable", y = "Relative occurrence", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "grey60"),
        legend.box.margin = margin(4, 6, 4, 6))

## ----strategies---------------------------------------------------------------
for (s in c("HCONCAT", "STANDARD", "SEEDED")) {
  ctrl <- trex_control(lloop_strategy = s, K = 10)
  s_sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 42, verbose = FALSE, control = ctrl)
  s_sel$select()
  cat(sprintf("%-10s L = %d  T_stop = %d  selected: %s\n", s, s_sel$L, s_sel$T_stop,
              paste(s_sel$selected_indices, collapse = " ")))
}

## ----dummies------------------------------------------------------------------
ctrl <- trex_control(dummy_distribution = trex_dummy_distribution("StudentT", df = 3), K = 10)
t_sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 42, verbose = FALSE, control = ctrl)
t_sel$select()
t_sel$selected_indices

## ----none---------------------------------------------------------------------
d <- trex_normalize_design(X, y, scaling_mode = "L2")
round(colMeans(d$X)[1:3], 12)
sqrt(colSums(d$X^2))[1:3]
ctrl <- trex_control(scaling_mode = "NONE")
n_sel <- TRexSelector$new(d$X, d$y, tFDR = 0.1, seed = 42, verbose = FALSE, control = ctrl)
n_sel$select()
identical(n_sel$selected_indices, sel$selected_indices)

## ----restore------------------------------------------------------------------
X_before <- X + 0
r_sel <- TRexSelector$new(X, y, tFDR = 0.1, seed = 42, verbose = FALSE)
max(abs(X - X_before))      # the working design
r_sel$restore()
max(abs(X - X_before))      # the original values are back

## ----mc-----------------------------------------------------------------------
tFDR_vec <- c(0.1, 0.2, 0.3)
MC <- 20
fdp <- tpp <- matrix(NA_real_, MC, length(tFDR_vec))
for (t in seq_along(tFDR_vec)) {
  for (mc in seq_len(MC)) {
    set.seed(1000 * t + mc)
    X_mc <- matrix(rnorm(n * p), n, p)
    y_mc <- as.numeric(X_mc %*% beta + rnorm(n))
    ctrl <- trex_control(K = 10)
    s <- TRexSelector$new(X_mc, y_mc, tFDR = tFDR_vec[t], seed = 1000 * t + mc,
                          verbose = FALSE, control = ctrl)
    s$select()
    fdp[mc, t] <- trex_fdp(s$selected_indices, true_support)
    tpp[mc, t] <- trex_tpp(s$selected_indices, true_support)
  }
}
data.frame(tFDR = tFDR_vec, FDR = colMeans(fdp), TPR = colMeans(tpp))

## ----mc-plot------------------------------------------------------------------
levels_fdr_tpr <- c("FDR", "TPR")
study <- rbind(data.frame(tFDR = tFDR_vec, value = colMeans(fdp), metric = "FDR"),
               data.frame(tFDR = tFDR_vec, value = colMeans(tpp), metric = "TPR"))
study$metric <- factor(study$metric, levels = levels_fdr_tpr)
diagonal <- data.frame(metric = factor("FDR", levels = levels_fdr_tpr), slope = 1, intercept = 0,
                       line = "FDR equal to the target")
limits <- data.frame(tFDR = 0.2, value = c(0, 0.5, 0, 1),
                     metric = factor(rep(levels_fdr_tpr, each = 2), levels = levels_fdr_tpr))
ggplot(study, aes(tFDR, value)) +
  geom_blank(data = limits) +
  geom_abline(data = diagonal, aes(slope = slope, intercept = intercept, linetype = line),
              colour = "#D1495B") +
  geom_line(aes(colour = "mean over 20 trials")) +
  geom_point(aes(colour = "mean over 20 trials"), size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = c("mean over 20 trials" = "#336C68")) +
  scale_linetype_manual(values = c("FDR equal to the target" = "dashed")) +
  labs(x = "Target FDR", y = NULL, colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "grey60"),
        legend.box.margin = margin(4, 6, 4, 6))

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

