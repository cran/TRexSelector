## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4, fig.retina = 2)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----data---------------------------------------------------------------------
set.seed(31)
n <- 100
p <- 50
L <- 50
true_support <- c(2, 9, 17, 33)
beta <- numeric(p)
beta[true_support] <- c(2, -2, 1.5, -1.5)
X <- matrix(rnorm(n * p), n, p)
D <- matrix(rnorm(n * L), n, L)
y <- as.numeric(X %*% beta + rnorm(n))

## ----tlars--------------------------------------------------------------------
solver <- TLARS_Solver$new(X, D, y)
solver$execute_step(T_stop = 3)
solver$get_num_steps()
solver$get_active_predictor_indices()
solver$get_active_dummy_indices()

## ----actives------------------------------------------------------------------
solver$get_actives()
solver$get_dummy_start_index()
solver$get_actions()[1:4]

## ----path---------------------------------------------------------------------
dim(solver$get_beta_path())
round(solver$get_beta()[true_support], 3)
path_diag <- data.frame(lambda = solver$get_lambda(), rss = solver$get_rss(),
                        r2 = solver$get_r2(), cp = solver$get_cp(), dof = solver$get_dof())
round(head(path_diag, 5), 3)
solver$restore()    # give X and D back before the next solver claims them

## ----path-figure, fig.width = 9, fig.height = 3.6, out.width = "100%"---------
library(ggplot2)

snapshot <- function(solver, label) {
  path <- solver$get_beta_path()                            # (p + L) x (steps + 1)
  entering <- abs(vapply(solver$get_actions(), function(a) a[1], numeric(1)))
  moved <- which(rowSums(path != 0) > 0)
  kind <- ifelse(moved > p, "dummy", ifelse(moved %in% true_support, "true predictor", "other predictor"))
  lines <- data.frame(panel = label, variable = rep(moved, each = ncol(path)),
                      step = rep(seq_len(ncol(path)) - 1, times = length(moved)),
                      coefficient = as.vector(t(path[moved, , drop = FALSE])),
                      kind = rep(kind, each = ncol(path)))
  marks <- data.frame(panel = label, step = seq_along(entering), dummy = entering > p,
                      label = ifelse(entering > p, paste0("D", entering - p), entering))
  list(lines = lines, marks = marks)
}

s <- TLARS_Solver$new(X, D, y)
s$execute_step(T_stop = 1)
first <- snapshot(s, "T_stop = 1")
s$execute_step(T_stop = 3)
second <- snapshot(s, "T_stop = 3, continued")
s$restore()

lines <- rbind(first$lines, second$lines)
marks <- rbind(first$marks, second$marks)
lines$panel <- factor(lines$panel, levels = c("T_stop = 1", "T_stop = 3, continued"))
marks$panel <- factor(marks$panel, levels = levels(lines$panel))
inclusions <- marks[marks$dummy, ]
inclusions$kind <- "dummy inclusion"
kinds <- c("true predictor", "other predictor", "dummy", "dummy inclusion")
ggplot(lines, aes(step, coefficient, group = variable, colour = kind, linetype = kind,
                  linewidth = kind)) +
  geom_vline(data = inclusions, aes(xintercept = step, colour = kind, linetype = kind,
                                    linewidth = kind), inherit.aes = FALSE) +
  geom_line() +
  geom_text(data = marks, aes(x = step, y = Inf, label = label), inherit.aes = FALSE,
            vjust = 1.4, size = 2.4, colour = "grey30") +
  facet_wrap(~ panel, scales = "free_x") +
  scale_colour_manual(values = c("true predictor" = "#336C68", "other predictor" = "grey45",
                                 "dummy" = "#D1495B", "dummy inclusion" = "#D1495B"),
                      breaks = kinds) +
  scale_linetype_manual(values = c("true predictor" = "solid", "other predictor" = "solid",
                                   "dummy" = "dashed", "dummy inclusion" = "dotted"),
                        breaks = kinds) +
  scale_linewidth_manual(values = c("true predictor" = 0.9, "other predictor" = 0.5,
                                    "dummy" = 0.6, "dummy inclusion" = 0.4), breaks = kinds) +
  labs(x = "Step", y = "Coefficient", colour = NULL, linetype = NULL, linewidth = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "white", colour = "grey60"))

## ----tenet--------------------------------------------------------------------
en <- TENET_Solver$new(X, D, y, options = tsolver_options(lambda2 = 0.5))
en$execute_step(T_stop = 3)
en$get_active_predictor_indices()
en$get_num_removals()
en$restore()

## ----writer-------------------------------------------------------------------
X_before <- X + 0
w <- TLASSO_Solver$new(X, D, y)
max(abs(X - X_before))      # the scaled columns
w$restore()
max(abs(X - X_before))      # the original values

## ----readers------------------------------------------------------------------
sx <- MLLpNormScaler$new(norm_type = 2, center = TRUE)
invisible(sx$fit_transform_inplace(X))
sd <- MLLpNormScaler$new(norm_type = 2, center = TRUE)
invisible(sd$fit_transform_inplace(D))
y_c <- y - mean(y)
reader_opts <- tsolver_options(normalize = FALSE, intercept = FALSE, precentered = TRUE)

## ----families-----------------------------------------------------------------
constructors <- list(
  TLARS = TLARS_Solver, TLASSO = TLASSO_Solver, TSTEPWISE = TSTEPWISE_Solver,
  TSTAGEWISE = TSTAGEWISE_Solver, TENET = TENET_Solver, TENETAug = TENETAug_Solver,
  TOMP = TOMP_Solver, TOOLS = TOOLS_Solver, TGP = TGP_Solver, TACGP = TACGP_Solver,
  TMP = TMP_Solver, TNCGMP = TNCGMP_Solver, TAFS = TAFS_Solver,
  TCDLASSO = TCDLASSO_Solver, TCDENET = TCDENET_Solver)
rows <- lapply(names(constructors), function(nm) {
  s <- constructors[[nm]]$new(X, D, y_c, options = reader_opts)
  s$execute_step(T_stop = 3)
  data.frame(solver = nm, steps = s$get_num_steps(),
             predictors = paste(s$get_active_predictor_indices(), collapse = " "))
})
do.call(rbind, rows)

## ----structured---------------------------------------------------------------
groups <- rep(1:5, each = 10)
ien <- TIENET_Solver$new(X, D, y_c, groups = groups,
                         options = tsolver_options(normalize = FALSE, intercept = FALSE,
                                                   precentered = TRUE, lambda2 = 0.5))
ien$execute_step(T_stop = 3)
ien$get_active_predictor_indices()
tik <- TCDTIKNET_Solver$new(X, D, y_c, tikhonov_K = diag(p),
                            options = tsolver_options(normalize = FALSE, intercept = FALSE,
                                                      precentered = TRUE, lambda2 = 0.5))
tik$execute_step(T_stop = 3)
tik$get_active_predictor_indices()

## ----save-load----------------------------------------------------------------
f <- tempfile(fileext = ".bin")
a <- TLARS_Solver$new(X, D, y_c, options = reader_opts)
a$execute_step(T_stop = 2)
a$save(f)
b <- TLARS_Solver$new(X, D, y_c, options = reader_opts)
b$load(f)
b$get_num_steps()
b$execute_step(T_stop = 4)
c <- TLARS_Solver$new(X, D, y_c, options = reader_opts)
c$execute_step(T_stop = 4)
identical(b$get_actives(), c$get_actives())
unlink(f)

## ----memory-------------------------------------------------------------------
a$get_memory_bytes()

## ----restore-data, include = FALSE--------------------------------------------
invisible(sx$inverse_transform_inplace(X))
invisible(sd$inverse_transform_inplace(D))

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

