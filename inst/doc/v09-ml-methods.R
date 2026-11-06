## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.width = 6, fig.height = 4)
options(width = 80)
library(TRexSelector)
old_threads <- trex_set_num_threads(min(2L, trex_get_max_threads()))

## ----data---------------------------------------------------------------------
set.seed(8)
n <- 120
p <- 30
X <- matrix(rnorm(n * p), n, p)
beta <- c(2, -1.5, 0, 0, 1, rep(0, p - 5))
y <- as.numeric(X %*% beta + rnorm(n))

## ----scalers------------------------------------------------------------------
z <- MLZScoreScaler$new()
z$fit(X)
round(z$get_centers()[1:3], 3)
round(z$get_scales()[1:3], 3)
X_before <- X + 0
invisible(z$transform_inplace(X))
round(colMeans(X)[1:3], 12)
round(apply(X, 2, sd)[1:3], 12)
invisible(z$inverse_transform_inplace(X))
max(abs(X - X_before))

## ----lp-----------------------------------------------------------------------
l2 <- MLLpNormScaler$new(norm_type = 2, center = TRUE)
invisible(l2$fit_transform_inplace(X))
sqrt(colSums(X^2))[1:3]
invisible(l2$inverse_transform_inplace(X))

## ----ridge--------------------------------------------------------------------
ridge <- MLRidge$new()
round(ridge$solve(X, y, lambda = 1)[1:5], 3)
cv <- MLRidgeCV$new()
cv$fit(X, y, n_folds = 5, n_lambda = 50)
c(lambda_min = cv$cv_min(), lambda_1se = cv$cv_1se())

## ----ridge-plot---------------------------------------------------------------
plot(cv$get_lambdas(), cv$get_cv_errors(), log = "x", type = "l", col = "#336C68", lwd = 2,
     xlab = "lambda", ylab = "CV error")
abline(v = c(cv$cv_min(), cv$cv_1se()), lty = c(1, 2), col = "red")

## ----enet---------------------------------------------------------------------
en <- MLElasticNet$new()
en$fit(X, y, alpha = 0.5, n_lambda = 50)
dim(en$get_coef())
length(en$get_lambdas())
en$converged()
B <- en$get_coef()               # the path can stop before n_lambda points
round(B[1:5, ncol(B)], 3)
round(en$predict(X[1:3, ])[, ncol(B)], 3)

## ----enet-cv------------------------------------------------------------------
encv <- MLElasticNetCV$new()
encv$fit(X, y, alpha = 0.5, n_folds = 5, n_lambda = 50)
c(lambda_min = encv$cv_min(), lambda_1se = encv$cv_1se())

## ----svd----------------------------------------------------------------------
svd_res <- MLSVD$new()$compute(X, M = 3)
names(svd_res)
dim(svd_res$U)
round(svd_res$S, 3)

## ----pca----------------------------------------------------------------------
X_before <- X + 0
pca <- MLPCA$new(X, M = 3, center = TRUE, normalize = TRUE)
pca$fit()
dim(pca$get_scores())
dim(pca$get_loadings())
round(pca$get_explained_variance(), 3)
round(pca$transform(X_before[1:2, ]), 3)
pca$restore()
max(abs(X - X_before))

## ----hac----------------------------------------------------------------------
set.seed(9)
centers <- rep(c(-3, 0, 3), each = 20)
Z <- centers + matrix(rnorm(60 * 4), 60, 4)
linkage <- ml_agglomerative_cluster(Z, method = ml_linkage_method$Ward,
                                    metric = ml_distance_metric$Euclidean)
dim(linkage)
labels <- ml_cut_tree(linkage, num_orig_objs = 60, num_clusters = 3)
table(labels, rep(1:3, each = 20))

## ----hac-vars-----------------------------------------------------------------
link_vars <- ml_agglomerative_cluster(t(X), method = ml_linkage_method$Single,
                                      metric = ml_distance_metric$Correlation)
head(ml_cut_tree(link_vars, num_orig_objs = p, num_clusters = 5))

## ----teardown, include = FALSE------------------------------------------------
trex_set_num_threads(old_threads)

