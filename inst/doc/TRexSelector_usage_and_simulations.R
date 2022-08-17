## ---- include = FALSE---------------------------------------------------------
# Store user's options()
old_options <- options()

library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.retina = 2,
  out.width = "85%",
  dpi = 96
  # pngquant = "--speed=1"
)
options(width = 80)

## ---- eval=FALSE--------------------------------------------------------------
#  # Install stable version from CRAN
#  install.packages("tlars")
#  
#  # Install development version from GitHub
#  install.packages("devtools")
#  devtools::install_github("jasinmachkour/tlars")

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("jasinmachkour/TRexSelector")

## ---- eval=FALSE--------------------------------------------------------------
#  library(TRexSelector)
#  help(package = "TRexSelector")
#  ?trex
#  ?random_experiments
#  ?lm_dummy
#  ?add_dummies
#  ?add_dummies_GVS
#  ?FDP
#  ?TPP
#  # etc.

## ---- eval=FALSE--------------------------------------------------------------
#  citation("TRexSelector")

## -----------------------------------------------------------------------------
library(TRexSelector)

# Setup
n <- 75 # number of observations
p <- 150 # number of variables
num_act <- 3 # number of true active variables
beta <- c(rep(1, times = num_act), rep(0, times = p - num_act)) # coefficient vector
true_actives <- which(beta > 0) # indices of true active variables
num_dummies <- p # number of dummy predictors (also referred to as dummies)

# Generate Gaussian data
set.seed(123)
X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
y <- X %*% beta + stats::rnorm(n)

## -----------------------------------------------------------------------------
# Seed
set.seed(1234)

# Numerical zero
eps <- .Machine$double.eps

# Variable selection via T-Rex
res <- trex(X = X, y = y, tFDR = 0.05, verbose = FALSE)
selected_var <- which(res$selected_var > eps)
paste0("True active variables: ", paste(as.character(true_actives), collapse = ", "))
paste0("Selected variables: ", paste(as.character(selected_var), collapse = ", "))

## -----------------------------------------------------------------------------
# Computations might take up to 10 minutes... Please wait... 

# Numerical zero
eps <- .Machine$double.eps

# Seed
set.seed(1234)

# Setup
n <- 100 # number of observations
p <- 150 # number of variables

# Parameters
num_act <- 10 # number of true active variables
beta <- rep(0, times = p) # coefficient vector (all zeros first)
beta[sample(seq(p), size = num_act, replace = FALSE)] <- 1 # coefficient vector (active variables with non-zero coefficients)
true_actives <- which(beta > 0) # indices of true active variables
tFDR_vec <- c(0.1, 0.15, 0.2, 0.25) # target FDR levels
MC <- 100 # number of Monte Carlo runs per stopping point

# Initialize results vectors
FDP <- matrix(NA, nrow = MC, ncol = length(tFDR_vec))
TPP <- matrix(NA, nrow = MC, ncol = length(tFDR_vec))

# Run simulations
for (t in seq_along(tFDR_vec)) {
  for (mc in seq(MC)) {
    
    # Generate Gaussian data
    X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    y <- X %*% beta + stats::rnorm(n)
    
    # Run T-Rex selector
    res <- trex(X = X, y = y, tFDR = tFDR_vec[t], verbose = FALSE)
    selected_var <- which(res$selected_var > eps)
    
    # Results
    FDP[mc, t] <- length(setdiff(selected_var, true_actives)) / max(1, length(selected_var))
    TPP[mc, t] <- length(intersect(selected_var, true_actives)) / max(1, length(true_actives))
  }
}

# Compute estimates of FDR and TPR by averaging FDP and TPP over MC Monte Carlo runs
FDR <- colMeans(FDP)
TPR <- colMeans(TPP)

## ----FDR_and_TPR, echo=FALSE, fig.align='center', message=FALSE, fig.width=12, fig.height=5, out.width = "95%"----
# Plot results
library(ggplot2)
library(patchwork)
tFDR_vec_percent <- 100 * tFDR_vec
plot_data <- data.frame(tFDR_vec = tFDR_vec_percent,
                        FDR = 100 * FDR,
                        TPR = 100 * TPR) # data frame containing data to be plotted (FDR and TPR in %)

# FDR vs. tFDR
FDR_vs_tFDR <-
  ggplot(plot_data, aes(x = tFDR_vec_percent, y = FDR)) +
    labs(x = "Target FDR",
         y = "FDR") +
    scale_x_continuous(breaks = tFDR_vec_percent, minor_breaks = c(), limits = c(tFDR_vec_percent[1], tFDR_vec_percent[length(tFDR_vec_percent)])) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), minor_breaks = c(), limits = c(0, 100)) +
    geom_line(size = 1.5, colour = "#336C68") +
    geom_abline(slope = 1, colour = "red", linetype = 2, size = 1) +
    geom_point(size = 2.5, colour = "#336C68") +
    theme_bw(base_size = 16) +
    theme(panel.background = element_rect(fill = "white", color = "black", size = 1)) +
    coord_fixed(ratio =  0.75 * (tFDR_vec_percent[length(tFDR_vec_percent)] - tFDR_vec_percent[1]) / (100 - 0))

# TPR vs. tFDR
TPR_vs_tFDR <- 
  ggplot(plot_data, aes(x = tFDR_vec_percent, y = TPR)) +
    labs(x = "Target FDR",
         y = "TPR") +
    scale_x_continuous(breaks = tFDR_vec_percent, minor_breaks = c(), limits = c(tFDR_vec_percent[1], tFDR_vec_percent[length(tFDR_vec_percent)])) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), minor_breaks = c(), limits = c(0, 100)) +
    geom_line(size = 1.5, colour = "#336C68") +
    geom_point(size = 2.5, colour = "#336C68") +
    theme_bw(base_size = 16) +
    theme(panel.background = element_rect(fill = "white", color = "black", size = 1)) +
    coord_fixed(ratio =  0.75 * (tFDR_vec_percent[length(tFDR_vec_percent)] - tFDR_vec_percent[1]) / (100 - 0))
FDR_vs_tFDR + TPR_vs_tFDR

## ----T-Rex_framework, echo=FALSE, fig.cap="Figure 1: Simplified overview of the T-Rex framework.", out.width = '85%'----
knitr::include_graphics("./figures/T-Rex_framework.png")

## ----EnlargedPredictorMatrix, echo=FALSE, fig.cap="Figure 2: The enlarged predictor matrix (predictor matrix with dummies).", out.width = '65%'----
knitr::include_graphics("./figures/predictor_matrix_with_dummies.png")

## ---- include = FALSE---------------------------------------------------------
# Reset user's options()
options(old_options)

