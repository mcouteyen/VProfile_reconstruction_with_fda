# ---- LOAD LIBRARY AND SOURCES ----
required_pkgs <- c("ggplot2", "dplyr", "fda")
invisible(lapply(required_pkgs, library, character.only = TRUE))

source("src/my_fda.R")
source("src/profile_reconstruction.R")
source("src/visualisation.R")

# ---- LOAD DATA ----

#data must have the following columns : "INDEX", "DEPTH"
# sparse_data <- ("INSERT DATA PATH")
#name of the variable to predict
my_var <- "TEMPERATUrRE"

## ---- USE EXAMPLE DATA ----
#
# uncomment to use example data
#
sparse_data  <- read.csv("data/example_data.csv")


# ---- SHINY APP : PROFILE VIZ ----
#
# uncomment to visualise your data smoothed 
# and test your smoothing parameters
#
# required_pkgs <- c("shiny", "bslib", "shinycssloaders")
# invisible(lapply(required_pkgs, library, character.only = TRUE))
# 
# source("src/my_app.R")
# shinyApp(ui, server)


# ---- INIT VARIABLES ----

#if a profiles starts after this depth it is considered truncated
start_limit_complete_profile <- 1.5 

#if a profiles ends before this depth it is considered truncated
end_limit_complete_profile <- 49.5 

#number of basis functions (K)
bspline_smooth_nbasis <- 100 

#smoothing parameter (lambda)
bspline_smooth_lambda <- 1e-3 

#number of principal components to use for reconstruction
n_pc <- 30

#sigma^2
#if small, model fits the observed data more closely
#if it is large, the fit is smoother
proj_pen_sigma <- .001

#depths to evaluate prediction
start <- 1
end <- 50
step <- .2


# ---- RUN RECONSTRUCTION----
source("src/run_pace_on_sparse_profiles.R")

# ---- VISUALISATION ----
# plot fpca plane
plot_scores_projected(fpca, scr)

# select which prediction to visualise in incomplete_indexes
index <- incomplete_indexes[8]

plot_one_profile(index, sparse_data, predicted, var_label = my_var, CI = TRUE) + theme(legend.position = "none")

# reconstructed incomplete profiles 
head(reconstructed_profiles)

# complete profiles and reconstructed incomplete profiles 
all_data <- rbind(complete_data, reconstructed_profiles)
head(all_data)

# reconstructed incomplete profiles with pointwise confidence intervals (CI) and confidence bands (CB)
head(predicted)
