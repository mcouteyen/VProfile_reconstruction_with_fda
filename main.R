# ---- load library and sources ----
required_pkgs <- c("ggplot2", "dplyr", "fda")
invisible(lapply(required_pkgs, library, character.only = TRUE))

source("src/my_fda.R")
source("src/profile_reconstruction.R")
source("src/visualisation.R")

# ---- load data ----

#data must have the following columns : "INDEX", "DEPTH"
# sparse_data <- ("")
#name of the variable to predict
my_var <- "TEMPERATURE"

# ---- check data ----


## ---- use example data ----
#
# uncomment to use example data
#
source("src/profile_tempering.R")
my_data <- read.csv("data/example_data.csv")
sparse_data <- shrink_random_profiles_with_seed(my_data)


# ---- profiles visualisation shiny ----
#
# uncomment to visualise your data smoothed 
#
# required_pkgs <- c("shiny", "bslib", "shinycssloaders")
# invisible(lapply(required_pkgs, library, character.only = TRUE))
# 
# source("src/my_app.R")
# shinyApp(ui, server)


# ---- init variables ----

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
proj_pen_sigma <- .01

#depths to evaluate prediction
start <- 1
end <- 50
step <- .2


# ---- run reconstruction ----
source("src/run_pace_on_sparse_profiles.R")

### ---- plots ----
# plot fpca plane
plot_scores_projected(fpca, scr)

# select reconstruction to visualise in incomplete_indexes
print(incomplete_indexes)

index <- 90
selected_indexes <- which(incomplete_indexes == index)

plot_one_profile(index, sparse_data, predicted, var_label = my_var, CI = TRUE) + theme(legend.position = "none")

# reconstructed incomplete profiles 
head(reconstructed_profiles)

# complete profiles and reconstructed incomplete profiles 
all_data <- rbind(complete_data, reconstructed_profiles)

# reconstructed incomplete profiles with pointwise confidence intervals (CI) and confidence bands (BI)
head(predicted)
