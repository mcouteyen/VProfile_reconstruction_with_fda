###############################################################################
# main.R
#
# Purpose:
# --------
# This script provides a complete workflow to reconstruct vertical profiles
# from sparse observations using functional data analysis (FDA).
# The method is described in the associated article.
#
# The script:
#   1. Loads required libraries and source files
#   2. Loads sparse vertical profile data
#   3. (Optional) Launches a Shiny app to explore smoothing parameters
#   4. Defines reconstruction and smoothing parameters
#   5. Runs the reconstruction algorithm
#   6. Visualizes reconstructed profiles 
#
# Author:
# -------
# M. Couteyen Carpaye
#
# Requirements:
# -------------
# - R >= 4.0
# - Packages: ggplot2, dplyr, fda
# - Folder structure as provided in the repository
#
###############################################################################


# ---- LOAD LIBRARIES AND SOURCE FILES ----------------------------------------

# Required R packages
required_pkgs <- c("ggplot2", "dplyr", "fda")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Source internal functions
source("src/my_fda.R")
source("src/profile_reconstruction.R")
source("src/visualisation.R")


# ---- LOAD DATA --------------------------------------------------------------

# Input data must contain at least the following columns:
#   - INDEX: profile identifier
#   - DEPTH: depth (positive values, increasing downward)
#   - one or more variables to be reconstructed (e.g. temperature, salinity)
#
# Example of custom data loading:
# sparse_data <- read.csv("INSERT_DATA_PATH")

# Name of the variable to reconstruct (must match a column name in the data)
my_var <- "TEMPERATURE"


# ---- USE EXAMPLE DATA --------------------------------------------------------
#
# Uncomment the following line to use the example dataset provided
# in the repository. This allows the user to test the full workflow
# without external data.
#
sparse_data <- read.csv("data/example_data.csv")


# ---- OPTIONAL: SHINY APP FOR PROFILE VISUALISATION --------------------------
#
# This Shiny application allows interactive exploration of smoothing
# parameters and visual inspection of smoothed profiles prior to
# reconstruction.
#
# Uncomment the block below to launch the application.
#
# required_pkgs <- c("shiny", "bslib", "shinycssloaders")
# invisible(lapply(required_pkgs, library, character.only = TRUE))
#
# source("src/my_app.R")
# shinyApp(ui, server)


# ---- INITIALISE RECONSTRUCTION PARAMETERS -----------------------------------

# Depth threshold (m) below which a profile is considered truncated
start_limit_complete_profile <- 1.5

# Depth threshold (m) above which a profile is considered truncated
end_limit_complete_profile <- 49.5

# Number of B-spline basis functions (K)
bspline_smooth_nbasis <- 100

# Smoothing parameter (lambda) for B-spline smoothing
bspline_smooth_lambda <- 1e-3

# Number of principal components used for reconstruction
n_pc <- 30

# Variance parameter (sigma^2) controlling the projection penalty
# - small values: closer fit to observed data
# - large values: smoother reconstructed profiles
proj_pen_sigma <- .1

# Depth grid used for prediction
start <- 1
end <- 50
step <- .2


# ---- RUN PROFILE RECONSTRUCTION ---------------------------------------------

# Execute the full reconstruction workflow on sparse profiles
source("src/run_pace_on_sparse_profiles.R")


# ---- VISUALISATION AND OUTPUT ------------------------------------------------

# Plot the FPCA score space (projected scores)
plot_scores_projected(fpca, scr)

# Select one incomplete profile to visualise
index <- incomplete_indexes[5]

# Plot observed and reconstructed profile with confidence intervals
plot_one_profile(
  index,
  sparse_data,
  reconstructed_profiles,
  var_label = my_var,
  CI = TRUE
) + theme(legend.position = "none")

# Reconstructed profiles with pointwise confidence intervals (CI)
# and confidence bands (CB)
head(reconstructed_profiles)

# Complete profiles and reconstructed incomplete profiles
head(all_data)
