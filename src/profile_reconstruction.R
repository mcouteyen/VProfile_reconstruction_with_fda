#' Check presence of required variables in a data frame
#'
#' Verifies that a data frame contains the mandatory columns required for
#' vertical profile reconstruction: INDEX, DEPTH, and the
#' variable of interest.
#'
#' @param data A data.frame containing profile observations.
#' @param my_var A character string giving the name of the variable to analyse.
#'
#' @return A character vector of missing column names. Returns an empty vector
#' if all required columns are present.
#'
#' @examples
check_variables_name <- function(data, my_var){
  required_cols <- c("INDEX", "DEPTH", my_var)
  
  return(setdiff(required_cols, names(data)))
}

#' Detect duplicated depth values within profiles
#'
#' Identifies duplicated combinations of INDEX and DEPTH, which
#' may affect functional smoothing and FPCA.
#'
#' @param data A data.frame containing at least INDEX and DEPTH.
#'
#' @return A data.frame listing duplicated DEPTH–INDEX pairs.
check_duplicates <- function(data){
  grouped_data <- data %>% 
    group_by(INDEX, DEPTH) %>% 
    summarise(count = n(), .groups = "keep") %>%
    as.data.frame()
  grouped_data <- grouped_data %>% filter(count > 1)
  return(grouped_data[, c("DEPTH", "INDEX")])
}

#' Identify rows containing missing values
#'
#' Returns the row indices for which at least one of the specified variables
#' contains missing values.
#'
#' @param data A data.frame.
#' @param var_label A character vector of column names to check.
#'
#' @return An integer vector of row indices with missing values.
all_na_index <- function(data, var_label){
  res <- c()
  for(var in var_label){
    indexes <- which(is.na(data[, var]))
    res <- c(res, indexes)
  }
  res <- unique(res)
  return(res)
}

#' Remove rows with missing values
#'
#' Removes rows containing missing values for the specified variables.
#'
#' @param data A data.frame.
#' @param var_label A character vector of column names.
#'
#' @return A data.frame with rows containing missing values removed.
remove_na_rows <- function(data, var_label){
  to_remove <- c()
  for(var in var_label)
    to_remove <- c(to_remove, which(is.na(data[, var_label])))
  if(length(to_remove > 0)){
    data <- data[-to_remove,]
  }
  return(data)
}

#' Retrieve complete vertical profiles
#'
#' Selects profiles that fully span a given argument range and contain
#' no missing values for the variables of interest.
#'
#' @param data A data.frame containing profile observations.
#' @param var_label Name of the variable of interest.
#' @param min_arg Minimum required argument value.
#' @param max_arg Maximum required argument value.
#' @param index Optional maximum profile index to retain.
#'
#' @return A data.frame containing only complete profiles.
retrieve_complete_profiles <- function(data, var_label, min_arg, max_arg, index = NA){
  profile_length <- data %>% group_by(INDEX) %>% summarise(MIN = min(DEPTH), MAX = max(DEPTH))
  nas <- all_na_index(data, var_label)
  if(length(nas > 0)){
    data <- data[-nas,]
  }
  if(is.na(index)){
    to_remove <- profile_length %>% filter(MIN > min_arg | MAX < max_arg)
  } else {
    to_remove <- profile_length %>% filter(MIN > min_arg | MAX < max_arg | INDEX > index)
  }
  
  # remove profiles that are too short
  complete_profiles <- data[-which(data$INDEX %in% to_remove$INDEX),]
  return(complete_profiles)
}

#' Create a functional data object from profile observations
#'
#' Converts discrete vertical profiles into a functional data object
#' using B-spline smoothing.
#'
#' @param data A data.frame containing profile observations.
#' @param x Name of the argument variable (e.g. depth).
#' @param all_y Name(s) of the variable(s) to smooth.
#' @param lambda Smoothing parameter.
#' @param K Number of B-spline basis functions.
#' @param min_prof Optional minimum argument value.
#' @param max_prof Optional maximum argument value.
#'
#' @return An object of class fd}.
create_fdobj <- function(data, x, all_y, lambda, K, min_prof = NA, max_prof = NA){
  indexes <- unique(data$INDEX)
  #create basis
  arg_prof <- unique(data[,x])
  if(is.na(min_prof)){
    min_prof <- min(arg_prof)
  }
  if(is.na(max_prof)){
    max_prof <- max(arg_prof)
  }
  basis <- create.bspline.basis(c(min_prof, max_prof), nbasis=K, norder=4)
  if(length(all_y) == 1){
    y <- all_y
    data <- data[!is.na(data[,y]),]
    coefs <- matrix(0, nrow = basis$nbasis, ncol = length(indexes))
    for(index in 1:length(indexes)){
      date_data <- data %>% filter(INDEX == indexes[index])
      date_data <- date_data %>% arrange(.data[[x]])
      date_argvals <- date_data[,x]
      date_y <- date_data[,y]
      date_smooth <- smooth.basisPar(argvals = date_argvals, y = date_y, fdobj = basis, lambda = lambda)
      coefs[,index] <- date_smooth$fd$coefs
    }
  } else {
    #compute coefs for each base function
    coefs <- array(NA, dim = c(basis$nbasis, length(indexes), 2))
    for(index in 1:length(indexes)){
      date_data <- data %>% filter(INDEX == indexes[index])
      for(var_index in 1:length(all_y)){
        y <- all_y[var_index]
        date_argvals <- date_data[!is.na(date_data[,y]),x]
        date_y <- date_data[!is.na(date_data[,y]),y]
        if(length(date_y) > 0){
          date_smooth <- smooth.basisPar(argvals = date_argvals, y = date_y, fdobj = basis, lambda = lambda)
          coefs[,index, var_index] <- date_smooth$fd$coefs
        }
      }
    }
  }
  res <- fd(coefs, basis)
  return(res)
}

#' Compute asymptotic variance for FPCA-based reconstruction
#'
#' Computes pointwise asymptotic variance used to derive confidence
#' intervals and confidence bands.
#'
#' @param fpca A FPCA object returned by pca.fd}.
#' @param n_pc Number of principal components used.
#' @param covestimate Functional covariance estimate.
#' @param proj_pen Projection penalty parameter.
#' @param step Resolution of the evaluation grid.
#'
#' @return A numeric vector of pointwise variances.
compute_asymptotic_sd <- function(fpca, n_pc, covestimate, proj_pen = .1, step= .1){
  range <- fpca$meanfd$basis$rangeval
  argvals <- seq(range[1], range[2], step)
  
  Covmat <- eval.bifd(argvals, argvals, covestimate)
  Sigma <- Covmat + diag(proj_pen, nrow = length(argvals))
  Phi <- t(eval.fd(argvals, fpca$harmonics))
  Lambda <- diag(fpca$values[1:n_pc])
  H <- Lambda %*% Phi
  Omega <- Lambda - H %*% solve(Sigma) %*% t(H)
  sd <- sapply(1:ncol(Phi), function(x) t(Phi[,x]) %*% Omega %*% Phi[,x])
  return(sd)
}

#' Reconstruct vertical profiles from FPCA scores
#'
#' Reconstructs full vertical profiles from FPCA scores, with optional
#' confidence intervals and confidence bands.
#'
#' @param fpca A FPCA object.
#' @param scores Matrix of FPCA scores.
#' @param profile_id Vector of profile identifiers.
#' @param n_pc Number of principal components used.
#' @param var Name of the reconstructed variable.
#' @param step Resolution of the output grid.
#' @param args_data Optional original argument values.
#' @param CI Logical; compute confidence intervals if TRUE.
#' @param sd Asymptotic variance.
#' @param alpha Significance level for confidence interval.
#'
#' @return A data.frame of reconstructed profiles.
reconstruct_profiles <- function(fpca, scores, profile_id, n_pc, var, step = .1, args_data = NULL, CI = FALSE, sd = NA, alpha = .05){
  n_profiles <- nrow(scores)
  range <- fpca$meanfd$basis$rangeval
  
  res <- data.frame()
  for(index in 1:n_profiles){
    if(is.null(args_data)){
      args <- seq(range[1], range[2], step)
    } else {
      args <- args_data[args_data$INDEX == profile_id[index], "DEPTH"]
    }
    profile <- fd(basisobj=fpca$meanfd$basis)
    profile$coefs <- fpca$meanfd$coefs
    for(pc in 1:n_pc){
      harmonic <- as.numeric(fpca$harmonics$coefs[,pc])
      profile$coefs <- profile$coefs + scores[index,pc] * harmonic
    }
    profile_df <- data.frame(INDEX = profile_id[index],
                             X = args,
                             VAR = eval.fd(args, profile)[,1])
    if(CI){
      profile_df <- profile_df %>% mutate(CI_up = VAR + qnorm(1-alpha/2)*sqrt(sd))
      profile_df <- profile_df %>% mutate(CI_low = VAR - qnorm(1-alpha/2)*sqrt(sd))
      profile_df <- profile_df %>% mutate(CB_up = VAR + sqrt(qchisq(1-alpha/2, df = n_pc)*sd))
      profile_df <- profile_df %>% mutate(CB_low = VAR - sqrt(qchisq(1-alpha/2, df = n_pc)*sd))
    }
    res <- rbind(res, profile_df)
  }
  colnames(res)[3] <- var
  colnames(res)[2] <- "DEPTH"
  return(res)
}

#' Reconstruct functional coefficients from FPCA scores
#'
#' @param fpca A FPCA object.
#' @param scores Matrix of FPCA scores.
#' @param n_pc Number of principal components.
#'
#' @return A matrix of reconstructed coefficients.
reconstruct_profiles_coefs <- function(fpca, scores, n_pc){
  n_profiles <- nrow(scores)
  basis <- fpca$meanfd$basis
  range <- basis$rangeval
  
  res <- array(dim = c(basis$nbasis, n_profiles))
  for(index in 1:n_profiles){
    res[, index] <- fpca$meanfd$coefs
    for(pc in 1:n_pc){
      harmonic <- as.numeric(fpca$harmonics$coefs[, pc])
      res[, index] <- res[, index] + scores[index,pc] * harmonic
    }
  }
  return(res)
}

#' Append reconstructed profiles to a functional data object
#'
#' @param fdobj A functional data object.
#' @param coefs Matrix or array of coefficients to append.
#'
#' @return A functional data object with appended profiles.
add_reps_to_fdobj <- function(fdobj, coefs){
  basis <- fdobj$basis
  dim_1 <- dim(fdobj$coefs)
  dim_2 <- dim(coefs)
  if(length(dim_1) > 2){
    new_coefs <- array(dim = c(dim_1[1], dim_1[2] + dim_2[2], dim_1[3]))
    for(var in 1:dim_1[3]){
      new_coefs[,,var] <- cbind(fdobj$coefs[,,var], coefs[,,var])
    }
  } else {
    new_coefs <- cbind(fdobj$coefs, coefs)
  }
  res <- fd(new_coefs, basis)
  return(res)
}

#' Project incomplete profiles onto FPCA space
#'
#' Projects sparse or truncated profiles using a PACE-based approach.
#'
#' @param incomplete_data A data.frame of incomplete profiles.
#' @param my_var Name of the variable to reconstruct.
#' @param fpca A FPCA object.
#' @param cov_fda Functional covariance estimate.
#' @param meanfd Functional mean.
#' @param sigma Measurement error variance.
#'
#' @return A matrix of FPCA scores.
project_incomplete_profiles <- function(incomplete_data, my_var, fpca, cov_fda, meanfd, sigma = NA){
  incomplete_indexes <- unique(incomplete_data$INDEX)
  n_args <- incomplete_data %>% group_by(INDEX) %>% summarise(count= n())
  max_args <- max(n_args$count)
  my_list <- list()
  time <- sort(unique(incomplete_data$DEPTH))
  for(index in 1:length(incomplete_indexes)){
    profile <- incomplete_data %>% filter(INDEX == incomplete_indexes[index])
    my_list[[index]] <- as.matrix(profile[, c("DEPTH", my_var)])
  }
  covestimate <- list(cov.estimate = cov_fda, meanfd = meanfd)
  scr <- scoresPACE_corrected(my_list, time, covestimate, fpca, sigma)
  return(scr)
}

#' Reconstruct complete profiles using FPCA
#'
#' Performs FPCA on complete profiles and reconstructs them, optionally
#' computing confidence intervals.
#'
#' @param complete_profiles A data.frame of complete profiles.
#' @param var_names Name(s) of the variable(s).
#' @param fdobj Functional data object.
#' @param n_pc Number of principal components.
#' @param CI Logical; compute confidence intervals if TRUE.
#' @param proj_pen Projection penalty.
#'
#' @return A data.frame of reconstructed profiles.
reconstruct_complete_profile_with_fpca <- function(complete_profiles, var_names, fdobj, n_pc, CI = FALSE, proj_pen = NA){
  complete_indexes <- unique(complete_profiles$INDEX)
  fpca <- pca.fd(fdobj, n_pc)
  cov_fda <- var.fd(fdobj)
  meanfd <- mean.fd(fdobj)
  scr <- fpca$scores
  if(CI){
    sd <- compute_asymptotic_sd(fpca, n_pc, cov_fda, proj_pen)
    res <- reconstruct_profiles(fpca, scr, complete_indexes, n_pc = n_pc, var = var_names, CI = TRUE, sd = sd)
  } else {
    res <- reconstruct_profiles(fpca, scr, complete_indexes, n_pc = n_pc, var = var_names)
  }
  return(res)
}
