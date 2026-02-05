
missing_cols <- check_variables_name(sparse_data, my_var)

if (length(missing_cols) > 0) {
  stop(
    sprintf(
      "MISSING COLUMN(S) : %s",
      paste(missing_cols, collapse = ", ")
    ),
    call. = FALSE
  )
}

duplicated_rows <- check_duplicates(sparse_data)

if(nrow(duplicated_rows) > 0){
  stop("DUPLICATED ROWS, check: duplicated_rows", call. = FALSE)
}

sparse_data <- remove_na_rows(sparse_data, my_var)

sparse_data <- sparse_data %>% arrange(INDEX, DEPTH)

complete_data <- retrieve_complete_profiles(sparse_data, 
                                            var_label = my_var, 
                                            min_arg = start_limit_complete_profile, 
                                            max_arg = end_limit_complete_profile)
incomplete_data <- setdiff(sparse_data, complete_data)

complete_indexes <- unique(complete_data$INDEX)
incomplete_indexes <- unique(incomplete_data$INDEX)


fdobj <- create_fdobj(complete_data, "DEPTH", my_var, lambda = bspline_smooth_lambda, K = bspline_smooth_nbasis)

fpca <- pca.fd(fdobj, n_pc)
cov_fda <- var.fd(fdobj)
meanfd <- mean.fd(fdobj)
scr <- project_incomplete_profiles(incomplete_data, my_var, fpca, cov_fda, meanfd, proj_pen_sigma)
scores_to_remove <- which(is.na(scr[,1]))

complete_scores <- cbind(data.frame(INDEX = unique(complete_data$INDEX)), 
                         as.data.frame(fpca$scores[,1:n_pc]))
incomplete_scores <- cbind(data.frame(INDEX = unique(incomplete_data$INDEX)), 
                           as.data.frame(scr[, 1:n_pc]))

if(length(scores_to_remove) > 0){
  incomplete_scores <- incomplete_scores[-scores_to_remove,]
  scr <- scr[-scores_to_remove,]
  incomplete_dates <- incomplete_dates[-scores_to_remove,]
}

reconstructed_coefs <- reconstruct_profiles_coefs(fpca, scr, n_pc)

new_fdobj <- add_reps_to_fdobj(fdobj, reconstructed_coefs)

sd <- compute_asymptotic_sd(fpca, n_pc, cov_fda, proj_pen_sigma)

reconstructed_profiles <- reconstruct_profiles(fpca, scr, incomplete_indexes, n_pc, my_var, CI = TRUE, sd = sd)

all_data <- rbind(complete_data, reconstructed_profiles[, c("INDEX", "DEPTH", my_var)])
