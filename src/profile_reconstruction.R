
all_na_index <- function(data, var_label){
  res <- c()
  for(var in var_label){
    indexes <- which(is.na(data[, var]))
    res <- c(res, indexes)
  }
  res <- unique(res)
  return(res)
}


remove_na_rows <- function(data, var_label){
  to_remove <- c()
  for(var in var_label)
    to_remove <- c(to_remove, which(is.na(data[, var_label])))
  if(length(to_remove > 0)){
    data <- data[-to_remove,]
  }
  return(data)
}


retrieve_complete_profiles <- function(data, arg_label = "DEPTH", var_label, min_arg, max_arg, index = NA){
  profile_length <- data %>% group_by(INDEX) %>% summarise(MIN = min(.data[[arg_label]]), MAX = max(.data[[arg_label]]))
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

plot_cov_fun <- function(cov_fda, start, end, step = .1){
  new_args <- seq(start, end, step)
  diag_cov <- evaldiag.bifd(new_args, cov_fda)
  res <- plot(new_args,diag_cov,type="l",
              xlab = "DEPTH", ylab="",
              cex.axis=1.25,cex.lab=1.25,cex.main=1.25,lwd=2)
  return(res)
}

plot_cor_fun <- function(fdobj, start, end, step = .1){
  new_args <- seq(start, end, step)
  cor_fda <- cor.fd(new_args, fdobj)
  require(fields)
  res <- image.plot(new_args,new_args,cor_fda,
                    xlab="DEPTH", ylab = "DEPTH",
                    cex.axis=1.25,cex.lab=1.25,cex.main=1.25,
                    legend.shrink=0.75,legend.line=-1.5)
  return(res)
}

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

reconstruct_profiles <- function(fpca, scores, profile_id, n_pc, var, arg = "DEPTH", step = .1, args_data = NULL, CI = FALSE, sd = NA, alpha = .05){
  n_profiles <- nrow(scores)
  range <- fpca$meanfd$basis$rangeval
  
  res <- data.frame()
  for(index in 1:n_profiles){
    if(is.null(args_data)){
      args <- seq(range[1], range[2], step)
    } else {
      args <- args_data[args_data$INDEX == profile_id[index], arg]
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
  colnames(res)[2] <- arg
  return(res)
}


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


plot_projected_scores <- function(fpca, new_scores){
  plot(fpca$scores[, 1], fpca$scores[, 2])
  points(new_scores[,1], new_scores[, 2], col = "red")
  plot(fpca$scores[, 3], fpca$scores[, 4])
  points(new_scores[,3], new_scores[, 4], col = "red")
}


retrieve_quantiles <- function(MSE_profiles){
  quant <- quantile(MSE_profiles$MSE, probs = c(.25, .50, .75))
  
  MSE_profiles <- MSE_profiles %>% mutate(diff_Q1 = abs(MSE - quant[1]),
                                          diff_Q2 = abs(MSE - quant[2]),
                                          diff_Q3 = abs(MSE - quant[3]))
  max <- MSE_profiles[which(MSE_profiles$MSE == max(MSE_profiles$MSE)),]$INDEX[1]
  min <- MSE_profiles[which(MSE_profiles$MSE == min(MSE_profiles$MSE)),]$INDEX[1]
  q1 <- MSE_profiles[which(MSE_profiles$diff_Q1 == min(MSE_profiles$diff_Q1)),]$INDEX[1]
  q2 <- MSE_profiles[which(MSE_profiles$diff_Q2 == min(MSE_profiles$diff_Q2)),]$INDEX[1]
  q3 <- MSE_profiles[which(MSE_profiles$diff_Q3 == min(MSE_profiles$diff_Q3)),]$INDEX[1]
  res <- data.frame(QUANTILE = c("MIN", "Q1", "MEDIAN", "Q3", "MAX"),
                    INDEX = c(min, q1, q2, q3, max))
  return(res)
}

retrieve_quantile_indexes <- function(predicted_data, observed_data, var){
  model_stat <- eval_pred(predicted_data, observed_data, var = var)
  MSE_profiles <- model_stat$MSE_profiles
  res <- retrieve_quantiles(MSE_profiles)
  return(res)
}

plot_profiles <- function(index_df, data_used, data_predicted, data_sup = NULL, var_label, CI = FALSE){
  indexes <- index_df$INDEX
  min_val <- min(data_used[, var_label])
  max_val <- max(data_used[, var_label])
  data_used <- data_used %>% filter(INDEX %in% indexes)
  data_predicted <- data_predicted %>% filter(INDEX %in% indexes)
  data_predicted$INDEX <- factor(data_predicted$INDEX,
                                 levels = indexes)
  facet_labels <- index_df$QUANTILE
  names(facet_labels) <- index_df$INDEX
  res <- ggplot() + 
    geom_point(data = data_used, aes(.data[[var_label]], DEPTH),
               alpha = .8,
               size = 2.2) +
    geom_path(data = data_predicted, 
              aes(.data[[var_label]], DEPTH), 
              lwd = .9,
              alpha = 1,
              col = "#E69F00") +
    
    scale_x_continuous() +
    scale_y_reverse() +
    ylab("DEPTH") +
    # xlab(var_label) +
    facet_wrap(.~INDEX,
               nrow = 1,
               labeller = labeller(INDEX = facet_labels),
               scales = "free_x") +
    theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 11),
          strip.background = element_rect(fill = "transparent"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
  if(!is.null(data_sup)){
    data_sup <- data_sup %>% filter(INDEX %in% indexes)
    res <- res + 
      geom_path(data = data_sup, aes(.data[[var_label]], DEPTH), 
                 col = my_palette[1],
                 alpha = .8,
                lwd = .9)
  }
  if(CI){
    var_col_number <- which(colnames(data_predicted) == var_label)
    BI_up <- colnames(data_predicted)[var_col_number + 1]
    BI_low <- colnames(data_predicted)[var_col_number + 2]
    res <- res + 
      geom_path(data = data_predicted, 
                aes(.data[[BI_up]], DEPTH), 
                alpha = .8,
                lty = "dotted",
                lwd = .9) +
      geom_path(data = data_predicted, 
                aes(.data[[BI_low]], DEPTH), 
                alpha = .8,
                lty = "dotted",
                lwd = .9)
  }
  
  return(res)
}


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

eval_pred <- function(predicted, observed, var, argvals = "DEPTH", subj = "INDEX"){
  observed <- observed[, c(subj, argvals, var)]
  predicted <- predicted[, c(subj, argvals, var)]
  res_data <- inner_join(observed, predicted, by = c(argvals, subj))
  colnames(res_data)[3:4] <- c("y", "yhat")
  res_data <- res_data %>% mutate(residuals = y-yhat)
  profile_stat <- res_data %>% group_by(.data[[subj]]) %>% summarise(MSE = 1/n()*sum(residuals^2))
  profile_stat <- profile_stat %>% arrange(desc(MSE))
  MSE <- 1/nrow(res_data)*sum(res_data$residuals^2)
  return(list(res_data = res_data, MSE_profiles = profile_stat, MSE = MSE))
}

eval_CI <- function(predicted, observed, var){
  merged_df <- inner_join(predicted, observed, by = c("INDEX", "DEPTH"))
  var_name <- paste0(var, ".y")
  pointwise <- merged_df %>% mutate(IN_CI = ifelse(.data[[var_name]] < CI_up & .data[[var_name]] > CI_low, TRUE, FALSE),
                                    IN_CB = ifelse(.data[[var_name]] < CB_up & .data[[var_name]] > CB_low, TRUE, FALSE))
  simultaneous_band <- pointwise %>% group_by(INDEX) %>% summarise(IN_CB = all(IN_CB))
  prop_CI <- sum(pointwise$IN_CI)/nrow(pointwise)
  prop_CB <- sum(simultaneous_band$IN_CB)/nrow(simultaneous_band)
  return(list(CI = prop_CI, CB = prop_CB))
}




