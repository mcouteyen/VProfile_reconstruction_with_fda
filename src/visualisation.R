#' Plot FPCA scores and projected incomplete profiles
#'
#' Produces a two-dimensional FPCA score plot (PC1 vs PC2), showing both
#' complete profiles used to compute the FPCA and incomplete profiles
#' projected onto the FPCA space.
#'
#' @param fpca A FPCA object returned by pca.fd.
#' @param scr A data.frame or matrix of FPCA scores for projected
#'   incomplete profiles.
#'
#' @return A ggplot object.
#'
#' @details
#' The percentage of variance explained by the first two principal components
#' is displayed on the axes labels.
plot_scores_projected <- function(fpca, scr){
  fpca_scores <- fpca$scores
  colnames <- c()
  for(n in 1:n_pc){
    colnames <- c(colnames, paste0("PC", n))
  }
  colnames(fpca_scores) <- colnames
  xlab <- paste0("PC1 ", round(fpca$varprop[1]*100, 2), "%")
  ylab <- paste0("PC2 ", round(fpca$varprop[2]*100, 2), "%")
  res <- ggplot() + 
    geom_point(data = fpca_scores, aes(PC1, PC2, col = "complete profiles"), 
               size = 2.7) +
    geom_point(data = scr, aes(PC1, PC2, col = "incomplete projected profiles"), 
               alpha = .5,
               size = 2.5) +
    xlab(xlab) +
    scale_color_manual(name = "", values = c("black", "grey")) + 
    ylab(ylab) +
    ggtitle(my_var) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "top")
  return(res)
}

#' Plot a reconstructed vertical profile
#'
#' Displays observed data points and reconstructed profiles for a given
#' profile index, with optional additional profiles and confidence intervals.
#'
#' @param index Profile identifier to plot.
#' @param data_used A data.frame containing the original observations.
#' @param data_predicted A data.frame containing reconstructed profiles.
#' @param data_sup Optional data.frame of supplementary profiles.
#' @param var_label Name of the variable to plot.
#' @param CI Logical; draw confidence intervals if TRUE.
#'
#' @return A ggplot object.
plot_one_profile <- function(index, data_used, data_predicted, data_sup = NULL, var_label, CI = FALSE){
  data_used <- data_used %>% filter(INDEX == index)
  data_predicted <- data_predicted %>% filter(INDEX == index)
  
  res <- ggplot() + 
    geom_point(data = data_used, aes(.data[[var_label]], DEPTH),
               alpha = .9,
               size = 2.2) +
    geom_path(data = data_predicted, 
              aes(.data[[var_label]], DEPTH), 
              color = 'blue',
              lwd = .8,
              alpha = .8) +
    scale_x_continuous(expand = c(0.07, 0.07)) +
    scale_y_reverse() + 
    ylab("DEPTH") +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10))
  if(!is.null(data_sup)){
    data_sup <- data_sup %>% filter(INDEX == index)
    res <- res + 
      geom_path(data = data_sup, aes(.data[[var_label]], DEPTH), 
                col = 'red',
                alpha = .4,
                lwd = .8)
  }
  if(CI){
    var_col_number <- which(colnames(data_predicted) == var_label)
    BI_up <- colnames(data_predicted)[var_col_number + 1]
    BI_low <- colnames(data_predicted)[var_col_number + 2]
    res <- res + 
      geom_path(data = data_predicted, 
                aes(.data[[BI_up]], DEPTH), 
                alpha = .8,
                lwd = .6,
                lty = "dotted") +
      geom_path(data = data_predicted, 
                aes(.data[[BI_low]], DEPTH), 
                alpha = .8,
                lwd = .6,
                lty = "dotted")
  }
  
  return(res)
}

#' Plot the diagonal of a functional covariance surface
#'
#' Plots the pointwise variance (diagonal of the covariance function)
#' as a function of the argument (e.g. depth).
#'
#' @param cov_fda A functional covariance object.
#' @param start Minimum argument value.
#' @param end Maximum argument value.
#' @param step Resolution of the evaluation grid.
#'
#' @return The output of the base plot function.
plot_cov_fun <- function(cov_fda, start, end, step = .1){
  new_args <- seq(start, end, step)
  diag_cov <- evaldiag.bifd(new_args, cov_fda)
  res <- plot(new_args,diag_cov,type="l",
              xlab = "DEPTH", ylab="",
              cex.axis=1.25,cex.lab=1.25,cex.main=1.25,lwd=2)
  return(res)
}

#' Plot the functional correlation surface
#'
#' Displays the correlation surface of a functional data object as
#' a two-dimensional image.
#'
#' @param fdobj A functional data object.
#' @param start Minimum argument value.
#' @param end Maximum argument value.
#' @param step Resolution of the evaluation grid.
#'
#' @return A correlation surface plot.
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
