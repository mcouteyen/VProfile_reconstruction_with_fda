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
    scale_x_continuous() +
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

