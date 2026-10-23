#' Construct a sparse data matrix or array from a list of observations
#'
#' This function builds a sparse matrix (or array for multivariate data)
#' filled with \code{NA} values from a list of irregularly sampled profiles.
#'
#' @details
#' This function is adapted from internal data-handling utilities used in the
#' \pkg{fda} package. It has been modified to better accommodate sparse and
#' multivariate profile data in the context of vertical oceanographic profiles.
#'
#' The original concepts and data structures are attributed to Ramsay et al.
#' (see the \pkg{fda} package documentation). This implementation represents
#' a user-side adaptation and extension.
#'
#' @param datalist A list. Each element is a matrix with at least two columns.
#'   The first column corresponds to the sampling index (e.g. depth or time),
#'   and the remaining columns correspond to observed variables.
#'
#' @return
#' If the data are univariate, a matrix of dimension
#' \code{length(unique(index)) x number_of_profiles}.
#' If multivariate, an array of dimension
#' \code{length(unique(index)) x number_of_profiles x number_of_variables}.
#'
#' @references
#' Ramsay, J. O., Wickham, H., Graves, S., & Hooker, G. (2023).
#' \emph{fda: Functional Data Analysis}. R package.
#'
#' @author
#' Adapted and extended by M. Couteyen Carpaye, 2025.
#'
#' @seealso \code{\link[fda]{smooth.basis}}, \code{\link[fda]{pca.fd}}
#'
#' @keywords internal
#'
#' @export
sparse.mat <- function(datalist){
  time = sort(unique(unlist(lapply(datalist, function(x) x[,1]))))
  data = matrix(NA, nrow = length(time), ncol = length(datalist))
  
  nvar = ncol(datalist[[1]]) - 1
  
  if(nvar == 1){
    data = matrix(NA, nrow = length(time), ncol = length(datalist))
    for(i in 1:length(datalist)){
      data[which((time %in% datalist[[i]][,1])),i] = datalist[[i]][,2]
    }
  }else{
    data = array(NA, dim = c(length(time),length(datalist),nvar))
    for(j in 1:nvar){
      for(i in 1:length(datalist)){
        data[which((time %in% datalist[[i]][,1])),i,j] = datalist[[i]][,(1 + j)]
      }
    }
  }
  return(data)
}

#' Compute PACE scores for sparse functional data (corrected version)
#'
#' This function computes Functional Principal Component scores for sparsely
#' observed profiles using a PACE-like approach.
#'
#' @details
#' This function is adapted from the internal logic of PACE score computation
#' routines used in the \pkg{fda} package. The original methodology is described
#' in Yao et al. (2005).
#'
#' Modifications were introduced to:
#' \itemize{
#'   \item handle sparse vertical profiles,
#'   \item explicitly account for truncated profiles,
#'   \item improve numerical stability,
#'   \item support multivariate functional data.
#' }
#'
#' The author of this script claims authorship only for these adaptations and
#' extensions; the underlying theoretical framework and original implementation
#' remain the intellectual property of the \pkg{fda} authors.
#'
#' @param data A list of matrices containing sparse observations.
#' @param time A numeric vector of observation points (e.g. depth).
#' @param covestimate A covariance estimate object produced by \pkg{fda}.
#' @param PC A functional PCA object from \pkg{fda}.
#' @param varest Optional noise variance. If \code{NA}, it is estimated from data.
#'
#' @return
#' A matrix of reconstructed FPCA scores for each profile.
#'
#' @references
#' Ramsay, J. O., & Silverman, B. W. (2005).
#' \emph{Functional Data Analysis}. Springer.
#' 
#' Yao, F., Mueller, H.G., Wang, J.L. (2005).
#'\emph{Functional data analysis for sparse longitudinal data}. J.American Statistical Association, 100, 577-590.
#'
#' Ramsay, J. O., Wickham, H., Graves, S., & Hooker, G. (2023).
#' \emph{fda: Functional Data Analysis}. R package.
#'
#' @author
#' Adapted and extended by M. Couteyen Carpaye, 2025.
#'
#' @seealso \code{\link[fda]{pca.fd}}, \code{\link[fda]{eval.fd}}
#'
#' @export
scoresPACE_corrected <- function(data, time, covestimate, PC, varest = NA){
  npc <- length(PC$varprop)
  coef <- PC$values[1:npc]
  nvar <- length(PC$meanfd$fdnames$funs)
  rangevals <- PC$meanfd$basis$rangeval
  tau <- rangevals[2] - rangevals[1]
  
  datamat <- sparse.mat(data)
  if(dim(datamat)[2] == 1){
    datalist <- data
  } else {
    datalist <- sparse.list(datamat, time)
  }
  indexes <- lapply(datalist, function(x) which(time %in% x[,1])) #sampling points for each subject
  phi <- lapply(indexes, function(x) eval.fd(time[x],PC$harmonics))
  
  if(nvar == 1){
    mean.point <- matrix(eval.fd(time, covestimate$meanfd),nrow = nrow(eval.fd(time, covestimate$meanfd)),ncol = nvar)
    data.var.mat <- (apply(datamat,2, function(x) x - mean.point))^2 #data given as matrix
    mu <- lapply(indexes,function(x) as.matrix(eval.fd(time[x],covestimate$meanfd),nrow = length(x),ncol = nvar))
    cent <- lapply(1:length(datalist), function(x) datalist[[x]][,-1]-mu[[x]])
    
    covdiag <- diag(eval.bifd(time,time,covestimate$cov.estimate))
    varest.mat <- as.matrix(apply(data.var.mat,2, function(x) x-covdiag))
    time_mat <- matrix(rep(time, ncol(varest.mat)), nrow = length(time), ncol = ncol(varest.mat))
    if(is.na(varest)){
      all_varest <- c()
      for(col in 1:ncol(varest.mat)){
        var_col <- varest.mat[, col]
        val <- var_col[!is.na(var_col)]
        time_col <- time_mat[!is.na(var_col), col]
        time_col <- time_col[which(time_col > (rangevals[1]+tau/4) & time_col < (rangevals[2]-tau/4))]
        val <- val[which(time_col > (rangevals[1]+tau/4) & time_col < (rangevals[2]-tau/4))]
        to_sum <- (val[-1] + val[-length(val)]) / 2 * diff(time_col)
        varest <- 2/tau * sum(to_sum)
        all_varest <- c(all_varest, varest)
      }
      
      varest <- mean(all_varest)
    } 
    diagon <- lapply(1:length(datalist), function(x) diag(varest, nrow = length(indexes[[x]]), ncol = length(indexes[[x]])))
    sigmai <- Map('+', lapply(indexes, function(x) eval.bifd(time[x],time[x],covestimate$cov.estimate)), diagon)
    
    scr <- t(do.call(cbind,lapply(1:length(datalist), function(x) coef*(t(phi[[x]])%*%solve(sigmai[[x]])%*%cent[[x]]) )))
    
  } else {
    
    phi <- lapply(phi, function(x) rbind(x[,,1], x[,,2]))
    
    mean.point <- matrix(eval.fd(time, covestimate$meanfd),nrow = length(time),ncol = nvar)
    mu <- lapply(indexes,function(x) matrix(eval.fd(time[x],covestimate$meanfd),nrow = length(x),ncol = nvar))
    cent <- lapply(1:length(datalist), function(x) datalist[[x]][,-1]-mu[[x]])
    cent <- lapply(cent, function(x) rbind(matrix(x[,1]), matrix(x[,2])))
    
    if(is.na(varest)){
      varest <- rep(0.1, nvar)
    } else if(length(varest) == 1){
      varest <- rep(varest, nvar)
    }
    
    # calculer covestimate à part 
    diagon = list()
    sigmai = list()
    cov_index <- c(1, 3)
    cross_cov_index <- c(2, 2)
    sigmai <- lapply(indexes, function(x) matrix(0, nrow = length(x)*nvar, ncol = length(x)*nvar))
    for(i in 1:nvar){
      for(j in 1:nvar){
        if(i == j){
          var.cov.estimate <- bifd(coef = covestimate$cov.estimate$coefs[,,,cov_index[i]], covestimate$cov.estimate$sbasis, covestimate$cov.estimate$tbasis, covestimate$cov.estimate$fdnames)
        } else if(i < j){
          var.cov.estimate <- bifd(coef = t(covestimate$cov.estimate$coefs[,,,cross_cov_index[i]]), covestimate$cov.estimate$sbasis, covestimate$cov.estimate$tbasis, covestimate$cov.estimate$fdnames)
        } else {
          var.cov.estimate <- bifd(coef = covestimate$cov.estimate$coefs[,,,cross_cov_index[i]], covestimate$cov.estimate$sbasis, covestimate$cov.estimate$tbasis, covestimate$cov.estimate$fdnames)
        }
        for(profile in 1:length(sigmai)){
          x <- indexes[[profile]]
          bloc_indexi <- 1:length(x) + (i - 1) * length(x)
          bloc_indexj <- 1:length(x) + (j - 1) * length(x)
          res <- eval.bifd(time[x],time[x],var.cov.estimate)
          if(i == j){
            diagon <- diag(varest[i], nrow = length(x), ncol = length(x))
            res <- res + diagon
          }
          sigmai[[profile]][bloc_indexi, bloc_indexj] <- res
        }
      }
    }
    
    scr <- t(do.call(cbind,lapply(1:length(datalist), function(x) coef*(t(phi[[x]])%*%solve(sigmai[[x]])%*%cent[[x]]) )))
    
  }
  
  return(scr)
  
}



