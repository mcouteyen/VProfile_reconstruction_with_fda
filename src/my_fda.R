library(fda)


sparse.mat <- function(datalist){
  #
  #Create a matrix of sparse data with NAs out of a list of sparse data
  #
  # Arguments:
  #
  # DATALIST .... A list object. Each element of the list is a matrix with ncol > 1. 
  #               The first column of each element corresponds to the point index per 
  #               observation. The following columns are the observations per variable.
  # 
  #Last Modified : 07 March 2025
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

mpca.fd <- function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj),
                   centerfns = TRUE)
{
  #  Carry out a functional PCA with regularization
  #  Arguments:
  #  FDOBJ      ... Functional data object
  #  NHARM     ... Number of principal components or harmonics to be kept
  #  HARMFDPAR ... Functional parameter object for the harmonics
  #  CENTERFNS ... If TRUE, the mean function is first subtracted from each 
  #                function.
  #
  #  Returns:  An object PCAFD of class "pca.fd" with these named entries:
  #  harmonics  ... A functional data object for the harmonics or eigenfunctions
  #  values     ... The complete set of eigenvalues
  #  scores     ... A matrix or array of scores on the principal components or 
  #                 harmonics
  #  varprop    ... A vector giving the proportion of variance explained
  #                 by each eigenfunction
  #  meanfd     ... A functional data object giving the mean function
  #
  
  #  Modified from pca.fd
  #  Last modified:  07 March 2025 
  
  #  Check FDOBJ
  
  if (!(is.fd(fdobj) || is.fdPar(fdobj)))  stop(
    "First argument is neither a functional data or a functional parameter object.")
  if (is.fdPar(fdobj)) fdobj <- fdobj$fd
  
  #  compute mean function and center if required
  
  meanfd <- mean.fd(fdobj)
  if (centerfns) {
    fdobj <- center.fd(fdobj)
  }
  
  #  get coefficient matrix and its dimensions
  
  coef  <- fdobj$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  nrep  <- coefd[2]
  coefnames <- dimnames(coef)
  if (nrep < 2) stop("PCA not possible without replications.")
  
  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  type     <- basisobj$type
  
  #  set up HARMBASIS
  
  harmbasis <- harmfdPar$fd$basis
  nhbasis   <- harmbasis$nbasis
  
  #  set up LFDOBJ and LAMBDA
  
  Lfdobj <- harmfdPar$Lfd
  lambda <- harmfdPar$lambda
  
  #  compute CTEMP whose cross product is needed
  
  if (ndim == 3) {
    nvar <- coefd[3]
    ctemp <- matrix(0, nvar * nbasis, nrep)
    for(j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      ctemp[index,  ] <- coef[,  , j]
    }
  } else {
    nvar  <- 1
    ctemp <- coef
  }
  X <- t(ctemp)
  
  #  set up cross product Lmat for harmonic basis,
  #  roughness penalty matrix Rmat, and
  #  penalized cross product matrix Lmat
  
  Wmat <- eval.penalty(harmbasis, 0)
  Wmat <- (Wmat + t(Wmat))/2
  
  #  compute the Choleski factor Mmat of Lmat
  
  Whalfmat    <- chol(Wmat)
  Whalfmatinv <- solve(Whalfmat)
  
  #  set up cross product and penalty matrices
  
  Vmat <- crossprod(X)/nrep
  
  #  set up matrix for eigenanalysis
  
  if(nvar == 1) {
    Cmat <- Whalfmat %*% Vmat %*% t(Whalfmat)
  } else {
    Cmat <- matrix(0,nvar*nhbasis,nvar*nhbasis)
    M <- matrix(0,nvar*nhbasis,nvar*nhbasis)
    W <- matrix(0,nvar*nhbasis,nvar*nhbasis)
    for (i in 1:nvar) {
      indexi <- 1:nbasis + (i - 1) * nbasis
      sigma.var <- sum(diag(Vmat[indexi,indexi] %*% Wmat))
      Mmat <- diag(1/sigma.var, nrow = length(indexi), ncol = length(indexi))
      M[indexi, indexi] <- Mmat
      W[indexi, indexi] <- Wmat
    }
    Mchol <- chol(M)
    Wchol <- chol(W)
    Cmat <- Mchol %*% Wchol %*% Vmat %*% t(Wchol) %*% Mchol
  }
  
  #  eigenalysis
  
  Cmat    <- (Cmat + t(Cmat))/2
  result  <- eigen(Cmat)
  eigvalc <- result$values
  eigvecc <- as.matrix(result$vectors[, 1:nharm])
  sumvecc <- apply(eigvecc, 2, sum)
  eigvecc[,sumvecc < 0] <-  - eigvecc[, sumvecc < 0]
  
  varprop <- eigvalc[1:nharm]/sum(eigvalc)
  
  #  set up harmfd
  
  if (nvar == 1) {
    harmcoef <- Whalfmatinv %*% eigvecc
  } else {
    harmcoef <- array(0, c(nbasis, nharm, nvar))
    beta <- solve(Mchol) %*% solve(Wchol) %*% eigvecc
    for (j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      harmcoef[,  , j] <- beta[index,]
    }
  }
  harmnames <- rep("", nharm)
  for(i in 1:nharm)
    harmnames[i] <- paste("PC", i, sep = "")
  if(length(coefd) == 2)
    harmnames <- list(coefnames[[1]], harmnames,"values")
  if(length(coefd) == 3)
    harmnames <- list(coefnames[[1]], harmnames, coefnames[[3]])
  harmfd   <- fd(harmcoef, harmbasis, harmnames)
  
  #  set up harmscr
  
  if (nvar == 1) {
    harmscr  <- inprod(fdobj, harmfd)
  } else {
    harmscr <- X %*% t(Wchol) %*% Mchol %*% eigvecc
  }
  
  #  set up the object pcafd of the pca.fd class containing the results
  
  pcafd        <- list(harmfd, eigvalc, harmscr, varprop, meanfd)
  class(pcafd) <- "pca.fd"
  names(pcafd) <- c("harmonics", "values", "scores", "varprop", "meanfd")
  
  return(pcafd)
}

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

cond.var.fd <- function (fdobj1, fdobj2 = fdobj1, mean1, mean2 = mean1, weights) {
  if (!(is.fd(fdobj1) || is.fdPar(fdobj1))) 
    stop("First argument is neither a functional data or a functional parameter object.")
  if (is.fdPar(fdobj1)) 
    fdobj1 <- fdobj1$fd
  if (!(is.fd(fdobj2) || is.fdPar(fdobj2))) 
    stop("Second argument is neither a functional data or a functional parameter object.")
  if (is.fdPar(fdobj2)) 
    fdobj2 <- fdobj2$fd
  fdobj1 <- wcenter.fd(fdobj1, mean1)
  fdobj2 <- wcenter.fd(fdobj2, mean2)
  coefx <- fdobj1$coefs
  coefy <- fdobj2$coefs
  coefdobj1 <- dim(coefx)
  coefdobj2 <- dim(coefy)
  basisx <- fdobj1$basis
  basisy <- fdobj2$basis
  nbasisx <- basisx$nbasis
  nbasisy <- basisy$nbasis
  Wmat <- diag(weights)
  if (coefdobj1[2] != coefdobj2[2]) 
    stop("Number of replications are not equal.")
  if (length(coefdobj1) == 2) {
    if (length(coefdobj2) == 2) {
      X <- Wmat %*% t(coefx)
      coefvar <- crossprod(X, t(coefy))
      coefnames <- list(dimnames(coefx)[[1]], dimnames(coefy)[[1]])
      varbifd <- bifd(coefvar, basisx, basisy, coefnames)
    }
    else stop("Both arguments must be univariate.")
  }
  else {
    nvar <- coefdobj1[3]
    npair <- nvar * (nvar + 1)/2
    coefvar <- array(0, c(nbasisx, nbasisx, 1, npair))
    varnames <- fdobj1$fdnames[[3]]
    m <- 0
    bivarnames <- vector("list", npair)
    for (i in 1:nvar) for (j in 1:i) {
      m <- m + 1
      X <- Wmat %*% t(coefx[, , i])
      coefvar[, , 1, m] <- crossprod(X, t(coefx[, , j]))
      bivarnames[m] <- paste(varnames[i], "vs", varnames[j])
    }
    bifdnames = fdobj1$fdnames
    bifdnames[[3]] <- bivarnames
    varbifd <- bifd(coefvar, basisx, basisx, bifdnames)
  }
  return(varbifd)
}

cond.mpca.fd <- function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj),
                         centerfns = TRUE, x, positions, h)
{
  #  Carry out a functional PCA with regularization
  #  Arguments:
  #  FDOBJ      ... Functional data object
  #  NHARM     ... Number of principal components or harmonics to be kept
  #  HARMFDPAR ... Functional parameter object for the harmonics
  #  CENTERFNS ... If TRUE, the mean function is first subtracted from each 
  #                function.
  #
  #  Returns:  An object PCAFD of class "pca.fd" with these named entries:
  #  harmonics  ... A functional data object for the harmonics or eigenfunctions
  #  values     ... The complete set of eigenvalues
  #  scores     ... A matrix or array of scores on the principal components or 
  #                 harmonics
  #  varprop    ... A vector giving the proportion of variance explained
  #                 by each eigenfunction
  #  meanfd     ... A functional data object giving the mean function
  #
  
  #  Modified from pca.fd
  #  Last modified:  07 March 2025 
  
  #  Check FDOBJ
  
  if (!(is.fd(fdobj) || is.fdPar(fdobj)))  stop(
    "First argument is neither a functional data or a functional parameter object.")
  if (is.fdPar(fdobj)) fdobj <- fdobj$fd
  
  # computes weights 
  weights <- compute_weights(x, positions, h)
  meanfd <- wmean.fd(fdobj, weights)
  if (centerfns) {
    fdobj <- wcenter.fd(fdobj, meanfd)
  }
  
  # #  compute mean function and center if required
  # 
  # meanfd <- mean.fd(fdobj)
  # if (centerfns) {
  #   fdobj <- center.fd(fdobj)
  # }
  
  #  get coefficient matrix and its dimensions
  
  coef  <- fdobj$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  nrep  <- coefd[2]
  coefnames <- dimnames(coef)
  if (nrep < 2) stop("PCA not possible without replications.")
  
  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  type     <- basisobj$type
  
  #  set up HARMBASIS
  
  harmbasis <- harmfdPar$fd$basis
  nhbasis   <- harmbasis$nbasis
  
  #  set up LFDOBJ and LAMBDA
  
  Lfdobj <- harmfdPar$Lfd
  lambda <- harmfdPar$lambda
  
  #  compute CTEMP whose cross product is needed
  
  if (ndim == 3) {
    nvar <- coefd[3]
    ctemp <- matrix(0, nvar * nbasis, nrep)
    for(j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      ctemp[index,  ] <- coef[,  , j]
    }
  } else {
    nvar  <- 1
    ctemp <- coef
  }
  X <- t(ctemp)
  
  #  set up cross product Lmat for harmonic basis,
  #  roughness penalty matrix Rmat, and
  #  penalized cross product matrix Lmat
  
  Wmat <- eval.penalty(harmbasis, 0)
  Wmat <- (Wmat + t(Wmat))/2
  
  #  compute the Choleski factor Mmat of Lmat
  
  Whalfmat    <- chol(Wmat)
  Whalfmatinv <- solve(Whalfmat)
  
  #  set up cross product and penalty matrices
  Omat <- diag(weights)
  Y <- Omat %*% X
  Vmat <- crossprod(Y, X)
  # Vmat <- crossprod(X)/nrep
  
  #  set up matrix for eigenanalysis
  
  if(nvar == 1) {
    Cmat <- Whalfmat %*% Vmat %*% t(Whalfmat)
  } else {
    Cmat <- matrix(0,nvar*nhbasis,nvar*nhbasis)
    M <- matrix(0,nvar*nhbasis,nvar*nhbasis)
    W <- matrix(0,nvar*nhbasis,nvar*nhbasis)
    for (i in 1:nvar) {
      indexi <- 1:nbasis + (i - 1) * nbasis
      sigma.var <- sum(diag(Vmat[indexi,indexi] %*% Wmat))
      Mmat <- diag(1/sigma.var, nrow = length(indexi), ncol = length(indexi))
      M[indexi, indexi] <- Mmat
      W[indexi, indexi] <- Wmat
    }
    Mchol <- chol(M)
    Wchol <- chol(W)
    Cmat <- Mchol %*% Wchol %*% Vmat %*% t(Wchol) %*% Mchol
  }
  
  #  eigenalysis
  
  Cmat    <- (Cmat + t(Cmat))/2
  result  <- eigen(Cmat)
  eigvalc <- result$values
  eigvecc <- as.matrix(result$vectors[, 1:nharm])
  sumvecc <- apply(eigvecc, 2, sum)
  eigvecc[,sumvecc < 0] <-  - eigvecc[, sumvecc < 0]
  
  varprop <- eigvalc[1:nharm]/sum(eigvalc)
  
  #  set up harmfd
  
  if (nvar == 1) {
    harmcoef <- Whalfmatinv %*% eigvecc
  } else {
    harmcoef <- array(0, c(nbasis, nharm, nvar))
    beta <- solve(Mchol) %*% solve(Wchol) %*% eigvecc
    for (j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      harmcoef[,  , j] <- beta[index,]
    }
  }
  harmnames <- rep("", nharm)
  for(i in 1:nharm)
    harmnames[i] <- paste("PC", i, sep = "")
  if(length(coefd) == 2)
    harmnames <- list(coefnames[[1]], harmnames,"values")
  if(length(coefd) == 3)
    harmnames <- list(coefnames[[1]], harmnames, coefnames[[3]])
  harmfd   <- fd(harmcoef, harmbasis, harmnames)
  
  #  set up harmscr
  
  if (nvar == 1) {
    harmscr  <- inprod(fdobj, harmfd)
  } else {
    harmscr <- X %*% t(Wchol) %*% Mchol %*% eigvecc
  }
  
  #  set up the object pcafd of the pca.fd class containing the results
  
  pcafd        <- list(harmfd, eigvalc, harmscr, varprop, meanfd)
  class(pcafd) <- "pca.fd"
  names(pcafd) <- c("harmonics", "values", "scores", "varprop", "meanfd")
  
  return(pcafd)
}

wmean.fd <- function(fdobj, weights){
  coefs <- fdobj$coefs
  coefs_dim <- dim(coefs)
  if(length(coefs_dim) == 2){
    mean_coefs <- apply(coefs, 1, function(x) weights*x)
    mean_coefs <- apply(mean_coefs, 2, sum)
  } else {
    nvar <- coefs_dim[3]
    mean_coefs <- array(0, c(coefs_dim[1], 1, nvar))
    for (j in 1:nvar) {
      coefs_weighted <- apply(coefs[, ,j], 1, function(x) weights*x)
      mean_coefs[, 1, j] <- apply(coefs_weighted, 2, sum)
    }
  }
  basis <- fdobj$basis
  fdnames <- fdobj$fdnames
  fdnames[[2]] <- "mean"
  fdnames[[3]] <- paste("mean", fdnames[[3]])
  meanfd <- fd(mean_coefs, basis, fdnames)
  return(meanfd)
}

wcenter.fd <- function(fdobj, wmean){
  coefs <- fdobj$coefs
  coefs_dim <- dim(coefs)
  if(length(coefs_dim) == 2){
    fdobj$coefs <- apply(coefs, 2, function(x) x-wmean$coefs)
  } else {
    nvar <- coefs_dim[3]
    for (j in 1:nvar) {
      fdobj$coefs[, , j] <- apply(coefs[,,j], 2, function(x) x-wmean$coefs[,,j])
    }
  }
  return(fdobj)
}


