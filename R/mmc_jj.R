
# Mixed Model Coexpression (MMC)
# (c) 2011-2020  Nick Furlotte
#
# This file may be used for your personal use.
# This file may be modified but the modified version must retain this copyright notice.
# This file or modified version must not be redistributed without prior permission from the author.
# This software is provided “AS IS”, without warranty of any kind.
# In no event shall the author be liable for any claim, damages or other liability
# in connection with the software or it's use.

# takes an expression matrix with genes on the rows and samples on the columns
# returns the covariance matrix for covariance between the columns
mmc.ncov <- function(exprs,useV='pairwise.complete'){
  exprs.norm = mmc.normMatrix(exprs)
  return(cov(exprs.norm,use=useV))
}

# takes an expression matrix with genes on the rows and samples on the columns
# produces a normalized matrix, where normalization is accomplished by subtracting rowmeans and dividing by row stdev for each row
mmc.normMatrix <- function(exprs){

  exprs.madj = exprs - matrix(rowMeans(exprs, na.rm=TRUE),nrow(exprs),ncol(exprs))
  exprs.sd = matrix(apply(exprs,1,sd,na.rm=TRUE),nrow(exprs),ncol(exprs))
  exprs.norm = exprs.madj / exprs.sd

  return(exprs.norm)
}


# ys - a matrix of gene expressions where the genes are rows and individuals columns
# K - the matrix generated with mmc.ncov - a matrix representing the correlation between samples
# The other parameters are for the variance component search and can generally be left alone
# faster is a boolean that tells the function to use the fast method or the slow method
#        in theory the absolute value from these two approaches will be the same.  However, they might be slighly different for small sample sizes.
#        The benefit of the slower method is that you will get signed coexpression values, while the fast method produces only the absolute values of the (coexpression)

mmc.cor <- function(ys, K, ngrids=100, llim=-10, ulim=10, esp=1e-10,faster=TRUE){

  X0 <- matrix(1,ncol(ys),1)

  N = nrow(ys)
  q0 <- ncol(X0)
  q1 <- q0 + 1

  dfs <- matrix(NA,nrow=N,ncol=N)
  stats <- matrix(NA,nrow=N,ncol=N)
  rs <- matrix(NA,nrow=N,ncol=N)
  vgs <- matrix(NA,nrow=N,ncol=N)
  ves <- matrix(NA,nrow=N,ncol=N)
  bs <- matrix(NA,nrow=N,ncol=N)
  ps <- matrix(NA,nrow=N,ncol=N)

  for(i in 1:N) {
    for(j in 1:N) {
      if(i == j){ rs[i,j] <- 1; next; }
      vids <- intersect(which(!is.na(ys[i,])),which(!is.na(ys[j,])))
      nv <- length(vids)

      X <- cbind(X0[vids,,drop=FALSE],ys[i,vids])
      if (det(crossprod(X,X)) == 0){ next; }
      eig.L <- emma.eigen.L.wo.Z(K[vids,vids])
      eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X)
      REMLE <- emma.REMLE(ys[j,vids],X,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
      U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),nv,nv,byrow=TRUE)
      if(sum(is.nan(U)) > 0){ next; }
      dfs[i,j] <- nv - q1

      yt <- crossprod(U,ys[j,vids])
      Xt <- crossprod(U,X)
      Xtt  = crossprod(Xt,Xt)

      if(rcond(Xtt) < 1e-7){
        next
      }

      iXX <- solve(Xtt)
      beta <- iXX%*%crossprod(Xt,yt)

      vgs[i,j] <- REMLE$vg
      ves[i,j] <- REMLE$ve
      bs[i,j] <- beta[q1]

      if(faster){
        stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
        rs[i,j] <- sqrt(stats[i,j]^2 / (stats[i,j]^2 + dfs[i,j]))
      }
      else{
        sigmaHat = vgs[i,j]*K[vids,vids] + ves[i,j]*diag(nv)
        sigmaHatInv = solve(sigmaHat)
        sigmaHatInvChol = chol(sigmaHatInv)

        y1adj = matrix(ys[j,vids] - mean(ys[j,vids]),ncol=1)
        y2adj = matrix(ys[i,vids] - mean(ys[i,vids]),ncol=1)

        y1adj.chol = sigmaHatInvChol %*% y1adj
        y2adj.chol = sigmaHatInvChol %*% y2adj

        var1 = t(y1adj.chol) %*% y1adj.chol
        var2 = t(y2adj.chol) %*% y2adj.chol

        rs[i,j] <- (t(y1adj.chol) %*% y2adj.chol) / sqrt(var1 * var2)
      }
    }
  }
  return (list(rs=rs,stats=stats,dfs=dfs,vgs=vgs,ves=ves,betas=bs))
}


# ys - a matrix of gene expressions where the genes are rows and individuals columns
# K - the matrix generated with mmc.ncov - a matrix representing the correlation between samples
mmc.corsym <- function(ys, K, ngrids=100, llim=-10, ulim=10, esp=1e-10){

  # Currently there are issues with using other covariates so I am only regressing out the mean
  # I'm looking into this issue
  #mmc.corsym <- function(ys, K, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10){
  #if ( is.null(X0) ) {
  #  X0 <- matrix(1,ncol(ys),1)
  #}

  X0 <- matrix(1,ncol(ys),1)
  X01 = X0
  X02 = X0

  N = nrow(ys)

  q01 <- ncol(X01)
  q11 <- q01 + 1
  q02 <- ncol(X02)
  q12 <- q02 + 1

  rs <- matrix(NA,nrow=N,ncol=N)

  for(i in 1:N){
    rs[i,i] <- 1
    if(i == N){ break }
    for(j in (i+1):N){
      vids <- intersect(which(!is.na(ys[i,])),which(!is.na(ys[j,])))
      nv <- length(vids)

      # Direction 1
      X1 <- cbind(X01[vids,,drop=FALSE],ys[i,vids])
      if (det(crossprod(X1,X1)) == 0){ next; }
      eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X1)
      REMLE <- emma.REMLE(ys[j,vids],X1,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)

      eig.L <- emma.eigen.L.wo.Z(K[vids,vids])
      U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),nv,nv,byrow=TRUE)
      if(sum(is.nan(U)) > 0){ next; }
      yt <- crossprod(U,ys[j,vids])
      Xt <- crossprod(U,X1)
      iXX <- solve(crossprod(Xt,Xt))
      beta1 <- iXX%*%crossprod(Xt,yt)

      vgs1 <- REMLE$vg
      ves1 <- REMLE$ve

      # Direction 2
      X2 <- cbind(X02[vids,,drop=FALSE],ys[j,vids])
      if (det(crossprod(X2,X2)) == 0){ next; }
      eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X2)
      REMLE <- emma.REMLE(ys[i,vids],X2,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)

      eig.L <- emma.eigen.L.wo.Z(K[vids,vids])
      U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),nv,nv,byrow=TRUE)
      if(sum(is.nan(U)) > 0){ next; }
      yt <- crossprod(U,ys[i,vids])
      Xt <- crossprod(U,X2)
      iXX <- solve(crossprod(Xt,Xt))
      beta2 <- iXX%*%crossprod(Xt,yt)

      vgs2 <- REMLE$vg
      ves2 <- REMLE$ve

      #cat(vgs1,",",vgs2,"\n")
      #cat(ves1,",",ves2,"\n\n")

      sigmaHat1 = vgs1*K[vids,vids] + ves1*diag(nv)
      sigmaHat2 = vgs2*K[vids,vids] + ves2*diag(nv)
      sigmaHat1Inv = solve(sigmaHat1)
      sigmaHat2Inv = solve(sigmaHat2)
      sigmaHat1InvChol = chol(sigmaHat1Inv)
      sigmaHat2InvChol = chol(sigmaHat2Inv)

      y1adj = matrix(ys[j,vids] - mean(ys[j,vids]),ncol=1)
      y2adj = matrix(ys[i,vids] - mean(ys[i,vids]),ncol=1)
      #y1adj = matrix(matrix(ys[j,vids],ncol=1) - (X1 %*% matrix(beta1,ncol=1)),ncol=1)
      #y2adj = matrix(matrix(ys[i,vids],ncol=1) - (X2 %*% matrix(beta2,ncol=1)),ncol=1)

      y1adj.chol = sigmaHat1InvChol %*% y1adj
      y2adj.chol = sigmaHat2InvChol %*% y2adj

      var1 = t(y1adj.chol) %*% y1adj.chol
      var2 = t(y2adj.chol) %*% y2adj.chol

      rs[i,j] <- (t(y1adj.chol) %*% y2adj.chol) / sqrt(var1 * var2)
      rs[j,i] <- rs[i,j]

    }
  }
  return(rs)
}


# cor <- function(x, use = "p"){
#   y=NULL;ngrids=100;llim=-10;ulim=10;esp=1e-10
#   x<-t(x)
#   K = mmc.ncov(x)
#   # Currently there are issues with using other covariates so I am only regressing out the mean
#   # I'm looking into this issue
#   #mmc.corsym <- function(x, K, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10){
#   #if ( is.null(X0) ) {
#   #  X0 <- matrix(1,ncol(x),1)
#   #}
#
#   X0 <- matrix(1,ncol(x),1)
#   X01 = X0
#   X02 = X0
#
#   N = nrow(x)
#
#   q01 <- ncol(X01)
#   q11 <- q01 + 1
#   q02 <- ncol(X02)
#   q12 <- q02 + 1
#
#   rs <- matrix(NA,nrow=N,ncol=N)
#
#   for(i in 1:N){
#     tt <- as.character(Sys.time())
#     cat(i,"/",nrow(x),": ",tt,"\n",sep = "")
#     rs[i,i] <- 1
#     if(i == N){ break }
#     for(j in (i+1):N){
#       # cat(j,"\n",sep = "")
#
#       vids <- intersect(which(!is.na(x[i,])),which(!is.na(x[j,])))
#       nv <- length(vids)
#
#       # Direction 1
#       X1 <- cbind(X01[vids,,drop=FALSE],x[i,vids])
#       if (det(crossprod(X1,X1)) == 0){ next; }
#       eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X1)
#       REMLE <- emma.REMLE(x[j,vids],X1,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
#
#       # eig.L <- emma.eigen.L.wo.Z(K[vids,vids])
#       # U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),nv,nv,byrow=TRUE)
#       # if(sum(is.nan(U)) > 0){ next; }
#       # yt <- crossprod(U,x[j,vids])
#       # Xt <- crossprod(U,X1)
#       # iXX <- solve(crossprod(Xt,Xt))
#       # beta1 <- iXX%*%crossprod(Xt,yt)
#
#       vgs1 <- REMLE$vg
#       ves1 <- REMLE$ve
#       # e <- lmmlite::eigen_rotation(K, x[j,vids], X1 )
#       # out <- lmmlite::fitLMM(e$Kva, e$y, e$X, tol = 1e-10)
#       # vgs1 <- out$sigmasq_g
#       # ves1 <- out$sigmasq_e
#
#       # Direction 2
#       X2 <- cbind(X02[vids,,drop=FALSE],x[j,vids])
#       if (det(crossprod(X2,X2)) == 0){ next; }
#       eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X2)
#       REMLE <- emma.REMLE(x[i,vids],X2,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
#
#       # eig.L <- emma.eigen.L.wo.Z(K[vids,vids])
#       # U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),nv,nv,byrow=TRUE)
#       # if(sum(is.nan(U)) > 0){ next; }
#       # yt <- crossprod(U,x[i,vids])
#       # Xt <- crossprod(U,X2)
#       # iXX <- solve(crossprod(Xt,Xt))
#       # beta2 <- iXX%*%crossprod(Xt,yt)
#
#       vgs2 <- REMLE$vg
#       ves2 <- REMLE$ve
#       # e <- lmmlite::eigen_rotation(K, x[i,vids], X2 )
#       # out <- lmmlite::fitLMM(e$Kva, e$y, e$X, tol = 1e-10)
#       # vgs1 <- out$sigmasq_g
#       # ves1 <- out$sigmasq_e
#
#
#
#       #cat(vgs1,",",vgs2,"\n")
#       #cat(ves1,",",ves2,"\n\n")
#
#       sigmaHat1 = vgs1*K[vids,vids] + ves1*diag(nv)
#       sigmaHat2 = vgs2*K[vids,vids] + ves2*diag(nv)
#       sigmaHat1Inv = solve(sigmaHat1)
#       sigmaHat2Inv = solve(sigmaHat2)
#       sigmaHat1InvChol = chol(sigmaHat1Inv)
#       sigmaHat2InvChol = chol(sigmaHat2Inv)
#
#       y1adj = matrix(x[j,vids] - mean(x[j,vids]),ncol=1)
#       y2adj = matrix(x[i,vids] - mean(x[i,vids]),ncol=1)
#       #y1adj = matrix(matrix(x[j,vids],ncol=1) - (X1 %*% matrix(beta1,ncol=1)),ncol=1)
#       #y2adj = matrix(matrix(x[i,vids],ncol=1) - (X2 %*% matrix(beta2,ncol=1)),ncol=1)
#
#       y1adj.chol = sigmaHat1InvChol %*% y1adj
#       y2adj.chol = sigmaHat2InvChol %*% y2adj
#
#       var1 = t(y1adj.chol) %*% y1adj.chol
#       var2 = t(y2adj.chol) %*% y2adj.chol
#
#       rs[i,j] <- (t(y1adj.chol) %*% y2adj.chol) / sqrt(var1 * var2)
#       rs[j,i] <- rs[i,j]
#
#     }
#   }
#   return(rs)
# }

cor <- function(x, use = "p"){
  cat("MMC cor\n")
  y=NULL;ngrids=100;llim=-10;ulim=10;esp=1e-10
  x<-t(x)
  K = mmc.ncov(x)
  # Currently there are issues with using other covariates so I am only regressing out the mean
  # I'm looking into this issue
  #mmc.corsym <- function(x, K, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10){
  #if ( is.null(X0) ) {
  #  X0 <- matrix(1,ncol(x),1)
  #}

  X0 <- matrix(1,ncol(x),1)
  X01 = X0
  X02 = X0

  N = nrow(x)

  q01 <- ncol(X01)
  q11 <- q01 + 1
  q02 <- ncol(X02)
  q12 <- q02 + 1

  rs <- matrix(NA,nrow=N,ncol=N)
  for(i in 1:N){
    tt <- as.character(Sys.time())
    cat(i,"/",nrow(x),": ",tt,"\n",sep = "")
    rs[i,i] <- 1

    if(i == N){ break }

    res <- parallel::mclapply((i+1):N, function(j){
      vids <- intersect(which(!is.na(x[i,])),which(!is.na(x[j,])))
      nv <- length(vids)

      # Direction 1
      X1 <- cbind(X01[vids,,drop=FALSE],x[i,vids])
      if (det(crossprod(X1,X1)) == 0){ next; }

      eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X1)
      REMLE <- emma.REMLE(x[j,vids],X1,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
      vgs1 <- REMLE$vg
      ves1 <- REMLE$ve

      # e <- lmmlite::eigen_rotation(K, x[j,vids], X1 )
      # out <- lmmlite::fitLMM(e$Kva, e$y, e$X, tol = 1e-10)
      # vgs1 <- out$sigmasq_g
      # ves1 <- out$sigmasq_e

      # Direction 2
      X2 <- cbind(X02[vids,,drop=FALSE],x[j,vids])
      if (det(crossprod(X2,X2)) == 0){ next; }
      eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X2)
      REMLE <- emma.REMLE(x[i,vids],X2,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
      vgs2 <- REMLE$vg
      ves2 <- REMLE$ve


      #cat(vgs1,",",vgs2,"\n")
      #cat(ves1,",",ves2,"\n\n")

      sigmaHat1 = vgs1*K[vids,vids] + ves1*diag(nv)
      sigmaHat2 = vgs2*K[vids,vids] + ves2*diag(nv)
      sigmaHat1Inv = solve(sigmaHat1)
      sigmaHat2Inv = solve(sigmaHat2)
      sigmaHat1InvChol = chol(sigmaHat1Inv)
      sigmaHat2InvChol = chol(sigmaHat2Inv)

      y1adj = matrix(x[j,vids] - mean(x[j,vids]),ncol=1)
      y2adj = matrix(x[i,vids] - mean(x[i,vids]),ncol=1)
      #y1adj = matrix(matrix(x[j,vids],ncol=1) - (X1 %*% matrix(beta1,ncol=1)),ncol=1)
      #y2adj = matrix(matrix(x[i,vids],ncol=1) - (X2 %*% matrix(beta2,ncol=1)),ncol=1)

      y1adj.chol = sigmaHat1InvChol %*% y1adj
      y2adj.chol = sigmaHat2InvChol %*% y2adj

      var1 = t(y1adj.chol) %*% y1adj.chol
      var2 = t(y2adj.chol) %*% y2adj.chol

      (t(y1adj.chol) %*% y2adj.chol) / sqrt(var1 * var2)
    }, mc.cores = 30) %>% unlist

    rs[(i+1):N,i] <- res
    rs[i,(i+1):N] <- res


    # for(j in (i+1):N){
    #
    #   vids <- intersect(which(!is.na(x[i,])),which(!is.na(x[j,])))
    #   nv <- length(vids)
    #
    #   # Direction 1
    #   X1 <- cbind(X01[vids,,drop=FALSE],x[i,vids])
    #   if (det(crossprod(X1,X1)) == 0){ next; }
    #
    #   eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X1)
    #   REMLE <- emma.REMLE(x[j,vids],X1,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
    #   vgs1 <- REMLE$vg
    #   ves1 <- REMLE$ve
    #
    #   # e <- lmmlite::eigen_rotation(K, x[j,vids], X1 )
    #   # out <- lmmlite::fitLMM(e$Kva, e$y, e$X, tol = 1e-10)
    #   # vgs1 <- out$sigmasq_g
    #   # ves1 <- out$sigmasq_e
    #
    #   # Direction 2
    #   X2 <- cbind(X02[vids,,drop=FALSE],x[j,vids])
    #   if (det(crossprod(X2,X2)) == 0){ next; }
    #   eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X2)
    #   REMLE <- emma.REMLE(x[i,vids],X2,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
    #   vgs2 <- REMLE$vg
    #   ves2 <- REMLE$ve
    #
    #
    #   #cat(vgs1,",",vgs2,"\n")
    #   #cat(ves1,",",ves2,"\n\n")
    #
    #   sigmaHat1 = vgs1*K[vids,vids] + ves1*diag(nv)
    #   sigmaHat2 = vgs2*K[vids,vids] + ves2*diag(nv)
    #   sigmaHat1Inv = solve(sigmaHat1)
    #   sigmaHat2Inv = solve(sigmaHat2)
    #   sigmaHat1InvChol = chol(sigmaHat1Inv)
    #   sigmaHat2InvChol = chol(sigmaHat2Inv)
    #
    #   y1adj = matrix(x[j,vids] - mean(x[j,vids]),ncol=1)
    #   y2adj = matrix(x[i,vids] - mean(x[i,vids]),ncol=1)
    #   #y1adj = matrix(matrix(x[j,vids],ncol=1) - (X1 %*% matrix(beta1,ncol=1)),ncol=1)
    #   #y2adj = matrix(matrix(x[i,vids],ncol=1) - (X2 %*% matrix(beta2,ncol=1)),ncol=1)
    #
    #   y1adj.chol = sigmaHat1InvChol %*% y1adj
    #   y2adj.chol = sigmaHat2InvChol %*% y2adj
    #
    #   var1 = t(y1adj.chol) %*% y1adj.chol
    #   var2 = t(y2adj.chol) %*% y2adj.chol
    #
    #   rs[i,j] <- (t(y1adj.chol) %*% y2adj.chol) / sqrt(var1 * var2)
    #   rs[j,i] <- rs[i,j]
    #   cat(j,":",rs[i,j],"\n",sep = "")
    # }

  }
  return(rs)
}



cor_lmmite <- function(x, use = "p"){
  y=NULL;ngrids=100;llim=-10;ulim=10;esp=1e-10
  x<-t(x)
  K = mmc.ncov(x)
  # Currently there are issues with using other covariates so I am only regressing out the mean
  # I'm looking into this issue
  #mmc.corsym <- function(x, K, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10){
  #if ( is.null(X0) ) {
  #  X0 <- matrix(1,ncol(x),1)
  #}

  X0 <- matrix(1,ncol(x),1)
  X01 = X0
  X02 = X0

  N = nrow(x)

  q01 <- ncol(X01)
  q11 <- q01 + 1
  q02 <- ncol(X02)
  q12 <- q02 + 1

  rs <- matrix(NA,nrow=N,ncol=N)

  for(i in 1:N){
    tt <- as.character(Sys.time())
    cat(i,"/",nrow(x),": ",tt,"\n",sep = "")
    rs[i,i] <- 1
    if(i == N){ break }
    for(j in (i+1):N){
      # cat(j,"\n",sep = "")

      vids <- intersect(which(!is.na(x[i,])),which(!is.na(x[j,])))
      nv <- length(vids)

      # Direction 1
      X1 <- cbind(X01[vids,,drop=FALSE],x[i,vids])
      if (det(crossprod(X1,X1)) == 0){ next; }
      # eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X1)
      # REMLE <- emma.REMLE(x[j,vids],X1,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
      # vgs1 <- REMLE$vg
      # ves1 <- REMLE$ve
      #
      e <- lmmlite::eigen_rotation(K, x[j,vids], X1 )
      out <- lmmlite::fitLMM(e$Kva, e$y, e$X, tol = 1e-10)
      vgs1 <- out$sigmasq_g
      ves1 <- out$sigmasq_e

      # Direction 2
      X2 <- cbind(X02[vids,,drop=FALSE],x[j,vids])
      if (det(crossprod(X2,X2)) == 0){ next; }
      # eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X2)
      # REMLE <- emma.REMLE(x[i,vids],X2,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.R1)
      # vgs2 <- REMLE$vg
      # ves2 <- REMLE$ve
      #
      e <- lmmlite::eigen_rotation(K, x[i,vids], X2 )
      out <- lmmlite::fitLMM(e$Kva, e$y, e$X, tol = 1e-10)
      vgs2 <- out$sigmasq_g
      ves2 <- out$sigmasq_e

      #cat(vgs1,",",vgs2,"\n")
      #cat(ves1,",",ves2,"\n\n")

      sigmaHat1 = vgs1*K[vids,vids] + ves1*diag(nv)
      sigmaHat2 = vgs2*K[vids,vids] + ves2*diag(nv)
      sigmaHat1Inv = solve(sigmaHat1)
      sigmaHat2Inv = solve(sigmaHat2)
      sigmaHat1InvChol = chol(sigmaHat1Inv)
      sigmaHat2InvChol = chol(sigmaHat2Inv)

      y1adj = matrix(x[j,vids] - mean(x[j,vids]),ncol=1)
      y2adj = matrix(x[i,vids] - mean(x[i,vids]),ncol=1)
      #y1adj = matrix(matrix(x[j,vids],ncol=1) - (X1 %*% matrix(beta1,ncol=1)),ncol=1)
      #y2adj = matrix(matrix(x[i,vids],ncol=1) - (X2 %*% matrix(beta2,ncol=1)),ncol=1)

      y1adj.chol = sigmaHat1InvChol %*% y1adj
      y2adj.chol = sigmaHat2InvChol %*% y2adj

      var1 = t(y1adj.chol) %*% y1adj.chol
      var2 = t(y2adj.chol) %*% y2adj.chol

      rs[i,j] <- (t(y1adj.chol) %*% y2adj.chol) / sqrt(var1 * var2)
      rs[j,i] <- rs[i,j]

    }
  }
  return(rs)
}

# The following set of functions are taken directly from code that is freely available at http://mouse.cs.ucla.edu/emma/install.html .
# This code (any function beginning with 'emma') is under any copyright or other licensing agreement stipulated at the above site.

emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
  t <- t1
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- lambda+delta
  return( 0.5*(nq*(sum(etasq/(ldelta*ldelta))+etas.2.sq/(delta*delta))/(sum(etasq/ldelta)+etas.2.sq/delta)-(sum(1/ldelta)+(n-t)/delta)) )
}

emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

  #  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if( is.nan(dLL[1]) ) { return (list(REML=0,delta=0,ve=0,vg=0)) }
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }

    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Lambdas)+(n-t)/delta))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }

    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
      }
    }
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/(n-q)
  }
  maxve <- maxva*maxdelta

  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.kinship <- function(snps, method="additive", use="all") {
  n0 <- sum(snps==0,na.rm=TRUE)
  nh <- sum(snps==0.5,na.rm=TRUE)
  n1 <- sum(snps==1,na.rm=TRUE)
  nNA <- sum(is.na(snps))

  stopifnot(n0+nh+n1+nNA == length(snps))

  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1

  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}
