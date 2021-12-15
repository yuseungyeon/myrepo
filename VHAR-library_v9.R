require(MASS)
require(glmnet)

VHAR.sim = function(N=500, A, Sigma){
  
  ## Transform A into VAR(22) coefficients
  k = dim(A)[1];
  burn = 500;
  
  D = A[,(1:k)];
  We = A[,((k+1):(2*k))];
  Mo = A[,-(1:(2*k))];
  
  AA= D + .2*We + 1/22*Mo;
  for(i in 1:4){
    AA = cbind(AA, .2*We + 1/22*Mo)
  }
  for(i in 1:17){
    AA = cbind(AA, 1/22*Mo)
  }
  
  #p = dim(AA)[2]/k;
  p=22
  
  inno = mvrnorm(n=N+burn, rep(0, k), Sigma);
  init = mvrnorm(n=p, rep(0, k), Sigma);
  init = matrix(init, nrow=p);
  
  # Find index for previous observations
  j=1;
  # ar term
  id = seq(from= j+p-1, to = j, by=-1);
  
  Y = matrix(0, (N+burn), k);
  for(r in 1:(N+burn)){
    Y[r,] = AA%*%as.vector(t(init[id,])) + inno[r,];
    init = rbind(init[-1,], Y[r,]);
  }
  
  return(Y[-(1:burn),]) # Final data is T*dim matrix
}



###########################################
## OLS estimation of VHAR
###########################################
VHAR.ols <- function(data){
  
  ## Input data is length*dim
  
  VHAR.design.matrix <- function(data){
    D <- t(data)[,-(1:21)]
    D <- D[,-(ncol(D))]
    
    W <- matrix(0, nrow = nrow(data)-22, ncol=ncol(data))
    W.data <- data[-(1:17),]
    for(i in 1:ncol(W)){
      temp <- na.omit( stats::filter(W.data[,i], rep(1/5, 5), sides=1) )
      W[,i] <- temp[-length(temp)]
    }
    W <- t(W)
    
    M <- matrix(0, nrow = nrow(data)-22, ncol=ncol(data))
    M.data <- data
    for(i in 1:ncol(M)){
      temp <- na.omit(stats::filter(M.data[,i], rep(1/22,22), sides=1) )
      M[,i] <- temp[-length(temp)]
    }
    M <- t(M)
  
    return(rbind(D,W,M))
  }

  VHAR.y.matrix <- function(data){
    # dimension of return matrix is [# of variables , 2067]
    return( t(data)[,-(1:22)] )
  }
  
  k = dim(data)[2]
  ######################################################################
  # dim(y) is [# of variables , n-22] and dim(x) is [# of variables * 3 , n-22]
  x <- VHAR.design.matrix(data)
  y <- VHAR.y.matrix(data)
  hatA <- y %*% t(x) %*% solve(x %*% t(x));
  return(list(hatA = hatA, D=hatA[,(1:k)], W = hatA[,((k+1):(2*k))], M = hatA[,-(1:(2*k))]))
}

###########################################
## Adaptive LASSO estimation of VHAR
###########################################

VHAR.adalasso <- function(data, nf=10, type="cv", gam=1, updateSigma =TRUE, debiasTF = FALSE){
  
  VHAR.design.matrix <- function(data){
    D <- t(data)[,-(1:21)]
    D <- D[,-(ncol(D))]
    
    W <- matrix(0, nrow = nrow(data)-22, ncol=ncol(data))
    W.data <- data[-(1:17),]
    for(i in 1:ncol(W)){
      temp <- na.omit( stats::filter(W.data[,i], rep(1/5, 5), sides=1) )
      W[,i] <- temp[-length(temp)]
    }
    W <- t(W)
    
    M <- matrix(0, nrow = nrow(data)-22, ncol=ncol(data))
    M.data <- data
    for(i in 1:ncol(M)){
      temp <- na.omit( stats::filter(M.data[,i], rep(1/22,22), sides=1) )
      M[,i] <- temp[-length(temp)]
    }
    M <- t(M)
    
    # return matrix dimension is [# of variables * 3 , n-22]
    return(rbind(D,W,M))
  }
  
  VHAR.y.matrix <- function(data){
    # dimension of return matrix is [# of variables , 2067]
    return( t(data)[,-(1:22)] )
  }

  ######################################
  # Sigma estimation in VAR
  ######################################
  VHAR.sigma = function(y, x, hA1){
    
    Resi = y - hA1%*%x;
    Sigma_z = Resi%*%t(Resi)/ncol(Resi);
    out = list();
    out$Resi = Resi;
    out$Sigma_z = Sigma_z;
    return(out);
  }
  
  x <- VHAR.design.matrix(data)
  y <- VHAR.y.matrix(data)
  k <- dim(y)[1]
  hatA0 = VHAR.ols(data)$hatA;
  
  if(updateSigma == TRUE){
    Sigma_z = VHAR.sigma(y, x, hatA0)$Sigma_z;
    sv=svd(Sigma_z)
    hal = sv$u%*%diag(1/sqrt(sv$d))%*%t(sv$v)
    
    newx = kronecker(t(x), diag(1, k));
    newy = as.matrix(as.vector(y));
    
    adjSig = kronecker(diag(1, dim(y)[2]), hal);
    newy = adjSig%*%newy;
    newx = adjSig%*%newx;
  } else {
    newx = kronecker(t(x), diag(1, k));
    newy = as.matrix(as.vector(y));
    
  }
#  # weight : initial estimate of alpha ( ridge estimate or lasso estimate )
#  ridge <- cv.glmnet(x=newx, y=newy, alpha=0, parallel=TRUE, intercept=FALSE)
#  weight <- 1/( abs( matrix( coef(ridge, s=ridge$lambda.min) ) ) )^gam
#  weight[which(weight==Inf)] <- 100000;

#  cvfit = cv.glmnet( x=newx, y=newy, alpha = 1, intercept=FALSE, type.measure = "mse", 
#                     nfolds=nf, penalty.factor = weight, parallel = TRUE )
#  cf.cv = coef(cvfit, s = "lambda.min")
 
  if(type == "cv"){
  pp = adalasso(newx, newy, k=nf, intercept=FALSE); 
  out=pp;
  A.lasso = matrix(pp$coefficients.lasso, nrow=k);
  A.adalasso = matrix(pp$coefficients.adalasso, nrow=k);
  } else {
  
  ## BIC selection of tuning parameters
#  lasso=ic.glmnet(newx,newy,crit = "bic")
#  A.lasso = matrix(lasso$coefficients[-1], nrow=k);
  ## == Adaptive LASSO == ##
  # Lasso as the first step model, intercept must be removed
  # to calculate the penalty factor.
#  first.step.coef=coef(lasso)[-1]
#  penalty.factor=abs(first.step.coef+1/sqrt(nrow(x)))^(-gam);
  penalty.factor = abs(as.vector(hatA0)+1/sqrt(nrow(x)))^(-gam);
  ada =ic.glmnet(newx,newy,crit="bic",penalty.factor=penalty.factor)
  A.adalasso = matrix(ada$coefficients[-1], nrow=k);
  out = ada;
#  
#  cc = glmnet(newx, newy, penalty.factor=penalty.factor, lambda = .005);
#  corrplot(matrix(cc$beta, nrow=k), method="color", is.corr=FALSE);
  }
  
#  out$A.lasso = A.lasso;
  out$A.adalasso = A.adalasso;
  out$A.ols = hatA0;
  
  if(debiasTF){
    
    const.ols = function(x,y, A.lasso){ 
      J = 1*(A.lasso != 0); p=3;
      vecA1 = c(J);
      A1 = matrix(vecA1, nrow=k);
      k2 = sum(vecA1);
      R = matrix(0, (k^2)*p, k2);
      
      for(i in 1:(k^2*p)){
        if(vecA1[i] == 1){ 
          id = sum(vecA1[1:i]);  
          R[i,id] = 1;
        }  
      }
      
      #  Sigma_update = VHAR.sigma(y, x, hatA0)$Sigma_z;
      Sigma_update = VHAR.sigma(y, x, A.lasso)$Sigma_z;
      sv=svd(Sigma_update)
      szInv = sv$u%*%diag(1/sv$d)%*%t(sv$v)
      varA1 = solve(t(R) %*% kronecker(x %*% t(x), szInv) %*% R);
      A1.db = matrix( R %*% varA1 %*% t(R) %*% kronecker(x, szInv) %*% as.vector(y), nrow=k);
      return(A1.db)
    }
    
  A.lasso.db = const.ols(x,y, A.lasso);
  A.adalasso.db = const.ols(x,y, A.adalasso);
  
  out$A.lasso.db = A.lasso.db;
  out$A.adalasso.db = A.adalasso.db;
  }
  
#  out$cv=cvfit;
#  out$hatA = matrix(cf.cv[-1], nrow=k, byrow=FALSE);
#  out$lam = cvfit$lambda.min;
#  out$mcv = min(cvfit$cvm);

  return(out)
}

# 
# ################################################
# ### forecasting
# ################################################
# 
# VHAR.input <- function(train){
#   
#   D <- t(train)[,-(1:21)]
#   D <- D[,(ncol(D))]
#   
#   W <- matrix(0, nrow = nrow(train)-21, ncol=ncol(train))
#   W.train <- train[-(1:17),]
#   for(i in 1:ncol(W)){
#     temp <- na.omit( stats::filter(W.train[,i], rep(1/5, 5), sides=1) )
#     W[,i] <- temp
#   }
#   W <- t(W)
#   W <- W[,(ncol(W))]
#   
#   M <- matrix(0, nrow = nrow(train)-21, ncol=ncol(train))
#   M.train <- train
#   for(i in 1:ncol(M)){
#     temp <- na.omit( stats::filter(M.train[,i], rep(1/22,22), sides=1) )
#     M[,i] <- temp
#   }
#   M <- t(M)
#   M <- M[,(ncol(M))]
#   
#   Yt <- as.matrix(c(D,W,M), ncol=1)
#   
#   return(Yt)
# }
# 


