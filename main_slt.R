
####the correlation function######
Tpm=function(p,rho){
  A=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i==j){A[i,j]=1}
      else{A[i,j]=rho^(abs(i-j))}
    }
    
  }
  return(A)
}
########The signal lasso Function lambda1*|X|+lambda2*|X-1|########
SIGNAL_l=function(Y,X,beta0,lambda1,lambda2,weights=1,constant=TRUE,iter_max,delta){
  #Y: the response vector
  #X: the covariate matrix
  #beta0: initial value of the estimated parameter
  #lambda1: the first penalty parameter
  #lambda2: the second penalty parameter
  #weights: if 1 is the signal lasso and a vector of length p is the adaptive lasso
  #iter_max: the maximal iterations
  #delta: control the accuracy 
  p=ncol(X)
  Re2=numeric(p)
  for (i in 1:p) {
    Re2[i]=t(X[,i])%*%X[,i]  
  } 
  lambda3=weights*lambda2
  ep1=(lambda1+lambda3)/Re2
  ep2=(lambda1-lambda3)/Re2
  ptm=proc.time()
  #beta1=numeric(p)
  #iter=1
  #   eps=Y-X%*%beta0
  # while (iter<=iter_max & sqrt(mean((beta1-beta0)^2))>delta) {
  #   beta1=beta0
  #   if(constant){mu0=mean(eps)}else{mu0=0}
  #   eps=eps-mu0
  #   for(i in 1:p){
  #      eps=eps+beta0[i]*X[,i]
  #      Re1=t(X[,i])%*%eps
  #      Re=Re1/Re2[i]
  #      if(Re<=0){
  #         beta0[i]=min(0,Re+ep1[i])
  #      }
  #     if(Re>0 & Re<=1+ep2[i]){
  #         beta0[i]=max(0,Re-ep2[i]) 
  #      }
  #      if(Re>1+ep2[i]){
  #         beta0[i]=max(1,Re-ep1[i]) 
  #      }
  #      eps=eps-beta0[i]*X[,i]
  #   } 
  #  iter=iter+1
  #}
  #mse=mean((Y-mu0-X%*%beta0)^2)
  #return(list(Mu=mu0,Beta=beta0,Mse=mse, Iter=iter,Times=proc.time()-ptm))
  #Those iterations code are replaced by the following Cpp Function
  con=1*constant
  fit.c=Signal_c(Y=Y, X=X, beta0=beta0, Re2=Re2, 
                 ep1=ep1, ep2=ep2,constant = con, iter_max=iter_max, p=p, delta=delta)
  mse=mean((Y-fit.c$Mu-X%*%fit.c$Beta)^2)
  return(list(Mu=fit.c$Mu, Beta=fit.c$Beta,Mse=mse,Iter=fit.c$iters,Times=proc.time()-ptm))
}

#######The CV function for signal lasso####
CV.SIGNAL=function(Y,X,beta0,nfolds,nlambda1,alpha=seq(0.3,0.7,length.out=5),weights=1,constant=TRUE,iter_max,delta){
  n=length(Y)
  #cvfit=cv.glmnet(x=X,y=Y,nfolds =nfolds,intercept=constant,type.measure = "mse")
  #lar1=glmnet(x=X,y=Y,lambda =cvfit$lambda.1se,alpha = 1,intercept = FALSE)
  #beta0=coef.glmnet(lar1)[-1]
  #if(length(nlambda)==1){
  #   lambda=cvfit$lambda[1:nlambda]
  #}
  #else{lambda=nlambda}
  
  folds=cv.folds(n,nfolds)
  
  re=NULL
  lambda3=NULL
  for (lambda2 in nlambda1) {
    for (j in alpha) {
      se=0
      lambda1=j*lambda2
      for (k in 1:nfolds) {
        Y1=Y[as.vector(unlist(folds[-k]))]
        X1=X[as.vector(unlist(folds[-k])),]
        Y2=Y[as.vector(unlist(folds[k]))]
        X2=X[as.vector(unlist(folds[k])),]
        fit=SIGNAL_l(Y=Y1,X=X1,beta0=beta0,lambda1,lambda2,weights=weights,constant=constant,iter_max,delta)
        se=se+mean((Y2-fit$Mu-X2%*%fit$Beta)^2)
      }
      re=c(re,se)
      lambda3=rbind(lambda3,c(lambda1,lambda2))
    }
  }
  index=which.min(re)
  return(list(lambda.1se=lambda3[index,],Re=re))   
}

###############The min penalty function   lambda*min(|X|,|X-1|)  ########
MDNR_l=function(Y,X,beta0,lambda,constant=TRUE,iter_max,delta){
  #Y: the response vector
  #X: the covariate matrix
  #beta0: initial value of the estimated parameter
  #lambda1: the first penalty parameter
  #lambda2: the second penalty parameter
  #iter_max: the maximal iterations
  #delta: control the accuracy
  p=ncol(X)
  Re2=numeric(p)
  for (i in 1:p) {
    Re2[i]=t(X[,i])%*%X[,i]  
  }
  ep=lambda/Re2
  ptm=proc.time()
  # beta1=numeric(p)
  # iter=1
  # while (iter<iter_max & sqrt(mean((beta1-beta0)^2))>delta) {
  #  beta1=beta0
  #  eps=Y-X%*%beta0
  #   if(constant){mu0=mean(eps)}else{mu0=0}
  #   eps=eps-mu0
  #  for(i in 1:p){
  #    eps=eps+beta0[i]*X[,i]
  #    Re1=t(X[,i])%*%eps
  #    Re=Re1/Re2[i]
  #   if(abs(Re)<=abs(Re-1)){
  #      beta0[i]=sign(Re)*max(0,abs(Re)-ep[i])
  #    }
  #    else{
  #      beta0[i]=1+sign(Re-1)*max(0,abs(Re-1)-ep[i])
  #     }
  #    eps=eps-beta0[i]*X[,i]
  # }
  #  iter=iter+1
  # }
  # mse=mean((Y-mu0-X%*%beta0)^2)
  # return(list(Mu=mu0,Beta=beta0,Mse=mse, Iter=iter,Times=proc.time()-ptm))
  #Those iterations code are replaced by the following Cpp Function
  con=1*constant
  fit.c=MDNR_c(Y=Y, X=X, beta0=beta0, Re2=Re2,
               ep=ep, constant=con,iter_max=iter_max, p=p, delta=delta)
  mse=mean((Y-fit.c$Mu-X%*%fit.c$Beta)^2)
  return(list(Mu=fit.c$Mu,Beta=fit.c$Beta, Mse=mse, Iter=fit.c$iters,Times=proc.time()-ptm))
}


###############The product penalty  function   lambda*|X|*|X-1|  ########
PNR_l=function(Y,X,beta0,lambda,constant=TRUE,iter_max,delta){
  #Y: the response vector
  #X: the covariate matrix
  #beta0: initial value of the estimated parameter
  #lambda1: the first penalty parameter
  #lambda2: the second penalty parameter
  #iter_max: the maximal iterations
  #delta: control the accuracy
  p=ncol(X)
  Re2=numeric(p)
  for (i in 1:p) {
    Re2[i]=t(X[,i])%*%X[,i]  
  }
  ptm=proc.time()
  # beta1=numeric(p)
  # iter=1
  # while (iter<iter_max & sqrt(mean((beta1-beta0)^2))>delta) {
  #  beta1=beta0
  #  eps=Y-X%*%beta0
  #   if(constant){
  #      mu0=mean(eps)
  #   } else{
  #         mu0=0
  #     }
  #   eps=eps-mu0
  #  for(i in 1:p){
  #    eps=eps+beta0[i]*X[,i]
  #    Re1=t(X[,i])%*%eps
  #    Re=Re1/Re2[i]
  #   if(Re<=0){
  #      beta0[i]=min(0,(Re1+lambda/2)/(Re2[i]+lambda))
  #    }
  #   if(Re>0 & Re<=0.5){
  #      if(Re2[i]>lambda){
  #      beta0[i]=max(0,(Re1-lambda/2)/(Re2[i]-lambda))
  #      } 
  #      else{
  #        beta0[i]=0
  #        }
  #    }
  #   if(Re>0.5 & Re<=1){
  #      if(Re2[i]>lambda){
  #        beta0[i]=min(1,(Re1-lambda/2)/(Re2[i]-lambda))
  #       } 
  #      else{
  #         beta0[i]=1
  #       }
  #    }
  #   if(Re>1){
  #      beta0[i]=max(1,(Re1+lambda/2)/(Re2[i]+lambda))
  #    }
  #   eps=eps-beta0[i]*X[,i]
  # }
  #  iter=iter+1
  # }
  # mse=mean((Y-mu0-X%*%beta0)^2)
  # return(list(Mu=mu0,Beta=beta0,Mse=mse, Iter=iter,Times=proc.time()-ptm))
  #Those iterations code are replaced by the following Cpp Function
  con=1*constant
  #lambda=weights*lambda
  fit.c=PNR_c(Y=Y, X=X, beta0=beta0, Re2=Re2,
              lambda=lambda, constant=con,iter_max=iter_max, p=p, delta=delta)
  mse=mean((Y-fit.c$Mu-X%*%fit.c$Beta)^2)
  return(list(Mu=fit.c$Mu,Beta=fit.c$Beta, Mse=mse, Iter=fit.c$iters,Times=proc.time()-ptm))
}
