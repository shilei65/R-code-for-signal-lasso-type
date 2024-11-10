AR.1=function(n,rho,sigma){
  x=numeric(n)
  x[1]=rnorm(1,sd=sigma)
  for(t in 2:n){
    x[t]=rho*x[t-1]+rnorm(1,sd=sigma)
  }
  return(x)
}

#计算误差相关指标的函数


error_compute<-function(beta, beta1){
  signal.1=1*I(beta1>0.9 & beta1<1.1)
  signal.0=1*I(beta1>-0.1 & beta1<0.1)
  TP=sum(signal.1*beta)
  FP=sum(signal.1*(1-beta))
  TN=sum(signal.0*(1-beta))
  FN=sum(signal.0*beta)
  UC=length(beta)-TP-FN-FP-TN
  UC1=sum(beta)-TP-FN
  UC0=sum(1-beta)-TN-FP
  return(list(TP = TP, FP = FP, TN = TN, FN = FN, UC = UC, UC1 = UC1, UC0 = UC0))
}


#######simulation function####
simulation_slt=function(n,beta,k,sigma,errortype="gaussian",disigntype="uncorrelated",printout=TRUE){
#n: number of samples
#beta: parameters  
#k: repeat times 
#sigma: scale parameters for different distributions
#errortype: include "gaussian", "exponential","gamma","AR(1)"
#disigntype: include "uncorrelated","correlated" 
#printout: print the iter times
p=length(beta)
p1=sum(I(beta!=0))

beta.l=TPR.l=TNR.l=MCC.l=FPR.l=SREL.l=SRNL.l=PCR.l=MCCs.l=G.l=UC.l=NULL  ##lasso
beta.al=TPR.al=TNR.al=MCC.al=FPR.al=SREL.al=SRNL.al=PCR.al=MCCs.al=G.al=UC.al=NULL  ##adaptive lasso
beta.sc=TPR.sc=TNR.sc=MCC.sc=FPR.sc=SREL.sc=SRNL.sc=PCR.sc=MCCs.sc=G.sc=UC.sc=NULL  ##SCAD lasso
beta.mc=TPR.mc=TNR.mc=MCC.mc=FPR.mc=SREL.mc=SRNL.mc=PCR.mc=MCCs.mc=G.mc=UC.mc=NULL  ##MCP lasso
beta.el=TPR.el=TNR.el=MCC.el=FPR.el=SREL.el=SRNL.el=PCR.el=MCCs.el=G.el=UC.el=NULL  ##elastic net
beta.s=TPR.s=TNR.s=MCC.s=FPR.s=SREL.s=SRNL.s=PCR.s=MCCs.s=G.s=UC.s=NULL  ##signal lasso
beta.as=TPR.as=TNR.as=MCC.as=FPR.as=SREL.as=SRNL.as=PCR.as=MCCs.as=G.as=UC.as=NULL ##adaptive signal lasso
beta.as1=TPR.as1=TNR.as1=MCC.as1=FPR.as1=SREL.as1=SRNL.as1=PCR.as1=MCCs.as1=G.as1=UC.as1=NULL ##adaptive signal lasso for fixed alpha
# beta.as2=TPR.as2=TNR.as2=MCC.as2=FPR.as2=SREL.as2=SRNL.as2=PCR.as2=F.as2=G.as2=UC.as2=NULL ##adaptive signal lasso
beta.p=TPR.p=TNR.p=MCC.p=FPR.p=SREL.p=SRNL.p=PCR.p=MCCs.p=G.p=UC.p=NULL ##product penalty 
beta.m=TPR.m=TNR.m=MCC.m=FPR.m=SREL.m=SRNL.m=PCR.m=MCCs.m=G.m=UC.m=NULL  ##min penalty


for(t in 1:k){
  if(disigntype=="uncorrelated"){
  X=matrix(rnorm(n*p),n,p)
  }
  if(disigntype=="correlated"){
    X=mvrnorm(n,mu=rep(0,p),Sigma = Tpm(p,0.5))  ##sigma(i,j)=rho^(|i-j|)
  }
  if(disigntype=="zerocolumn"){
    X=matrix(rnorm(n*p),n,p)  
    #X[,(p-1):p]<-X[,1:2] ## let m columns are the same
    X[,p]<-X[,p]-X[,p] ## let last columns equal to 0
    X[,p-1]<-X[,p-1]-X[,p-1]
  }
  if(disigntype=="illness"){
    X=matrix(rnorm(n*p),n,p)  
    X[,p-2]<-X[,p] ## three columns are the same
    X[,p-1]<-X[,p]
  }
  if(errortype=="gaussian"){
  Y=X%*%beta+rnorm(n,sd=sigma)
  }
  if(errortype=="exponential"){
    Y=X%*%beta+rexp(n,rate=1/sigma)
  }
  if(errortype=="gamma"){
    Y=X%*%beta+rgamma(n,shape = 4,scale=sigma/2)
  }
  if(errortype=="AR1"){
    Y=X%*%beta+AR.1(n,rho=0.8,sigma)
  }
  
  if(errortype=="mixed"){
    Y=X%*%beta+0.5*rnorm(n,mean=-0.2,sd=sigma/2)+0.5*rnorm(n,mean=0.2,sd=sigma)
  }
  
  #### if p1>p2 #####
  #Y=X%*%rep(1,length(beta))-Y
  #beta=1-beta
  
  ###lasso fit#####
  cvfit0=cv.glmnet(x=X,y=Y,nfolds=5,intercept=TRUE,type.measure = "mse")
  # fit.l=glmnet(x=X,y=Y,lambda =cvfit0$lambda.1se,alpha = 1,intercept = TRUE)
  # beta1=coef.glmnet(fit.l)[-1]
  fit.l=msgps(X=X,y=as.vector(Y), penalty="enet", alpha=0, lambda=cvfit0$lambda.1se, intercept=TRUE)
  beta1=coef(fit.l)[,3][-1]
  
  beta.l=rbind(beta.l, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.l<-c(SREL.l, TP/sum(beta))
  SRNL.l<-c(SRNL.l, TN/sum(1-beta))         
  TPR.l<-c(TPR.l, TP/(TP+FN+0.001))
  FPR.l<-c(FPR.l, FP/(TN+FP+0.001))
  MCC.l<-c(MCC.l, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.l<-c(PCR.l, TP/(TP+FP))
  # F.l<-c(F.l, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.l<-c(G.l, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.l<-c(UC.l, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.l<-c(MCCs.l, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  ####adaptive lasso####
  
  fit.al=msgps(X=X,y=as.vector(Y), penalty="alasso", gamma=1, lambda=cvfit0$lambda.1se, intercept=TRUE)
  beta1=coef(fit.al)[,3][-1]
  beta.al=rbind(beta.al, beta1)
  beta0=beta1
  # signal.1=1*I(beta1>0.95 & beta1<1.05)
  # signal.0=1*I(beta1>-0.05 & beta1<0.05)
  # TPR.al=c(TPR.al,sum(signal.1[1:p1])/p1)
  # TNR.al=c(TNR.al,sum(signal.0[(p1+1):p])/(p-p1))
  # PCR.al=c(PCR.al,sum(signal.1[1:p1])/sum(signal.1))
  # 
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.al<-c(SREL.al, TP/sum(beta))
  SRNL.al<-c(SRNL.al, TN/sum(1-beta))         
  TPR.al<-c(TPR.al, TP/(TP+FN+0.001))
  FPR.al<-c(FPR.al, FP/(TN+FP+0.001))
  MCC.al<-c(MCC.al, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.al<-c(PCR.al, TP/(TP+FP))
  # F.al<-c(F.al, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.al<-c(G.al, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.al<-c(UC.al, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.al<-c(MCCs.al, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  
  ####SCAD lasso####
  cvfit=cv.ncvreg(X=X,y=Y,nfolds=5,penalty="SCAD",intercept = TRUE)
  fit.sc=cvfit$fit
  beta1=fit.sc$beta[,cvfit$min][-1]
  beta.sc=rbind(beta.sc, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.sc<-c(SREL.sc, TP/sum(beta))
  SRNL.sc<-c(SRNL.sc, TN/sum(1-beta))         
  TPR.sc=c(TPR.sc,TP/(TP+FN+0.001))
  FPR.sc<-c(FPR.sc, FP/(TN+FP+0.001))
  MCC.sc<-c(MCC.sc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.sc<-c(PCR.sc, TP/(TP+FP))

  UC.sc<-c(UC.sc, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.sc<-c(MCCs.sc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  ####MCP lasso####
  cvfit=cv.ncvreg(X=X,y=Y,nfolds=5,penalty="MCP",intercept = TRUE)
  fit.mc=cvfit$fit
  beta1=fit.mc$beta[,cvfit$min][-1]
  beta.mc=rbind(beta.mc, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.mc<-c(SREL.mc, TP/sum(beta))
  SRNL.mc<-c(SRNL.mc, TN/sum(1-beta))         
  TPR.mc=c(TPR.mc,TP/(TP+FN))
  FPR.mc<-c(FPR.mc, FP/(TN+FP))
  MCC.mc<-c(MCC.mc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.mc<-c(PCR.mc, TP/(TP+FP))
  
  UC.mc<-c(UC.mc, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.mc<-c(MCCs.mc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  ####elasticnet lasso####
  lam=cvm=NULL
  alpha0=seq(0.1,0.9,0.1)
  for(alpha1 in seq(0.1,0.9,0.1)){
  cvfit=cv.glmnet(x=X,y=Y,nfolds=5,intercept=TRUE,type.measure = "mse",alpha=alpha1)
  lam=c(lam,cvfit$lambda.1se)
  cvm=c(cvm, min(cvfit$cvm))
  }
  fit.el=glmnet(x=X,y=Y,lambda =lam[which.min(cvm)],alpha =alpha0[which.min(cvm)],intercept = TRUE)
  beta1=coef(fit.el)[-1]
  beta.el=rbind(beta.el, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.el<-c(SREL.el, TP/sum(beta))
  SRNL.el<-c(SRNL.el, TN/sum(1-beta))         
  TPR.el<-c(TPR.el, TP/(TP+FN))
  FPR.el<-c(FPR.el, FP/(TN+FP))
  MCC.el<-c(MCC.el, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.el<-c(PCR.el, TP/(TP+FP))
  
  
  UC.el<-c(UC.el, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.el<-c(MCCs.el, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  ####signal lasso####
  #fit=CV.SIGNAL(Y=Y,X=X,beta0 = beta0,nfolds=5,nlambda1=c(1,0.1,0.01,0.001,0.0001),alpha=c(seq(0.1,0.5,0.05)), 
  #              weights = 1,constant=TRUE,iter_max = 2000,delta=1e-7)
  
  fit=CV.SIGNAL(Y=Y,X=X,beta0 = beta0,nfolds=5,nlambda1=cvfit0$lambda[1:10],alpha=c(seq(0.1,0.9,0.1),1/seq(0.1,0.9,0.1)), 
                weights = 1,constant=TRUE,iter_max = 2000,delta=1e-7)
  fit.s=SIGNAL_l(Y=Y,X=X,beta0 = beta0,lambda1=fit$lambda.1se[1],
                 lambda2 =fit$lambda.1se[2],weights=1,iter_max = 2000,delta=1e-7)
  
  beta1=fit.s$Beta
  beta.s=rbind(beta.s, beta1)
  # signal.1=1*I(beta1>0.95 & beta1<1.05)
  # signal.0=1*I(beta1>-0.05 & beta1<0.05)
  # TPR.s=c(TPR.s,sum(signal.1[1:p1])/p1)
  # TNR.s=c(TNR.s,sum(signal.0[(p1+1):p])/(p-p1))
  # PCR.s=c(PCR.s,sum(signal.1[1:p1])/sum(signal.1))
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.s<-c(SREL.s, TP/sum(beta))
  SRNL.s<-c(SRNL.s, TN/sum(1-beta))         
  TPR.s<-c(TPR.s, TP/(TP+FN+0.0001))
  FPR.s<-c(FPR.s, FP/(TN+FP+0.0001))
  MCC.s<-c(MCC.s, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.s<-c(PCR.s, TP/(TP+FP))
  # F.s<-c(F.s, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.s<-c(G.s, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.s<-c(UC.s, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.s<-c(MCCs.s, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  
  ####adaptive signal lasso####
  fit=glmnet(x=X,y=Y,lambda =cvfit0$lambda.1se,alpha =0,intercept = TRUE)
  w=abs(coef.glmnet(fit)[-1])
  fit=CV.SIGNAL(Y=Y,X=X,beta0 = beta0,nfolds=5,nlambda1=2000,alpha=c(seq(0,0.8,0.01)), 
                weights =w,constant=TRUE,iter_max = 2000,delta=1e-7)
  fit.as=SIGNAL_l(Y=Y,X=X,beta0 = beta0,lambda1=fit$lambda.1se[1],
                 lambda2 =fit$lambda.1se[2],weights=w,iter_max = 2000,delta=1e-7)
  beta1=fit.as$Beta
  beta.as=rbind(beta.as, beta1)
  # 
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.as<-c(SREL.as, TP/sum(beta))
  SRNL.as<-c(SRNL.as, TN/sum(1-beta))         
  TPR.as<-c(TPR.as, TP/(TP+FN+0.001))
  FPR.as<-c(FPR.as, FP/(TN+FP+0.001))
  MCC.as<-c(MCC.as, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.as<-c(PCR.as, TP/(TP+FP))
  # F.as<-c(F.as, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.as<-c(G.as, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.as<-c(UC.as, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.as<-c(MCCs.as, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  ####adaptive signal lasso for fixed alpha####
  fit=glmnet(x=X,y=Y,lambda =cvfit0$lambda.1se,alpha =0,intercept = TRUE)
  w=abs(coef.glmnet(fit)[-1])
  # fit=CV.SIGNAL(Y=Y,X=X,beta0 = beta0,nfolds=5,nlambda1=1000,alpha=0.5, 
  #               weights =w,constant=TRUE,iter_max = 2000,delta=1e-7)
  fit.as1=SIGNAL_l(Y=Y,X=X,beta0 = beta0,lambda1=1000*0.5,
                  lambda2 =2000,weights=w,iter_max = 2000,delta=1e-7)
  beta1=fit.as1$Beta
  beta.as1=rbind(beta.as1, beta1)

  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.as1<-c(SREL.as1, TP/sum(beta))
  SRNL.as1<-c(SRNL.as1, TN/sum(1-beta))         
  TPR.as1<-c(TPR.as1, TP/(TP+FN))
  FPR.as1<-c(FPR.as1, FP/(TN+FP))
  MCC.as1<-c(MCC.as1, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.as1<-c(PCR.as1, TP/(TP+FP))
 
  UC.as1<-c(UC.as1, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.as1<-c(MCCs.as1, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  
  #### product penalty lasso####
  fit.p=PNR_l(Y,X,beta0,lambda=2000,constant=TRUE,iter_max=2000,delta=1e-7)
  beta1=fit.p$Beta
  beta.p=rbind(beta.p, beta1)
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  
  SREL.p<-c(SREL.p, TP/sum(beta))
  SRNL.p<-c(SRNL.p, TN/sum(1-beta))         
  TPR.p<-c(TPR.p, TP/(TP+FN+0.0001))
  FPR.p<-c(FPR.p, FP/(TN+FP+0.0001))
  MCC.p<-c(MCC.p, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))+.001))
  PCR.p<-c(PCR.p, TP/(TP+FP))
  # F.p<-c(F.p, 2*((TP/(TP+FP))*(TP/sum(beta)))/((TP/(TP+FP+.001))+(TP/sum(beta))))
  # G.p<-c(G.p, sqrt((TN/sum(1-beta))*(TP/sum(beta))))
  UC.p<-c(UC.p, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.p<-c(MCCs.p, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  #### min penalty lasso####
  fit.m=MDNR_l(Y,X,beta0,lambda=2000,constant=TRUE,iter_max=2000,delta=1e-7)
  beta1=fit.m$Beta
  beta.m=rbind(beta.m, beta1)
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.m<-c(SREL.m, TP/sum(beta))
  SRNL.m<-c(SRNL.m, TN/sum(1-beta))         
  TPR.m<-c(TPR.m, TP/(TP+FN+0.0001))
  FPR.m<-c(FPR.m, FP/(TN+FP+0.0001))
  MCC.m<-c(MCC.m, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))+.001))
  PCR.m<-c(PCR.m, TP/(TP+FP))
  # F.m<-c(F.m, 2*((TP/(TP+FP))*(TP/sum(beta)))/((TP/(TP+FP+.001))+(TP/sum(beta))))
  # G.m<-c(G.m, sqrt((TN/sum(1-beta))*(TP/sum(beta))))
  # 
  UC.m<-c(UC.m, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.m<-c(MCCs.m, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  
  
  
  if(printout){
  print(t)
    }
}
MCC.l[which(is.na(MCC.l))]=0
MCC.al[which(is.na(MCC.al))]=0
MCC.sc[which(is.na(MCC.sc))]=0
MCC.mc[which(is.na(MCC.mc))]=0
MCC.el[which(is.na(MCC.el))]=0
MCC.s[which(is.na(MCC.s))]=0
MCC.as[which(is.na(MCC.as))]=0
MCC.as1[which(is.na(MCC.as1))]=0
#MCC.as2[which(is.na(MCC.as2))]=0
MCC.p[which(is.na(MCC.p))]=0
MCC.m[which(is.na(MCC.m))]=0

PCR.l[which(is.na(PCR.l))]=0
PCR.al[which(is.na(PCR.al))]=0
PCR.sc[which(is.na(PCR.sc))]=0
PCR.mc[which(is.na(PCR.mc))]=0
PCR.el[which(is.na(PCR.el))]=0
PCR.s[which(is.na(PCR.s))]=0
PCR.as[which(is.na(PCR.as))]=0
PCR.as1[which(is.na(PCR.as1))]=0
#PCR.as2[which(is.na(PCR.as2))]=0
PCR.p[which(is.na(PCR.p))]=0
PCR.m[which(is.na(PCR.m))]=0

SREL.l[which(is.na(SREL.l))]=0
SREL.al[which(is.na(SREL.al))]=0
SREL.sc[which(is.na(SREL.sc))]=0
SREL.mc[which(is.na(SREL.mc))]=0
SREL.el[which(is.na(SREL.el))]=0
SREL.s[which(is.na(SREL.s))]=0
SREL.as[which(is.na(SREL.as))]=0
SREL.as1[which(is.na(SREL.as1))]=0
SREL.p[which(is.na(SREL.p))]=0
SREL.m[which(is.na(SREL.m))]=0




TPR.l[which(is.na(TPR.l))]=0
TPR.al[which(is.na(TPR.al))]=0
TPR.sc[which(is.na(TPR.sc))]=0
TPR.mc[which(is.na(TPR.mc))]=0
TPR.el[which(is.na(TPR.el))]=0
TPR.s[which(is.na(TPR.s))]=0
TPR.as[which(is.na(TPR.as))]=0
TPR.as1[which(is.na(TPR.as1))]=0
TPR.p[which(is.na(TPR.p))]=0
TPR.m[which(is.na(TPR.m))]=0
#TPR.as2[which(is.na(TPR.as2))]=0

FPR.l[which(is.na(FPR.l))]=0
FPR.al[which(is.na(FPR.al))]=0
FPR.sc[which(is.na(FPR.sc))]=0
FPR.mc[which(is.na(FPR.mc))]=0
FPR.el[which(is.na(FPR.el))]=0
FPR.s[which(is.na(FPR.s))]=0
FPR.as[which(is.na(FPR.as))]=0
FPR.as1[which(is.na(FPR.as1))]=0
FPR.p[which(is.na(FPR.p))]=0
FPR.m[which(is.na(FPR.m))]=0
#FPR.as2[which(is.na(FPR.as2))]=0

# G.l[which(is.na(G.l))]=0
# G.al[which(is.na(G.al))]=0
# G.sc[which(is.na(G.sc))]=0
# G.mc[which(is.na(G.mc))]=0
# G.el[which(is.na(G.el))]=0
# G.s[which(is.na(G.s))]=0
# G.as[which(is.na(G.as))]=0
# G.as1[which(is.na(G.as1))]=0
# G.as2[which(is.na(G.as2))]=0
# 
# F.l[which(is.na(F.l))]=0
# F.al[which(is.na(F.al))]=0
# F.sc[which(is.na(F.sc))]=0
# F.mc[which(is.na(F.mc))]=0
# F.el[which(is.na(F.el))]=0
# F.s[which(is.na(F.s))]=0
# F.as[which(is.na(F.as))]=0
# F.as1[which(is.na(F.as1))]=0
# F.as2[which(is.na(F.as2))]=0

MCCs.l[which(is.na(MCCs.l))]=0
MCCs.al[which(is.na(MCCs.al))]=0
MCCs.sc[which(is.na(MCCs.sc))]=0
MCCs.mc[which(is.na(MCCs.mc))]=0
MCCs.el[which(is.na(MCCs.el))]=0
MCCs.s[which(is.na(MCCs.s))]=0
MCCs.as[which(is.na(MCCs.as))]=0
MCCs.as1[which(is.na(MCCs.as1))]=0
#MCC.as2[which(is.na(MCC.as2))]=0
MCCs.p[which(is.na(MCCs.p))]=0
MCCs.m[which(is.na(MCCs.m))]=0

#####results####
Mse=c(mean((colMeans(beta.l)-beta)^2)+mean(apply(beta.l, 2, var)),
mean((colMeans(beta.al)-beta)^2)+mean(apply(beta.al, 2, var)),
mean((colMeans(beta.sc)-beta)^2)+mean(apply(beta.sc, 2, var)),
mean((colMeans(beta.mc)-beta)^2)+mean(apply(beta.mc, 2, var)),
mean((colMeans(beta.el)-beta)^2)+mean(apply(beta.el, 2, var)),
mean((colMeans(beta.s)-beta)^2)+mean(apply(beta.s, 2, var)),
mean((colMeans(beta.as)-beta)^2)+mean(apply(beta.as, 2, var)),
mean((colMeans(beta.as1)-beta)^2)+mean(apply(beta.as1, 2, var)),
mean((colMeans(beta.p)-beta)^2)+mean(apply(beta.p, 2, var)),
mean((colMeans(beta.m)-beta)^2)+mean(apply(beta.m, 2, var)))

TPR=c(mean(TPR.l),mean(TPR.al),mean(TPR.sc),mean(TPR.mc),
      mean(TPR.el),mean(TPR.s),mean(TPR.as),mean(TPR.as1),mean(TPR.p),mean(TPR.m))
MCC=c(mean(MCC.l),mean(MCC.al),mean(MCC.sc),mean(MCC.mc),
      mean(MCC.el),mean(MCC.s),mean(MCC.as),mean(MCC.as1),mean(MCC.p),mean(MCC.m))
PCR=c(mean(PCR.l),mean(PCR.al),mean(PCR.sc),mean(PCR.mc),
      mean(PCR.el),mean(PCR.s),mean(PCR.as),mean(PCR.as1),mean(PCR.p),mean(PCR.m))

# F1=c(mean(F.l),mean(F.al),mean(F.sc),mean(F.mc),
#      mean(F.el),mean(F.s),mean(F.as),mean(F.as1),mean(F.as2))
# G_mean=c(mean(G.l),mean(G.al),mean(G.sc),mean(G.mc),
#          mean(G.el),mean(G.s),mean(G.as),mean(G.as1),mean(G.as2))

SREL=c(mean(SREL.l),mean(SREL.al),mean(SREL.sc),mean(SREL.mc),
      mean(SREL.el),mean(SREL.s),mean(SREL.as),mean(SREL.as1),mean(SREL.p),mean(SREL.m))
SRNL=c(mean(SRNL.l),mean(SRNL.al),mean(SRNL.sc),mean(SRNL.mc),
       mean(SRNL.el),mean(SRNL.s),mean(SRNL.as),mean(SRNL.as1),mean(SRNL.p),mean(SRNL.m))
FPR=c(mean(FPR.l),mean(FPR.al),mean(FPR.sc),mean(FPR.mc),
      mean(FPR.el),mean(FPR.s),mean(FPR.as),mean(FPR.as1),mean(FPR.p),mean(FPR.m))

UCR=c(mean(UC.l),mean(UC.al),mean(UC.sc),mean(UC.mc),
      mean(UC.el),mean(UC.s),mean(UC.as),mean(UC.as1),mean(UC.p),mean(UC.m))

MCCs=c(mean(MCCs.l),mean(MCCs.al),mean(MCCs.sc),mean(MCCs.mc),
       mean(MCCs.el),mean(MCCs.s),mean(MCCs.as),mean(MCCs.as1),mean(MCCs.p),mean(MCCs.m))

out=rbind(Mse,SREL,SRNL,TPR,PCR,FPR,MCC,MCCs,UCR)
colnames(out)=c("lasso","alasso","SCAD","MCP",
                "elastic","signal","asignal","asignal_1", "product", "min")
return(t(out))
}


