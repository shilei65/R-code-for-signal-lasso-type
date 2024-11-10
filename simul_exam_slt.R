

setwd("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ signal-lasso-type")
Rcpp::sourceCpp("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ signal-lasso-type/main_slt.cpp")


source("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ signal-lasso-type/main_slt.R")
source("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ signal-lasso-type/simulation_slt.R")


library(Rcpp)
library(lars)
library(msgps)
library(MASS)
library(Matrix)
library(glmnet)
library(ncvreg)


###caseI 1 2 3 gaussian####
n=150
p=60
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(1,2)){
print(sigma)  
print(simulation_asl(n,beta,k=30,sigma,errortype="AR1",disigntype="correlated",printout = FALSE))
}




n=50
p=150
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(1,2)){
  print(sigma)  
  print(simulation_asl(n,beta,k=50,sigma,errortype="gaussian",disigntype="illness",printout = FALSE))
}



n=50
p=150
p1=120
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(0.4,1,2)){
  print(sigma)  
  print(simulation_asl(n,beta,k=100,sigma,errortype="AR1",disigntype="correlated",printout = FALSE))
}



# ###caseI 4 5 6 gaussian####
# n=100
# p=150
# p1=10
# beta=c(rep(1,p1),rep(0,p-p1))
# set.seed(1)
# for(sigma in c(0.4,1,2)){
#   print(sigma)
# print(simulation(n,beta,k=500,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
# 
# }


###caseI 7 8 9 gaussian####
n=100
p=30
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(2)){
  print(sigma) 
  print(simulation_asl(n,beta,k=200,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
}

##########
n=100
p=30
p1=20
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(0.4,1,2)){
  print(sigma) 
  print(simulation_asl(n,beta,k=100,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
}



# ###caseI 10 11 12 gaussian####
# n=50
# p=500
# p1=8
# beta=c(rep(1,p1),rep(0,p-p1))
# set.seed(1)
# for(sigma in c(0.4,1,2)){
#   print(sigma) 
#   print(simulation(n,beta,k=500,sigma,errortype="gaussian",disigntype="uncorrelated",printout = FALSE))
# }



###caseI 13 14 15 gaussian####
n=50
p=150
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(2)){
  print(sigma)  
  print(simulation_asl(n,beta,k=500,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
}


# ###caseI 16 17 18 gaussian####
# n=100
# p=150
# p1=50
# beta=c(rep(1,p1),rep(0,p-p1))
# set.seed(1)
# for(sigma in c(0.4,1,2)){
#   print(sigma) 
#   print(simulation(n,beta,k=500,sigma,errortype="gaussian",disigntype="uncorrelated",printout = FALSE))
#   
# }

###caseI 19 20 21 gaussian####
n=50
p=150
p1=20
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(2)){
  print(sigma) 
  print(simulation_asl(n,beta,k=200,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
} 

# ###caseI 22 23 24 gaussian####
n=50
p=150
p1=100
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(0.4,1,2)){
  print(sigma) 
  print(simulation_asl(n,beta,k=100,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
} 
