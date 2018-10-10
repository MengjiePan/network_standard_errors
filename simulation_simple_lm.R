#This script considers a one predictor model y=beta0+beta1X1 where X1 is the 
#bernoulli: Xij is 1(Xi=Xj), where Xi and Xj follow bernoulli distribution with pi, which is correlated with gi
#normal: Xij follow normal distribution with mean zero and sd correlated with gi and gj
#diff: absolute difference between X1i and X1j, X1i and X1j follow standard normal with different sd.

library(gtools)
library(mvtnorm)
source('/homes/mpan1/network_se/code_au18/function_file.R')
main=function(Iter,n,ng,beta,sigma_epsilon,sigma_a_vec, sigma_b_vec,sigma_z_vec,sigma_gamma_m,rou,d,niter,covariate_options,covariate_params){
  CI_record=simulation(n,ng,beta,sigma_epsilon,sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,rou,d,niter,covariate_options,covariate_params)
  #return(CI_record)
  save(CI_record,file=paste0('/homes/mpan1/network_se/results_au18/',covariate_options,'/CI_record_n_',n,'_ng_',ng,'_sd_',covariate_params[1],'_',covariate_params[2],'_Iter_',Iter,'.RData'))
}


args = commandArgs(trailingOnly=TRUE)
Iter=as.integer(args[1])
n=as.integer(args[2])
ng=as.integer(args[3])
covariate_params=c(as.numeric(args[4]),as.numeric(args[5]))
covariate_options=args[6]
print(n)
print(ng)
K=2
beta=rep(1,2)
sigma_epsilon=0.866
sigma_a=0.957
sigma_b=0.677
sigma_gamma=0.677
sigma_z=0.677
rou=0.5
d=2

sigma_a_vec=c(sigma_a*2,sigma_a/2)
sigma_b_vec=c(sigma_b*2,sigma_b/2)
sigma_z_vec=c(sigma_z*2,sigma_z/2)
#should change this to sqrt
sigma_gamma_m=sqrt(sigma_z_vec)%*%t(sqrt(sigma_z_vec))

main(Iter,n,ng,beta,sigma_epsilon,sigma_a_vec, sigma_b_vec,sigma_z_vec,sigma_gamma_m,rou,d,niter=1000,covariate_options,covariate_params)

