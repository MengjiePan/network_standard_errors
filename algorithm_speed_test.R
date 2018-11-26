#.libPaths( c( .libPaths(), "~/Rlibs") )
library(tictoc)
library(mvtnorm)
library(gtools)
source('/Users/Mengjie/Desktop/network_se/code_au18/utility_file.R')
source('/Users/Mengjie/Desktop/network_se/code_au18/simulation_file.R')

n=80
ng=n/2

covariate_params=c(0.66,0.33)
covariate_options="normal"
K=2
beta=rep(1,2)
sigma_epsilon=0.866/4
sigma_a=0.957/4
sigma_b=0.677/4
sigma_gamma=0.677/4
sigma_z=0.677/4
rou=0.5
d=2

sigma_a_vec=c(sigma_a*2,sigma_a/2)
sigma_b_vec=c(sigma_b*2,sigma_b/2)
sigma_z_vec=c(sigma_z*2,sigma_z/2)
sigma_gamma_m=sqrt(sigma_z_vec)%*%t(sqrt(sigma_z_vec))

start_time <- Sys.time()
tmp=simulation(n,ng,beta,sigma_epsilon,sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,rou,d,niter=1,covariate_options,covariate_params)
end_time <- Sys.time()
end_time - start_time
