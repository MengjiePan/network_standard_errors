# simulate errors using the error generative model
simulate_errors=function(n,g,sigma_epsilon,sigma_a_vec, sigma_b_vec,sigma_z_vec,sigma_gamma_m,rou,d){
  a=rep(NA,n)
  b=rep(NA,n)
  
  for (i in 1:n){
    tmp=rmvnorm(1, mean=rep(0,2), sigma=matrix(c(sigma_a_vec[g[i]]^2, rou*sigma_a_vec[g[i]]*sigma_b_vec[g[i]],rou*sigma_a_vec[g[i]]*sigma_b_vec[g[i]],sigma_b_vec[g[i]]^2),nrow=2,ncol=2))
    a[i]=tmp[1]
    b[i]=tmp[2]
  }
  
  z=matrix(rnorm(d*n,mean=0,sd=rep(sigma_z_vec[g],d)),nrow=n,ncol=d,byrow=F)
  gamma=matrix(NA,nrow=n,ncol=n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      gamma[i,j]=rnorm(1,mean=0,sd=sigma_gamma_m[g[i],g[j]])
      gamma[j,i]=gamma[i,j]
    }
  }
  
  xi=matrix(NA,nrow=n,ncol=n)
  for (i in 1:n){
    for (j in (1:n)[-i]){
      xi[i,j]=a[i]+b[j]+sum(z[i,]*z[j,])+gamma[i,j]+rnorm(1,mean=0,sd=sigma_epsilon)
    }
  }
  return(xi)
}

# randomly assign ng nodes to block 1 and the rest to block 2
simulate_g=function(n,ng){
  g_ind=sample(1:n,size=ng,replace=F)
  g=numeric(n)+1
  g[g_ind]=2
  return(g)
}

# simulate X_{ij}=|X_i-X_j|, and X_i and X_j follow normal with zero mean
simulate_X_absdiff=function(n,sd_options,g){
  X_vec=rnorm(n=n,mean=0,sd=sd_options[g])
  X=outer(X_vec, X_vec, FUN = function(x1,x2){return (abs(x1-x2))})
  return(X)
}

# simulate X_{ij} follow normal with zero mean
simulate_X_normal=function(n,sd_options,g){
  sigma_matrix=sqrt(sd_options[g] %*% t(sd_options[g] ))
  X=rnorm(n=n^2,mean=0,sd=c(sigma_matrix))
  X=matrix(X,nrow=n,ncol=n,byrow=F)
  return(X)
}

# simulate X_{ij}=1_{X_i==X_j}
simulate_X_bernoulli=function(n,probs,g){
  X_vec=rbinom(n=n,size=1,prob=probs[g])
  X=outer(X_vec, X_vec, FUN = function(x1,x2){return (as.numeric(x1==x2))})
  return(X)
}

# get observations y
simulate_y=function(xi,X,beta){
  n=dim(xi)[1]
  y=matrix(beta[1],nrow=n,ncol=n)+beta[2]*X+xi
  return(y)
}

# run one simulation of X and g, with niter simulated errors
# output niter confidence intervals for differenr estimators
simulation=function(n,ng,beta,sigma_epsilon,sigma_a_vec, sigma_b_vec,sigma_z_vec,sigma_gamma_m,rou,d,niter,covariate_options,covariate_params){
  tic("Total time to do steps that only need to be done once for a given X and g")
  tic("getting exch node set, given network size n")
  nodeSetExch6=node.set.exch(n)
  nodeSetExch5=nodeSetExch6
  nodeSetExch5[[5]]=rbind(nodeSetExch5[[5]],nodeSetExch5[[6]])
  toc()
  p=length(beta)
  g=simulate_g(n,ng)
  K=max(g)
  if (covariate_options=='bernoulli'){
    X=simulate_X_bernoulli(n,covariate_params,g)
  }else if (covariate_options=='normal'){
    X=simulate_X_normal(n,covariate_params,g)
  }else if (covariate_options=='absdiff'){
    X=simulate_X_absdiff(n,covariate_params,g)
  }else{
    stop("covariate option not supported")
  }
  X_m=cbind(rep(1,n*(n-1)),c(X)[-seq(1,n^2,n+1)])
  tic("get oracle omega")
  sigmaE=get_sigmaE_true(n,g,sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,K,d,rou,nodeSetExch6)
  toc()
  beta_covariance_trueparams=get_beta_covariance(X_m,sigmaE)
  tic("get block-exchangeable node set with block oracle")
  nodeSetTrue=node.set(g,nodeSetExch6)
  toc()
  
  CI_trueparams_record=NULL
  CI_estparams_trueblock_record=NULL
  CI_estparams_estblock_record=NULL
  CI_exch_record=NULL
  CI_dyadic_record=NULL
  toc()
  tic("Total time to complete one iteration of simulating errors and calculating confidence intervals using different estimators (repeated 1000 times given X and g)")
  for (iter in 1:niter){
    xi=simulate_errors(n,g,sigma_epsilon,sigma_a_vec, sigma_b_vec,sigma_z_vec,sigma_gamma_m,rou,d)    
    y=simulate_y(xi,X,beta)
    mod=get_linear_model(y,X)
    residuals=mod$residuals
    coefs=mod$coefs
    tic("get block-exchangeable omega with block oracle")
    est_sigmaE_trueblock=get_sigmaE_est(n,nodeSetTrue,residuals)
    toc()
    beta_covariance_estparams_trueblock=get_beta_covariance(X_m,est_sigmaE_trueblock)
    beta_covariance_estparams_trueblock=make.positive.var(beta_covariance_estparams_trueblock)$V
    
    tic("get estimated blocks")
    est_blocks=get_est_blocks(n,residuals,nodeSetExch6,K)
    toc()
    tic("get block-exchangeable node set with estimated blocks")
    nodeSetEst=node.set(est_blocks,nodeSetExch6)
    toc()
    tic("get block-exchangeable omega with estimated blocks")
    est_sigmaE_estblock=get_sigmaE_est(n,nodeSetEst,residuals)
    toc()
    beta_covariance_estparams_estblock=get_beta_covariance(X_m,est_sigmaE_estblock)
    beta_covariance_estparams_estblock=make.positive.var(beta_covariance_estparams_estblock)$V
    
    tic("get exch omega")
    est_sigmaE_exch=get_sigmaE_est(n,nodeSetExch5,residuals)
    toc()
    beta_covariance_exch=get_beta_covariance(X_m,est_sigmaE_exch)
    beta_covariance_exch=make.positive.var(beta_covariance_exch)$V
    
    tic("get DC omega")
    est_sigmaE_dyadic=get_sigmaE_dyadic(n,nodeSetExch5,residuals)
    toc()
    beta_covariance_dyadic=get_beta_covariance(X_m,est_sigmaE_dyadic)
    beta_covariance_dyadic=make.positive.var(beta_covariance_dyadic)$V
    
    CI_trueparams=rep(NA,p*2)
    CI_estparams_trueblock=rep(NA,p*2)
    CI_estparams_estblock=rep(NA,p*2)
    CI_exch=rep(NA,p*2)
    CI_dyadic=rep(NA,p*2)
    for (beta_index in 1:p){
      CI_trueparams[2*beta_index-1]=coefs[beta_index]-qnorm(.975)*sqrt(diag(beta_covariance_trueparams))[beta_index]
      CI_trueparams[2*beta_index]=coefs[beta_index]+qnorm(.975)*sqrt(diag(beta_covariance_trueparams))[beta_index]
      CI_estparams_trueblock[2*beta_index-1]=coefs[beta_index]-qnorm(.975)*sqrt(diag(beta_covariance_estparams_trueblock))[beta_index]
      CI_estparams_trueblock[2*beta_index]=coefs[beta_index]+qnorm(.975)*sqrt(diag(beta_covariance_estparams_trueblock))[beta_index]
      CI_estparams_estblock[2*beta_index-1]=coefs[beta_index]-qnorm(.975)*sqrt(diag(beta_covariance_estparams_estblock))[beta_index]
      CI_estparams_estblock[2*beta_index]=coefs[beta_index]+qnorm(.975)*sqrt(diag(beta_covariance_estparams_estblock))[beta_index]
      CI_exch[2*beta_index-1]=coefs[beta_index]-qnorm(.975)*sqrt(diag(beta_covariance_exch))[beta_index]
      CI_exch[2*beta_index]=coefs[beta_index]+qnorm(.975)*sqrt(diag(beta_covariance_exch))[beta_index]
      CI_dyadic[2*beta_index-1]=coefs[beta_index]-qnorm(.975)*sqrt(diag(beta_covariance_dyadic))[beta_index]
      CI_dyadic[2*beta_index]=coefs[beta_index]+qnorm(.975)*sqrt(diag(beta_covariance_dyadic))[beta_index]
    }
    CI_trueparams_record=rbind(CI_trueparams_record,CI_trueparams)
    CI_estparams_trueblock_record=rbind(CI_estparams_trueblock_record,CI_estparams_trueblock)
    CI_estparams_estblock_record=rbind(CI_estparams_estblock_record,CI_estparams_estblock)
    CI_exch_record=rbind(CI_exch_record,CI_exch)
    CI_dyadic_record=rbind(CI_dyadic_record,CI_dyadic)
  }
  toc()
  ans=list(X=X,ng=ng,g=g,CI_trueparams_record=CI_trueparams_record,CI_estparams_trueblock_record=CI_estparams_trueblock_record,CI_estparams_estblock_record=CI_estparams_estblock_record,CI_exch_record=CI_exch_record,CI_dyadic_record=CI_dyadic_record)
  return(ans)
}

