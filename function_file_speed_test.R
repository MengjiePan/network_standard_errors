dyad <- function(i.in, j.in, n.tot, directed=T) 
{
  # transform n indices to dyad indices
  if (directed==T){ 
    dyad.result <- ((i.in-1) + (j.in-1)*(n.tot-1) + (j.in > i.in)) * (i.in != j.in)
  } else if (directed==F) {
    #     if (j.in>i.in){stop('Only lower triangular matrix for undirected network')}
    ij <- cbind(i.in, j.in)
    ij <- t(apply(ij, 1, sort, decreasing=T))
    j.1 <- ij[,2]    # want invariance to ordering
    i.1 <- ij[,1]    # want invariance to ordering
    dyad.result <- ( i.1 - 1 -.5*j.1^2 + (n.tot - 1/2)*j.1 - n.tot + 1  ) * (i.1 != j.1)
  } else {stop('Logical T/F only for directed input')}
  return(dyad.result)
}


node.set <- function(g)
{  
  # Generate node sets of various overlapping dyad pairs
  # Return list of each set of nodes, with null for first set (diagonal)
  
  n.tot=length(g)
  nodes.1 <- cbind(rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)],
                   rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)])
  
  nodes.2 <- nodes.3 <- nodes.4 <- nodes.5 <-nodes.6 <-  c()
  
  for (i in 1:n.tot){ 
    
    # ij,ji
    if (i<n.tot){
      c1 <- rep(i,(n.tot-i))
      #       c2 <- (1:n.tot)[-i]
      c2 <- ((i+1):n.tot)
      
      nodes.2 <- rbind(nodes.2,cbind(c1,c2,c2,c1))
    }
    
    # ij,il  ;  ij,kj
    c1 <- rep(i,(n.tot - 1)*(n.tot-2)/2)
    c2 <- rep( (1:(n.tot-1))[-i], times=(n.tot-2):(1 - 1*(i==n.tot) ) )
    c3 <- c1
    c4.mat <- outer(1:n.tot,rep(1,n.tot))[-i,-i]
    c4 <- c4.mat[lower.tri(diag(n.tot-1))]
    
    nodes.3 <- rbind(nodes.3, cbind(c1,c2,c3,c4))
    nodes.4 <- rbind(nodes.4, cbind(c2,c1,c4,c3))
    
    # ij,jl and ij,ki
    nodes.5 <- rbind(nodes.5, cbind(c1,c2,c4,c3))  
    nodes.6 <- rbind(nodes.6, cbind(c2,c1,c3,c4))  
  }
  
  K=max(g)
  nodes.1.list=list()
  for (i in 1:K){
    for (j in 1:K){
      nodes.1.list=c(nodes.1.list,list(nodes.1[which(g[nodes.1[,1]]==i & g[nodes.1[,2]]==j),]))
    }
  }
  
  nodes.2.list=list()
  for (i in 1:K){
    for (j in i:K){
      if(j==i){nodes.2.list=c(nodes.2.list,list(nodes.2[which(g[nodes.2[,1]]==i & g[nodes.2[,2]]==j),]))}
      else{
        nodes.2.list=c(nodes.2.list,list(rbind(nodes.2[which(g[nodes.2[,1]]==i & g[nodes.2[,2]]==j),],
                                               nodes.2[which(g[nodes.2[,1]]==j & g[nodes.2[,2]]==i),])))}
    }
  }
  
  ##fix for 1 by 4 matrix
  for (i in 1:length(nodes.2.list)){
    if (length(nodes.2.list[[i]])==4){
      nodes.2.list[[i]]=matrix(nodes.2.list[[i]],ncol=4,nrow=1,byrow = F)
    }
  }
  
  nodes.3.list=list()
  for (i in 1:K){
    for (j in 1:K){
      for (k in j:K){
        if(j==k){nodes.3.list=c(nodes.3.list,list(nodes.3[which(g[nodes.3[,1]]==i & g[nodes.3[,2]]==j & g[nodes.3[,4]]==k),]))}
        else{
          nodes.3.list=c(nodes.3.list,list(rbind(nodes.3[which(g[nodes.3[,1]]==i & g[nodes.3[,2]]==j & g[nodes.3[,4]]==k),],
                                                 nodes.3[which(g[nodes.3[,1]]==i & g[nodes.3[,2]]==k & g[nodes.3[,4]]==j),])))}
      }
    }
  }
  
  nodes.4.list=list()
  for (i in 1:K){
    for (j in 1:K){
      for (k in j:K){
        if(j==k){nodes.4.list=c(nodes.4.list,list(nodes.4[which(g[nodes.4[,2]]==i & g[nodes.4[,1]]==j & g[nodes.4[,3]]==k),]))}
        else{
          nodes.4.list=c(nodes.4.list,list(rbind(nodes.4[which(g[nodes.4[,2]]==i & g[nodes.4[,1]]==j & g[nodes.4[,3]]==k),],
                                                 nodes.4[which(g[nodes.4[,2]]==i & g[nodes.4[,1]]==k & g[nodes.4[,3]]==j),])))}
      }
    }
  }
  
  nodes.5.list=list()
  for (i in 1:K){
    for (j in 1:K){
      for (k in j:K){
        if(j==k){
          nodes.5.list=c(nodes.5.list,list(rbind(nodes.5[which(g[nodes.5[,1]]==i & g[nodes.5[,2]]==j & g[nodes.5[,3]]==k),],
                                                 nodes.6[which(g[nodes.6[,2]]==i & g[nodes.6[,1]]==j & g[nodes.6[,4]]==k),])))          
        }else{
          nodes.5.list=c(nodes.5.list,list(rbind(nodes.5[which(g[nodes.5[,1]]==i & g[nodes.5[,2]]==j & g[nodes.5[,3]]==k),],
                                                 nodes.5[which(g[nodes.5[,1]]==i & g[nodes.5[,2]]==k & g[nodes.5[,3]]==j),],
                                                 nodes.6[which(g[nodes.6[,2]]==i & g[nodes.6[,1]]==j & g[nodes.6[,4]]==k),],
                                                 nodes.6[which(g[nodes.6[,2]]==i & g[nodes.6[,1]]==k & g[nodes.6[,4]]==j),])))}
      }
    }
  }
  
  return(c(nodes.1.list,nodes.2.list,nodes.3.list,nodes.4.list,nodes.5.list))
}



node.set.exch <- function(n.tot)
{  
  # Generate node sets of various overlapping dyad pairs
  # Return list of each set of nodes, with null for first set (diagonal)
  
  nodes.1 <- cbind(rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)],
                   rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)])
  
  nodes.2 <- nodes.3 <- nodes.4 <- nodes.5 <-nodes.6 <-  c()
  
  for (i in 1:n.tot){ 
    
    # ij,ji
    c1 <- rep(i,(n.tot-1))
    #       c2 <- (1:n.tot)[-i]
    c2 <- c(1:n.tot)[-i]
    
    nodes.2 <- rbind(nodes.2,cbind(c1,c2,c2,c1))
    
    
    # ij,il  ;  ij,kj
    c1 <- rep(i,(n.tot - 1)*(n.tot-2)/2)
    c2 <- rep( (1:(n.tot-1))[-i], times=(n.tot-2):(1 - 1*(i==n.tot) ) )
    c3 <- c1
    c4.mat <- outer(1:n.tot,rep(1,n.tot))[-i,-i]
    c4 <- c4.mat[lower.tri(diag(n.tot-1))]
    
    nodes.3 <- rbind(nodes.3, cbind(c1,c2,c3,c4))
    nodes.4 <- rbind(nodes.4, cbind(c2,c1,c4,c3))
    
    # ij,jl and ij,ki
    nodes.5 <- rbind(nodes.5, cbind(c1,c2,c4,c3))  
    nodes.6 <- rbind(nodes.6,cbind(c2,c1,c3,c4))  
  }
  return(list(nodes.1,nodes.2,nodes.3,nodes.4,nodes.5,nodes.6))
}

getTrueSigma=function(nodeSet,true_sigma,n){
  sigmaE=matrix(0,nrow=n*(n-1),ncol=n*(n-1))
  for (i in 1:length(nodeSet)){
    d1 <- dyad(nodeSet[[i]][,1], nodeSet[[i]][,2], n)
    d2 <- dyad(nodeSet[[i]][,3], nodeSet[[i]][,4], n)
    sigmaE[cbind(d1,d2)]=true_sigma[i]
    sigmaE[cbind(d2,d1)]=true_sigma[i]
  }
  return(sigmaE)
}

getTrueSigmaVec=function(sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,K,d,rou){
  ans=NULL
  for (i in 1:K){
    for (j in 1:K){
      ans=c(ans,sigma_a_vec[i]^2+sigma_b_vec[j]^2+d*sigma_z_vec[i]^2*sigma_z_vec[j]^2+sigma_gamma_m[i,j]^2+sigma_epsilon^2)
    }
  }
  for (i in 1:K){
    for (j in i:K){
      ans=c(ans,rou*sigma_a_vec[i]*sigma_b_vec[i]+rou*sigma_a_vec[j]*sigma_b_vec[j]+d*sigma_z_vec[i]^2*sigma_z_vec[j]^2+sigma_gamma_m[i,j]^2)
    }
  }
  ans=c(ans,rep(sigma_a_vec^2,each=K*(K+1)/2))
  ans=c(ans,rep(sigma_b_vec^2,each=K*(K+1)/2))
  ans=c(ans,rep(rou*sigma_a_vec*sigma_b_vec,each=K*(K+1)/2))
  return(ans)
}


make.positive.var <- function(V.test)
{
  eig.V <- eigen(as.matrix(V.test), symmetric=T)
  eig.flag <- sum(eig.V$values < 0)                
  if (eig.flag>0) { 
    eig.V$values[which(eig.V$values<0)] <- 0        
    V.corrected <- eig.V$vectors %*% diag(eig.V$values) %*% t(eig.V$vectors)
  } else { 
    V.corrected <- V.test
  }
  return(list(V=V.corrected, flag=1*(eig.flag>0) ))
}

cross.prod.one=function(resCrossProd,nodeSetExch,node.ind){
  ans=vector("list", 6) 
  ans[[1]]=resCrossProd[[1]][c(which(nodeSetExch[[1]][,1]==node.ind),which(nodeSetExch[[1]][,2]==node.ind))]
  ans[[2]]=resCrossProd[[2]][which(nodeSetExch[[2]][,1]==node.ind)]
  ans[[3]]=resCrossProd[[3]][which(nodeSetExch[[3]][,1]==node.ind)]
  ans[[4]]=resCrossProd[[4]][which(nodeSetExch[[4]][,2]==node.ind)]
  ans[[5]]=resCrossProd[[5]][which(nodeSetExch[[5]][,1]==node.ind)]
  ans[[6]]=resCrossProd[[6]][which(nodeSetExch[[6]][,2]==node.ind)]
  return(ans)
}

make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j]=A[j,i]  <- S[i,j]
      }
    }
  }
  return (A)
}

spectralClustering=function(W,D,K){
  L=D-W
  evL <- eigen(L, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-K+1):ncol(evL$vectors)]
  km <- kmeans(Z, centers=K, nstart=5)
  return(km$cluster)
}

spectralClustering_norm=function(W,D,K){
  n=dim(W)[1]
  L=diag(n)-solve(D)%*%W
  evL <- eigen(L, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-K+1):ncol(evL$vectors)]
  km <- kmeans(Z, centers=K, nstart=5)
  return(km$cluster)
}


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

simulate_g=function(n,ng){
  g_ind=sample(1:n,size=ng,replace=F)
  g=numeric(n)+1
  g[g_ind]=2
  return(g)
}

simulate_X_absdiff=function(n,sd_options,g){
  X_vec=rnorm(n=n,mean=0,sd=sd_options[g])
  X=outer(X_vec, X_vec, FUN = function(x1,x2){return (abs(x1-x2))})
  return(X)
}

simulate_X_normal=function(n,sd_options,g){
  sigma_matrix=sqrt(sd_options[g] %*% t(sd_options[g] ))
  X=rnorm(n=n^2,mean=0,sd=c(sigma_matrix))
  X=matrix(X,nrow=n,ncol=n,byrow=F)
  return(X)
}

simulate_X_bernoulli=function(n,probs,g){
  X_vec=rbinom(n=n,size=1,prob=probs[g])
  X=outer(X_vec, X_vec, FUN = function(x1,x2){return (as.numeric(x1==x2))})
  return(X)
}


simulate_y=function(xi,X,beta){
  n=dim(xi)[1]
  y=matrix(beta[1],nrow=n,ncol=n)+beta[2]*X+xi
  return(y)
}

get_linear_model=function(y,X){
  n=dim(y)[1]
  y_tmp=c(y)[-seq(1,n^2,n+1)]
  X_tmp=c(X)[-seq(1,n^2,n+1)]
  mod=lm(y_tmp~X_tmp)
  residuals=y-mod$coefficients[1]-mod$coefficients[2]*X
  return(list(residuals=residuals,coefs=mod$coefficients))
}

get_sigmaE_true=function(n,g,sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,K,d,rou){
  true_sigma=getTrueSigmaVec(sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,K,d,rou)
  sigmaE=getTrueSigma(node.set(g),true_sigma,n)
  return(sigmaE)
}

get_sigmaE_est=function(n,nodeSet,residuals){
  est_sigmaE=matrix(0,nrow=n*(n-1),ncol=n*(n-1))
  for (i in 1:length(nodeSet)){
    d1 <- dyad(nodeSet[[i]][,1], nodeSet[[i]][,2], n)
    d2 <- dyad(nodeSet[[i]][,3], nodeSet[[i]][,4], n)
    est_sigma=mean(residuals[cbind(nodeSet[[i]][,1],nodeSet[[i]][,2])]*residuals[cbind(nodeSet[[i]][,3],nodeSet[[i]][,4])])
    est_sigmaE[cbind(d1,d2)]=est_sigma
    est_sigmaE[cbind(d2,d1)]=est_sigma
  }
  return(est_sigmaE)
}

get_est_blocks=function(n,residuals,nodeSetExch6,K,normalized=F){
  tic("get residual products for each dyad")
  #get residual products for each dyad
  #6 dyads because the fifth dyad is separated into two cases (ij,jk) and (jk,ij)
  resCrossProd=vector("list", 6) 
  for (i in 1:6){
    resCrossProd[[i]]=residuals[nodeSetExch6[[i]][,c(1,2)]]*residuals[nodeSetExch6[[i]][,c(3,4)]]
  }
  toc()
  tic("get residuals products for each actor i")
  #get residuals products for each actor i
  resCrossProdInd=vector("list", n) 
  for (i in 1:n){
    resCrossProdInd[[i]]=cross.prod.one(resCrossProd,nodeSetExch6,i)
  }
  toc()
  tic("get similarity matrix")
  #loop through eachpair of actors i and j, and calculate 5 ks statistic and take the average
  similarity.ks.mean=matrix(NA,nrow=n,ncol=n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      ks.record=numeric(5)
      for (k in 1:4){
        ks.record[k]=ks.test(resCrossProdInd[[i]][[k]],resCrossProdInd[[j]][[k]])$statistic
      }
      ks.record[5]=ks.test(c(resCrossProdInd[[i]][[5]],resCrossProdInd[[i]][[6]]),c(resCrossProdInd[[j]][[5]],resCrossProdInd[[j]][[6]]))$statistic
      similarity.ks.mean[i,j] =similarity.ks.mean[j,i] <- 1-mean(ks.record)
    }
  }
  diag(similarity.ks.mean)=1
  toc()
  tic("perform spectral clustering")
  #perform spectral clustering
  S=make.affinity(similarity.ks.mean,n*0.2)
  D=diag(rowSums(S))
  if (normalized){
    est_blocks=spectralClustering(S,D,K)
  }else{
    est_blocks=spectralClustering_norm(S,D,K)
  }
  toc()
  return(est_blocks)
}



get_sigmaE_dyadic=function(n,nodeSetExch5,residuals){
  est_sigmaE_dyadic=matrix(0,nrow=n*(n-1),ncol=n*(n-1))
  nodeSetDyad=rbind(nodeSetExch5[[1]],nodeSetExch5[[2]],nodeSetExch5[[3]],nodeSetExch5[[4]],nodeSetExch5[[5]])
  d1 <- dyad(nodeSetDyad[,1], nodeSetDyad[,2], n)
  d2 <- dyad(nodeSetDyad[,3], nodeSetDyad[,4], n)
  est_sigmaE_dyadic[cbind(d2,d1)]<-est_sigmaE_dyadic[cbind(d1,d2)]<-residuals[nodeSetDyad[,1:2]]*residuals[nodeSetDyad[,3:4]]
  return(est_sigmaE_dyadic)
}

get_beta_covariance=function(X,sigmaE){
  return(solve(crossprod(X))%*%t(X)%*%sigmaE%*%X%*%solve(crossprod(X)))
}

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
  sigmaE=get_sigmaE_true(n,g,sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,K,d,rou)
  toc()
  beta_covariance_trueparams=get_beta_covariance(X_m,sigmaE)
  tic("get block-exchangeable node set with block oracle")
  nodeSetTrue=node.set(g)
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
    nodeSetEst=node.set(est_blocks)
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

