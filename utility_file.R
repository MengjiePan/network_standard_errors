# transform n indices to dyad indices
# borrowed from old code
dyad <- function(i.in, j.in, n.tot, directed=T) 
{
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

# Generate node sets of various overlapping dyad pairs, with exchangeability assumption
# Return list of each set of nodes, returns 6 lists instead of 5 to make the last case separate
# so that it is easier to get block-exch node set from exch node-set
# modified from old code
node.set.exch <- function(n.tot)
{  
  
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
  return(list(nodes.1,nodes.2,nodes.3,nodes.4,nodes.5,nodes.6))
}

# Generate node sets of various overlapping dyad pairs, with block-exchangeability assumption
# Return list of each set of nodes
# take exch node set and block membership as input
node.set <- function(g,nodeSetExch)
{  
  nodes.1=nodeSetExch[[1]]
  nodes.2=nodeSetExch[[2]]
  nodes.3=nodeSetExch[[3]]
  nodes.4=nodeSetExch[[4]]
  nodes.5=nodeSetExch[[5]]
  nodes.6=nodeSetExch[[6]]
  
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


# make the estimated error covariance matrix positive definite
# by zerong out negative eigenvalues
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

# for a specific node, returns a list of length 6
# where each item is residual products involving that specific node for each dyad
# the fifth and sixth item in the list will be combined later
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

#get ks statistic without calculating p-value etc.
get_ks_stat_manual=function(x,y){
  n.x=length(x)
  n.y=length(y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  return( max(abs(z)))
}

# get a similarity value for a pair of nodes i and j
# the function is to be used with outer function to replace the loop
get_avg_ks=Vectorize(function(i,j,resCrossProdInd,get_ks_stat_manual){
  ans=0
  if(i>j){
    ks.record=numeric(5)
    for (k in 1:4){
      ks.record[k]=get_ks_stat_manual(resCrossProdInd[[i]][[k]],resCrossProdInd[[j]][[k]])
    }
    ks.record[5]=get_ks_stat_manual(c(resCrossProdInd[[i]][[5]],resCrossProdInd[[i]][[6]]),c(resCrossProdInd[[j]][[5]],resCrossProdInd[[j]][[6]]))
    similarity.ks.mean=1-mean(ks.record)
    ans=similarity.ks.mean
  }
  return (ans)
}, vectorize.args=c("i", "j"))

#transform the similarity matrix into a similarity graph
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

#perform unnormalized spectral clustering
spectralClustering=function(W,D,K){
  L=D-W
  evL <- eigen(L, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-K+1):ncol(evL$vectors)]
  km <- kmeans(Z, centers=K, nstart=5)
  return(km$cluster)
}

#perform normalized spectral clustering
spectralClustering_norm=function(W,D,K){
  n=dim(W)[1]
  L=diag(n)-solve(D)%*%W
  evL <- eigen(L, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-K+1):ncol(evL$vectors)]
  km <- kmeans(Z, centers=K, nstart=5)
  return(km$cluster)
}

# get estimated block membership
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
  #cl <- makeCluster(1)
  #similarity.ks.mean=parSapply(cl,1:n,1:n,FUN=get_avg_ks,resCrossProdInd=resCrossProdInd,get_ks_stat_manual=get_ks_stat_manual)
  similarity.ks.mean=outer(1:n,1:n,FUN=get_avg_ks,resCrossProdInd=resCrossProdInd,get_ks_stat_manual=get_ks_stat_manual)
  #stopCluster(cl)
  similarity.ks.mean=t(similarity.ks.mean)+similarity.ks.mean
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

# get linear model for simple linear regression
# returns coefficients and a residual matrix
get_linear_model=function(y,X){
  n=dim(y)[1]
  y_tmp=c(y)[-seq(1,n^2,n+1)]
  X_tmp=c(X)[-seq(1,n^2,n+1)]
  mod=lm(y_tmp~X_tmp)
  residuals=y-mod$coefficients[1]-mod$coefficients[2]*X
  return(list(residuals=residuals,coefs=mod$coefficients))
}

# take parameters in error generative model as input
# return a vector of true parameters in error covariance matrix as output
# the output is used as input in function getTrueSigma
# used to calculate confidence intervals for Omega oracle
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

# take a vector of true parameters in error covariance matrix, and block-exch node set as input
# return Sigma (error covariance matrix) oracle
# used to calculate confidence intervals for Omega oracle
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

# combine the above functions to get Omega oracle
get_sigmaE_true=function(n,g,sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,K,d,rou,nodeSetExch){
  true_sigma=getTrueSigmaVec(sigma_a_vec,sigma_b_vec,sigma_z_vec,sigma_gamma_m,K,d,rou)
  sigmaE=getTrueSigma(node.set(g,nodeSetExch),true_sigma,n)
  return(sigmaE)
}

# get estimated Omega with node set as input
# used to get estimated Omega with exch, block-exch (block oracle and block estimated)
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

#  get estimated Omega for dyadic clustering
get_sigmaE_dyadic=function(n,nodeSetExch5,residuals){
  est_sigmaE_dyadic=matrix(0,nrow=n*(n-1),ncol=n*(n-1))
  nodeSetDyad=rbind(nodeSetExch5[[1]],nodeSetExch5[[2]],nodeSetExch5[[3]],nodeSetExch5[[4]],nodeSetExch5[[5]])
  d1 <- dyad(nodeSetDyad[,1], nodeSetDyad[,2], n)
  d2 <- dyad(nodeSetDyad[,3], nodeSetDyad[,4], n)
  est_sigmaE_dyadic[cbind(d2,d1)]<-est_sigmaE_dyadic[cbind(d1,d2)]<-residuals[nodeSetDyad[,1:2]]*residuals[nodeSetDyad[,3:4]]
  return(est_sigmaE_dyadic)
}

# get covariance matrix for beta hat with Omega matrix as input
get_beta_covariance=function(X,sigmaE){
  return(solve(crossprod(X))%*%t(X)%*%sigmaE%*%X%*%solve(crossprod(X)))
}


