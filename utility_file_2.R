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
  n=length(g)
  last_nodes_len=(n-1)*(n-2)/2
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
  nodes.4.list=list()
  nodes.5.list=list()
  for (i in 1:K){
    sub_i=c(sapply(which(g==i),FUN=function(ind){((ind-1)*last_nodes_len+1):(ind*last_nodes_len)}))
    nodes.3_sub=nodes.3[sub_i,]
    nodes.4_sub=nodes.4[sub_i,]
    nodes.5_sub=nodes.5[sub_i,]
    nodes.6_sub=nodes.6[sub_i,]
    for (j in 1:K){
      for (k in j:K){
        if(j==k){
          nodes.3.list=c(nodes.3.list,list(nodes.3_sub[which( g[nodes.3_sub[,2]]==j & g[nodes.3_sub[,4]]==k),]))
          nodes.4.list=c(nodes.4.list,list(nodes.4_sub[which(g[nodes.4_sub[,1]]==j & g[nodes.4_sub[,3]]==k),]))
          nodes.5.list=c(nodes.5.list,list(rbind(nodes.5_sub[which( g[nodes.5_sub[,2]]==j & g[nodes.5_sub[,3]]==k),],
                                                 nodes.6_sub[which( g[nodes.6_sub[,1]]==j & g[nodes.6_sub[,4]]==k),])))          
        }
        else{
          nodes.3.list=c(nodes.3.list,list(rbind(nodes.3_sub[which( g[nodes.3_sub[,2]]==j & g[nodes.3_sub[,4]]==k),],
                                                 nodes.3_sub[which(g[nodes.3_sub[,2]]==k & g[nodes.3_sub[,4]]==j),])))
          nodes.4.list=c(nodes.4.list,list(rbind(nodes.4_sub[which( g[nodes.4_sub[,1]]==j & g[nodes.4_sub[,3]]==k),],
                                                 nodes.4_sub[which( g[nodes.4_sub[,1]]==k & g[nodes.4_sub[,3]]==j),])))
          nodes.5.list=c(nodes.5.list,list(rbind(nodes.5_sub[which( g[nodes.5_sub[,2]]==j & g[nodes.5_sub[,3]]==k),],
                                                 nodes.5_sub[which( g[nodes.5_sub[,2]]==k & g[nodes.5_sub[,3]]==j),],
                                                 nodes.6_sub[which( g[nodes.6_sub[,1]]==j & g[nodes.6_sub[,4]]==k),],
                                                 nodes.6_sub[which( g[nodes.6_sub[,1]]==k & g[nodes.6_sub[,4]]==j),])))
        }
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
# cross.prod.one=function(resCrossProd,nodeSetExch,node.ind){
#   ans=vector("list", 6) 
#   ans[[1]]=resCrossProd[[1]][c(which(nodeSetExch[[1]][,1]==node.ind),which(nodeSetExch[[1]][,2]==node.ind))]
#   ans[[2]]=resCrossProd[[2]][c(which(nodeSetExch[[2]][,1]==node.ind),which(nodeSetExch[[2]][,2]==node.ind))]
#   ans[[3]]=resCrossProd[[3]][which(nodeSetExch[[3]][,1]==node.ind)]
#   ans[[4]]=resCrossProd[[4]][which(nodeSetExch[[4]][,2]==node.ind)]
#   ans[[5]]=resCrossProd[[5]][which(nodeSetExch[[5]][,1]==node.ind)]
#   ans[[6]]=resCrossProd[[6]][which(nodeSetExch[[6]][,2]==node.ind)]
#   return(ans)
# }

#get ks statistic without calculating p-value etc.
get_ks_stat_manual=function(x,y){
  #x=quantile(x,seq(0.01,0.99,0.01))
  #y=quantile(y,seq(0.01,0.99,0.01))
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
      ks.record[k]=get_ks_stat_manual(resCrossProdInd[[k]][i,],resCrossProdInd[[k]][j,])
    }
    ks.record[5]=get_ks_stat_manual(c(resCrossProdInd[[5]][i,],resCrossProdInd[[6]][i,]),c(resCrossProdInd[[5]][j,],resCrossProdInd[[6]][j,]))
    similarity.ks.mean=1-mean(ks.record)
    ans=similarity.ks.mean
  }
  return (ans)
}, vectorize.args=c("i", "j"))

# get_avg_ks=function(i,j,resCrossProdInd,get_ks_stat_manual){
#   ans=0
#   #if(i>j){
#     ks.record=numeric(5)
#     for (k in 1:4){
#       ks.record[k]=get_ks_stat_manual(resCrossProdInd[[k]][i,],resCrossProdInd[[k]][j,])
#     }
#     ks.record[5]=get_ks_stat_manual(c(resCrossProdInd[[5]][i,],resCrossProdInd[[6]][i,]),c(resCrossProdInd[[5]][j,],resCrossProdInd[[6]][j,]))
#     similarity.ks.mean=1-mean(ks.record)
#     ans=similarity.ks.mean
#   #}
#   return (ans)
# }

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
  resCrossProdInd=vector("list", 6) 
  # resCrossProdInd[[1]]=matrix(NA,nrow=n,ncol=2*(n-1))
  # resCrossProdInd[[2]]=matrix(NA,nrow=n,ncol=n-1)
  # resCrossProdInd[[3]]=matrix(NA,nrow=n,ncol=(n-1)*(n-2)/2)
  # resCrossProdInd[[4]]=matrix(NA,nrow=n,ncol=(n-1)*(n-2)/2)
  # resCrossProdInd[[5]]=matrix(NA,nrow=n,ncol=(n-1)*(n-2)/2)
  # resCrossProdInd[[6]]=matrix(NA,nrow=n,ncol=(n-1)*(n-2)/2)
  quantile_prob_len=200
  quantile_prob=seq(0.001,0.999,(0.999-0.001)/(quantile_prob_len-1))
  resCrossProdInd[[1]]=matrix(NA,nrow=n,ncol=2*(n-1))
  resCrossProdInd[[2]]=matrix(NA,nrow=n,ncol=n-1)
  resCrossProdInd[[3]]=matrix(NA,nrow=n,ncol=quantile_prob_len)
  resCrossProdInd[[4]]=matrix(NA,nrow=n,ncol=quantile_prob_len)
  resCrossProdInd[[5]]=matrix(NA,nrow=n,ncol=quantile_prob_len)
  resCrossProdInd[[6]]=matrix(NA,nrow=n,ncol=quantile_prob_len)
  last_four_len=(n-1)*(n-2)/2
  for (node.ind in 1:n){
    last_four_indices=((node.ind-1)*last_four_len+1):((node.ind)*last_four_len)
    resCrossProdInd[[1]][node.ind,]=resCrossProd[[1]][c(which(nodeSetExch6[[1]][,1]==node.ind),which(nodeSetExch6[[1]][,2]==node.ind))]
    resCrossProdInd[[2]][node.ind,]=resCrossProd[[2]][c(which(nodeSetExch6[[2]][,1]==node.ind),which(nodeSetExch6[[2]][,2]==node.ind))]
    resCrossProdInd[[3]][node.ind,]=quantile(resCrossProd[[3]][last_four_indices],quantile_prob)
    resCrossProdInd[[4]][node.ind,]=quantile(resCrossProd[[4]][last_four_indices],quantile_prob)
    resCrossProdInd[[5]][node.ind,]=quantile(resCrossProd[[5]][last_four_indices],quantile_prob)
    resCrossProdInd[[6]][node.ind,]=quantile(resCrossProd[[6]][last_four_indices],quantile_prob)
  }
  
  toc()
  tic("get similarity matrix")
  #cl <- makeCluster(1)
  #similarity.ks.mean=parSapply(cl,1:n,1:n,FUN=get_avg_ks,resCrossProdInd=resCrossProdInd,get_ks_stat_manual=get_ks_stat_manual)
  similarity.ks.mean=outer(1:n,1:n,FUN=get_avg_ks,resCrossProdInd=resCrossProdInd,get_ks_stat_manual=get_ks_stat_manual)
  #similarity.ks.mean.vec=mapply(FUN=get_avg_ks,i=nodeSetExch6[[2]][,1],j=nodeSetExch6[[2]][,2],MoreArgs=list(resCrossProdInd=resCrossProdInd,get_ks_stat_manual=get_ks_stat_manual))
  #stopCluster(cl)
  #similarity.ks.mean=matrix(0,nrow=n,ncol=n)
  #similarity.ks.mean[nodeSetExch6[[2]][,c(1,2)]]=similarity.ks.mean.vec
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


# get covariance matrix for beta hat with Omega matrix as input
get_beta_covariance=function(X,sigmaE){
  return(chol2inv(chol(crossprod(X)))%*%t(X)%*%sigmaE%*%X%*%chol2inv(chol(crossprod(X))))
}

vec.net <- function(A, directed=T)
{
  # convert matrix A to vector
  # unfold along columns, omitting diagonal
  # return n(n-1) vector or n(n-1)/2 vector depending on directed option input
  # Uses lower triangular unfolding for undirected case
  
  n <- nrow(A)
  if (directed==T){ 
    vec.out <- as.vector(A)[-seq(1,n^2,n+1)] 
  } else if (directed==F) {
    vec.out <- A[lower.tri(A)]
  } else {stop('Logical input only for directed input')}
  return(vec.out)
}

# Matricize a network vector (without diagonal)
mat.net <- function(V, directed=T)
{
  # convert vector V to a matrix
  # build along columns, omitting diagonal
  # return n(n-1) vector or n(n-1)/2 vector depending on directed option input
  # Uses lower triangular unfolding for undirected case
  
  if (directed==T){   
    d <- length(as.vector(V))    
    n <- floor((1+sqrt(1+4*d))/2)
    
    A <- matrix(1:n^2,n,n)
    ind <- vec.net(A, directed)            # indices in matrix
    Mat.out <- matrix(0,n,n)
    Mat.out[ind] <- V
  } else if (directed==F){
    d <- length(as.vector(V))    
    n <- floor((1+sqrt(1+8*d))/2)
    
    Mat.out <- matrix(0,n,n)
    Mat.out[lower.tri(Mat.out)] <- V
    Mat.out <- Mat.out + t(Mat.out)
  } else {stop('Logical only for directed input')}
  return(Mat.out)
}

sandwich.var <- function(bread, meat) 
{
  # Calculate estimate of variance using input bread and meat matrices
  V.out <- bread %*% meat %*% bread
  
  return(V.out)
}

meat.S <- function(node.list,X.1,e.1,K=1)
{
  # Build meat for structured/exchangeable sandwich standard error estimator
  # Takes in design matrix X.1, residuals e.1, and list of overlapping nodes 
  
  # sizes
  d <- nrow(X.1)
  n <- floor((1+sqrt(1+4*d))/2)
  p <- ncol(X.1)
  e.mat <- mat.net(e.1)
  
  # Estimates
  phi <- rep(0,length(node.list))
  meat.1 <- matrix(0,p,p)
  for (i in 1:K^2){
    nodes.tmp <- node.list[[i]]
    phi[i] <- mean(e.mat[ nodes.tmp[,1:2]]*e.mat[nodes.tmp[,3:4]])
    d1 <- dyad(nodes.tmp[,1], nodes.tmp[,2], n)
    d2 <- dyad(nodes.tmp[,3], nodes.tmp[,4], n)
    meat.1 <- meat.1 + crossprod(X.1[d1,], X.1[d2,])*phi[i]  
  }
  
  meat.2 <- matrix(0,p,p)
  for (i in (K^2+1):length(node.list)){
    nodes.tmp <- node.list[[i]]
    phi[i] <- mean(e.mat[ nodes.tmp[,1:2]]*e.mat[nodes.tmp[,3:4]])
    d1 <- dyad(nodes.tmp[,1], nodes.tmp[,2], n)
    d2 <- dyad(nodes.tmp[,3], nodes.tmp[,4], n)
    meat.2 <- meat.2 + crossprod(X.1[d1,], X.1[d2,])*phi[i]  + crossprod(X.1[d2,], X.1[d1,])*phi[i]
  }
  
  # params <- c(phi,0)
  # S.list <- Sigma.ind(n, directed=TRUE)    # list of indicator matrices, 6 parameter  
  # C = Reduce("+", lapply(1:length(S.list), function(x) params[x]*S.list[[x]]))
  # if(rankMatrix(C,method = "qr")[1]<nrow(C)){sing_flag<-TRUE}else{sing_flag<-FALSE}
  # if(det(C)<0 && sing_flag==FALSE){negdef_flag<-TRUE}else{negdef_flag<-FALSE}
  
  meat.out <- meat.2 + meat.1
  # return(list(M=meat.out,sing_flag = sing_flag, negdef_flag = negdef_flag))
  return(list(M=meat.out,sing_flag = NA, negdef_flag = NA))
  
}

meat.U <- function(node.list,X.1,e.1)
{
  # Build meat for unstructured/DC sandwich standard error estimator
  # Takes in design matrix X.1, residuals e.1, and list of overlapping nodes 
  
  # sizes
  d <- nrow(X.1)
  n <- floor((1+sqrt(1+4*d))/2)
  p <- ncol(X.1)
  
  # Estimates
  X.e <- X.1*tcrossprod(e.1, rep(1,p))
  meat.1 <- crossprod(X.e)
  
  meat.2 <- matrix(0,p,p)
  for (i in 2:length(node.list)){
    nodes.tmp <- node.list[[i]]
    d1 <- dyad(nodes.tmp[,1], nodes.tmp[,2], n)
    d2 <- dyad(nodes.tmp[,3], nodes.tmp[,4], n)
    meat.2 <- meat.2 + crossprod(X.e[d1,], X.e[d2,]) + crossprod(X.e[d2,], X.e[d1,])
  }
  
  # S.list <- Sigma.ind(n, directed=TRUE)    # list of indicator matrices, 6 parameter  
  # C_DC = outer(c(e.1),c(e.1),'*')
  # C_DC <- C_DC*(1-S.list[[6]])
  # if(rankMatrix(C_DC,method = "qr")[1]<nrow(C_DC)){sing_flag<-TRUE}else{sing_flag<-FALSE}
  # if(det(C_DC)<0 && sing_flag==FALSE){negdef_flag<-TRUE}else{negdef_flag<-FALSE}
  
  meat.out <- meat.2 + meat.1
  # return(list(M=meat.out,sing_flag=sing_flag,negdef_flag=negdef_flag))
  return(list(M=meat.out,sing_flag=NA,negdef_flag=NA))
}


