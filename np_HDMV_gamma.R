

rm(list=ls())
gamma.sim.data <- function(mu1, mu2, Gamma,p,n1,n2,k,D1,D2,D3,B,Nsim=100,trace=TRUE,seed=100)
{
  set.seed(seed)
  require(MASS)
  #population mu1 here =0
  
  I<-diag(p)
  alpha <- 4
  beta <- 2
  
  MU1<-matrix(rep(mu1,n1),nrow = p)
  
  #population mu2 =0
  MU2<-matrix(rep(mu2,n2),nrow = p)
  
  I<-diag(p)
  data.sim <- as.list(numeric(Nsim))
  
  
  for (i in 1:Nsim){
    
    #population X1
    U1 <- matrix((rgamma(n1*p,alpha,beta)-alpha/beta)/sqrt(alpha/beta^2),nrow=p,ncol=n1,byrow = T)
    X1 <- t(Gamma%*%U1+MU1)  ##n1*p
    
    #population X2
    U2 <- matrix((rgamma(n2*p,alpha,beta)-alpha/beta)/sqrt(alpha/beta^2),nrow=p,ncol=n2,byrow = T)
    X2 <- t(Gamma%*%U2+MU2)  ##n2*p
    
    data.sim[[i]] <- list(id=i,X1=X1,X2=X2)
  }
  
  attr(data.sim,"k") <- k
  attr(data.sim,"D1") <- D1
  attr(data.sim,"D2") <- D2
  attr(data.sim,"D3") <- D3
  attr(data.sim,"B") <- B
  attr(data.sim,"seed") <- seed
  data.sim
}


gamma.sim.single <- function(datasingle,D1,D2,D3,B,trace)
{ 
  trD = function(x, D){
    # x is n*p matrix
    nn <- nrow(x)
    p <- ncol(x)
    x_sumcol <- matrix(colSums(x), nn, p, byrow = TRUE)
    q1 <- x_sumcol - x
    q2 <- x%*%D
    q3 <- q2%*%t(x)
    diagq3 <- diag(q3)
    
    A1 <- sum(colSums(q2 * q1))
    A2 <- (sum(colSums(q3 * q3)) - (q4 <- sum(diagq3^2)))
    A3 <- (sum(colSums(t(q3) * (q5 <- q1 %*% t(q2)))) -  sum(diagq3 * (q6 <- diag(q5))))
    A4 <- (sum(diagq3 * colSums(t(x) * (D%*%t(x_sumcol)))) - q4)
    A5 <- sum(q6 * colSums(t(q1) * (D%*%t(q1))))
    
    temp1 <- (nn - 1) * nn
    temp2 <- (nn - 2) * temp1
    temp3 <- (nn - 3) * temp2
    
    trD_est <- ((q7 <- sum(diagq3))/nn - A1/temp1)
    tr2D_est <- (A2/temp1 - (A3 - A4) * 2/temp2 + (A5 - 2 * A3 + 3 * A4 - A1 * q7)/temp3)
    
    traceD <- c(trD_est = trD_est, tr2D_est = tr2D_est)
    
    traceD
  }
  
  HDmv.same <- function(X1, X2, D){
    
    ## X1 and X2 are n1*p and n2*p matrix
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    tau <- (n1 + n2)/(n1 * n2)
    X1bar <- colMeans(X1)
    X2bar <- colMeans(X2)
    
    ##' 
    ##' proposed projection test for D
    ##' 
    
    Tnew <- sum(colSums((X1bar - X2bar) *D) * (X1bar - X2bar))
    temp1 <- trD(X1,D)
    temp2 <-  trD(X2,D)
    trD_sigma <- (temp1[1] + temp2[1])/2
    trD_sigma2 <- (temp1[2] + temp2[2])/2
    Tn <- (Tnew - tau * trD_sigma) / (tau * sqrt(2 * trD_sigma2))
    pnorm(Tn,lower.tail = FALSE)
    # Tn
  }
  statistics <- function(X1,X2,n1,n2,D1,D2,D3,BS_trS0){
    n <- n1+n2-2
    tau <- (n1+n2)/(n1*n2)
    X <- rbind(X1, X2)
    group <- rep(c(1,2),times=c(n1,n2))
    X1bar <- colMeans(X[group==1,]) 
    X2bar <- colMeans(X[group==2,]) 
    
    Xbardiff <- X1bar-X2bar
    
    
    #proposed projection test
    Tpt1 <- sum(colSums(Xbardiff *D1) * Xbardiff)
    Tpt2 <- sum(colSums(Xbardiff *D2) * Xbardiff)
    Tpt3 <- sum(colSums(Xbardiff *D3) * Xbardiff)
    
    
    
    #Chen and Qin (2010)
    tcq1 <- (X1%*%t(X1))
    diag(tcq1) <- 0
    tcq2 <- (X2%*%t(X2))
    diag(tcq2) <- 0
    tcq3 <- X1%*%t(X2)
    Tcq <- sum(tcq1)/(n1*(n1-1))+sum(tcq2)/(n2*(n2-1))-(2/(n1*n2))*sum(tcq3)
    
    
    # Bai and Saranadasa (1996)
    BS_A1 <- (n1-1)*var(X1)
    BS_A2 <- (n2-1)*var(X2)
    BS_S <- (BS_A1+BS_A2)/n
    BS_trS <- sum(diag(BS_S))
    Tbs0 <- sum(Xbardiff^2) - tau*BS_trS0
    Tbs <-  sum(Xbardiff^2)-tau*BS_trS
    
    c(Tbs0,Tbs,Tcq,Tpt1,Tpt2,Tpt3)
    
  }
  
  if (trace) print(datasingle$id)
  
  X1 <- datasingle$X1 ###n1*p
  X2 <- datasingle$X2  ##n2*p
  X <- rbind(X1, X2)
  
  ppt10 <- as.numeric(HDmv.same(X1,X2,D1))
  ppt20 <- as.numeric(HDmv.same(X1,X2,D2))
  ppt30 <- as.numeric(HDmv.same(X1,X2,D3))
  
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1+n2-2
  group <- rep(c(1,2),times=c(n1,n2))
  
  BS_A1 <- (n1-1)*var(X1)
  BS_A2 <- (n2-1)*var(X2)
  BS_S <- (BS_A1+BS_A2)/n
  BS_trS0 <- sum(diag(BS_S))
  
  # observed statistics
  statistics.obs <- statistics(X1,X2,n1,n2,D1,D2,D3,BS_trS0)
  
  #permutaiton
  statistics.star <- ( sapply(1:B, function(x){
    set.seed(x)
    ind = sample(x= (1:(n1+n2)), size= n1+n2, replace = FALSE)
    Xs = X[ind,]
    Xs1 <- Xs[group==1,]
    Xs2 <- Xs[group==2,]
    statistics(Xs1,Xs2,n1,n2,D1,D2,D3,BS_trS0)
  }) )
  
  pbs0 <- mean(drop(statistics.star[1,])>=drop(statistics.obs[1]))
  pbs <- mean(drop(statistics.star[2,])>=drop(statistics.obs[2]))
  pcq <- mean(drop(statistics.star[3,])>=drop(statistics.obs[3]))
  ppt1 <- mean(drop(statistics.star[4,])>=drop(statistics.obs[4]))
  ppt2 <- mean(drop(statistics.star[5,])>=drop(statistics.obs[5]))
  ppt3 <- mean(drop(statistics.star[6,])>=drop(statistics.obs[6]))

  #### result
  data.frame(PT1=ppt10,PT2=ppt20,PT3=ppt30,npBS0 = pbs0,npBS = pbs, npCQ=pcq, npPT1 = ppt1, npPT2 = ppt2,npPT3 = ppt3)
  
}



gamma.sim.all <- function(dataall,cores=1,trace=FALSE)
{
  require(plyr)
  #require(doMC)
  #registerDoMC(cores)
  require(foreach)
  require(doParallel)
  registerDoParallel(cores=cores)
  
  # cl <- makeCluster(cores)
  # registerDoParallel(cores=cl)
  
  
  D1 <- attributes(dataall)$D1
  D2 <- attributes(dataall)$D2
  D3 <- attributes(dataall)$D3
  B <- attributes(dataall)$B
  
  
  out <- ldply(.data=dataall,.fun=gamma.sim.single,
               D1=D1,D2=D2,D3=D3,B=B,trace=trace,.parallel=TRUE)
  out
}


##############################
###' update 2024-9-28
###' example running the code
#############################
Nsim <- 10000
n1 <- 5
n2 <- 5
p <- 400
n <- n1+n2-2 
B <- factorial(n1+n2)/factorial(n1)/factorial(n2)
alpha1 <- c(1,rep(0,p-1))
alpha2 <- c(rep(sqrt(p/2)/(p/2),p/2),rep(0,p-p/2))
alpha3 <- c(rep(sqrt(p)/p,p))


k<-sqrt(p/log(p))
I<-diag(p)
D1 <- I+k*(alpha1%*%t(alpha1))
D2 <- I+k*(alpha2%*%t(alpha2))
D3 <- I+k*(alpha3%*%t(alpha3))

#mu
mu1 <- rep(0,p)
mu2 <- rep(0,p)
###sigma with MA model covariance matrix structure
Sigma <- matrix(rep(0,p*p),nrow=p)
rho <- 0
Sigma <- rho^abs(col(Sigma)-row(Sigma))

### Sigma with dense covariance matrices
rho <- 0.5
Sigma <- (1-rho) * as.matrix(diag(rep(1,p))) + rho * (rep(1,p)%*% t(rep(1,p)))

##gamma
decomp <- eigen(Sigma)
value <- decomp$values
vectors <- decomp$vectors
Gamma <- vectors%*%diag(sqrt(value))%*%t(vectors)

datasim <- gamma.sim.data(mu1=mu1, mu2=mu2, Gamma = Gamma,
                          p=p,n1=n1,n2=n2,
                          k=k,D1=D1,D2=D2,D3=D3,
                          B= B,Nsim=Nsim,seed=100)


system.time({np.gamma.mean.vector <- gamma.sim.all(dataall = datasim,cores = 5,
                                                   trace = TRUE)})


(filename=paste("np_gamma_null_n1",n1,"_n2",n2,"_p",p,"_rho2",rho,"_B",B,"_Nsim",Nsim,"_dense.RData",sep =""))
save(np.gamma.mean.vector=np.gamma.mean.vector,file=filename)


### empirical type I error

size.plot <- function(out,level=0.05){
  ###############################
  ##corrected significance level
  #############################
  
  BS_s <-  mean( out$npBS < level)
  CQ_s <- mean(out$npCQ < level)
  NPT1_s <- mean(out$npPT1 < level)
  NPT2_s <- mean(out$npPT2 < level)
  NPT3_s <- mean(out$npPT3 < level)
  PT1_s <- mean(out$PT1 < level)
  PT2_s <- mean(out$PT2 < level)
  PT3_s <- mean(out$PT3 < level)
  
  
  
  BS_c <-  quantile( out$npBS, probs = level)
  CQ_c <- quantile(out$npCQ, probs = level)
  NPT1_c <- quantile(out$npPT1, probs = level)
  NPT2_c <- quantile(out$npPT2, probs = level)
  NPT3_c <- quantile(out$npPT3, probs = level)
  PT1_c <- quantile(out$PT1, probs = level)
  PT2_c <- quantile(out$PT2, probs = level)
  PT3_c <- quantile(out$PT3, probs = level)
  
  
  
  #result
  corrected.alpha <- cbind(NPT1_c, NPT2_c,NPT3_c,BS_c,CQ_c,PT1_c, PT2_c,PT3_c)
  size <- cbind(NPT1_s, NPT2_s,NPT3_s,BS_s,CQ_s,PT1_s, PT2_s,PT3_s)
  
  rownames(corrected.alpha) <- c("corrected.alpha")
  rownames(size) <- c("size")
  colnames(size) <- c("npPT1","npPT2","npPT3","npBS","npCQ","PT1","PT2","PT3")
  ###combine the result
  rbind( size,corrected.alpha)
}

size.plot(out = np.gamma.mean.vector)
