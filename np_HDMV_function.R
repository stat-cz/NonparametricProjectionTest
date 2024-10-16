

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
  
  #Tnew <- t(X1bar - X2bar) %*% D %*% (X1bar - X2bar)
  Tnew <- sum(colSums((X1bar - X2bar) *D) * (X1bar - X2bar))
  temp1 <- trD(X1,D)
  temp2 <-  trD(X2,D)
  trD_sigma <- (temp1[1] + temp2[1])/2
  trD_sigma2 <- (temp1[2] + temp2[2])/2
  Tn <- (Tnew - tau * trD_sigma) / (tau * sqrt(2 * trD_sigma2))
  #Tn_P <- pnorm(Tn, lower.tail = FALSE) * (Tn >= 0)  + pnorm(Tn, lower.tail = TRUE) * (Tn < 0) 
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
  # data.frame(npBS0 = Tbs0,npBS = Tbs, npCQ=Tcq, npPT1 = Tpt1, npPT2 = Tpt2,npPT3 = Tpt3)
  
}


np.HDmv.same <- function(X1,X2,B=10000,alpha.level=0.05){
  ## X1 and X2 are n1*p and n2*p matrix
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1+n2-2
  p <- ncol(X1)
  
  alpha1 <- c(1,rep(0,p-1))
  alpha2 <- c(rep(sqrt(ceiling(p/2))/ceiling(p/2),ceiling(p/2)),rep(0,p-ceiling(p/2)))
  alpha3 <- c(rep(sqrt(p)/p,p))
  k<-sqrt(p/log(p))
  I<-diag(p)
  D1<-I+k*(alpha1%*%t(alpha1))
  D2<-I+k*(alpha2%*%t(alpha2))
  D3<-I+k*(alpha3%*%t(alpha3))
  
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
  
  statistics.obs <- statistics(X1,X2,n1,n2,D1,D2,D3,BS_trS0)
  
  # set.seed(1)
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

  
  if(ppt10 < alpha.level){
    results1 <- list( P_value=ppt10,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results1 <- list(P_value=ppt10,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  if(ppt20 < alpha.level){
    results2 <- list( P_value=ppt20,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results2 <- list( P_value=ppt20,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  if(ppt30 < alpha.level){
    results3 <- list( P_value=ppt30,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results3 <- list( P_value=ppt30,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  
  if(pbs0<alpha.level){
    results_npBS0 <- list( P_value=pbs0,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results_npBS0 <- list( P_value=pbs0,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  if(pbs<alpha.level){
    results_npBS <- list( P_value=pbs,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results_npBS <- list( P_value=pbs,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  if(pcq < alpha.level){
    results_npCQ <- list( P_value=pcq,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results_npCQ <- list( P_value=pcq,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  
  if(ppt1 < alpha.level){
    results1_np <- list( P_value=ppt1,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results1_np <- list(P_value=ppt1,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  if(ppt2 < alpha.level){
    results2_np <- list( P_value=ppt2,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results2_np <- list( P_value=ppt2,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  if(ppt3 < alpha.level){
    results3_np <- list( P_value=ppt3,symbol="<", nominal_level=alpha.level,test_results="Reject null hypothesis ")
  }else{
    results3_np <- list( P_value=ppt3,symbol=">", nominal_level=alpha.level,test_results="Do not reject null hypothesis ")
  }
  
  return(rbind(results1,results2,results3,results_npBS0,results_npBS,results_npCQ,results1_np,results2_np,results3_np))
  
}


