
rm(list = ls())

setwd(".../stat-cz/NonparametricProjectionTest")
source("np_HDMV_function.R")

##############################
### example1 : simulation data
###############################
require(mvtnorm)
X1 <- rmvnorm(20,mean=rep(0,20),sigma=diag(20))
X2 <- rmvnorm(20,mean=rep(0.5,20),sigma=diag(20))

(results.same <- np.HDmv.same(X1,X2))


### zhizaoye

##########################################################################################
################################  zhizaoye      ########################################
###the first(may) sample X1(p*n)
zhizaoye_May <- read.csv("zhizaoyemay.csv",header = T) #????ҵ
zhizaoye_May_r <- zhizaoye_May$returns##############################
zhizaoye_May_p <- length(unique(zhizaoye_May$company))
zhizaoye_May_n <- length(unique(zhizaoye_May$time))

X1 <- matrix(zhizaoye_May_r,nrow = zhizaoye_May_p,ncol = zhizaoye_May_n,byrow = T) ##X1=p*n
zhizaoye_May_names <- unique(zhizaoye_May$company)
zhizaoye_May_time <- unique(zhizaoye_May$time)
rownames(X1) <- zhizaoye_May_names
colnames(X1) <- zhizaoye_May_time



###the other sample X2(p*n)
zhizaoye_others <- read.csv("zhizaoyeother.csv",header = T) #????ҵ
zhizaoye_others_r <- zhizaoye_others$returns
zhizaoye_others_p <- length(unique(zhizaoye_others$company))
zhizaoye_others_n <- length(unique(zhizaoye_others$time))

X2 <- matrix(zhizaoye_others_r,nrow = zhizaoye_others_p,ncol = zhizaoye_others_n,byrow = T) ##X2=p*n
zhizaoye_others_names <- unique(zhizaoye_others$company)
zhizaoye_others_time <- unique(zhizaoye_others$time)
rownames(X2) <- zhizaoye_others_names
colnames(X2) <- zhizaoye_others_time


(results.same <- np.HDmv.same(t(X1),t(X2)))
