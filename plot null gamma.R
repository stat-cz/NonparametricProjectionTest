


rm(list=ls())
dir <- ".../stat-cz/NonparametricProjectionTest"
save_path <- ".../stat-cz/NonparametricProjectionTest/figure"
setwd(dir)

##' vary p and fix rho and sample size
n <- n1 <-  n2 <- 5
p <- c(20, 50, 100, 150, 200, 400)
Nsim <- 10000
B <- factorial(n1+n2) / factorial(n1) / factorial(n2)

## MA model
rho <- rho2 <- 0
(pattern=paste("np_gamma_null_n1",n1,"_n2",n2,"_p",p,"_rho2",rho2,"_B",B,"_Nsim",Nsim,".RData",sep =""))

## dense 
rho <- rho2 <- 0.5
(pattern=paste("np_gamma_null_n1",n1,"_n2",n2,"_p",p,"_rho2",rho2,"_B",B,"_Nsim",Nsim,"_dense.RData",sep =""))





collect_statistics <- function(dir,pattern,n1,n2,pi,rho,Nsim=100000) 
{
  
  ###statistics
  
  files <- dir(dir, pattern = pattern, full.names = TRUE)
  load(files)
  
  
  ##p-value
  
  PT1 <- np.gamma.mean.vector$PT1
  PT2 <- np.gamma.mean.vector$PT2
  PT3 <- np.gamma.mean.vector$PT3
  npbs <- np.gamma.mean.vector$npBS
  npcq <- np.gamma.mean.vector$npCQ
  nppt1 <- np.gamma.mean.vector$npPT1
  nppt2 <- np.gamma.mean.vector$npPT2
  nppt3 <- np.gamma.mean.vector$npPT3
  
  
  stat_P <-  c(PT1,PT2,PT3, npbs, npcq, nppt1, nppt2, nppt3)
  
  ####ppoints(p-value)
  PT1_p <- ppoints(PT1)
  PT2_p <- ppoints(PT2)
  PT3_p <- ppoints(PT3)
  npbs_p <- ppoints(npbs)
  npcq_p <- ppoints(npcq)
  nppt1_p <- ppoints(nppt1)
  nppt2_p <- ppoints(nppt2)
  nppt3_p <- ppoints(nppt3)
  
  stat_P_p <-  c(PT1_p,PT2_p,PT3_p,npbs_p,npcq_p,nppt1_p,nppt2_p,nppt3_p)
  
  ####sort(p-value)
  
  PT1_s <- sort(PT1)
  PT2_s <- sort(PT2)
  PT3_s <- sort(PT3)
  npbs_s <- sort(npbs)
  npcq_s <- sort(npcq)
  nppt1_s <- sort(nppt1)
  nppt2_s <- sort(nppt2)
  nppt3_s <- sort(nppt3)
  
  stat_P_s <-  c(PT1_s,PT2_s,PT3_s,npbs_s,npcq_s,nppt1_s,nppt2_s,nppt3_s)
  
  number <- length(PT1)
  
  method <- as.factor(rep(c("PT1","PT2","PT3","npBS","npCQ","npPT1","npPT2","npPT3"),each=number))
  methods <- as.factor(rep(c("PT1","PT2","PT3","npBS","npCQ","npPT1","npPT2","npPT3"),each=number))
  

  
  method.n <- length(levels(method))
  #setting <- as.factor(rep(paste("c==",c," ","gamma==",gamma," ","kappa==",kappa,sep =""),3000))
  setting <- as.factor(rep(paste("p==",pi,sep =""),(method.n*number)))
  setting_rho <- as.factor(rep(paste("rho==",rho,sep =""),(method.n*number)))
  setting_n <- as.factor(rep(paste("n==",n1,sep =""),(method.n*number)))
  
  color <- as.factor(rep(c(LETTERS[1:method.n]),each=number))
  df <- data.frame(
    statistic_P = stat_P,
    statistic_P_p = stat_P_p,
    statistic_P_s = stat_P_s,
    method = method,
    methods = methods,
    setting = setting,
    setting_rho = setting_rho,
    setting_n = setting_n,
    color = color)
}

##'
##' vary p and fix rho and sample size
##'
##'

df <- NULL
for (i in 1:length(p)) 
{
  stats <- collect_statistics(dir,
                              pattern = pattern[i],n1=n1,n2=n2,
                              pi=p[i],rho=rho,Nsim=Nsim)
  df <- rbind(df,stats)
}



####density plot
require(ggplot2)
# windows()
df$methods = factor(df$methods, levels=c("npPT1","npPT2","npPT3","npBS","npCQ","PT1","PT2","PT3"), ordered = T)


###p-value###

require("RColorBrewer")
WAcol <- brewer.pal(6, "Set2")
require(ggsci)
###point
pvalue <- ggplot(subset(df,  methods %in% c("npPT1","npPT2","npPT3","npBS","npCQ","PT1","PT2","PT3")),
                 aes(x=statistic_P_p,y=statistic_P_s,color = color))+
  scale_color_aaas()+
  geom_abline(intercept = 0, slope = 1,size=0.5,color="darkgray")+
  geom_point(size=0.5)+
  scale_x_continuous(breaks=seq(0, 1, 0.5)) +
  scale_y_continuous(breaks=seq(0, 1, 0.5)) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  labs(x = "Uniform quantiles", y = expression(paste("p", "-value ")))+
  coord_cartesian(y = c(0,1),x=c(0,1))


##' 
##' various p and fixed sample size n=5 with MA model Sigma
##' 
pvalue.plot <- pvalue+facet_grid(setting~methods,scales = "free",labeller = label_parsed)
pvalue.plot
setwd(save_path)
ggsave(file=paste("pvalue_point_gamma_non_mvtest_n1",n1,"_n2",n2,".png",sep=""),plot=pvalue.plot,width=10,height=10)
ggsave(file=paste("pvalue_point_gamma_non_mvtest_n1",n1,"_n2",n2,".pdf",sep=""),plot=pvalue.plot,width=10,height=10)

