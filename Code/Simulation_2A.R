#Load necessary libraries
library(genius)
library(AER)

#Set random number gen
set.seed(12345)

N<-100000

Z_number<-6

GZ_number<-6

sampleGZ<-1:Z_number

iterations<-100

it.GxEB<-matrix(rep(0, len=Z_number*iterations), nrow = Z_number)
it.GxEse<-matrix(rep(0, len=Z_number*iterations), nrow = Z_number)
it.GxEF<-matrix(rep(0, len=Z_number*iterations), nrow = Z_number)
it.GxEFp<-matrix(rep(0, len=Z_number*iterations), nrow = Z_number)

#Generate empty vectors for estimates
Obs_beta<-rep(0,iterations)
GENIUS_beta<-rep(0,iterations)
GENIUS_se<-rep(0,iterations)
bpstat<-rep(0,iterations)
bp_pval<-rep(0,iterations)

for(i in 1:iterations){
  
  print(i)
  
  gamma3<-c(0,0.018,0.0265,0.041,0.0588,0.08242)
  
  #Generate single continuous instrument G
  G<-matrix(rep(0, len=N), nrow = N)
  G[,1]<-rnorm(N,0,1)
  
  #Generate interaction covariates Z
  Z<-matrix(rep(0, len=N*Z_number), nrow = N)
  
  for(k in 1:Z_number){
    Z[,k]<-rnorm(N,0,1)
  }
  
  #Generate GZ interaction variable
  GZ<-matrix(rep(0, len=N*Z_number), nrow = N)
  
  for(k in 1:Z_number){
    GZ[,k]<-G[,1] * Z[,k]
  }
  
  #Generate single continuous confounder U
  U<-rnorm(N,0,1)
  
  #Generate single continuous exposure X
  X<-matrix(rep(0, len=N), nrow = N)
  X[,1]<- 1 + 1*G[,1] + 1*U + rnorm(N,0,1)
  for(k in 1:Z_number){
    X[,1]<- X[,1] + 1*Z[,k] + gamma3[k]*GZ[,k]
  }
  
  #Generate single continuous outcome Y
  Y<- 1 + 1*G[,1] + 1*X + 1*U + rnorm(N,0,1)
  for(k in 1:Z_number){
    Y[,1]<- Y[,1] + 0*Z[,k] + 0*GZ[,k]
  }
  
  #OLS estimate using only exposure
  Obs_beta[i]<-summary(lm(Y~X))$coef[2,1]
  
  #Specify MR-GENIUS model (interaction agnostic)
  GENIUS_1<-genius_addY(c(Y), c(X), c(G), formula = c(X) ~ c(G), alpha = 0.05, lower = -10,
                        upper = 10)
  
  #MR-GENIUS estimate
  GENIUS_beta[i]<-GENIUS_1$beta.est
  
  #MR-GENIUS standard error
  GENIUS_se[i]<-sqrt(GENIUS_1$beta.var)
  
  #Het.test for MR-GENIUS
  bpstat[i]<-bptest(X~G)$statistic
  bp_pval[i]<-bptest(X~G)$p.value
  
  #GxE results matrix
  GxEB<-matrix(rep(0, len=Z_number*4), nrow = Z_number)
  for(j in 1:Z_number){
    
    t.model<-summary(ivreg(Y~X+G+Z[,j]|G+Z[,j]+GZ[,j]),diagnostics=T)
    
    GxEB[j,1]<-t.model$coef[2,1]
    GxEB[j,2]<-t.model$coef[2,2]
    GxEB[j,3]<-t.model$diagnostics[1,3]
    GxEB[j,4]<-t.model$diagnostics[1,4]
  }
  
  for(k in 1:Z_number){
    
    it.GxEB[,i]<-GxEB[,1]
    it.GxEse[,i]<-GxEB[,2]
    it.GxEF[,i]<-GxEB[,3]
    it.GxEFp[,i]<-GxEB[,4]
    
  }
  
}

#MR-GENIUS estimate
mean(GENIUS_beta)

#MR-GENIUS standard error
mean(GENIUS_se)

#Het.test for MR-GENIUS
mean(bpstat)
mean(bp_pval)

#Avg.GxE values

res<-matrix(rep(0, len=Z_number*4), nrow = Z_number)

for(k in 1:Z_number){
  
  res[k,1]<-mean(it.GxEB[k,])
  res[k,2]<-mean(it.GxEse[k,])
  res[k,3]<-mean(it.GxEF[k,])
  res[k,4]<-mean(it.GxEFp[k,])
  
}



res<-data.frame(res)
names(res)<-c("Mean_GxEbeta","Mean_GxEse","Mean_GxEF","Mean_GxEFpvalue")

#Plot results

res$ge_lci<-res$Mean_GxEbeta-(1.96*res$Mean_GxEse)

res$ge_uci<-res$Mean_GxEbeta+(1.96*res$Mean_GxEse)

res$index<-1:6

#Plot results

sim2plot1<-ggplot(res)+geom_point(aes(x=index, y=Mean_GxEbeta),size=1)+
  coord_flip()+geom_errorbar(aes(ymin=ge_lci,ymax=ge_uci,x=index),width=0,orientation = "x")+
  scale_y_continuous(limits = c(-1,3.1),name = "Mean effect estimate")+geom_hline(yintercept=1,linetype = "dashed")+theme_bw()+
  scale_x_continuous(limits = c(1,7),breaks=rev(1:6),labels= rev(round(res$Mean_GxEF)),name = "Mean F-statistic")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_segment(aes(x = index[1], xend = index[1], y = Mean_GxEbeta[1], yend = 3.1), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = index[1], xend = index[1], y = Mean_GxEbeta[1], yend = -1), arrow = arrow(length = unit(0.3, "cm")))

#png("Figure2B2.png",width=3,height=3,units="in",res=300)
#sim2plot1
#dev.off()

mean(GENIUS_beta)

mean(GENIUS_beta-(1.96*GENIUS_se))
mean(GENIUS_beta+(1.96*GENIUS_se))





