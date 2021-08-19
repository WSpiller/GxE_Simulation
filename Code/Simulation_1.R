#Load necessary libraries
library(genius)
library(AER)
library(ggplot2)
library(ggrepel)

#Set random number gen
set.seed(12345)

N<-100000

Z_number<-100

GZ_number<-10

sampleGZ<-sample(1:Z_number,GZ_number, replace = FALSE)

iterations<-1000

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
  
  gamma3<-rnorm(Z_number,0,0.01)
  
  for(k in sampleGZ){
    gamma3[k]<-rnorm(1,2,2)
    while(abs(gamma3[k]) < 1){
      gamma3[k]<-rnorm(1,2,2)
    }
  }
  
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
  X[,1]<- 1 + G[,1] + 1*U + rnorm(N,0,10)
  for(k in 1:Z_number){
    X[,1]<- X[,1] + Z[,k] + gamma3[k]*GZ[,k]
  }
  
  #Generate single continuous outcome Y
  Y<- 1 + 1*G[,1] + 1*X + 1*U + rnorm(N,0,1)
  for(k in 1:Z_number){
    Y[,1]<-Y[,1] + 0*Z[,k] + 0*GZ[,k]
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

#MR-GENIUS CI

mean(GENIUS_beta) - (1.96*mean(GENIUS_se))
mean(GENIUS_beta) + (1.96*mean(GENIUS_se))

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

sigvec<-rep(0,Z_number)
for(i in 1:Z_number){
  if(res$Mean_GxEFpvalue[i]<0.05/Z_number){
    sigvec[i]<-1
  }
}

res$sig<-as.factor(sigvec)

res$labels<-NULL

for(i in 1:length(res[,1])){
  
  if(res$sig[i]=="1"){
    
    res$labels[i]<-row.names(res)[i]
  }else{
    res$labels[i]<-""
  }
}

B<-ggplot(res,aes(x=c(1:Z_number),y=-log10(Mean_GxEFpvalue)))+
  geom_point(aes(color=sig),size=1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Interaction covariate")+ylab(expression(-log[10](P)))+
  geom_segment(aes(x = 0, xend = Z_number, y = -log10(0.05/Z_number), yend = -log10(0.05/Z_number)))+
  scale_x_continuous(limits = c(0,(Z_number+1)),expand=c(0,0))+
  scale_y_continuous(limits = c(0,max(-log10(res$Mean_GxEFpvalue)+10)),expand=c(0,0))+
  theme(legend.title=element_blank())+theme(legend.position = "none")+
  geom_text_repel(data=res, aes(label=labels))

B

png("Figure2A.png",width=4,height=3,units="in",res=200)
B
dev.off()

#mean values for MR-GENIUS and across valid instruments

res[res$sig ==1,]
mean(res[res$sig ==1,1])
mean(res[res$sig ==1,1]) - 1.96*mean(res[res$sig ==1,2])
mean(res[res$sig ==1,1]) + 1.96*mean(res[res$sig ==1,2])
mean(res[res$sig ==1,3])

mean(GENIUS_beta)
mean(GENIUS_beta) - 1.96 * mean(GENIUS_se)
mean(GENIUS_beta) + 1.96 * mean(GENIUS_se)

mean(bpstat)
mean(bp_pval)
