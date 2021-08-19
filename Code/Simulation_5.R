#Load necessary libraries
library(genius)
library(AER)
library(ggplot2)
library(ggrepel)

G_beta<-NULL
G_se<-NULL
medFvec<-NULL
GxE_beta_all<-NULL
GxE_se_all<-NULL
GxE_beta_valid<-NULL
GxE_se_valid<-NULL
bp_out<-NULL
bpb_out<-NULL

for(m in c(0,1,5,10,50,99)){
  
  print(m)
  
  #Set random number gen
  set.seed(12345)
  
  N<-10000
  
  No.G2invalid<-m
  
  Z_number<-100
  
  GZ_number<-100
  
  sampleGZ<-sample(1:Z_number,GZ_number, replace = FALSE)
  
  sample_invalid<-sample(1:Z_number,No.G2invalid, replace = FALSE)
  
  Zvec<-1:Z_number
  
  validZ<-Zvec[!Zvec %in% sample_invalid]
  
  validZ<-validZ[1]
  
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
    
    gamma3<-rep(0,100)
    
    for(k in sampleGZ){
      while(abs(gamma3[k]) < 1){
        gamma3[k]<-abs(rnorm(1,2,2))
      }
    }
    
    gamma3[validZ]<-max(gamma3)*2
    
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
    Y<- 1 + G[,1] + 1*X + 1*U + rnorm(N,0,1)
    for(k in 1:Z_number){
      Y[,1]<-Y[,1] + 0*Z[,k] + 0*GZ[,k]
    }
    
    for(k in sample_invalid){
      Y<- Y + 1*GZ[,k]
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
  
  
  mean(GENIUS_se)
  
  #Het.test for MR-GENIUS
  mean(bpstat)
  mean(bp_pval)
  
  #Avg.GxE values
  
  res<-matrix(rep(0, len=Z_number*4), nrow = Z_number)
  
  for(k in 1:Z_number){
    
    res[k,1]<-median(it.GxEB[k,])
    res[k,2]<-median(it.GxEse[k,])
    res[k,3]<-median(it.GxEF[k,])
    res[k,4]<-median(it.GxEFp[k,])
    
  }
  
  res<-data.frame(res)
  names(res)<-c("Mean_GxEbeta","Mean_GxEse","Mean_GxEF","Mean_GxEFpvalue")
  
  G_beta<-c(G_beta,mean(GENIUS_beta))
  G_se<-c(G_se,mean(GENIUS_se))
  medFvec<-c(medFvec,mean(res[,3]))
  GxE_beta_all<-c(GxE_beta_all,mean(res[,1]))
  GxE_se_all<-c(GxE_se_all,mean(res[,2]))
  bpb_out<-c(bpb_out,mean(bp_pval))
  bp_out<-c(bp_out,mean(bpstat))
  GxE_beta_valid<-c(GxE_beta_valid,res[validZ,1])
  GxE_se_valid<-c(GxE_se_valid,res[validZ,2])
  
}

all_out<-data.frame(G_beta,G_se,GxE_beta_all,GxE_se_all,GxE_beta_valid,GxE_se_valid,medFvec,bpb_out,bp_out)

all_out$G_beta - 1.96 * all_out$G_se
all_out$G_beta + 1.96 * all_out$G_se

all_out$GxE_beta_all - 1.96 * all_out$GxE_se_all
all_out$GxE_beta_all + 1.96 * all_out$GxE_se_all

all_out$GxE_beta_valid - 1.96 * all_out$GxE_se_valid
all_out$GxE_beta_valid + 1.96 * all_out$GxE_se_valid
