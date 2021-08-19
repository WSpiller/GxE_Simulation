#Load necessary libraries
library(genius)
library(AER)
library(ggplot2)
library(ggrepel)

#Set random number gen
set.seed(12345)

N<-10000

G_number<-100

proportions<-c(0,1,5,10,50,100)

iterations<-1000

resmat<-matrix(rep(0, len=11*(length(proportions))), nrow = length(proportions))

for(m in 1:length(proportions)){
  
  No.invalid<-proportions[m]
  
  print(No.invalid)
  
  invalid_GZs<-sample(1:G_number,No.invalid, replace = FALSE)
  
  #Generate empty vectors for estimates
  Obs_beta<-rep(0,iterations)
  GENIUS_beta<-rep(0,iterations)
  GENIUS_se<-rep(0,iterations)
  GxE_beta<-rep(0,iterations)
  GxE_se<-rep(0,iterations)
  bpstat<-rep(0,iterations)
  bp_pval<-rep(0,iterations)
  F_stat<-rep(0,iterations)
  F_p<-rep(0,iterations)
  Sargan_stat<-rep(0,iterations)
  Sargan_p<-rep(0,iterations)
  
  for(i in 1:iterations){
    
    print(i)
    
    #Generate interaction covariates Z
    G<-matrix(rep(0, len=N*G_number), nrow = N)
    
    for(k in 1:G_number){
      G[,k]<-rnorm(N,0,1)
    }
    
    #Single interaction Z
    
    #Generate single continuous instrument G
    Z<-matrix(rep(0, len=N), nrow = N)
    Z[,1]<-rnorm(N,0,1)
    
    #Generate GZ interaction variable
    GZ<-matrix(rep(0, len=N*G_number), nrow = N)
    
    for(k in 1:G_number){
      GZ[,k]<-G[,k] * Z[,1]
    }
    
    #Generate single continuous confounder U
    U<-rnorm(N,0,1)
    
    #Generate single continuous exposure X
    X<-matrix(rep(0, len=N), nrow = N)
    X[,1]<- 1 + Z[,1] + 1*U + rnorm(N,0,10)
    for(k in 1:G_number){
      X[,1]<- X[,1] + G[,k] + 1*GZ[,k]
    }
    
    beta4<-rep(0,G_number)
    
    for(j in invalid_GZs){
      beta4[j]<-1
    }
    
    #Generate single continuous outcome Y
    Y<- 1 + G[,1] + 1*X + 1*U + rnorm(N,0,1)
    for(k in 1:G_number){
      Y[,1]<-Y[,1] + 0*G[,k] + beta4[k]*GZ[,k]
    }
    
    #Create PRS of G
    PRS<-rep(0,N)
    
    for(k in 1:G_number){
      PRS<-PRS + G[,k]
    }
    
    
    #Create PRS of G
    PRSGZ<-rep(0,N)
    
    for(k in 1:G_number){
      PRSGZ<-PRS * Z
    }
    
    #OLS estimate using only exposure
    Obs_beta[i]<-summary(lm(Y~X))$coef[2,1]
    
    #Specify MR-GENIUS model (interaction agnostic)
    GENIUS_1<-genius_addY(c(Y), c(X), c(PRS), formula = c(X) ~ c(PRS), alpha = 0.05, lower = -10,
                          upper = 10)
    
    #MR-GENIUS estimate
    GENIUS_beta[i]<-GENIUS_1$beta.est
    
    #MR-GENIUS standard error
    GENIUS_se[i]<-sqrt(GENIUS_1$beta.var)
    
    #Het.test for MR-GENIUS
    bpstat[i]<-bptest(X~PRS)$statistic
    bp_pval[i]<-bptest(X~PRS)$p.value
    
    t.model<-summary(ivreg(Y~X+Z+PRS|PRS+Z+PRSGZ),diagnostics=T)
    
    t.modelS<-summary(ivreg(Y~X+Z+G|G+Z+GZ),diagnostics=T)
    
    GxE_beta[i]<-t.model$coef[2,1]
    GxE_se[i]<-t.model$coef[2,2]
    
    F_stat[i]<-t.modelS$diagnostics[1,3]
    F_p[i]<-t.modelS$diagnostics[1,4]
    
    Sargan_stat[i]<-t.modelS$diagnostics[3,3]
    Sargan_p[i]<-t.modelS$diagnostics[3,4]
    
  }
  
  resmat[m,1]<-mean(Obs_beta)
  resmat[m,2]<-median(GENIUS_beta)
  resmat[m,3]<-median(GENIUS_se)
  resmat[m,4]<-mean(GxE_beta)
  resmat[m,5]<-mean(GxE_se)
  resmat[m,6]<-mean(F_stat)
  resmat[m,7]<-mean(F_p)
  resmat[m,8]<-mean(bpstat)
  resmat[m,9]<-mean(bp_pval)
  resmat[m,10]<-mean(Sargan_stat)
  resmat[m,11]<-mean(Sargan_p)
  
}

results<-data.frame(resmat)
names(results)<-c("Obs_beta","GENIUS_beta","GENIUS_se","GxE_beta","GxE_se",
                  "F_stat","F_p","bpstat","bp_pval","Sargan_stat","Sargan_p")


results$GENIUS_beta - 1.96 * results$GENIUS_se
results$GENIUS_beta + 1.96 * results$GENIUS_se

results$GxE_beta - 1.96 * results$GxE_se
results$GxE_beta + 1.96 * results$GxE_se

