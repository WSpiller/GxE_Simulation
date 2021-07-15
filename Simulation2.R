#Load necessary libraries
library(genius)
library(AER)
library(ggplot2)
library(ggrepel)

#Set random number gen
set.seed(12345)

N<-100000

iterations<-100

#Generate empty vectors for estimates
Obs_beta<-rep(0,iterations)
GENIUS_beta<-rep(0,iterations)
GENIUS_se<-rep(0,iterations)
bphet<-rep(0,iterations)
GxE_beta1<-rep(0,iterations)
GxE_beta2<-rep(0,iterations)
GxE_beta3<-rep(0,iterations)
GxE_beta4<-rep(0,iterations)
GxE_beta5<-rep(0,iterations)
GxE_beta6<-rep(0,iterations)

GxE_se1<-rep(0,iterations)
GxE_se2<-rep(0,iterations)
GxE_se3<-rep(0,iterations)
GxE_se4<-rep(0,iterations)
GxE_se5<-rep(0,iterations)
GxE_se6<-rep(0,iterations)

F_stats1<-rep(0,iterations)
F_stats2<-rep(0,iterations)
F_stats3<-rep(0,iterations)
F_stats4<-rep(0,iterations)
F_stats5<-rep(0,iterations)
F_stats6<-rep(0,iterations)

pF_stats1<-rep(0,iterations)
pF_stats2<-rep(0,iterations)
pF_stats3<-rep(0,iterations)
pF_stats4<-rep(0,iterations)
pF_stats5<-rep(0,iterations)
pF_stats6<-rep(0,iterations)

for(i in 1:iterations){

#Generate single  continuous instrument G
G<-matrix(rep(0, len=N), nrow = N)
for(j in 1:N){
  G[j,]<-rnorm(1,0,1)
}

#Generate single  continuous instrument G
Z<-matrix(rep(0, len=N*6), nrow = N)
for(k in 1:6){
  for(j in 1:N){
    Z[j,k]<-rnorm(1,0,1)
  }
}

#Generate single  continuous instrument G
GZ<-matrix(rep(0, len=N*6), nrow = N)
for(k in 1:6){
  for(j in 1:N){
    GZ[j,k]<-G[j,]*Z[j,k]
  }
}

gamma3<-c(0.01,0.642,1.0365,1.54,2.237,3.15)

#Generate single  continuous instrument G
X<-matrix(rep(0, len=N), nrow = N)
for(j in 1:N){
  X[j,]<-1 + G[j,] + rnorm(1,0,100)
  for(k in 1:6){

      X[j,]<-X[j,] + Z[j,k] + gamma3[k]*GZ[j,k]
    }
}

#Generate confounder
U<- 1 + rnorm(N,0,1)
X<- X + 7*U

#Generate outcome
Y<- 1 + 0*G + 1*X + 7*U + rnorm(N,0,1)

#OLS estimate using only exposure
tObs_beta<-summary(lm(Y~X))$coef[2,1]

#Specify MR-GENIUS model (interaction agnostic)
GENIUS_1<-genius_addY(c(Y), c(X), c(G), formula = c(X) ~ c(G), alpha = 0.05, lower = -10,
                      upper = 10)

#Het.test for MR-GENIUS
tbphet<-bptest(X~G)$statistic

#MR-GENIUS estimate
GENIUS_b<-GENIUS_1$beta.est

#MR-GENIUS standard error
GENIUS_tse<-sqrt(GENIUS_1$beta.var)

#GxE results matrix
GxEB<-matrix(rep(0, len=6*4), nrow = 6)
for(j in 1:6){
  
  t.model<-summary(ivreg(Y~X+G+Z[,j]|G+Z[,j]+GZ[,j]),diagnostics=T)
  
  GxEB[j,1]<-t.model$coef[2,1]
  GxEB[j,2]<-t.model$coef[2,2]
  GxEB[j,3]<-t.model$diagnostics[1,3]
  GxEB[j,4]<-t.model$diagnostics[1,4]
  
  F_stats1[i]<-GxEB[1,3]
  F_stats2[i]<-GxEB[2,3]
  F_stats3[i]<-GxEB[3,3]
  F_stats4[i]<-GxEB[4,3]
  F_stats5[i]<-GxEB[5,3]
  F_stats6[i]<-GxEB[6,3]
  
  GxE_beta1[i]<-GxEB[1,1]
  GxE_beta2[i]<-GxEB[2,1]
  GxE_beta3[i]<-GxEB[3,1]
  GxE_beta4[i]<-GxEB[4,1]
  GxE_beta5[i]<-GxEB[5,1]
  GxE_beta6[i]<-GxEB[6,1]
  
  GxE_se1[i]<-GxEB[1,2]
  GxE_se2[i]<-GxEB[2,2]
  GxE_se3[i]<-GxEB[3,2]
  GxE_se4[i]<-GxEB[4,2]
  GxE_se5[i]<-GxEB[5,2]
  GxE_se6[i]<-GxEB[6,2]
  
  Obs_beta[i]<-tObs_beta
  GENIUS_beta[i]<-GENIUS_b
  GENIUS_se[i]<-GENIUS_tse
  bphet[i]<-tbphet
}

}

meanFs<-c(mean(F_stats1),mean(F_stats2),mean(F_stats3),
          mean(F_stats4),mean(F_stats5),mean(F_stats6))

mean_ge_beta<-c(mean(GxE_beta1),mean(GxE_beta2),mean(GxE_beta3),
          mean(GxE_beta4),mean(GxE_beta5),mean(GxE_beta6))

mean_ge_se<-c(mean(GxE_se1),mean(GxE_se2),mean(GxE_se3),
                mean(GxE_se4),mean(GxE_se5),mean(GxE_se6))

ge_lci<-mean_ge_beta-(1.96*mean_ge_se)
  
ge_uci<-mean_ge_beta+(1.96*mean_ge_se)

res_data<-data.frame(mean_ge_beta,mean_ge_se,ge_lci,ge_uci,meanFs)

write.csv(res_data,"results.csv",row.names=F)

res_data<-read.csv("results.csv",header=T)

res_data$index<-1:6

#Plot results

sim2plot1<-ggplot(res_data)+geom_point(aes(x=index, y=mean_ge_beta),size=0.8)+
  coord_flip()+geom_errorbar(aes(ymin=ge_lci,ymax=ge_uci,x=index),width=0,orientation = "x")+
  scale_y_continuous(limits = c(0.25,1.75),name = "Mean effect estimate")+geom_hline(yintercept=1,linetype = "dashed")+theme_bw()+
  scale_x_continuous(limits = c(1,7),breaks=rev(1:6),labels= rev(round(res_data$meanFs)),name = "Mean F-statistic")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_segment(aes(x = index[1], xend = index[1], y = 0.25, yend = 1.75), arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = index[1], xend = index[1], y = 1.75, yend = 0.25), arrow = arrow(length = unit(0.3, "cm")))

png("Figure2B.png",width=3,height=3,units="in",res=200)
sim2plot1
dev.off()

