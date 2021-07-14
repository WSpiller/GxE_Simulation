#Load necessary libraries
library(genius)
library(AER)
library(ggplot2)
library(ggrepel)

#Set random number gen
set.seed(12345)

N<-100000

Non_zeroGZ<-sample(1:100,10, replace = FALSE)

#Generate single  continuous instrument G
G<-matrix(rep(0, len=N), nrow = N)
for(j in 1:N){
  G[j,]<-rnorm(1,0,1)
}

#Generate single  continuous instrument G
Z<-matrix(rep(0, len=N*100), nrow = N)
for(k in 1:100){
  for(j in 1:N){
    Z[j,k]<-rnorm(1,0,1)
  }
}

#Generate single  continuous instrument G
GZ<-matrix(rep(0, len=N*100), nrow = N)
for(k in 1:100){
  for(j in 1:N){
    GZ[j,k]<-G[j,]*Z[j,k]
  }
}

#Generate single  continuous instrument G
X<-matrix(rep(0, len=N), nrow = N)
for(j in 1:N){
  X[j,]<-1 + G[j,] + rnorm(1,0,10)
  for(k in 1:100){
    if(k %in% Non_zeroGZ){
      gamma3<-rnorm(1,0,10)
      if(gamma3 < 2){
        gamma3<-gamma3 + (5*sign(gamma3))
      }
      X[j,]<-X[j,] + Z[j,k] + gamma3*GZ[j,k]
    }else{
      gamma3<-rnorm(1,0,0.001)
      X[j,]<-X[j,] + Z[j,k] + gamma3*GZ[j,k]
    }
  }
}

#Generate confounder
U<- 1 + rnorm(N,0,1)
X<- X + 7*U

#Generate outcome
Y<- 1 + 0*G + 1*X + 7*U + rnorm(N,0,1)

#OLS estimate using only exposure
Obs_beta<-summary(lm(Y~X))$coef[2,1]

#Specify MR-GENIUS model (interaction agnostic)
GENIUS_1<-genius_addY(c(Y), c(X), c(G), formula = c(X) ~ c(G), alpha = 0.05, lower = -10,
                      upper = 10)

#Het.test for MR-GENIUS
bphet<-bptest(X~G)$statistic

#MR-GENIUS estimate
GENIUS_beta<-GENIUS_1$beta.est

#MR-GENIUS standard error
GENIUS_se<-sqrt(GENIUS_1$beta.var)

#GxE results matrix
GxEB<-matrix(rep(0, len=100*4), nrow = 100)
for(j in 1:100){
  
  t.model<-summary(ivreg(Y~X+G+Z[,j]|G+Z[,j]+GZ[,j]),diagnostics=T)
  
  GxEB[j,1]<-t.model$coef[2,1]
  GxEB[j,2]<-t.model$coef[2,2]
  GxEB[j,3]<-t.model$diagnostics[1,3]
  GxEB[j,4]<-t.model$diagnostics[1,4]
}

#Generate Manhattan-style plot

plot_data<-data.frame(GxEB)
names(plot_data)<-c("beta","se","Fstat","Fpvalue")

sigvec<-rep(0,100)
for(i in 1:100){
  if(plot_data$Fpvalue[i]<0.05/100){
    sigvec[i]<-1
  }
}

plot_data$sig<-as.factor(sigvec)

plot_data$labels<-NULL

for(i in 1:length(plot_data[,1])){
  
  if(plot_data$sig[i]=="1"){
    
    plot_data$labels[i]<-row.names(plot_data)[i]
  }else{
    plot_data$labels[i]<-""
  }
}

B<-ggplot(plot_data,aes(x=c(1:100),y=-log10(Fpvalue)))+
  geom_point(aes(color=sig),size=1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Interaction covariate")+ylab(expression(-log[10](P)))+
  geom_segment(aes(x = 0, xend = 100, y = -log10(0.05/100), yend = -log10(0.05/100)))+
  scale_x_continuous(limits = c(0,101),expand=c(0,0))+
  scale_y_continuous(limits = c(0,max(-log10(plot_data$Fpvalue)+10)),expand=c(0,0))+
  theme(legend.title=element_blank())+theme(legend.position = "none")+
  geom_text_repel(data=plot_data, aes(label=labels),ylim = c(-log10(0.05/100),100))

#non-zero GZ for reference

Non_zeroGZ

png(filename = "F_man.png",
    width = 1600, height = 1200, units = "px", res=200,
    bg = "white")

B

dev.off()


