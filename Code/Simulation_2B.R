res1<-read.csv("F1_res.csv",header=T)
res2<-read.csv("F5_res.csv",header=T)
res3<-read.csv("F10_res.csv",header=T)
res4<-read.csv("F25_res.csv",header=T)
res5<-read.csv("F50_res.csv",header=T)
res6<-read.csv("F100_res.csv",header=T)

res<-rbind(res1,res2,res3,res4,res5,res6)
