R
library(data.table)
library(plotrix)
library(scales)

#Tinnitus explains 
#N=36921 cases and 63507 contorls 90.3676365 prevalence)
res2 <- fread('tinnitus_OR.csv',data.table=F)
names(res2) <- c("Project","Decile","Estimate","Std. Error")

res2[which(res2$Project=="UKBMVPtoR4"),]$Decile <- 1:5 + 0.075
res2[which(res2$Project=="UKBtoR3"),]$Decile <- 1:5
res2[which(res2$Project=="MVPtoUKB"),]$Decile <- 1:5 - 0.075

res2[which(res2$Project=="UKBMVPtoAAM"),]$Decile <- 1:5 + 0.075*2
res2[which(res2$Project=="UKBMVPtoHNA"),]$Decile <- 1:5 + 0.075*3


res2$OR <- exp(res2$Estimate)
res2$LCI <- exp(res2$Estimate - 1.96*res2[,"Std. Error"])
res2$UCI <- exp(res2$Estimate + 1.96*res2[,"Std. Error"])

res2$color <- NA
caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)
 
 green=rgb(147,205, 221,max=255)
 purple=rgb(179,162,199,max=255)


res2[which(res2$Project=="UKBMVPtoR4"),]$color <-adam_blue
res2[which(res2$Project=="UKBtoR3"),]$color <-caro_orange
res2[which(res2$Project=="MVPtoUKB"),]$color <- beauty_red
res2[which(res2$Project=="UKBMVPtoHNA"),]$color <-purple
res2[which(res2$Project=="UKBMVPtoAAM"),]$color <-green



#this will be fore the second data

#Nagelkerke is 0.0096 on the observed scale - 0.0157 assuming the sample prevalence matches the population prevalence.

cis2 <- rbind(res2)


pdf('tinnitus_allgroups_ethnic.pdf',7,8)
par(mar=c(5, 4, 4, 2) + 0.5)
plotCI(x=cis2$Decile,y=cis2$OR,li=cis2$LCI,ui=cis2$UCI,lwd=2,ylim=c(1,1.75),pch=19,cex.axis=1.25,xlab="PRS Quintile",ylab="Quintile Odds Ratio (95% CI)",main="",cex.lab=1.45,col=cis2$color,scol=alpha(cis2$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8)

legend("topleft",col=c("white",beauty_red,caro_orange,adam_blue,green,purple),legend=c("Training -> Target", "MVP -> UKB", "UKB -> MVP", "Meta-analysis -> MVP R4", "Meta-analysis -> MVP AAM", "Meta-analysis -> MVP HNA"),bty="n",pch=19,cex=1.5)

axis(1,at=c(1:5),cex.axis=1.25)
axis(2,at=c(1,1.25,1.5,1.75), labels=c("1","1.25","1.5","1.75"),cex.axis=1.25)


dev.off()

0.3676365


Sample R code:
0.005201	5.89E-16

#To user: Put in these values:
K=0.125 #Population prevalence
P=0.36  #Sample prevalence
h2=0.005201 #h2 on the observed scale
seh2=5.89E-16 #standard error of h2 on the observed scale

#Do calculations
zv <- dnorm(qnorm(K))
h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
var_h2_liab <- ( seh2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2) ^2

#Report liability scale h2snp and se 
h2_liab #0.0797
sqrt(var_h2_liab) #0.0159
