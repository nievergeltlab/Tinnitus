R
library(data.table)
library(plotrix)
library(scales)


caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)
adam_blue2 <-  rgb(60,160,255,maxColorValue=255)

aqua_ppt <- rgb(0,176,240,max=255)
red_ptt <- red
orange_ppt <- rgb(255,192,0,max=255)


goterm <- fread('lui_pvalues.txt',data.table=F)
names(goterm)[1] <- 'Name'
goterm$color <- adam_blue2


goterm[which(goterm$Name == "Outer"),]$color <-  rgb(174,252,254,maxColorValue=255)
goterm[which(goterm$Name == "Inner"),]$color <-  rgb(174,252,254,maxColorValue=255)
goterm[which(goterm$Name == "Deiter"),]$color <-  rgb(250,252,179,maxColorValue=255)
goterm[which(goterm$Name == "Pillar"),]$color <-  rgb(250,252,179,maxColorValue=255)
goterm[which(goterm$Name == "Melanocytes"),]$color <- "darkgrey" #rgb(10,10,10,maxColorValue=255)


row.names(goterm) <- goterm$Name


goterm <- goterm[order(goterm$p, decreasing=TRUE),]


dm <- as.data.frame(goterm[,c(2)])
row.names(dm) <- goterm$Name

dm <- -log10(dm)


#dm2 <-

#dm <- dm[nrow(dm):1,]

dm[6:34,] <- 1

row.names(dm)[34] <- "Inner.border.phalangeal..Hensen.cells"

#pdf('PTSD_goterm.pdf',7,8)
jpeg('lui.jpeg',2000,3000,quality=100,antialias='none')
par(mar=c(15, 75, 5, 5) + 0.5,lwd=5,cex.axis=5)
myplot <- barplot(height=as.matrix(t(dm)),col=goterm$color,las=1,horiz=TRUE,cex.axis=2,cex.lab=5,xaxt='n',beside=T,xlim=c(0,3),space=c(rep(0.25,6),rep(0,28)))
abline(v=-log10(0.05/5),col='black',lwd=7,lty=2)
abline(v=-log10(0.05),col='black',lwd=7,lty=1)

axis(1,c(0,1,2,3),cex.axis=5,cex.lab=5,lwd=7,lwd.ticks=7,mgp=c(0,5,0))


dev.off()

