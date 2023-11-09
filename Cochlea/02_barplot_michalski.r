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


goterm <- fread('michalski_pvalues.txt',data.table=F)
goterm$color <- adam_blue2

goterm[which(goterm$Type == "Circulating"),]$color <- rgb(223,90,81,maxColorValue=255)
goterm[which(goterm$Type == "GliaCells"),]$color <- rgb(107,150,255,maxColorValue=255)
goterm[which(goterm$Type == "HairCells"),]$color <-  rgb(174,252,254,maxColorValue=255)
goterm[which(goterm$Type == "LateralWall"),]$color <- "grey"
goterm[which(goterm$Type == "Neurons"),]$color <- rgb(88,103,248,maxColorValue=255)
goterm[which(goterm$Type == "SupportingCells"),]$color <-  rgb(250,252,179,maxColorValue=255)
goterm[which(goterm$Type == "SurroundingStructures"),]$color <- rgb(202,182,219,maxColorValue=255)

"grey"

row.names(goterm) <- goterm$Name


goterm <- goterm[order(goterm$Type, goterm$MAGMA, decreasing=TRUE),]


dm <- as.data.frame(goterm[,c(2)])
row.names(dm) <- goterm$Name

dm <- -log10(dm)


#dm2 <-

#dm <- dm[nrow(dm):1,]

spacings <- c(rep(0,7),0.5,rep(0,2),0.5,0.5,rep(0,6),0.5,0,0.5,rep(0,3),0.5,rep(0,8),0.5)
colorset <-


#pdf('PTSD_goterm.pdf',7,8)
jpeg('michalski.jpeg',2000,3000,quality=100,antialias='none')
par(mar=c(15, 75, 5, 5) + 0.5,lwd=5,cex.axis=5)
myplot <- barplot(height=as.matrix(t(dm)),col=goterm$color,las=1,horiz=TRUE,cex.axis=2,cex.lab=5,xaxt='n',beside=T,xlim=c(0,3),space=spacings[length(spacings):1])
abline(v=-log10(0.05/34),col='black',lwd=7,lty=2)
abline(v=-log10(0.05),col='black',lwd=7,lty=1)

axis(1,c(0,1,2,3),cex.axis=5,cex.lab=5,lwd=7,lwd.ticks=7,mgp=c(0,5,0))


dev.off()


     C1 C2 C3 C4 C5
Supporting
LDSC


##Lui 
goterm <- fread('lui_pvalues.txt',data.table=F)
goterm$color <- adam_blue2

#goterm[which(goterm$Term == "Cellular Component"),]$color <- beauty_red
#goterm[which(goterm$Term == "Molecular Function"),]$color <- aqua_ppt
#goterm[which(goterm$Term == "Biological Process"),]$color <-  orange_ppt # rgb(127,253,166,maxColorValue=255) 

row.names(goterm) <- goterm$Name


goterm <- goterm[order(goterm$Name,decreasing=TRUE),]

dm <- goterm[,c(2,3)]


dm <- -log10(dm)




#dm2 <-

#dm <- dm[nrow(dm):1,]

#pdf('PTSD_goterm.pdf',7,8)
jpeg('lui.jpeg',2000,3000,quality=100,antialias='none')
par(mar=c(5, 45, 5, 5) + 0.5,lwd=5,cex.axis=5)
myplot <- barplot(height=as.matrix(t(dm)),col=c(adam_blue,orange_ppt),las=1,horiz=TRUE,cex.axis=2,cex.lab=5,xaxt='n',beside=T,xlim=c(0,7))
abline(v=-log10(0.05/5),col='black',lwd=7,lty=2)

axis(1,c(0,1,2,3,4,5,6,7),cex.axis=5,cex.lab=5,lwd=7,lwd.ticks=7,mgp=c(0,5,0))


dev.off()



goterm <- fread('michalski_magma_transpose.txt',data.table=F)
goterm$color <- adam_blue2
row.names(goterm) <- goterm$Name

dm <- as.data.frame(goterm[,-c(1,dim(goterm)[2])])
row.names(dm) <- goterm$Name

dm <- -log10(dm)


#dm2 <-

#dm <- dm[nrow(dm):1,]

#pdf('PTSD_goterm.pdf',7,8)
jpeg('michalski2.jpeg',2000,3000,quality=100,antialias='none')
par(mar=c(5, 75, 5, 5) + 0.5,lwd=5,cex.axis=5)
myplot <- barplot(height=as.matrix(t(dm)),col=c(adam_blue),las=1,horiz=TRUE,cex.axis=2,cex.lab=5,xaxt='n',beside=T,xlim=c(0,3))
abline(v=-log10(0.05/35),col='black',lwd=7,lty=2)
axis(1,c(0,1,2,3),cex.axis=5,cex.lab=5,lwd=7,lwd.ticks=7,mgp=c(0,5,0))

dev.off()


     C1 C2 C3 C4 C5
Supporting
LDSC


##Lui 
goterm <- fread('lui_pvalues.txt',data.table=F)
goterm$color <- adam_blue2

#goterm[which(goterm$Term == "Cellular Component"),]$color <- beauty_red
#goterm[which(goterm$Term == "Molecular Function"),]$color <- aqua_ppt
#goterm[which(goterm$Term == "Biological Process"),]$color <-  orange_ppt # rgb(127,253,166,maxColorValue=255) 

row.names(goterm) <- goterm$Name


goterm <- goterm[order(goterm$Name,decreasing=TRUE),]

dm <- goterm[,c(2,3)]


dm <- -log10(dm)




#dm2 <-

#dm <- dm[nrow(dm):1,]

#pdf('PTSD_goterm.pdf',7,8)
jpeg('lui.jpeg',2000,3000,quality=100,antialias='none')
par(mar=c(5, 45, 5, 5) + 0.5,lwd=5,cex.axis=5)
myplot <- barplot(height=as.matrix(t(dm)),col=c(adam_blue,orange_ppt),las=1,horiz=TRUE,cex.axis=2,cex.lab=5,xaxt='n',beside=T,xlim=c(0,7))
abline(v=-log10(0.05/5),col='black',lwd=7,lty=2)
axis(1,c(0,1,2,3,4,5,6,7),cex.axis=5,cex.lab=5,lwd=7,lwd.ticks=7,mgp=c(0,5,0))


dev.off()






,xlab="-log10 p-value",

names.arg=goterm$Name