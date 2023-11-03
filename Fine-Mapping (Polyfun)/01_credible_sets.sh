#Read polyfun output files (concatenated files)
d1 <- read.table('finemapLD.inform.TRANS.complete.txt',stringsAsFactors=F,header=T)

library(plyr)
library(data.table)

range2 <- function(x)
{
r1 <- range(x)
return(c(r1[1], r1[2]))
}

#Count how many SNPs were tested (i.e. this will include the SNPs in the credible set and those that do not make it to the credible set)
d2 <- read.table('finemapLD.inform.EUR.complete.allsnps.txt',stringsAsFactors=F,header=T)
credsets3 <- ddply(d2, ~ Locus, colwise(length,"BP"))

write.table(credsets3,file="all_set_nsnps.txt",row.names=F)

#Calculate width of credible set
credsets1 <- ddply(d1, ~ Locus, colwise(min,"BP"))
credsets2 <- ddply(d1, ~ Locus, colwise(max,"BP"))

write.table(cbind(credsets1,credsets2),file="credible_set_windows.txt",row.names=F)

#Count how many SNPs are in credible set
credsets2 <- ddply(d1, ~ Locus, colwise(length,"BP"))

write.table(credsets2,file="credible_set_nsnps.txt",row.names=F)



#Merge FUMA snps.txt with the SNP results from Polyfun (i.e. to annotate SNPs in credible sets)
d3 <- fread('eur_snps.txt',data.table=F)
d3$SNP <- d3$rsID

dm <- merge(d3,d2,by="SNP",suffixes=c("_fuma","_finemap"))

credsets4 <- ddply(dm, ~ Locus, colwise(length,"BP"))

write.table(credsets4,file="eur_set_nsnps_fuma.txt")
