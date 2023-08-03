 
#PRS of MVP + EHR -> UKBB

#Format sumstats
 zcat /mnt/ukbb/adam/tinnitus_gwas/MVP_tinnitus/eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz | awk '{ if (NR==1) {$1="SNP";$4="A1";$5="A2";BETA="BETA";$12="P"}; if(NR>1) BETA=log($9); print $1,toupper($4),toupper($5),BETA,$12}' > eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.sst

#Get UKBB SNP list
 cat /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_*.bim > ukbb_snps.bim

#Run PRS-CS

#UKB

 python /mnt/ukbb/adam/ptsd/cnv_aug_2021/PRScs-master/PRScs.py --ref_dir=/mnt/ukbb/adam/ptsd/prs_csx/ldblk_1kg_eur --bim_prefix=ukbb_snps \
 --sst_file=eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.sst --n_gwas=308879   \
 --out_dir=ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25 --phi=1e-2  




#Combine all outputs
cat ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr*.txt | awk '{print $2,$4,$6}' > ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt
awk '{print $1}' ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt  >  ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist

for chr in {1..22}
do
cat ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt | awk '{print $2,$4,$6}' > ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt.plink
done


#extract the list of SNPs relevant to the PRS
 for chr in {1..22}
 do
  ../plink2 --bed /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_"$chr".bed \
 --bim /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_"$chr".bim \
 --fam /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_"$chr".fam \
 --extract   ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist  \
 --make-bed --out ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr"  
 done
 
#Get PRS and un average it.


for chr in {1..22}
do
mv ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim.bk
done

for chr in {1..22}
do
Rscript rename_bim.r ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim.bk ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim
done


for chr in {1..22}
do
 #will this break if the chromosome is not shown?
# ../plink2 --bfile ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr"   --score ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt.plink --out ukbb_mvpehr_pred_25/ukbb_$chr
 #sed 's/#//g' ukbb_mvpehr_pred_25/ukbb_"$chr".sscore | awk -v CHR=$chr '{if(NR==1) TSCORE="TSCORE"CHR; if(NR>1) TSCORE=$3*$5; print $1,$2, TSCORE}' >  ukbb_mvpehr_pred_25/ukbb_"$chr".sscore.fixed
 
 ../plink --bfile ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr"   --score ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt.plink sum --out ukbb_mvpehr_pred_25/ukbb_$chr
done

for chr in {1..22}
do
awk -v chr=$chr '{if(NR==1) $6="SCORESUM"chr; print}' ukbb_mvpehr_pred_25/ukbb_"$chr".profile >  ukbb_mvpehr_pred_25/ukbb_"$chr".profile.fixed
done

R
library(plyr)
library(data.table)

#read each file into a data frame with the same name
for (i in 1:22)
{
	assign(
		paste('ukbb_',i,sep=''), fread(paste('ukbb_mvpehr_pred_25/ukbb_',i,'.profile.fixed',sep=''), header=T, data.table=F)[,c(1,2,6)]
		) 
}

#parse the text list of data frame names as a list of data frames
data_list <- eval( 
			parse( 
				text=paste(
					"list(", paste(grep("ukbb",ls(),value=TRUE), collapse=','), ")" 
					)
				)
			)


#combine all data frames by id_visit (won't work for subjects missing this variable!!!!!)
datA <- join_all(data_list,by=c("FID","IID"), type="left", match="first")
#sum over all score columns

#datA$SCORE <- apply(datA[,-c(1,2)],1,sum) #plink2
datA$SCORE <- apply(datA[,-c(1:2)],1,sum) #plink1
write.table(datA,'ukbb_mvpehr_pred_25/ukbb_mvpehr_pred.profile',quote=F,row.names=F)

#take the association between this and PTSD..

R
nagelkerke2 <- function(lintercept,lfull,n)
{
 #specify the log likelihood of the intercept model, the full model, and the n obs
 numer <- 1-(exp(2*(lintercept - lfull)/n))
 denom <- 1-exp(2*lintercept/n)
 nagelkerke=numer/denom
 return(as.numeric(nagelkerke))
}




library(data.table)
library(ordinal)
prs <- fread('ukbb_mvpehr_pred_25/ukbb_mvpehr_pred.profile',data.table=F)
pcs <- fread('/mnt/ukbb/adam/tinnitus_gwas/phenotype/UKB_tinnitus_eur_related_april2_2019.pheno',data.table=F)


d2 <- merge(prs,pcs,by=c("FID","IID"),suffixes=c("","_na"))

d2$SCALE <- scale(d2$SCORE)

#ologit on scores
m1 <- clm(as.factor(f.4803.max_coding3) ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) +SCALE,data=subset(d2))
m1a <- clm(as.factor(f.4803.max_coding3) ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) ,data=subset(d2))
m1n <- clm(as.factor(f.4803.max_coding3) ~ 1,data=subset(d2))

nagelkerke2(logLik(m1n),logLik(m1),172995) - nagelkerke2(logLik(m1n),logLik(m1a),172995)
# 0.007173332

#binary on dx (dx is > 0)
library(fmsb)
lm1 <- glm(f.4803.max_coding2 ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) +SCALE,data=subset(d2),family='binomial')
lm1a <- glm(f.4803.max_coding2 ~f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) ,data=subset(d2),family='binomial')




NagelkerkeR2(lm1)$R2 - NagelkerkeR2(lm1a)$R2

#0.007357539

 #coding 1 (middle people removed)
 Xlm1 <- glm(f.4803.max_coding1 ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) +SCALE,data=subset(d2),family='binomial')


#Quantile odds ratios?

d2$prsquantil <- cut(d2$SCORE,breaks=quantile(d2$SCORE,seq(0,1,0.2)))
Jm1 <- as.data.frame(summary(glm(f.4803.max_coding2  ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+prsquantil,data=subset(d2),family='binomial'))$coefficients[8:11,1:2])

#Export regression coefficients to file for meta analysis
 m1xbetas <- Jm1
 names(m1xbetas) <- c("beta","se")

 study1='ukbb'
 study2='ukbb'
 m1xbetas2 <- t(m1xbetas[,1])
 row.names(m1xbetas2) <- c(paste(study1,study2,sep="_"))

 m1xbetas3 <- t(m1xbetas[,2])
 row.names(m1xbetas2) <- c(paste(study1,study2,sep="_"))
 row.names(m1xbetas3) <- c(paste(study1,study2,sep="_"))

 write.table(m1xbetas2, file=paste("ukbb_ukbb_mvpehrpred.results.beta",sep=''),quote=F,row.names=T,col.names=F)
 write.table(m1xbetas3, file=paste("ukbb_ukbb_mvpehrpred.results.se",sep=''),quote=F,row.names=T,col.names=F)


m2 <- as.data.frame(t(matrix(c(0,0.00000000001))))

names(m2) <- names(Jm1)

ms <- rbind(m2,Jm1)
rownames(ms) <- paste("Q",c(1:5),sep="")

names(ms)[2] <- "se"

ms$OR <- exp(ms$Estimate)
ms$LCI = exp(ms$Estimate -1.96*ms$se)
ms$UCI = exp(ms$Estimate +1.96*ms$se)
ms$Decile <-c(1:5)

petrolblue <- rgb(62,81,192,maxColorValue=255 )

ms$color <- petrolblue

res2 <- ms

library(Hmisc)
library(fmsb)
library(data.table)
library(plotrix)
library(lmtest)

res2$color="blue"
res2 <- subset(res2,select=c(Decile,OR,LCI,UCI,color))


pdf('prs_decile_ukb_pred_mvp.pdf',7,7)

plotCI(x=res2$Decile,y=res2$OR,li=res2$LCI,ui=res2$UCI,lwd=2,ylim=c(1,1.6),pch=19,cex.axis=1.25,xlab="PRS Decile",ylab="Quintile Odds Ratio (95% CI)",main="",cex.lab=1.45,col=res2$color,scol=alpha(res2$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8)
legend("topleft",col=c("white","blue","red"),legend=c("MVP -> UKB"),bty="n",pch=19,cex=1.5)
axis(1,at=c(1:10),cex.axis=1.25)
axis(2,at=c(1,1.25,1.5,1.75), labels=c("1","1.25","1.5","1.75"),cex.axis=1.25)

dev.off()


#Compare PRSes


R
library(data.table)
d1 <- fread('/mnt/ukbb/adam/tinnitus_gwas/output_prsice/ukbb_tinnitus_oct5_2021.best',data.table=F)
c1 <- fread('/mnt/ukbb/adam/tinnitus_gwas/phenotype/UKB_tinnitus_eur_related_may2_2019.pheno',data.table=F)

dm <- merge(d1,c1,by=c("FID","IID"))

dm$prs_center <- scale(dm$PRS )

dm2 <- merge(dm,d2,by="FID")

png('prs_comparison.png')
 plot(dm2$SCALE,dm2$prs_center)
dev.off()
