#We have an average score for all chunks. We need a score summed across all chunks to get a full PRS.

#Motivation: PLINK gives average scores, but not total scores.  Because different N SNPs included in each chunk, summing averages would be wrong.We need to get an actual allelic score and not an average.

#Make a directory for these adjusted outputs 
mkdir  /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/adjusted_prs


#Call into raw PRS directory and calculate PRS
cd /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/

#do not include the X chromosome

for file in $(ls | grep sscore | grep -v chrX)
do
 cat $file | sed 's/#//g' | awk  '{if (NR==1) TSCORE="TSCORE"; if(NR>1) TSCORE=$3*$5; print $1,$2,TSCORE}' > adjusted_prs/$file
done


#To get a total score, stack all N files then just use dplyr to sum within subjects

cat adjusted_prs/* | awk '{if(NR==1 || $1!="FID") print $1,$3}' | sort -k1b,1 > prs_f3_to_mvp_allsubjects.presum.prs

grep -v FID prs_f3_to_mvp_allsubjects.presum.prs > prs_f3_to_mvp_allsubjects.presum.prs2
#Or split this file by subject in awk then calclate the sum in awk, then add all sum files together (if R doesnt work)

#Compute total PRS using chunk PRSes 

R 
library(data.table)
library(plyr)
d1 <- fread('prs_f3_to_mvp_allsubjects.presum.prs2',data.table=F)
names(d1) <-c ("FID","TSCORE")

#Sum within all subjects
d1_scored <- ddply(d1, ~FID, colwise(sum,"TSCORE"))

write.table(d1_scored, file='/group/research/mvp039/AM/tinnitus_paper2_prs/prs_f3_to_mvp_allsubjects.prs',row.names=F)


#Now actually do analysis


R 
library(data.table)
library(plyr)

prs <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/prs_f3_to_mvp_allsubjects.prs',stringsAsFactors=F,header=T)

peur <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/tinnitus_holdout_eur_nov01_2022.pheno',stringsAsFactors=F,header=T)
#paam <- fread('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/PCLsum_aam_oct24_2022.pheno',data.table=F)

peur$race <- "eur"
#paam$race <- "aam"

pceur <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/eur_pca.eigenpc.txt',stringsAsFactors=F,header=T)
#pcaam <- fread('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/aam_pca.eigenpc.txt',data.table=F)

deur0 <- merge(peur,pceur,by=c("FID","IID"))
deur1 <- merge(deur0,prs,by="FID")

#aamcount <- dim(paam)[1] #Only do this many europeans
deur <- deur1#[sample(1:dim(deur1)[1], aamcount,replace=F),]


#daam0 <- merge(paam,pcaam,by=c("FID","IID"))
#daam <- merge(daam0,prs,by="FID")

deur$prsquantile <- cut(deur$TSCORE,breaks=quantile(deur$TSCORE,seq(0,1,0.2)))
#daam$prsquantile <- cut(daam$TSCORE,breaks=quantile(daam$TSCORE,seq(0,1,0.2)))

deur$prsscale <- scale(deur$TSCORE,center=TRUE,scale=TRUE)
m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+prsscale,data=deur,family='binomial'))$coefficients[7:10,1:2]

m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+prsquantile,data=deur,family='binomial'))$coefficients[7:10,1:2]


#m1aam <- summary(glm(I(PCL_sum >= 45) ~ PC1+PC2+PC3+PC4+PC5+prsquantile,data=daam,family='binomial'))$coefficients[7:10,1:2]

library(fmsb)
#NagelkerkeR2(glm(I(PCL_sum >= 45) ~ PC1+PC2+PC3+PC4+PC5+TSCORE,data=daam,family='binomial'))
NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5+TSCORE,data=deur,family='binomial'))

NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5,data=deur,family='binomial'))


save.image('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/mvp_tinprs_v1.R')


#Africans

R 
library(data.table)
library(plyr)

prs <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/prs_f3_to_mvp_allsubjects.prs',stringsAsFactors=F,header=T)

peur <- read.table('/scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_redo/tinnitus18_v2_aam_nov01_2022.pheno',stringsAsFactors=F,header=T)

peur$race <- "aam"

pceur <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/aam_pca.eigenpc.txt',stringsAsFactors=F,header=T)

deur0 <- merge(peur,pceur,by=c("FID","IID"))
deur1 <- merge(deur0,prs,by="FID")


deur <- deur1 #[sample(1:dim(deur1)[1], aamcount,replace=F),]


deur$prsquantile <- cut(deur$TSCORE,breaks=quantile(deur$TSCORE,seq(0,1,0.2)))

deur$prsscale <- scale(deur$TSCORE,center=TRUE,scale=TRUE)
m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+prsscale,data=deur,family='binomial'))$coefficients[7:10,1:2]

m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+prsquantile,data=deur,family='binomial'))$coefficients[7:10,1:2]



library(fmsb)

NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5+TSCORE,data=deur,family='binomial'))

NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5,data=deur,family='binomial'))


save.image('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/mvp_tinprs_v1_aam.R')

#hispanics
R 
library(data.table)
library(plyr)

prs <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/prs_f3_to_mvp_allsubjects.prs',stringsAsFactors=F,header=T)

peur <- read.table('/scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_redo/tinnitus18_v2_his_nov01_2022.pheno',stringsAsFactors=F,header=T)

peur$race <- "eur"

pceur <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/his_pca.eigenpc.txt',stringsAsFactors=F,header=T)

deur0 <- merge(peur,pceur,by=c("FID","IID"))
deur1 <- merge(deur0,prs,by="FID")


deur <- deur1 #[sample(1:dim(deur1)[1], aamcount,replace=F),]


deur$prsquantile <- cut(deur$TSCORE,breaks=quantile(deur$TSCORE,seq(0,1,0.2)))

deur$prsscale <- scale(deur$TSCORE,center=TRUE,scale=TRUE)
m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+prsscale,data=deur,family='binomial'))$coefficients[7:10,1:2]

m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+prsquantile,data=deur,family='binomial'))$coefficients[7:10,1:2]



library(fmsb)

NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5+TSCORE,data=deur,family='binomial'))

NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5,data=deur,family='binomial'))


save.image('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/mvp_tinprs_v1_his.R')



