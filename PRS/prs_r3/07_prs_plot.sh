#We have an average score for all chunks. We need a score summed across all chunks to get a full PRS.

#Motivation: PLINK gives average scores, but not total scores.  Because different N SNPs included in each chunk, summing averages would be wrong.We need to get an actual allelic score and not an average.

#Make a directory for these adjusted outputs 
mkdir  /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs_mvpr3/adjusted_prs


#Call into raw PRS directory and calculate PRS
cd /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs_mvpr3/

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

write.table(d1_scored, file='/group/research/mvp039/AM/tinnitus_paper2_prs_mvpr3/prs_f3_to_mvp_allsubjects.prs',row.names=F)


#Now actually do analysis


R 
library(data.table)
library(plyr)

prs <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs_mvpr3/prs_f3_to_mvp_allsubjects.prs',stringsAsFactors=F,header=T)

peur <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs_mvpr3/tinnitus18_v2_eur_nov01_2022.pheno',stringsAsFactors=F,header=T)
#paam <- fread('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs_mvpr3/PCLsum_aam_oct24_2022.pheno',data.table=F)

peur$race <- "eur"

pceur <- read.table('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs_mvpr3/eur_pca.eigenpc.txt',stringsAsFactors=F,header=T)

deur0 <- merge(peur,pceur,by=c("FID","IID"))
deur1 <- merge(deur0,prs,by="FID")

deur <- deur1#[sample(1:dim(deur1)[1], aamcount,replace=F),]


#quantilized for plotting, scaled for inference
deur$prsscale <- scale(deur$TSCORE,center=TRUE,scale=TRUE)

deur$prsquantile <- cut(deur$TSCORE,breaks=quantile(deur$TSCORE,seq(0,1,0.2)))

m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+prsquantile,data=deur,family='binomial'))$coefficients[7:10,1:2]

m1eur <- summary(glm(broad1~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+prsscale,data=deur,family='binomial'))$coefficients[7:10,1:2]



library(fmsb)

NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+TSCORE,data=deur,family='binomial'))
NagelkerkeR2(glm(broad1 ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=deur,family='binomial'))


save.image('/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs_mvpr3/mvp_tinprs_v1.R')




