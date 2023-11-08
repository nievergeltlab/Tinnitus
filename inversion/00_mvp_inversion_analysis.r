library(data.table)
pheno <- fread('tinnitus_eur_jul8_2021_royce.pheno',data.table=F)
names(pheno)[1] <- c("FID")
ancpcs <- fread('eur_pca.eigenpc.txt',data.table=F)
da <- fread("Inversion_calling_8p22_MVPRel3full.txt",data.table=F) #inversion status
names(da)[1] <- c("FID")

d1 <- merge(pheno,ancpcs,by=c("FID"))
dm <- merge(d1,da,by="FID")

#additive coding
dm$inv2 <- as.numeric(as.factor(dm$inversion)) -1

#Compare inversion status
j1 <- summary(glm(any_tinnitus - 1 ~ as.factor(inversion) + PC1 +PC2 + PC3 +PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dm,family="binomial"))
summary(glm(any_tinnitus - 1 ~ as.factor(inversion) ,data=dm,family="binomial"))

#additive coding
j1A <- summary(glm(any_tinnitus - 1 ~ inv2 ,data=dm,family="binomial"))

#Compare only non heterozygotes
#j2 <- summary(glm(any_tinnitus - 1 ~ inv2 + PC1 +PC2 + PC3 +PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=subset(dm,inversion != "IV|NI"),family="binomial"))


capture.output(j1A,file="/home/home1/vhasdcmaihoa/mvp010shared/sums_oct7_2021/eur_inversion_glm_additive_v2.txt")

capture.output(table(dm$inversion,dm$any_tinnitus),file="/home/home1/vhasdcmaihoa/mvp010shared/sums_oct7_2021/eur_inversion_freqs.txt")




#basically the same inferences, the extra parameter is not better than an additive coding
m1 <- glm(any_tinnitus - 1 ~ inversion + PC1 +PC2 + PC3 +PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dm,family="binomial")
m0 <- glm(any_tinnitus - 1 ~  + PC1 +PC2 + PC3 +PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dm,family="binomial")
library(lmtest)

lrtest(m1,m0)

#Export phenotype data stratified by inversion type then do analysis. How do handle heterozygotes?
write.table(subset(dm,inv2==0,select=c("FID","IID","any_tinnitus")),file="anytinnitus_eur_jul8_2021_iviv.pheno",quote=F,row.names=F)
write.table(subset(dm,inv2==1,select=c("FID","IID","any_tinnitus")),file="anytinnitus_eur_jul8_2021_ivnv.pheno",quote=F,row.names=F)
write.table(subset(dm,inv2==2,select=c("FID","IID","any_tinnitus")),file="anytinnitus_eur_jul8_2021_nvnv.pheno",quote=F,row.names=F)

write.table(subset(dm,inv2==0,select=c(FID,IID,inv2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)),file="anytinnitus_eur_jul8_2021_iv.covar",quote=F,row.names=F) 

#At this point run GWAS stratified by inversion status

#Concatenate inversion status stratified  GWAS results
cat /scratch/scratch3/mvp010/vhasdcmaihoa/any_tinnitus/eur_jul8_2021_ivivonly_chr8.8.1.20000000.any_tinnitus.glm.logistic.hybrid | awk '{if (NR==1 || ($1 != "CHROM" && $9 >= 0.01 && $9 <= 0.99 && $12 > 0.6 && $2 >= 6000000 && $2 <= 16000000)) print}' | sort -g -k 21 > /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/eur_jul8_2021_ivivonlynopcs_chr8.8.1.20000000.any_tinnitus.glm.logistic.hybrid.maf01.invregion

cat /scratch/scratch3/mvp010/vhasdcmaihoa/any_tinnitus/eur_jul8_2021_ivnvonly_chr8.8.1.20000000.any_tinnitus.glm.logistic.hybrid | awk '{if (NR==1 || ($1 != "CHROM" && $9 >= 0.01 && $9 <= 0.99 && $12 > 0.6 && $2 >= 6000000 && $2 <= 16000000)) print}' | sort -g -k 21 > /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/eur_jul8_2021_ivnvonlynopcs_chr8.8.1.20000000.any_tinnitus.glm.logistic.hybrid.maf01.invregion

cat /scratch/scratch3/mvp010/vhasdcmaihoa/any_tinnitus/eur_jul8_2021_nvnvonly_chr8.8.1.20000000.any_tinnitus.glm.logistic.hybrid | awk '{if (NR==1 || ($1 != "CHROM" && $9 >= 0.01 && $9 <= 0.99 && $12 > 0.6 && $2 >= 6000000 && $2 <= 16000000)) print}' | sort -g -k 21 > /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/eur_jul8_2021_nvnvonlynopcs_chr8.8.1.20000000.any_tinnitus.glm.logistic.hybrid.maf01.invregion


