R

library(fread)

options(stringsAsFactors=F)

famfile <- read.table('20180823.GenotypeData.Release3full.fam')
names(famfile)[1] <- "mvp_010_id"

ancestry <- read.csv('MVP_010_HARE.csv',header=T)
names(ancestry)[1] <- "mvp_010_id"

pheno <- read.csv('MVP010_tinnitus_18_2.csv',header=T)
names(pheno)[1] <- "mvp_010_id"

relateds <- read.table('20180823.GenotypeData.Release3full.kin0.txt',header=T)

 

#PHENOTYPE CONSTRUCTION


#Broad definition: anyone with any indication of case status (presence of ICD or self report) is a case. All others are controls
pheno$any_tinnitus <- 1
pheno[which ( !is.na(pheno$tin_icd) | pheno$hvtin == 1),]$any_tinnitus <- 2 


#Confirm: Controls do not have ICD codes
 table(pheno$any_tinnitus,pheno$tin_icd,useNA="always")
#Confirm: Controls do not SR tinnitus
 table(pheno$any_tinnitus,pheno$hvtin,useNA="always")

#Confirm: N cases from broad tinnitus 
 Ncounts <- as.matrix(table(pheno$hvtin,pheno$tin_icd,useNA="always"))

#Confirm: Everyone who is ICD + or HVtin = 1 adds up to the number of cases based on my any_tinnitus variable
 sum(Ncounts[,-dim(Ncounts)[2]] ) + Ncounts[2,dim(Ncounts)[2]] 

#Confirm: The negation of this (the remaining 2 cell entries) adds up to the number of controls
#Note that this also implies that there are 113879 subjects who are NA for both in the broad definition of controls! 
 Ncounts[1,dim(Ncounts)[2]] + Ncounts[3,dim(Ncounts)[2]] 

#Construct PLINK versions of phenotypes
 pheno$hvtin_use <- pheno$hvtin + 1

 pheno$tin_icd_use <- 1
 pheno[which(!is.na(pheno$tin_icd)),]$tin_icd_use <- 2
 
 #Confirm match
 table(pheno$hvtin_use,pheno$hvtin,useNA="always")
 table(pheno$tin_icd_use,pheno$tin_icd,useNA="always")


#aug27,2021: hearing variables
 pheno$hvhear_use <- pheno$hvhear + 1

 pheno$hear_icd_use <- 1
 pheno[which(!is.na(pheno$hear_icd)),]$hear_icd_use <- 2
 

 pheno$any_hear <- 1
 pheno[which ( !is.na(pheno$hear_icd) | pheno$hvhear == 1),]$any_hear <- 2


 
 #Confirm match
 table(pheno$hvtin_use,pheno$hvtin,useNA="always")
 table(pheno$tin_icd_use,pheno$tin_icd,useNA="always")


#I need to remove related people, preferring to retain cases. There are a few tinnitus phenotypes. I will remove based on the broadest possible
pheno_s <- subset(pheno,select=c("mvp_010_id","any_tinnitus"))
pheno_s2 <- pheno_s
names(pheno_s)[1] <- "ID1"
names(pheno_s2)[1] <- "ID2"

r1 <- merge(pheno_s,relateds,by="ID1",suffixes=c("","_id1"))
r2 <- merge(pheno_s2,r1,by="ID2",suffixes=c("","_id2"))

r3 <- subset(r2,select=c("ID1","ID2","Kinship","any_tinnitus","any_tinnitus_id2"))


r3$prune <- NA
for (i in 1:dim(r3)[1])
{
 if(is.na(r3[i,]$any_tinnitus))
 {
  r3[i,]$prune <- r3[i,]$ID1
 } else if (is.na(r3[i,]$any_tinnitus_id2))
  {
   r3[i,]$prune <- r3[i,]$ID2
} else {
 minsub <- which.min(c(r3[i,]$any_tinnitus, r3[i,]$any_tinnitus_id2))
 r3[i,]$prune <- r3[i,minsub] #this works as long as column 1 is ID 1 and column 2 is id2, otherwise need some indexing
 rm(minsub)}
}

pheno2 <- subset(pheno, mvp_010_id %in% famfile$mvp_010_id & !(mvp_010_id %in% r3$prune)) #prune relateds  and to only people with genotypes in general


#Export phenotypes
pheno_eur <- subset(pheno2,hare=="EUR",select=c(mvp_010_id,mvp_010_id,any_tinnitus,hvtin_use,tin_icd_use))
pheno_aam <- subset(pheno2,hare=="AFR",select=c(mvp_010_id,mvp_010_id,any_tinnitus,hvtin_use,tin_icd_use))
pheno_his <- subset(pheno2,hare=="HIS",select=c(mvp_010_id,mvp_010_id,any_tinnitus,hvtin_use,tin_icd_use))
pheno_asn <- subset(pheno2,hare=="ASN",select=c(mvp_010_id,mvp_010_id,any_tinnitus,hvtin_use,tin_icd_use))

write.table(pheno_eur,file="tinnitus_eur_jul8_2021.pheno",quote=F,row.names=F)
write.table(pheno_aam,file="tinnitus_aam_jul8_2021.pheno",quote=F,row.names=F)
write.table(pheno_his,file="tinnitus_his_jul8_2021.pheno",quote=F,row.names=F)
write.table(pheno_asn,file="tinnitus_asn_jul8_2021.pheno",quote=F,row.names=F)

#At this point you should perform PCA

#hearing covariate

d1 <- read.table('eur_pca.eigenpc.txt',header=T) #PCAs
pheno <- read.csv('/data/data1/mvp010/Vinci_Data/20210707/MVP010_tinnitus_18_2.csv',header=T)
names(pheno)[1] <- "FID"

 pheno$hvhear_use <- pheno$hvhear + 1

 pheno$hear_icd_use <- 1
 pheno[which(!is.na(pheno$hear_icd)),]$hear_icd_use <- 2
 

 pheno$any_hear <- 1
 pheno[which ( !is.na(pheno$hear_icd) | pheno$hvhear == 1),]$any_hear <- 2
 
dm <- merge(d1,pheno,by="FID")

write.table(subset(dm,select=c(FID,IID,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,hvhear)),file="eur_pca.eigenpc.txt_hvhear",row.names=F,quote=F)
write.table(subset(dm,select=c(FID,IID,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,any_hear)),file="eur_pca.eigenpc.txt_anyhear",row.names=F,quote=F)
