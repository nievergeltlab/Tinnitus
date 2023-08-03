
##Quality control - do PCs on this region call inversion status?
for chr in 8
do

 ./plink2 --bfile /mnt/ukbb/adam/ptsd/bgeneur/ukb_unrel_$chr --indep-pairwise 1000 50 0.05 --rm-dup --out temp/UKB_ptsd_eur_unrelated_$chr
 ./plink2 --bfile /mnt/ukbb/adam/ptsd/bgeneur/ukb_unrel_$chr --indep-pairwise 1000 50 0.05 --exclude  temp/UKB_ptsd_eur_unrelated_$chr.rmdup.mismatch --out temp/UKB_ptsd_eur_unrelated_$chr
 ./plink2 --bed /mnt/ukbb/preimputation_genotypes/ukb_chr"$chr"_v2.bed --bim  /mnt/ukbb/preimputation_genotypes/ukb_chr"$chr"_v2.bim --fam ukb41209_cal_chr1_v2_s488292_nophen.fam --keep UKB_ptsd_eur_unrelated_m0_july7_2020.pheno --extract temp/UKB_ptsd_eur_unrelated_"$chr".prune.in --exclude range exclusion_regions_hg19.txt --make-bed --out temp/UKB_ptsd_eur_unrelated_"$chr"

done

 ./plink2 --bfile /mnt/ukbb/adam/ptsd/bgeneur/ukb_unrel_$chr --indep-pairwise 1000 50 0.05  --maf 0.05 --rm-dup --out temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinonly
 ./plink2 --bfile /mnt/ukbb/adam/ptsd/bgeneur/ukb_unrel_"$chr" --indep-pairwise 1000 50 0.05 --chr 8 --geno 0.01 --from-bp 8000000 --to-bp 13000000 --maf 0.05 --keep UKB_ptsd_eur_unrelated_m0_july7_2020.pheno  --exclude  temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinonly.rmdup.mismatch --out temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinonly
 ./plink2 --bed /mnt/ukbb/preimputation_genotypes/ukb_chr"$chr"_v2.bed  --bim  /mnt/ukbb/preimputation_genotypes/ukb_chr"$chr"_v2.bim --fam ukb41209_cal_chr1_v2_s488292_nophen.fam --keep UKB_ptsd_eur_unrelated_m0_july7_2020.pheno  --extract temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinonly.prune.in  --make-bed --out temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinvonly
 
 
 ./flashpca_x86-64 --bfile temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinvonly --suffix ukbb_unrel_eur_flashpca_chr8_noexclusionsinvonly
 ./plink2 --bfile  temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinvonly --glm allow-no-covars cols=chrom,pos,ref,alt,ax,a1freq,machr2,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p --keep UKB_ptsd_eur_unrelated_m0_july7_2020.pheno --pheno pcsukbb_unrel_eur_flashpca_chr8_noexclusionsinvonly --out pcsukbb_unrel_eur_flashpca_chr8_noexclusionsinvonly.test
 ./plink --bfile  temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinvonly --recodeA --out temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinvonly 
 
#rs12541254 for PC1, rs4256558 for for PC2


#Determine if inversion PCs associated with outcomes

R
library(data.table)
d1 <- fread('pcsukbb_unrel_eur_flashpca_chr8_noexclusionsinvonly',data.table=F)
d2 <- fread('temp/UKB_ptsd_eur_unrelated_8_noexclusionsinvonly.raw',data.table=F)
d3 <- d2[,c(1,2,7:dim(d2)[2])]
snpcols <- names(d3)[-c(1,2)]

dm <- merge(d1,d2,by=c("FID","IID"))

for (snp in snpcols)
{
mf <- as.formula(paste("PC1 ~", snp))
print(summary(lm(mf,data=dm)))
}

summary(lm(PC1 ~ as.factor(f.22000.0.0),data=dm)
summary(lm(PC1 ~ as.factor(f.21003.0.0),data=dm))
summary(lm(PC1 ~ as.factor(f.22000.0.0),data=dm))
summary(lm(PC2 ~ as.factor(f.21003.0.0),data=dm))
summary(lm(PC2 ~ as.factor(f.22000.0.0),data=dm))

cor(dm[,-c(1,2)])

dm$f2 <- as.factor(interaction(dm$rs7004456_C ,dm$rs4240658))
png('pcsukbb_unrel_eur_flashpca_chr8_noexclusionsinvonly.png')

plot(dm$PC1,dm$PC2,col=as.numeric(dm$f2)+5)
dev.off()



plot(d1$PC3,d1$PC4)
plot(d1$PC5,d1$PC6)
plot(d1$PC7,d1$PC8)
plot(d1$PC9,d1$PC10)


dev.off()
#Not sure what it is reflecting - cpuild be non eur, coudl just be british -non british eur!
#regardless, there arent so many, so this should be relatively fine for analysis.
d2 <- d1[which(d1$PC1 < 0.03 & d1$PC3 <= 0.03),]
pdf('pcsukbb_unrel_eur_flashpca_filtered.pdf',7,7)
plot(d2$PC1,d2$PC2)
plot(d2$PC2,d2$PC3)
plot(d2$PC4,d2$PC5)
dev.off()


##Start calling inversion status
phen <- fread('UKB_ptsd_eur_unrelated_m0_july7_2020.pheno',data.table=F)
dm <- merge(d1,phen,by=c("FID","IID"))
write.table(dm,file="UKB_ptsd_eur_unrelated_m0_may13_2021.pheno",quote=F,row.names=F)


    
    
 ./plink2 --bfile /mnt/ukbb/adam/ptsd/bgeneur/ukb_unrel_$chr --indep-pairwise 1000 50 0.05  --maf 0.05 --rm-dup --out temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinonly
 
 ./plink2 --bfile /mnt/ukbb/adam/ptsd/bgeneur/ukb_unrel_"$chr" --indep-pairwise 1000 50 0.1 --chr 8 --from-bp 7934925 --to-bp 11824441  --geno 0.01  --maf 0.05 --exclude  temp/UKB_ptsd_eur_unrelated_"$chr"_noexclusionsinonly.rmdup.mismatch --out temp/UKB_ptsd_eur_unrelated_"$chr"_8p23prune

 head -n 5000 /mnt/ukbb/adam/tinnitus_gwas/phenotype/UKB_tinnitus_eur_related_may2_2019.pheno > /mnt/ukbb/adam/tinnitus_gwas/phenotype/UKB_tinnitus_eur_related_may2_2019.pheno.first5k
 
 
  ./plink2 --bfile /mnt/ukbb/preimputation_genotypes/ukb_chr"$chr"_v2  --keep  /mnt/ukbb/adam/tinnitus_gwas/phenotype/UKB_tinnitus_eur_related_may2_2019.pheno --chr 8 --from-bp 7934925 --to-bp 11824441  --geno 0.01  --maf 0.05  --make-bed --out temp/UKB_8p23
 
R 
library(invClust)
library(snpStats)
d1 <- fread('/mnt/ukbb/adam/tinnitus_gwas/phenotype/UKB_tinnitus_eur_related_may2_2019.pheno',data.table=F)

for(batch in unique(d1$f.22000.0.0))
{
 print(batch)
 d1s <- subset(d1,f.22000.0.0== batch,select=c(FID,IID))
 write.table(d1s,file='temp/d1s.subjects',quote=F,row.names=F)
 system('./plink2 --bfile temp/UKB_8p23 --keep temp/d1s.subjects --make-bed --out temp/UKB_8p23s')
 
 geno.data <- read.plink('temp/UKB_8p23s',na.strings = c("0", "-9","NA"))
 geno<-geno.data$genotypes
 roi<-data.frame(chr=8,LBP=7934925, RBP=11824441, reg= "inv1")
 annot.read<-geno.data$map
 annot<-annot.read[,c(1,2,4)]
 identical(annot[,2],colnames(geno))
 
 invcall<-invClust(roi=roi, wh = 1, geno=geno, annot=annot, dim=2)
 
 outname=batch
 
 pdf(paste('chr8p23_calls/chr8p23inversion',outname,'.pdf',sep=''),7,7)
 plot(invcall) 
 dev.off()
 
 invUnc<-as.data.frame(invcall["genotypes"])
 
 #only call 
 invUnc$callmax <- apply(invUnc[,c(1:3)],1,which.max)
 invUnc$maxprob <- apply(invUnc[,c(1:3)],1,max)

 invUnc[invUnc$maxprob <= 0.9,]$callmax <- NA #if max prob is too low, set to NA
 
 write.table(invUnc,file=paste('chr8p23_',outname,'.txt',sep=''),quote=F)
 #just in case..
 rm(geno.data)
 rm(d1s)
 }
 
#coding is 1 = NI/NI 2 = NI/I 3 = I/I 
 cat chr8p23_calls/*.txt | sort -g -k 4 | awk '{if (NR ==1) print "FID","callmax"; if (NR>1 && $4 != "callmax") print $1, $5}'  > chr8p23_inversion.txt
 
 

#hard calls
inv<-invGenotypes(invcall)

#export data

#probabilities

#leave the bad calls as na