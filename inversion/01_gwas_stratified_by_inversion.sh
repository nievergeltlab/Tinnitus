##SNP Association analysis adjusted for inversion region (chr 8 region only)


#Stratify phenotype by inversion
R

library(data.table)
d1 <- fread('phenotype/UKB_tinnitus_eur_related_may2_2019.pheno',data.table=F)
calls <- fread('/mnt/ukbb/adam/ptsd/chr8p23_inversion.txt',data.table=F)

dm <- merge(d1,calls,by="FID")

dm[which(dm$f.22000.0.0 == -9),]$ f.22000.0.0 <- 9999


#write.table(dm, file='UKB_tinnitus_eur_related_may2_2019_inversion.pheno',quote=F,row.names=F)

summary(lm(f.4803.max_coding3 ~ callmax 
             ,data=dm))
             
table(subset(dm,!is.na(f.4803.max_coding3))$callmax)


##Association analysis stratified by inversion status

  BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_related_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_related_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_related_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr8_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=f.4803.max_coding3_allcov_related_chr8inv.bed.stats.gz \
    --statsFileBgenSnps=f.4803.max_coding3_allcov_related_chr8inv.bgen.stats.gz \
    --phenoFile=phenotype/UKB_tinnitus_eur_related_may2_2019.pheno \
    --phenoCol=f.4803.max_coding3 \
    --covarFile=UKB_tinnitus_eur_related_may2_2019_inversion.pheno --covarMaxLevels=107  --qCovarCol=f.21003.max --qCovarCol=callmax  --qCovarCol=f.31.0.0 --qCovarCol=f.22009.0.{1:6} --covarCol=f.54.0.0 --covarCol=f.22000.0.0 
    
    
 
 zcat f.4803.max_coding3_allcov_related_chr8inv.bgen.stats.gz | awk '{if (NR ==1 || ($3>= 7500000 && $3<= 13500000)) print}' | sort -g -k 14 > f.4803.max_coding3_allcov_related_chr8inv_regiononly.txt
 
 
 awk '{if (NR == 1 || $18 == "1") print}' UKB_tinnitus_eur_related_may2_2019_inversion.pheno > UKB_tinnitus_eur_related_may2_2019_inversion_iviv.pheno
 awk '{if (NR == 1 || $18 == "2") print}' UKB_tinnitus_eur_related_may2_2019_inversion.pheno > UKB_tinnitus_eur_related_may2_2019_inversion_ivnv.pheno
  awk '{if (NR == 1 || $18 == "3") print}' UKB_tinnitus_eur_related_may2_2019_inversion.pheno > UKB_tinnitus_eur_related_may2_2019_inversion_nvnv.pheno
  
 
  BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_related_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_related_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_related_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr8_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=f.4803.max_coding3_allcov_related_chr8inv_iviv.bed.stats.gz \
    --statsFileBgenSnps=f.4803.max_coding3_allcov_related_chr8inv_iviv.bgen.stats.gz \
    --phenoFile=UKB_tinnitus_eur_related_may2_2019_inversion_iviv.pheno \
    --phenoCol=f.4803.max_coding3 \
    --covarFile=UKB_tinnitus_eur_related_may2_2019_inversion_iviv.pheno --covarMaxLevels=107  --qCovarCol=f.21003.max  --qCovarCol=f.31.0.0 --qCovarCol=f.22009.0.{1:6} --covarCol=f.54.0.0 --covarCol=f.22000.0.0 
    
  BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_related_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_related_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_related_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr8_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=f.4803.max_coding3_allcov_related_chr8inv_ivnv.bed.stats.gz \
    --statsFileBgenSnps=f.4803.max_coding3_allcov_related_chr8inv_ivnv.bgen.stats.gz \
    --phenoFile=UKB_tinnitus_eur_related_may2_2019_inversion_ivnv.pheno \
    --phenoCol=f.4803.max_coding3 \
    --covarFile=UKB_tinnitus_eur_related_may2_2019_inversion_ivnv.pheno --covarMaxLevels=107  --qCovarCol=f.21003.max  --qCovarCol=f.31.0.0 --qCovarCol=f.22009.0.{1:6} --covarCol=f.54.0.0 --covarCol=f.22000.0.0 
        
     BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_related_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_related_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_related_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr8_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=f.4803.max_coding3_allcov_related_chr8inv_nvnv.bed.stats.gz \
    --statsFileBgenSnps=f.4803.max_coding3_allcov_related_chr8inv_nvnv.bgen.stats.gz \
    --phenoFile=UKB_tinnitus_eur_related_may2_2019_inversion_nvnv.pheno \
    --phenoCol=f.4803.max_coding3 \
    --covarFile=UKB_tinnitus_eur_related_may2_2019_inversion_nvnv.pheno --covarMaxLevels=107  --qCovarCol=f.21003.max  --qCovarCol=f.31.0.0 --qCovarCol=f.22009.0.{1:6} --covarCol=f.54.0.0 --covarCol=f.22000.0.0 
     