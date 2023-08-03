### Tinnitus UKBB GWAS ###

#Note 1: These PLINK files are unrelated europeans. 
 #For analysis of related europeans and non europeans, use files in ../ancestry_unzip
  #In that case, default fam file pheno might not work, replace pheno with 1 value so bolt runs, 
  #i.e. cat ../ancestry_unzip/ukb_chr1_v2.fam | awk '{print $1,$2,$3,$4,$5,"1"}' > ukb_chr1_v2x.fam 
  #Then use this .fam file for the --fam line

#Note 2 (IMPORTANT):

#Bolt will perform an analysis of the genotyped SNPS and save it to it's own file. 
#You need to specify output names for this file (statsFile) AND for the output file based on imputed data  (statsFileBgenSNps)!

#Note 3: If you use --bed/bim/fam files that are subset of the total subjects, need to add --noBgenIDcheck to force analysis


#Phenotype files:
#Unrelated subjects phenotype file is called: UKB_tinnitus_eur_unrelated_april2_2019.pheno
#Related subjects phenotype file is called: UKB_tinnitus_eur_related_april2_2019.pheno

#Phenotypes
#f.4803.max_coding1 = Neale lab tinnitus (Yes, now or most all of the time vs no tinnitus. Lower levels of tinnitus are REMOVED)
#f.4803.max_coding2 = Clifford binary coding (All forms of tinnitus endorsement vs no tinnitus)
#f.4803.max_coding3 = Clifford continuous trait tinnitus coding 

#Covariates
#f.2247.max_coding1 = Hearing loss at time of highest recorded tinnitus
#f.21003.max = Age at time of highest recorded tinnitus
#f.31.0.0 = Sex 
#f.22009.0.1, f.22009.0.2, f.22009.0.3, f.22009.0.4, f.22009.0.5 = UKB Principal components 1-5
#f.22000.0.0 = genotyping batch (MUST be run as a categorical covariate)
#f.54.0.0 assessment center  (MUST be run as a categorical covariate)

#Covariate usage:
#Notice that categorical and quant covariates have slightly different flags
#See bolt manual

#Main GWAS
  BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_unrelated_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_unrelated_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_unrelated_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{1..22}_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=male_female_h2_rg/f.4803.max_coding1_nocov_related.bed.stats.gz \
    --statsFileBgenSnps=male_female_h2_rg/f.4803.max_coding1_allcov_related.bgen.stats.gz \
    --phenoFile=phenotype/UKB_tinnitus_eur_related_may2_2019.pheno \
    --phenoCol=f.4803.max_coding1 --covarFile=UKB_tinnitus_eur_unrelated_may2_2019.pheno --covarMaxLevels=107  --qCovarCol=f.2247.max_coding1 --qCovarCol=f.21003.max --qCovarCol=f.31.0.0 --qCovarCol=f.22009.0.{1:6} --covarCol=f.54.0.0 --covarCol=f.22000.0.0 \
 
 
  #Alternative: hearing covariate
  BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_unrelated_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_unrelated_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_unrelated_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{1..22}_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=f.4803.max_coding3_hearcov_related.bed.stats.gz \
    --statsFileBgenSnps=f.4803.max_coding3_hearcov_related.bgen.stats.gz \
    --phenoFile=phenotype/UKB_tinnitus_eur_related_may2_2019.pheno \
    --phenoCol=f.4803.max_coding3 \
    --covarFile=phenotype/UKB_tinnitus_eur_related_may2_2019.pheno  --covarMaxLevels=107  --qCovarCol=f.2247.max_coding1 --qCovarCol=f.21003.max --qCovarCol=f.31.0.0 --qCovarCol=f.22009.0.{1:6} --covarCol=f.54.0.0 --covarCol=f.22000.0.0 
  
 


#Tinnitus without subjects with/without hearing loss

 #Stratify phenotype based on hearing loss
 awk '{if(NR==1 || $6 == 0) print}' UKB_tinnitus_eur_related_april2_2019.pheno > UKB_tinnitus_eur_related_april2_2019_nohearingloss.pheno
 awk '{if(NR==1 || $6 == 1) print}' UKB_tinnitus_eur_related_april2_2019.pheno > UKB_tinnitus_eur_related_april2_2019_hearingloss.pheno
  
 #GWAS in subjects without hearing loss
  BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_unrelated_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_unrelated_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_unrelated_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{1..22}_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=f.4803.max_coding3_nocov_related_nohearingloss.bed.stats.gz \
    --statsFileBgenSnps=f.4803.max_coding3_nocov_related_nohearingloss.bgen.stats.gz \
    --phenoFile=UKB_tinnitus_eur_related_april2_2019_nohearingloss.pheno \
    --phenoCol=f.4803.max_coding3 \

 #GWAS in subjects with hearing loss
  BOLT-LMM_v2.3.2/bolt \
    --bed=pca/UKB_tinnitus_eur_unrelated_{1..22}.bed \
    --bim=pca/UKB_tinnitus_eur_unrelated_{1..22}.bim \
    --fam=pca/UKB_tinnitus_eur_unrelated_1.fam  \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --geneticMapFile=BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --lmmForceNonInf \
    --numThreads=6 \
    --bgenFile=/mnt/ukbb/adam/bgen/ukb_imp_chr{1..22}_v3.bgen \
    --bgenMinMAF=1e-3 \
    --bgenMinINFO=0.3 \
    --sampleFile=ukb40951_imp_chr1_v3_s487324.sample \
    --verboseStats \
    --modelSnps all.pruned \
    --remove bolt.in_plink_but_not_imputed.FID_IID.968.txt \
    --noBgenIDcheck \
    --statsFile=f.4803.max_coding3_nocov_related_hearingloss.bed.stats.gz \
    --statsFileBgenSnps=f.4803.max_coding3_nocov_related_hearingloss.bgen.stats.gz \
    --phenoFile=UKB_tinnitus_eur_related_april2_2019_hearingloss.pheno \
    --phenoCol=f.4803.max_coding3


#Prune analysis to relevant results, make summary data for fuma
 zcat f.4803.max_coding3_nocov_unrelated.bgen.stats.gz |  awk '{ if (NR==1) $14="P"; if (NR ==1 || ($7 >= 0.01 && $7 <= 0.99 && $8 > 0.6)) print $1,$2,$3,$5,$6,$7,$8,$11,$12,$14}'  | sort -g -k 10 | gzip > f.4803.max_coding3_nocov_unrelated.bgen.stats.fuma.gz
 cut -d " " -f 5 UKB_tinnitus_eur_unrelated_feb27_2019.pheno  | grep -v NA | wc -l # This gets N. YOU MUST WRITE IN CORRECT N BELOW!
 zcat f.4803.max_coding3_nocov_unrelated.bgen.stats.fuma.gz | awk 'BEGIN{OFS="\t"}{if (NR == 1) {N="N"; $4="A1";$5="A2";$6="FRQ";} else N=155821; print $1,$4,$5,$6,$7,$8,$9,$10,N}' > f.4803.max_coding3_nocov_unrelated.bgen.stats.ldsc

