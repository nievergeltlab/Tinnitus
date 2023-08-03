#!/bin/bash

 # module load 2020
 # module load Boost/1.72.0-GCCcore-9.3.0-no_mpi
 # module load Tk/8.6.10-GCCcore-9.3.0
 
 #When you start LISA:
 LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 
 
 cd  /home/maihofer/freeze3_gwas
 
##Convert sumstats

sumstatspath=/home/maihofer/mixer/python_convert/sumstats.py

 # #Freeze 3
 # python $sumstatspath csv --sumstats results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz --auto  --frq Freq1 --pval P-value --a1 Allele1 --a2 Allele2  --head 10 --out mixer/eur_ptsd_pcs_v4_aug3_2021.fuma.gz.mixer1
 # python $sumstatspath zscore --sumstats mixer/eur_ptsd_pcs_v4_aug3_2021.fuma.gz.mixer1 | \
 # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/eur_ptsd_pcs_v4_aug3_2021.fuma.gz.mixer2
 # gzip mixer/eur_ptsd_pcs_v4_aug3_2021.fuma.gz.mixer2

 # #F3 subsets
 # #PTSD EHR only
  # python $sumstatspath csv --sumstats results_filtered/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz --auto  --frq Freq1 --pval P-value --a1 Allele1 --a2 Allele2  --head 10 --out mixer/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.mixer1
  # python $sumstatspath zscore --sumstats mixer/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.mixer1 | \
  # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.mixer2
  # gzip mixer/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.mixer2

 # #Freeze 3, no EHR, no MVP
  # python $sumstatspath csv --sumstats results_filtered/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz --auto  --frq Freq1 --pval P-value --a1 Allele1 --a2 Allele2  --head 10 --out mixer/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.mixer1
  # python $sumstatspath zscore --sumstats mixer/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.mixer1 | \
  # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.mixer2
  # gzip mixer/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.mixer2
  
  
 #MVP
  # python $sumstatspath csv --sumstats sumstats/TotalPCL_MVP_eur.gz --auto --snp rsid --frq Freq1 --pval P --a1 Allele1 --a2 Allele2 --n-val 186689  --head 10 --out mixer/TotalPCL_MVP_eur.gz.mixer1
  # python $sumstatspath zscore --sumstats mixer/TotalPCL_MVP_eur.gz.mixer1 | \
  # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/TotalPCL_MVP_eur.gz.mixer2
  # gzip mixer/TotalPCL_MVP_eur.gz.mixer2

 #Case/control 
   # python $sumstatspath csv --sumstats results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz --auto  --frq Freq1 --pval P-value --a1 Allele1 --a2 Allele2  --n neff --head 10 --out mixer/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.mixer1
   # python $sumstatspath zscore --sumstats mixer/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.mixer1 | \
   # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.mixer2
   # gzip mixer/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.mixer2
  
 #MDD
  #zcat mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz  | tail -n+76 | sed 's/#//g' > pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2
  
   # python $sumstatspath csv --sumstats pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2 --auto  --chr CHROM --bp POS --snp ID --frq FCON --pval PVAL --a1 A1 --a2 A2  --n NEFF --head 10 --out mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz.mixer1
   # python $sumstatspath zscore --sumstats mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz.mixer1 | \
   # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz.mixer2
   # gzip mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz.mixer2
  
  #MDD3
  #zcat pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz  | tail -n+76 | sed 's/#//g' > pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz2
  
   # python $sumstatspath csv --sumstats pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz2 --auto  --chr CHROM --bp POS --snp ID --frq FCON --pval PVAL --a1 A1 --a2 A2  --n NEFF --head 10 --out mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz.mixer1
   # python $sumstatspath zscore --sumstats mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz.mixer1 | \
   # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz.mixer2
   # gzip mixer/pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz.mixer2
  

 #SCZ 3 
  #zcat PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz  | tail -n+74 | sed 's/#//g' > PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz2
  
   # python $sumstatspath csv --sumstats PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz2 --auto  --chr CHROM --bp POS --snp ID --frq FCON --pval PVAL --a1 A1 --a2 A2  --n NEFF --head 10 --out mixer/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz.mixer1
   # python $sumstatspath zscore --sumstats mixer/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz.mixer1 | \
   # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz.mixer2
   # gzip mixer/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz.mixer2
  
  
  
  
 #BIP 3 
   #zcat pgc-bip2021-all.vcf.tsv.gz  | tail -n+73 | sed 's/#//g' | awk '{if(NR>1) $13=$13*2; print}' > pgc-bip2021-all.vcf.tsv.gz2
  
   python $sumstatspath csv --sumstats pgc-bip2021-all.vcf.tsv.gz2 --auto  --chr CHROM --bp POS --snp ID --frq FCON --pval PVAL --a1 A1 --a2 A2  --n NEFFDIV2 --head 10 --out mixer/pgc-bip2021-all.vcf.tsv.gz.mixer1
   python $sumstatspath zscore --sumstats mixer/pgc-bip2021-all.vcf.tsv.gz.mixer1 | \
   python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/pgc-bip2021-all.vcf.tsv.gz.mixer2
   gzip mixer/pgc-bip2021-all.vcf.tsv.gz.mixer2
  
  
  
  
  #IBS
    # python $sumstatspath csv --sumstats 81_GCST90016564_buildGRCh37.tsv --auto  --chr CHROM --bp POS --snp SNP --frq MAF --pval P --a1 A1 --a2 A2  --n samplesize --head 10 --out mixer/81_GCST90016564_buildGRCh37.tsv.mixer1
   # python $sumstatspath zscore --sumstats mixer/81_GCST90016564_buildGRCh37.tsv.mixer1 | \
   # python $sumstatspath qc  --max-or 1e37 --max-or 1e37  --exclude-ranges 6:26000000-34000000 --out mixer/81_GCST90016564_buildGRCh37.tsv.mixer2
   # gzip mixer/81_GCST90016564_buildGRCh37.tsv.mixer2
  
  
 
 
  
#sbatch --time=2:05:00 --error errandout/convert.e --output errandout/convert.o 00_prepare_files.sh --export=ALL 
