#!/bin/bash

 
 #module load 2019
 #module load Boost.Python/1.67.0-intel-2019b-Python-3.6.6
 #module load Tk/8.6.8-GCCcore-8.3.0 
 
 #When you start LISA:
 #LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 
 # study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
 # study2=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz
 
 # study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
 # study2=TotalPCL_MVP_eur.gz

 # study1=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz
 # study2=TotalPCL_MVP_eur.gz

# study1=eur_ptsd_pcs_v4_aug3_2021.fuma.gz 
# study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
# study1=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz 
# study1=TotalPCL_MVP_eur.gz
# study1=eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz

#study1=pgc-bip2021-all.vcf.tsv.gz 
#study2=PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz 
#study2=pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz

# study=pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz
  
#study1=


#study2=pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz
#study2=ckqny.scz2snpres.af
#study2=daner_PGC_BIP32b_mds7a_0416a.gz
#study2=daner_PGC_MDD_noUKB_no23andMe.txt
#study2=PGC_UKB_depression_genome-wide.txt_2
  
cd  /home/maihofer/freeze3_gwas/mixer

##Convert sumstats

#Studies:



studycomb="$study1"_"$study2"
python3 /home/maihofer/mixer/precimed/mixer.py fit2 \
      --trait1-file "$study1".mixer2.gz \
      --trait2-file "$study2".mixer2.gz \
      --trait1-params-file "$study1".rep${SLURM_ARRAY_TASK_ID}.json \
      --trait2-params-file "$study2".rep${SLURM_ARRAY_TASK_ID}.json \
      --out "$studycomb".rep${SLURM_ARRAY_TASK_ID} \
      --extract  /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file  /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so \
      

python3 /home/maihofer/mixer/precimed/mixer.py  test2 \
      --trait1-file "$study1".mixer2.gz \
      --trait2-file "$study2".mixer2.gz \
      --load-params-file "$studycomb".rep${SLURM_ARRAY_TASK_ID}.json \
      --out "$studycomb".test.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so 
      
       # sbatch --array=1-20  --time=16:05:00 --error errandout/bivariate_"$study1"_"$study2"_%a.e --output errandout/bivariate_"$study1"_"$study2"_%a.o --export=ALL,study1=$study1,study2=$study2 03_mixer_bivariate.sh 
