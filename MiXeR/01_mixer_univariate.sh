#!/bin/bash


 
 #module load 2019
 #module load Boost.Python/1.67.0-intel-2019b-Python-3.6.6
 #module load Tk/8.6.8-GCCcore-8.3.0 
 
 #When you start LISA:
 #LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 
 cd  /home/maihofer/freeze3_gwas/mixer/
 
#Limited SNP subset analysis
if [ ! -f "$study".rep${SLURM_ARRAY_TASK_ID}.json ] 
then
python3 /home/maihofer/mixer/precimed/mixer.py fit1 \
      --trait1-file "$study".mixer2.gz \
      --out "$study".rep${SLURM_ARRAY_TASK_ID} \
      --extract /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so 
fi


#Full SNp analysis
if [ ! -f "$study".mixer_results.test.rep${SLURM_ARRAY_TASK_ID}.json ]
then

python3 /home/maihofer/mixer/precimed/mixer.py test1 \
      --trait1-file "$study".mixer2.gz \
      --load-params-file "$study".rep${SLURM_ARRAY_TASK_ID}.json \
      --out "$study".mixer_results.test.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so 
fi

# study=eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz
#   study=pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz
  

#Example code for running 1 study

# study=eur_ptsd_pcs_v4_aug3_2021.fuma.gz.mixer2
# study=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.mixer2
# study=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.mixer2
# study=TotalPCL_MVP_eur.gz.mixer2
#study=81_GCST90016564_buildGRCh37.tsv


# sbatch --array=1-20 --time=5:05:00 --error errandout/mixer_13_"$study"_%a.e --output errandout/mixer_13_"$study"_%a.o --export=ALL,study="$study" 01_mixer_univariate.sh 

#example code for running many 

# for study in pgc-bip2021-all.vcf.tsv.gz PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz 
# do
 # sbatch --array=1-20 --time=5:05:00 --error errandout/mixer_13_"$study"_%a.e --output errandout/mixer_13_"$study"_%a.o --export=ALL,study="$study" 01_mixer_univariate.sh 
# done
