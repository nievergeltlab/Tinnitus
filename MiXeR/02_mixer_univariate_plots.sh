#!/bin/bash

 
 #module load 2019
 #module load Boost.Python/1.67.0-intel-2019b-Python-3.6.6
 #module load Tk/8.6.8-GCCcore-8.3.0 
 
 #When you start LISA:
 #LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 
 cd  /home/maihofer/freeze3_gwas/mixer
 
# study=eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz
#   study=pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz
  
 for study in pgc-bip2021-all.vcf.tsv.gz PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz  #  eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz # eur_ptsd_pcs_v4_aug3_2021.fuma.gz eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz TotalPCL_MVP_eur.gz
 do

#Generate plots using just the subset snps 
python3 /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$study".rep@.json --out "$study".mixer_results.fit
python3 /home/maihofer/mixer/precimed/mixer_figures.py one --json "$study".mixer_results.fit.json --out "$study".mixer_results.fit.plots --statistic mean std
 
#Generate plots based on whole model 
python3 /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$study".mixer_results.test.rep@.json --out "$study".mixer_results.fitB
python3 /home/maihofer/mixer/precimed/mixer_figures.py one --json "$study".mixer_results.fitB.json --out "$study".mixer_results2.test.plots --statistic mean std

done

# studies
# eur_broad_mar172020_allchr.hearing.maf01.resultsa.fuma.gz 
# f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz 
# HD_EA_gwas_sumstats.txt.gz 
# ukbbcoding3_related_nocov_mvp_eur_broad_nov17_2020.tbl.gz 
# ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz 
# eur_broad_mar172020_allchr.broad_ehr_mixer.maf01.resultsa.gz
study=PGC_UKB_depression_genome-wide.txt_2

#can just run all at once
 # sbatch --time=1:05:00 --error errandout/mixeruniplots_13_%a.e --output errandout/mixeruniplots_13_%a.o --export=ALL 02_mixer_univariate_plots.sh 


