#Perform genetic correlations between a trait and all existing PGC psychiatric traits (as of july 2021)

#need to make a file where the possible alleles are given - will crash otherwise, because sometimes (rare) the same rs has different coded alleles 
#Basically just make a file consisting of the SNP, A1, A2 from your main dataset where you want to see what traits are genetically correlated to it.

LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz | awk '{if(NR==1) {$4="A1"; $5="A2"}; print $1,$4,$5}' | LC_ALL=C sort -k1b,1)  > f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles

#Reformat main summary statistics

#UKB TIn
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz    | LC_ALL=C sort -u -k1b,1 ) > f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge  --N 175995 --a1 ALLELE1 --a2 ALLELE0 --out f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.munge.gz
#MVP TIN
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat MVP_tinnitus/eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz   | awk '{if (NR==1) {$5="A2"; $8="N";$10="SE";} print}' | LC_ALL=C sort -u -k1b,1 ) > MVP_tinnitus/eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats MVP_tinnitus/eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.premunge --N 308879 --a2 AX --merge-alleles UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles --out MVP_tinnitus/eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
#META TIN
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz | awk '{if (NR == 1) $3="SNP"; print}'  | LC_ALL=C sort -u -k3b,3 ) > metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge  --merge-alleles UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles  --out metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size

#HEARING
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz | awk '{if (NR == 1) $3="SNP"; print}'  | LC_ALL=C sort -u -k3b,3 ) > metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.premunge  --merge-alleles UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles  --out metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size


#MVP HEARING
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat MVP_hearing/eur_broad_mar172020_allchr.hearing.maf01.resultsa.fuma.gz   | awk '{if (NR==1) {$5="A2"; $8="N";$10="SE";} print}' | LC_ALL=C sort -u -k1b,1 ) > MVP_hearing/eur_broad_mar172020_allchr.hearing.maf01.resultsa.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats MVP_hearing/eur_broad_mar172020_allchr.hearing.maf01.resultsa.fuma.gz.premunge --N 223797 --a2 AX --merge-alleles UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles --out MVP_hearing/eur_broad_mar172020_allchr.hearing.maf01.resultsa.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size

#UKB hearing
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat hearing_difficulty_wells_ukbb/HD_EA_gwas_sumstats.txt.gz.linprobit.fuma.gz    | LC_ALL=C sort -u -k1b,1 ) > hearing_difficulty_wells_ukbb/HD_EA_gwas_sumstats.txt.gz.linprobit.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats hearing_difficulty_wells_ukbb/HD_EA_gwas_sumstats.txt.gz.linprobit.fuma.gz.premunge  --merge-alleles UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles --N 227152 --a1 ALLELE1 --a2 ALLELE0 --out hearing_difficulty_wells_ukbb/HD_EA_gwas_sumstats.txt.gz.linprobit.fuma.gz.munge.gz

#Reformat all other summary statistics for ldsc

LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/adhd_eur_jun2017.gz   |  LC_ALL=C sort -u -k2b,2 ) > mendelian_randomization_inputs/adhd_eur_jun2017.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/adhd_eur_jun2017.gz.premunge --N 53293 --merge-alleles UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles --out mendelian_randomization_inputs/adhd_eur_jun2017.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/pgc_alcdep.eur_discovery.aug2018_release_FIXED.txt.gz  |  LC_ALL=C sort -u -k2b,2 | cut -d " " -f1-8 ) > pgc_alcdep.eur_discovery.aug2018_release.txt.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/pgc_alcdep.eur_discovery.aug2018_release.txt.gz.premunge --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/pgc_alcdep.eur_discovery.aug2018_release.txt.gz.munge.gz
 
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/pgcAN2.2019-07.vcf.tsv.gz | grep -v "#" | awk '{if (NR==1) $3="SNP"; print}'  |  LC_ALL=C sort -u -k3b,3 ) > pgcAN2.2019-07.vcf.tsv.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/pgcAN2.2019-07.vcf.tsv.gz.premunge --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --a1 REF --a2 ALT --N 72517 --out mendelian_randomization_inputs/pgcAN2.2019-07.vcf.tsv.gz.munge.gz
 
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/anxiety.meta.full.cc.tbl.gz | awk '{if (NR==1) {$1="SNP"; $10="N"}; print}'  |  LC_ALL=C sort -u -k1b,1 ) > anxiety.meta.full.cc.tbl.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/anxiety.meta.full.cc.tbl.gz.premunge --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/anxiety.meta.full.cc.tbl.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/iPSYCH-PGC_ASD_Nov2017.gz   |  LC_ALL=C sort -u -k2b,2 ) > mendelian_randomization_inputs/iPSYCH-PGC_ASD_Nov2017.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/iPSYCH-PGC_ASD_Nov2017.gz.premunge --N 46351 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/iPSYCH-PGC_ASD_Nov2017.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/daner_PGC_BIP32b_mds7a_0416a.gz   |  LC_ALL=C sort -u -k2b,2 ) > mendelian_randomization_inputs/daner_PGC_BIP32b_mds7a_0416a.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/daner_PGC_BIP32b_mds7a_0416a.gz.premunge --N 51710 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/daner_PGC_BIP32b_mds7a_0416a.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/Cannabis_ICC_UKB_het.txt.gz   |  LC_ALL=C sort -u -k2b,2 ) > mendelian_randomization_inputs/Cannabis_ICC_UKB_het.txt.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/Cannabis_ICC_UKB_het.txt.gz.premunge --ignore Z --N 162082 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/Cannabis_ICC_UKB_het.txt.gz.munge.gz
 
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat mendelian_randomization_inputs/PGC_UKB_depression_genome-wide.txt  | awk '{if (NR==1) $1="SNP"; print}'  |  LC_ALL=C sort -u -k1b,1 ) > mendelian_randomization_inputs/PGC_UKB_depression_genome-wide.txt.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/PGC_UKB_depression_genome-wide.txt.premunge --signed-sumstats LogOR,0 --N 500199 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/PGC_UKB_depression_genome-wide.txt.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/ocd_aug2017.gz   |  LC_ALL=C sort -u -k2b,2 ) > mendelian_randomization_inputs/ocd_aug2017.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/ocd_aug2017.gz.premunge --N 9725 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/ocd_aug2017.gz.munge.gz
 
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz | awk '{if (NR==1) $3="SNP"; print}'  |  LC_ALL=C sort -u -k3b,3 ) > mendelian_randomization_inputs/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.premunge --N 55374 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.munge.gz
 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/PGC3_SCZ_wave3_public.v2.tsv.gz   |  LC_ALL=C sort -u -k2b,2 ) > mendelian_randomization_inputs/PGC3_SCZ_wave3_public.v2.tsv.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/PGC3_SCZ_wave3_public.v2.tsv.gz.premunge --N 161405 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/PGC3_SCZ_wave3_public.v2.tsv.gz.munge.gz
 
 LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/SA_in_MDD_BIP_SCZ_2019.gz   |  LC_ALL=C sort -u -k2b,2 ) > mendelian_randomization_inputs/SA_in_MDD_BIP_SCZ_2019.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/SA_in_MDD_BIP_SCZ_2019.gz.premunge --N 23801 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/SA_in_MDD_BIP_SCZ_2019.gz.munge.gz

LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/TS_Oct2018.gz  |  LC_ALL=C sort -u -k1b,1 ) > mendelian_randomization_inputs/TS_Oct2018.gz.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/TS_Oct2018.gz.premunge --N 14307 --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out mendelian_randomization_inputs/TS_Oct2018.gz.munge.gz
 
 
  
 
#Hearing

LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat mendelian_randomization_inputs/GERA-EUR-ARHI.mtcojo  |  LC_ALL=C sort -u -k1b,1 ) > mendelian_randomization_inputs/GERA-EUR-ARHI.mtcojo.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/GERA-EUR-ARHI.mtcojo.premunge --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out GERA-EUR-ARHI.mtcojo.munge.gz
 

#Insomnia
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat mendelian_randomization_inputs/Insomnia_sumstats_Jansenetal.txt.gz  | awk '{if (NR == 1 || (length($5) == 1 && length($6) == 1)) print}' |  LC_ALL=C sort -u -k1b,1 ) > mendelian_randomization_inputs/Insomnia_sumstats_Jansenetal.txt.gz.premunge
zcat  metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz | awk '{print $1,$2,$3}' >  metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz.alleles
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/Insomnia_sumstats_Jansenetal.txt.gz.premunge --merge-alleles  UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles --out mendelian_randomization_inputs/Insomnia_sumstats_Jansenetal.txt.gz.munge.gz
 


#Hearing2

LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat mendelian_randomization_inputs/GERA-EUR-ARHI.mtcojo  |  LC_ALL=C sort -u -k1b,1 ) > mendelian_randomization_inputs/GERA-EUR-ARHI.mtcojo.premunge
/mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats mendelian_randomization_inputs/GERA-EUR-ARHI.mtcojo.premunge --merge-alleles f.4803.max_coding3_no_covar_related.bgen.stats.fuma.alleles --out GERA-EUR-ARHI.mtcojo.munge.gz
 



hearing_difficulty_wells_ukbb/HD_EA_gwas_sumstats.txt.gz.linprobit.fuma.gz.munge.gz


#Get genetic correlations
for analysis in  adhd_eur_jun2017.gz.munge.gz.sumstats.gz Insomnia_sumstats_Jansenetal.txt.gz.munge.gz.sumstats.gz  GERA-EUR-ARHI.mtcojo.munge.gz.sumstats.gz  SA_in_MDD_BIP_SCZ_2019.gz.munge.gz.sumstats.gz  adhd_eur_jun2017.gz.munge.gz.sumstats.gz pgc_alcdep.eur_discovery.aug2018_release.txt.gz.munge.gz.sumstats.gz pgcAN2.2019-07.vcf.tsv.gz.munge.gz.sumstats.gz anxiety.meta.full.cc.tbl.gz.munge.gz.sumstats.gz iPSYCH-PGC_ASD_Nov2017.gz.munge.gz.sumstats.gz daner_PGC_BIP32b_mds7a_0416a.gz.munge.gz.sumstats.gz Cannabis_ICC_UKB_het.txt.gz.munge.gz.sumstats.gz PGC_UKB_depression_genome-wide.txt.munge.gz.sumstats.gz ocd_aug2017.gz.munge.gz.sumstats.gz eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.munge.gz.sumstats.gz PGC3_SCZ_wave3_public.v2.tsv.gz.munge.gz.sumstats.gz  TS_Oct2018.gz.munge.gz.sumstats.gz 
do
 outname=$(echo $analysis | cut -d "." -f1)
 echo $outname
  /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg  hearing_difficulty_wells_ukbb/HD_EA_gwas_sumstats.txt.gz.linprobit.fuma.gz.munge.gz.sumstats.gz,mendelian_randomization_inputs/$analysis \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out rgs/ukbhearing_$outname
 done
  
#Send rg outputs to a single file
for analysis in Insomnia_sumstats_Jansenetal.txt.gz.munge.gz.sumstats.gz GERA-EUR-ARHI.mtcojo.munge.gz.sumstats.gz  SA_in_MDD_BIP_SCZ_2019.gz.munge.gz.sumstats.gz  adhd_eur_jun2017.gz.munge.gz.sumstats.gz pgc_alcdep.eur_discovery.aug2018_release.txt.gz.munge.gz.sumstats.gz pgcAN2.2019-07.vcf.tsv.gz.munge.gz.sumstats.gz anxiety.meta.full.cc.tbl.gz.munge.gz.sumstats.gz iPSYCH-PGC_ASD_Nov2017.gz.munge.gz.sumstats.gz daner_PGC_BIP32b_mds7a_0416a.gz.munge.gz.sumstats.gz Cannabis_ICC_UKB_het.txt.gz.munge.gz.sumstats.gz PGC_UKB_depression_genome-wide.txt.munge.gz.sumstats.gz ocd_aug2017.gz.munge.gz.sumstats.gz eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz.munge.gz.sumstats.gz PGC3_SCZ_wave3_public.v2.tsv.gz.munge.gz.sumstats.gz  TS_Oct2018.gz.munge.gz.sumstats.gz
do
outname=$(echo $analysis | cut -d "." -f1)
 echo $outname
  grep "Summary of Genetic Correlation Results" -A2  rgs/ukbhearing_$outname.log | tail -n1 >> ukbhearing_rgs.txt
   done
   
   
   #rg to meta hearing
   
  python2   /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg  UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.munge.gz.sumstats.gz,metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out rgs/ukbtinnitus_metahearing
 
  python2    /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg  MVP_tinnitus/eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.munge.gz.sumstats.gz,metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out rgs/mvptinnitus_metahearing
 
  python2    /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg  metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz,metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out rgs/metatinnitus_metahearing
 
 
 
 #

  
  