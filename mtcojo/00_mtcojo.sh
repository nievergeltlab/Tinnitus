

Tinnitus adjusting for hearing

#Use 10k ukbb samples as reference individuals
ls 1000_unrel/* | grep maf | grep bed | sed 's/.bed//g' > mtcojo_ref_data.txt

#reformat data to cojo format
 #tinnitus_dat$BETA = tinnitus_dat$Zscore * sed / sqrt(2*tinnitus_dat$Weight*tinnitus_dat$Freq1*(1-tinnitus_dat$Freq1))
 #tinnitus_dat$SE <- tinnitus_dat$B / tinnitus_dat$Zscore 
zcat metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz | awk 'OFS="\t" {beta="beta";se="se"; if(NR>1) {beta=3 * $8 / sqrt(2*$6*(1-$6)*$7); se = beta/$8} ; if (NR==1) {$3="SNP"; $4="A1";$5="A2";$6="freq";$7="N";$9="p"; beta="beta"; se="se"}; print $3,toupper($4),toupper($5),$6,beta,se,$9,$7} ' | gzip > metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.cojo.gz
zcat metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz | sed 's/ -0.000 / -0.0001 /g' | sed 's/ 0.000 / 0.0001 /g' | awk 'OFS="\t" {beta="beta";se="se"; if(NR>1) {beta=3 * $8 / sqrt(2*$6*(1-$6)*$7); se = beta/$8} ; if (NR==1) {$3="SNP"; $4="A1";$5="A2";$6="freq";$7="N";$9="p"; beta="beta"; se="se"}; print $3,toupper($4),toupper($5),$6,beta,se,$9,$7} ' | gzip > metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.cojo.gz


#Remove hearing from tinnitus
echo "tinn metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.cojo.gz
hear metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.cojo.gz" > mtcojo_summary_data.list

gcta_1.93.2beta/gcta64  --mbfile mtcojo_ref_data.txt --mtcojo-file mtcojo_summary_data.list --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ --w-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ --out tinnitus_hearing_mtcojo_result

#Remove tinnitus from hearing
echo "hear metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.cojo.gz
tinn metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.cojo.gz" > mtcojo_summary_data_hear.list

gcta_1.93.2beta/gcta64  --mbfile mtcojo_ref_data.txt --mtcojo-file mtcojo_summary_data_hear.list --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ --w-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ --out hearing_tinnitus_mtcojo_result




#genetic correlation to adjusted results

tinnitus_hearing_mtcojo_result.mtcojo.cma

#THis 
code is also in rgs to psychiatric traits file:
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat tinnitus_hearing_mtcojo_result.mtcojo.cma | awk '{print $1,$2,$3}' | LC_ALL=C sort -k1b,1)  > tinnitus_hearing_mtcojo_result.mtcojo.cma.alleles

#Reformat main summary statistics
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat tinnitus_hearing_mtcojo_result.mtcojo.cma  | awk '{if (NR ==1) {$5="bbb";$6="sss";$7="ppp"};print}'  | LC_ALL=C sort -u -k1b,1 ) > tinnitus_hearing_mtcojo_result.mtcojo.cma.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats tinnitus_hearing_mtcojo_result.mtcojo.cma.premunge --p bC_pval --signed-sumstats bC,0  --out tinnitus_hearing_mtcojo_result.mtcojo.cma.munge.gz

LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz | awk '{if (NR==1) {$3="SNP"} print}' | LC_ALL=C sort -u -k3b,3 ) > metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.premunge  --merge-alleles tinnitus_hearing_mtcojo_result.mtcojo.cma.alleles  --out metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size

zcat metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz


        #Estimate rgs
   ldsc-master/ldsc.py \
  --rg tinnitus_hearing_mtcojo_result.mtcojo.cma.munge.gz.sumstats.gz,metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz \
  --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
  --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
  --out  mtcojo_vs_directadj_gcov
  
  #very high rg between methods
  zcat tinnitus_hearing_mtcojo_result.mtcojo.cma.gz | awk '{print $1,$2,$3,$4,$8,$9,$10,$11,$12}' | gzip > tinnitus_hearing_mtcojo_result.mtcojo.cma.gz.fuma.gz
  
  gzip hearing_tinnitus_mtcojo_result.mtcojo.cma
   zcat hearing_tinnitus_mtcojo_result.mtcojo.cma | awk '{print $1,$2,$3,$4,$8,$9,$10,$11,$12}' | gzip > hearing_tinnitus_mtcojo_result.mtcojo.cma.gz.fuma.gz
   
  
  #rg with hearing, after adjustment
   LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz | awk '{if (NR == 1) $3="SNP"; print}'  | LC_ALL=C sort -u -k3b,3 ) > metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.premunge
 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.premunge  --merge-alleles UKB_tinnitus/f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz.premunge.alleles  --out metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size

  ldsc-master/ldsc.py \
  --rg metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz,metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.munge.gz.sumstats.gz \
  --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
  --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
  --out  tinnitus_hearcov_vs_hearing_meta
  
  ldsc-master/ldsc.py \
  --rg metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz,metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz.munge.gz.sumstats.gz \
  --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
  --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
  --out  tinnitus_nohearcov_vs_hearing_meta
  
  #correlation to original results
      ldsc-master/ldsc.py \
  --rg metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz, metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz \
  --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
  --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
  --out  tinnitus_hearcov_vs_tinnitus_nocov
  
  
    ldsc-master/ldsc.py \
  --rg metal_results/ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz,metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz \
  --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
  --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
  --out  tinnitus_hearcov_vs_tinnitus_meta
  
  
  
  #calculate if for the significant snps, the effect size changed significantly. Calculate the covariance between z-scores, which should just be the cov of the z scores of the merged data?
  #well its not standardized... but what is the standardization factor here??
  