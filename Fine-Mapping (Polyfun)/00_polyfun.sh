
#install polyfun
 git clone https://github.com/omerwe/polyfun
 cd polyfun
 sudo /usr/local/bin/anaconda3/condabin/conda env create -f polyfun.yml
 python test_polyfun.py #seems to run

#Get existing functional annotations
 #wget https://data.broadinstitute.org/alkesgroup/LDSCORE/baselineLF_v2.2.UKB.polyfun.tar.gz
 

#Start polyfun
conda activate polyfun
 
 
#Munge Sumstats

#Premunge - capitalize A1 and A2, name MAF and N columns, change X to 23.

#The X chromosome is not present in polyfun.. just remove it.
zcat ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz | awk '{if(NR==1) { $3="SNP"; $4="A1";$5="A2";$6="MAF";$7="N"; $8="Z";$9="P"}; if (NR>1) {$4=toupper($4); $5=toupper($5)};  print}' | grep -v X > ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge

python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun

#Calculate priors - using approach 1 of pre-computed prior stuff
 python /home/genetics/polyfun/extract_snpvar.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior
 
 
 
#Remove the not foudn SNPs - this is a bit dangerous
R
library(data.table)
d1 <- fread('ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge',data.table=F)
d2 <- fread('zcat ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior.miss.gz', data.table=F)

d1a <- subset(d1,!(SNP %in% d2$SNP))

write.table(d1a,file='ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge',quote=F,row.names=F)

python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun



#Calculate priors - using approach 1 of pre-computed prior stuff
 python /home/genetics/polyfun/extract_snpvar.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior
 
mkdir output

#We should set max-num causal to 1, per author recommendations as we don't ahve good LD, but this generates NULL results! 
dos2unix tinnitus_finemaplist_withrefs.csv

#Need SNPVAR column
IFS=$'\n'
for snpset in $(cat tinnitus_finemaplist_withrefs.csv | head -n2 )
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')


  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".gz
  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".npz
  
    # python /home/genetics/polyfun/finemapper.py \
    # --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
    # --chr $chr \
    # --n 470336 \
    # --start  $start 	  \
    # --end $stop   \
    # --method susie \
    # --max-num-causal 1 \
    # --out output/finemapLD_1.inform.EUR."$snp"."$chr"."$start"."$stop".gz
    
    python /home/genetics/polyfun/finemapper.py \
        --ld "$fileld" \
    --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
    --chr $chr \
    --n 470336 \
    --start  $start 	  \
    --end $stop   \
    --method susie \
    --max-num-causal 2 \
    --out output/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz
        
    # python /home/genetics/polyfun/finemapper.py \
        # --ld "$fileld" \
    # --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
    # --chr $chr \
    # --n 470336 \
    # --start  $start 	  \
    # --end $stop   \
    # --method susie \
    # --max-num-causal 3 \
    # --out output/finemapLD_3.inform.EUR."$snp"."$chr"."$start"."$stop".gz
            
      # python /home/genetics/polyfun/finemapper.py \
          # --ld "$fileld" \
    # --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
    # --chr $chr \
    # --n 470336 \
    # --start  $start 	  \
    # --end $stop   \
    # --method susie \
    # --max-num-causal 4 \
    # --out output/finemapLD_4.inform.EUR."$snp"."$chr"."$start"."$stop".gz
       

    # python /home/genetics/polyfun/finemapper.py \
        # --ld "$fileld" \
    # --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
    # --chr $chr \
    # --n 470336 \
    # --start  $start 	  \
    # --end $stop   \
    # --method susie \
    # --max-num-causal 5 \
    # --out output/finemapLD_5.inform.EUR."$snp"."$chr"."$start"."$stop".gz
      
        
      rm "$fileld".gz "$fileld".npz
    
      
done


IFS=$'\n'
for snpset in $(cat tinnitus_finemaplist_withrefs.csv | head -n2)
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')

#Get markers in credible set
zcat output/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{if (NR==1 || $15=="1") print leadsnp, $0}' > output_filtered/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".credible

#all markers
zcat output/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{print leadsnp, $0}' > output_filtered/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".allsnps

done


#Creible set
cat output_filtered/*.credible | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}' > finemapLD.inform.EUR.complete.txt

#all markers
cat output_filtered/*.allsnps | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}'   > finemapLD.inform.EUR.complete.allsnps.txt




### transancestry analysis

#Premunge - capitalize A1 and A2, name MAF and N columns, change X to 23.


zcat ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz | awk '{if(NR==1) {$1="CHR"; $2="POS"; $3="SNP"; $4="A1";$5="A2";$6="MAF";$7="N"; $8="Z";$9="P"}; if (NR>1) {$4=toupper($4); $5=toupper($5)};  print  $1"_"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9}' |  LC_ALL=C sort -u -k1,1 | cut -d " " -f2- | sort -g -k 9 > ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.premunge


python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.premunge --out ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.polyfun

python /home/genetics/polyfun/extract_snpvar.py --sumstats ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.polyfun --out ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.polyfun_prior
 
 
 
#Remove the not foudn SNPs - this is a bit dangerous
R
library(data.table)
d1 <- fread('ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.premunge',data.table=F)
d2 <- fread('zcat ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.polyfun_prior.miss.gz', data.table=F)

d1a <- subset(d1,!(SNP %in% d2$SNP))

write.table(d1a,file='ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.premunge',quote=F,row.names=F)

python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.premunge --out ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.polyfun

#Calculate priors - using approach 1 of pre-computed prior stuff
 python /home/genetics/polyfun/extract_snpvar.py --sumstats ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.polyfun --out ukbbcoding3_plus_mvp_transethnic_broad_jul8_2021_withhet1.tbl2.fuma.gz.polyfun_prior
 
 
 
#We should set max-num causal to 1, per author recommendations as we don't ahve good LD, but this generates NULL results! 
dos2unix trans_finemaplist_withrefs.csv

#Need SNPVAR column
IFS=$'\n'
for snpset in $(cat cat trans_finemaplist_withrefs.csv | tail -n 16 | sed 's/"//g' | grep rs823154 )
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')


  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".gz
  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".npz
  

    python /home/genetics/polyfun/finemapper.py \
        --ld "$fileld" \
    --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
    --chr $chr \
    --n 484874 \
    --start  $start 	  \
    --end $stop   \
    --method susie \
    --max-num-causal 2 \
    --out output/finemapLD_2.inform.TRANS."$snp"."$chr"."$start"."$stop".gz
     
   
      rm "$fileld".gz "$fileld".npz
    
      
done


IFS=$'\n'
for snpset in $(cat trans_finemaplist_withrefs.csv | sed 's/"//g' | grep rs823154 )
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')

#Get markers in credible set
zcat output/finemapLD_2.inform.TRANS."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{if (NR==1 || $15=="1") print leadsnp, $0}' > output_filtered/finemapLD.inform.TRANS."$snp"."$chr"."$start"."$stop".credible

#all markers
zcat output/finemapLD_2.inform.TRANS."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{print leadsnp, $0}' > output_filtered/finemapLD.inform.TRANS."$snp"."$chr"."$start"."$stop".allsnps

done


#Creible set
cat output_filtered/*.TRANS.*.credible | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}' > finemapLD.inform.TRANS.complete.txt

#all markers
cat output_filtered/*.TRANS.*.allsnps | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}'   > finemapLD.inform.TRANS.complete.allsnps.txt

