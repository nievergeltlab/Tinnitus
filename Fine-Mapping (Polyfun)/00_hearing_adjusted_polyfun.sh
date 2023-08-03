
#Start polyfun
conda activate polyfun
 
 
#Munge Sumstats

#Premunge - capitalize A1 and A2, name MAF and N columns, change X to 23.

#The X chromosome is not present in polyfun.. just remove it.
zcat ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz | awk '{if(NR==1) { $3="SNP"; $4="A1";$5="A2";$6="MAF";$7="N"; $8="Z";$9="P"}; if (NR>1) {$4=toupper($4); $5=toupper($5)};  print}' | grep -v X > ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.premunge

python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.premunge --out ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun

#Calculate priors - using approach 1 of pre-computed prior stuff
 python /home/genetics/polyfun/extract_snpvar.py --allow-missing --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun --out ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior
 
 
 
# #Remove the not foudn SNPs - this is a bit dangerous
# R
# library(data.table)
# d1 <- fread('ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.premunge',data.table=F)
# d2 <- fread('zcat ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior.miss.gz', data.table=F)

# d1a <- subset(d1,!(SNP %in% d2$SNP))

# write.table(d1a,file='ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.premunge',quote=F,row.names=F)

# python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.premunge --out ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun



#Calculate priors - using approach 1 of pre-computed prior stuff
 python /home/genetics/polyfun/extract_snpvar.py --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun --out ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior
 
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
    # --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior \
    # --chr $chr \
    # --n 470336 \
    # --start  $start 	  \
    # --end $stop   \
    # --method susie \
    # --max-num-causal 1 \
    # --out output/finemapLD_1.inform.EUR."$snp"."$chr"."$start"."$stop".gz
    
    python /home/genetics/polyfun/finemapper.py \
        --ld "$fileld" \
    --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior \
    --chr $chr \
    --n 470336 \
    --start  $start 	  \
    --end $stop   \
    --method susie \
    --max-num-causal 2 \
    --out output_hearadj/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz
        
    # python /home/genetics/polyfun/finemapper.py \
        # --ld "$fileld" \
    # --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior \
    # --chr $chr \
    # --n 470336 \
    # --start  $start 	  \
    # --end $stop   \
    # --method susie \
    # --max-num-causal 3 \
    # --out output/finemapLD_3.inform.EUR."$snp"."$chr"."$start"."$stop".gz
            
      # python /home/genetics/polyfun/finemapper.py \
          # --ld "$fileld" \
    # --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior \
    # --chr $chr \
    # --n 470336 \
    # --start  $start 	  \
    # --end $stop   \
    # --method susie \
    # --max-num-causal 4 \
    # --out output/finemapLD_4.inform.EUR."$snp"."$chr"."$start"."$stop".gz
       

    # python /home/genetics/polyfun/finemapper.py \
        # --ld "$fileld" \
    # --sumstats ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz.polyfun_prior \
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
zcat output_hearadj/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{if (NR==1 || $15!="0") print leadsnp, $0}' > output_filtered_hearadj/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".credible

#all markers
zcat output_hearadj/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{print leadsnp, $0}' > output_filtered_hearadj/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".allsnps

done

rs3809161
rs9795522


#Creible set
cat output_filtered_hearadj/*.credible | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}' > finemapLD.inform.EUR.complete_hearingadj.txt

#all markers
cat output_filtered_hearadj/*.allsnps | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}'   > finemapLD.inform.EUR.complete_hearingadj.allsnps.txt

