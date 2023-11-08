#BSUB -q plink_gpfs #mvp010
#BSUB -M 15000
#BSUB -R "rusage[mem=15000]"
#BSUB -W 6:00
#BSUB -n 1
#BSUB -o errandout/tineur.o%J_%I
#BSUB -e errandout/tineur.e%J_%I


cd /scratch/scratch2/mvp010/vhasdcmaihoa/

jobnum=${LSB_JOBINDEX}

echo "job number is $jobnum $LSB_JOBINDEX ${LSB_JOBINDEX}"

if [ $jobnum -le 9 ]
then
 jobnum=$(echo "00"$jobnum) #add 0 padding
fi

if [ $jobnum -ge 10 ] && [ $jobnum -le 99 ]
then
 jobnum=$(echo "0"$jobnum) #add 0 padding
fi


IFS=$'\n'
 for files in $(cat /home/home1/vhasdcmaihoa/genotype_filenames/genotype_files_"$jobnum") #be sure to rename the 00 file to something else
do
 
 outfilename=$(echo $files | awk 'BEGIN{FS="/"}{ print $NF}' | awk 'BEGIN{FS="eagle."}{ print $2}' | awk 'BEGIN{FS=".minimimac"}{print $1}' ) 

anc=eur
awk '{print $3,$7}' "$files".pvar | awk 'BEGIN{FS=";";OFS="\t"}{print $1,$3}' | sed 's/AF=//g' | sed 's/R2=//g' |awk '{if($2 > 0.001 && $2 <0.999 && $3 >= 0.6 && $3 != "")print $1}'  > /scratch/scratch2/mvp010/vhasdcmaihoa/european_broad_mar172020_"$outfilename".snplist

for anc in  eur # asn # eur 
do
#The covariate file is split up by inversion status

/group/tools/plink2a/plink2a_dev_20190429 --pfile  $files \
--pheno /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/anytinnitus_"$anc"_jul8_2021.pheno \
--glm cols=chrom,pos,ref,alt,ax,a1freq,a1freqcc,alt1,firth,test,machr2,nobs,orbeta,se,ci,tz,p hide-covar firth-fallback \
--ci 0.95 --memory 15000 --threads 2 \
--extract /scratch/scratch2/mvp010/vhasdcmaihoa/european_broad_mar172020_"$outfilename".snplist \
--covar-variance-standardize \
--pheno-name any_tinnitus \
--covar /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/anytinnitus_eur_jul8_2021_iv.covar \
--out /scratch/scratch3/mvp010/vhasdcmaihoa/any_tinnitus/"$anc"_jul8_2021_ivcov_"$outfilename" 


/group/tools/plink2a/plink2a_dev_20190429 --pfile  $files \
--pheno /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/anytinnitus_"$anc"_jul8_2021_nvnv.pheno \
--glm cols=chrom,pos,ref,alt,ax,a1freq,a1freqcc,alt1,firth,test,machr2,nobs,orbeta,se,ci,tz,p hide-covar firth-fallback \
--ci 0.95 --memory 15000 --threads 2 \
--extract /scratch/scratch2/mvp010/vhasdcmaihoa/european_broad_mar172020_"$outfilename".snplist \
--covar-variance-standardize \
--pheno-name any_tinnitus \
--out /scratch/scratch3/mvp010/vhasdcmaihoa/any_tinnitus/"$anc"_jul8_2021_nvnvonlynopcs_"$outfilename" 

/group/tools/plink2a/plink2a_dev_20190429 --pfile  $files \
--pheno /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/anytinnitus_"$anc"_jul8_2021_ivnv.pheno \
--glm cols=chrom,pos,ref,alt,ax,a1freq,a1freqcc,alt1,firth,test,machr2,nobs,orbeta,se,ci,tz,p hide-covar firth-fallback \
--ci 0.95 --memory 15000 --threads 2 \
--extract /scratch/scratch2/mvp010/vhasdcmaihoa/european_broad_mar172020_"$outfilename".snplist \
--covar-variance-standardize \
--pheno-name any_tinnitus \
--out /scratch/scratch3/mvp010/vhasdcmaihoa/any_tinnitus/"$anc"_jul8_2021_ivnvonlynopcs_"$outfilename" 


done

done


#bsub  -J  "ar[74]" -G mvp010 < 02_tinnitus_gwas_mar17_2020_chr8.sh




#rerun failed jobs..


