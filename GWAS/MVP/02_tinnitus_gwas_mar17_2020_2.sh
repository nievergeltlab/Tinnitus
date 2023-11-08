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
 for files in $(cat genotype_files_"$jobnum") #Increment over split genotype files (since this is done as an array, this loop only contains one file)
do
 
 outfilename=$(echo $files | awk 'BEGIN{FS="/"}{ print $NF}' | awk 'BEGIN{FS="eagle."}{ print $2}' | awk 'BEGIN{FS=".minimimac"}{print $1}' ) 


for anc in   eur aam his
do
#Note: june3 2021, replaced eiegenpc cov with hearing, this is for analysis conditioned on hearing. change back for replicating tinnitus gwas
/group/tools/plink2a/plink2a_dev_20190429 --pfile  $files \
--pheno /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/tinnitus_"$anc"_jul8_2021.pheno \
--covar /home/home1/vhasdcmaihoa/mvp010shared/AM/royce_gwases/tinnitus_v2/"$anc"_pca.eigenpc.txt \
--glm cols=chrom,pos,ref,alt,ax,a1freq,a1freqcc,alt1,firth,test,machr2,nobs,orbeta,se,ci,tz,p hide-covar firth-fallback \
--ci 0.95 --memory 15000 --threads 2 \
--extract /scratch/scratch2/mvp010/vhasdcmaihoa/european_broad_mar172020_"$outfilename".snplist \
--covar-variance-standardize \
--pheno-name any_tinnitus \
--out /scratch/scratch3/mvp010/vhasdcmaihoa/any_tinnitus/"$anc"_jul8_2021_"$outfilename" 
done

done


#bsub  -J  "ar[1-154]" -G mvp010 < 02_tinnitus_gwas_mar17_2020_2.sh



