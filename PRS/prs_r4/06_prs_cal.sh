
#BSUB -q short #mvp039
#BSUB -M 20000
#BSUB -R "rusage[mem=20000]"
#BSUB -W 06:00
#BSUB -o /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/errandout/prscal.o%J_%I
#BSUB -e /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/errandout/prscal.e%J_%I

scratchpath=/scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs


#Concatenate PRS results

#cat /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/prs_f3_to_mvp_pst_eff_a1_b0.5_phiauto_chr*.txt |  awk '{print $2,$4,$6}' > prs_f3_to_mvp__pst_eff_a1_b0.5_phiauto_allchr.txt 

#IDs need to be linked back to MVP IDs...
#LC_ALL=C join <(awk '{print $1,$2}' "$mvpprspath"//tinprs.txt.prscs.sorted.linked) <(LC_ALL=C sort -k1b,1 prs_f3_to_mvp__pst_eff_a1_b0.5_phiauto_allchr.txt ) | awk '{print $2,$3,$4}' >  prs_f3_to_mvp__pst_eff_a1_b0.5_phiauto_allchr.txt.mvp

module load plink2a

cd /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs

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
 for files in $(cat /group/research/mvp039/genotype_filenames/x"$jobnum") #be sure to rename the 00 file to something else
do
 
 outfilename=$(echo $files | awk 'BEGIN{FS="/"}{print $NF}'| sed 's/.dose//g') 


group=all #males #females

for anc in eur # asn #eur #aam his
do

#Make sure that the phenotype and covariate file names, the phenotype name, match with what you are doing
#Make sure that the output file path directory exists
if [ $group == "all" ]
then
plink2a --pfile  $files --memory 20000 --threads 6  \
--score /home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs/prs_f3_to_mvp__pst_eff_a1_b0.5_phiauto_allchr.txt.mvp \
--out /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/prs_f3_to_mvp_"$outfilename" 

fi



done

done

#Job running scripts. DO NOT UNCOMMENT IN THE FILE ITSELF!!!!

#Test job
##bsub  -J  "prscal_[1]" -G mvp039 < 06_prs_cal.sh

#Run all jobs



