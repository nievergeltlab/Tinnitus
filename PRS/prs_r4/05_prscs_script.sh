#BSUB -q short #mvp010
#BSUB -M 20000
#BSUB -R "rusage[mem=20000]"
#BSUB -W 06:00
#BSUB -o /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/errandout/prscs.o%J_%I
#BSUB -e /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/errandout/prscs.e%J_%I


#Load PRScs and PLINK
module load PRScs/0878d30
module load plink
#mkdir /scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/errandout

job=${LSB_JOBINDEX} #22 for testing
rd=/group/Public_data/PRS-CS/ldblk_1kg_eur
bfile=/scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/mvp_all.bim.sorted.prscsin
sst=/scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/tinprs.txt.prscs.sorted.linked.header1
od=/scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs/prs_f3_to_mvp

PRScs.py --ref_dir $rd --bim_prefix $bfile --chrom $job --n_gwas 481874 --sst $sst --out_dir $od

#bsub  -J  "prscs_[1-21]" -G mvp039 < 05_prscs_script.sh

#"$scratchpath"//mvp_all.bim.sorted.prscsin.bim
