#BSUB -J euro_pca
#BSUB -q short
#BSUB -o errandout/pca.out
#BSUB -e errandout/pcae.out



for race in eur his aam 
do


#get 10,000 random subjects of European ancestry to thin markers
sort -R tinnitus_"$race"_jul8_2021.pheno | head -n 10000 | awk '{print $1,$2}'  > 10k_random_"$race"

plink2a_dev_20190429 --bfile 20180823.GenotypeData.Release3full \
--keep 10k_random_"$race" --memory 8000 --make-bed --out  10k_random_ldprune_"$race"_1

#Calculate PCA
plink2a_dev_20190429 --bfile 10k_random_ldprune_"$race"_1 --maf 0.01 \
 --indep-pairwise 1000 50 0.05 --memory 32000 --out 10k_random_ldprune_"$race" 

plink2a_dev_20190429 --bfile 20180823.GenotypeData.Release3full  \
--keep tinnitus_"$race"_jul8_2021.pheno --extract 10k_random_ldprune_"$race".prune.in \
--make-bed --memory 5000 --out broad_tinnitus_"$race"_ldprune

flashpca2 --seed 17 --bfile broad_tinnitus_"$race"_ldprune \
--outvec "$race"_pca.eigenveg.txt --outval  "$race"_pca.eigenval.txt \
--outpc "$race"_pca.eigenpc.txt

done

#bsub -G mvp010 -J euro_pca -M 8000 -W 06:00 < 01_european_pca_mvp010.sh

