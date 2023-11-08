#UPDATED CODE (used for age analysis):

race=eur

#Doble check: make sure all jobs completed 
grep "Successfully completed" errandout/gwas_age_eur.o902902_* | wc -l

#Get number of SNPs tested
cat european_broadtinnitus_hearonly_mar172020_*.snplist | wc -l 

#Make file of GWAS results
cat "$race"_broad_withage_mar172020_chr*.broad_tinnitus.glm.logistic.hybrid | sed 's/#//g' | awk '{if (NR==1 || ($1 != "CHROM" )) print}' > "$race"_broad_withage_mar172020_allchr.broad_tinnitus.maf001.resultsa

#Chek line length of results. Should match SNPs tested (or virtually match)
wc -l "$race"_broad_withage_mar172020_allchr.broad_tinnitus.maf001.resultsa

#Filter to info > 0.6 and maf > 1%
cat "$race"_broad_withage_mar172020_allchr.broad_tinnitus.maf001.resultsa | awk '{if (NR==1 || ($1 != "CHROM" && $9 >= 0.01 && $9 <= 0.99 && $12 > 0.6)) print}' | sort -g -k 21 > "$race"_broad_withage_mar172020_allchr.broad_tinnitus.maf01.resultsa


#maf 01 output should have ~ 9 million lines for europenas

wc -l "$race"_broad_withage_mar172020_allchr.broad_tinnitus.maf01.resultsa

cp "$race"_broad_withage_mar172020_allchr.broad_tinnitus.maf001.resultsa /group/research/mvp010/sums_nov18_2020/.
cp "$race"_broad_withage_mar172020_allchr.broad_tinnitus.maf01.resultsa /group/research/mvp010/sums_nov18_2020/.


