#Convert results to SNP locations: SNP ID, CHR, BP
zcat metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz | awk '{print $3,$1,$2}' > metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.locations
zcat metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz | awk '{if(NR==1) {$3="SNP"; $7="N"; $9="P"}; print $3,$7,$9}' > metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.pvalues

#The B37 data that comes with MAGMA
./magma --annotate --snp-loc metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.locations --gene-loc NCBI37.3.gene.loc  --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.locations.magma

#Gene based analysis
./magma --bfile /home/genetics/reference_panels/g1000_eur \
--pval metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.pvalues ncol=2 \
--gene-annot ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.locations.magma.genes.annot \
--out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based
