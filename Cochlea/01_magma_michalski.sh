library(data.table)
 
#Gene expression data from Michaelski, as processed by Yanning.
 d1 <- fread('P20_pseudobulk-mean-rawcounts.csv',data.table=F) 
  nameset <- d1[,1]
  d1 <- as.data.frame(t(d1[,-1]))
  names(d1) <- nameset

  log2f <- function(X)
  {
   log2(as.numeric(X) + 1)
  }
  d1 <- as.data.frame(apply(d1,c(1,2),log2f))
  
#Get mean across each tissue type
 d1$genome_avg <- apply(d1,1,mean,na.rm=T)
 
 d1c <- d1
 d1c$mousegeneuc <- row.names(d1)
 
 
#Load mouse to human annotations, in entrez form
 library(dplyr)
 mouse_human_genes = fread("HOM_MouseHumanSequence.rpt",data.table=F)
 names(mouse_human_genes)[1] <- "DB.Class.Key"
#Convert to wide format
 mousers <- reshape(mouse_human_genes, idvar="DB.Class.Key", timevar="Common Organism Name",direction='wide')
 mousers2 <- subset(mousers, select=c("EntrezGene ID.human","Symbol.mouse, laboratory"))
 names(mousers2) <- c("humangene","mousegeneuc")

 d2 <-  merge(d1c,mousers2,by="mousegeneuc")
 
#Remove duplicates (strict, both removed)
 d3 <-  d2[!(duplicated(d2$humangene) | duplicated(d2$humangene, fromLast=TRUE)), ]

 row.names(d3) <- d3$humangene

 d3 <- d3[,!(names(d3) %in% c("humangene","mousegeneuc"))]
 
##MAGMA export

#Add human names back in to starting position
 dmagma <- as.data.frame(cbind(row.names(d3),d3),check.names=TRUE)
 
 names(dmagma)[1] <- "Entrez"
 names(dmagma) <- make.names(names(dmagma))
 
 
#Export data for MAGMA
 write.table(dmagma,"mouseexpressionmichalski_human_mapped.txt",quote=F,row.names=F)

###MAGMA analysis

#Setp one: do not condition on covariates
/mnt/ukbb/software/magma --gene-results ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based.genes.raw \
--gene-covar mouseexpressionmichalski_human_mapped.txt  \
--out michalski_unconditioned

cat ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based.genes.raw.gsa1.gsa.out

#Condition on mean of other cell types
/mnt/ukbb/software/magma --gene-results ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based.genes.raw \
--gene-covar mouseexpressionmichalski_human_mapped.txt  \
--model condition-hide=genome_avg  \
--out michalski_genome_avg


##LDSC export

#Calculate t statistics for each gene and each tissue
colcount <- ncol(d3)

tvalcalc <- function(dframe, colcount, input_col,Xval)
{
 X <- Xval
 Y <- as.numeric(dframe)
 tval <- summary(lm(Y ~ X))$coefficients[2,3]
 return(tval)
}

d3_tv <- as.data.frame(matrix(nrow=nrow(d3),ncol=ncol(d3)))
colnames(d3_tv) <- names(d3)
row.names(d3_tv) <- row.names(d3)

for (input_col in 1:colcount)
{
 print(input_col)
 Xval <- sign(1:colcount == input_col)
 d3_tv[,input_col] <- apply(d3,1,tvalcalc,colcount=colcount,input_col=input_col,Xval=Xval)
}

#For LDSC, need to:
#1 Get rid of genes that are not expressed
#2 Replace NA values
#NAs is tautological with non-expression. All of the rows with NA are NA everywhere

#Remove NA rows
 row.has.na <- apply(d3_tv, 1, function(x){any(is.na(x))})
 d3_tv2 <- d3_tv[!row.has.na,]

 d3_tv2$genome_avg <- apply (d3_tv2[,-(names(d3_tv3) %in% "genome_avg")],1,median) #background is the median test statistic

 d3_tv3  <- round(d3_tv2,4)


#Take top 10% of expressed genes for each
for (input_col in 1:colcount)
{
 qt <- quantile(d3_tv3[,input_col],0.9)
 pickgenes <- which(d3_tv3[,input_col] >= qt)
 write.table(row.names(d3_tv3[pickgenes,]),file=paste("michalski_",make.names(names(d3_tv)[input_col]),"_genes.txt",sep=''),row.names=F,col.names=F,quote=F)
}

#For genome avg, do all 
write.table(row.names(d3_tv3),file=paste("michalski_genome_avg_genes.txt",sep=''),row.names=F,col.names=F,quote=F)


##Make LDSC files https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial

#Entrez coordinate file (only need to make once
 #cat ../glist-hg19_entrezid.txt | tail -n +2 |  sort -n -k1 -k2 | grep -v NA | grep -v CHR  | awk '{print $4,"chr"$1,$2,$3}' | cat geneheader.txt - > entrezid.coords

#Make LD files
conda activate  ldsc
for celltype in $( ls | grep _genes.txt | tail -n+33 | head -n 12 | grep -v Type.II ) #   
do
 echo $celltype
 for chr in   {1..22}
 do
  echo $chr
  outname=$(echo $celltype | sed 's/_genes.txt//g' | sed 's/michaelski_/michaelski./g')
  /mnt/ukbb/software/ldsc/make_annot.py --gene-set-file "$celltype" --gene-coord-file entrezid.coords --windowsize 100 --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC."$chr".bim --annot-file  michaelski_ldsc/"$outname"."$chr".gz
  /mnt/ukbb/software/ldsc/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC."$chr" --ld-wind-kb 100 --annot michaelski_ldsc/"$outname"."$chr".gz --thin-annot --out michaelski_ldsc/"$outname"."$chr" 
 done
done



#Get positions of all genes

ls | grep michal | grep genes | grep -v avg | sed 's/_genes.txt//g'  | awk '{print $1, "michaelski_ldsc/"$1".,michaelski_ldsc/michalski_genome_avg."}' > michalski.ldcts


split -l 12  michalski.ldcts --additional-suffix _ldctsmichal

xaa_ldctsmichal
xab_ldctsmichal
xac_ldctsmichal
##LDSC analysis

 /mnt/ukbb/software/ldsc/ldsc.py \
    --h2-cts /mnt/ukbb/adam/tinnitus_gwas/paper2/metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz  \
    --ref-ld-chr /mnt/ukbb/software/ldsc/baseline_v1.2/baseline. \
    --out ldsc_cts/first10_$celltype \
    --ref-ld-chr-cts xaa_ldctsmichal \
    --w-ld-chr /mnt/ukbb/software/ldsc/weights_hm3_no_hla/weights.
	
 /mnt/ukbb/software/ldsc/ldsc.py \
    --h2-cts /mnt/ukbb/adam/tinnitus_gwas/paper2/metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz  \
    --ref-ld-chr /mnt/ukbb/software/ldsc/baseline_v1.2/baseline. \
    --out ldsc_cts/second10_$celltype \
    --ref-ld-chr-cts xab_ldctsmichal \
    --w-ld-chr /mnt/ukbb/software/ldsc/weights_hm3_no_hla/weights.
	

 /mnt/ukbb/software/ldsc/ldsc.py \
    --h2-cts /mnt/ukbb/adam/tinnitus_gwas/paper2/metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz  \
    --ref-ld-chr /mnt/ukbb/software/ldsc/baseline_v1.2/baseline. \
    --out ldsc_cts/third10_$celltype \
    --ref-ld-chr-cts xac_ldctsmichal \
    --w-ld-chr /mnt/ukbb/software/ldsc/weights_hm3_no_hla/weights.
	



#Get top 10% genes for each category, then add 100kb buffers



