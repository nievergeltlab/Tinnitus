#Map mouse genes to human genes using biomart
 library(data.table)
 
#Gene expression data from Liu_Lu_2018_supplement.xlsx
 d1 <- fread('mouse_expression_lui.csv',data.table=F) 
 
#Create a column for merging
 d1$mousegene <- d1[,"Gene name"]
 
#For log2 means, take log2 + 1 for each, then take mean.
 l2mean <- function(x,...)
 {
  mean(log2(x+1),...)
 }


#Get sum across all categories
 pillars_set <- grep("Pillars",names(d1))
  d1$pillars_mean <- apply(d1[,pillars_set],1,l2mean,na.rm=T)
 deiters_set <- grep("Deiters",names(d1))
  d1$deiters_mean <- apply(d1[,deiters_set],1,l2mean,na.rm=T)
 ohc_set <- grep("OHC",names(d1))
  d1$ohc_mean <- apply(d1[,ohc_set],1,l2mean,na.rm=T)
 ihc_set <- grep("IHC",names(d1))
  d1$ihc_mean <- apply(d1[,ihc_set],1,l2mean,na.rm=T)
 melanocytes_set <- grep("Melanocytes",names(d1))
  d1$melanocytes_mean <- sapply(d1[,melanocytes_set],l2mean,na.rm=T)  #Only 1 melanocyte read
 liver_set <- grep("Liver",names(d1))
  d1$liver_mean <- apply(d1[,liver_set],1,l2mean,na.rm=T)
 rbc_set <- grep("RBC",names(d1))
  d1$rbc_mean <- apply(d1[,rbc_set],1,l2mean,na.rm=T)
  
#Get mean across each tissue type
 d1$genome_avg <- apply(d1[,grep("_mean",names(d1))],1,mean,na.rm=T)
 
 d1ear <- d1[,-which(names(d1) %in% c("rbc_mean","liver_mean"))]
 d1$genome_avg_ear <- apply(d1ear[,grep("_mean",names(d1ear))],1,mean,na.rm=T)
  
 
 d1c <- d1
 d1c$mousegeneuc <- d1c[,"Feature ID"]
 
 
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
 
 d3 <- subset(d3,select=c(pillars_mean,deiters_mean,ohc_mean,ihc_mean,melanocytes_mean)) # ,liver_mean,rbc_mean,genome_avg,genome_avg_ear
 
##MAGMA export

#Add human names back in to starting position
 dmagma <- as.data.frame(cbind(row.names(d3),d3),check.names=TRUE)
 
 names(dmagma)[1] <- "Entrez"
 names(dmagma) <- make.names(names(dmagma))
 
 
#Export data for MAGMA
 write.table(dmagma,"mouseexpressionlui_human_mapped.txt",quote=F,row.names=F)

  #Condition on mean of other cell types
/mnt/ukbb/software/magma --gene-results ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based.genes.raw \
--gene-covar mouseexpressionlui_human_mapped.txt  \
--model condition-hide=genome_avg_ear  \
--out michalski_genome_avg_ear


  


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


 d3_tv3  <- round(d3_tv2,4)


#Take top 10% of expressed genes for each
for (input_col in 1:colcount)
{
 qt <- quantile(d3_tv3[,input_col],0.9)
 pickgenes <- which(d3_tv3[,input_col] >= qt)
 write.table(row.names(d3_tv3[pickgenes,]),file=paste("lui_",make.names(names(d3_tv)[input_col]),"_genes_earonly.txt",sep=''),row.names=F,col.names=F,quote=F)
}

#For genome avg, do all 
write.table(row.names(d3_tv3),file=paste("lui_genome_avg_genes.txt",sep=''),row.names=F,col.names=F,quote=F)


##Make LDSC files https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial

#Entrez coordinate file (only need to make once
 #cat ../glist-hg19_entrezid.txt | tail -n +2 |  sort -n -k1 -k2 | grep -v NA | grep -v CHR  | awk '{print $4,"chr"$1,$2,$3}' | cat geneheader.txt - > entrezid.coords

#Make LD files
conda activate  ldsc
for celltype in $(ls   | grep lui | grep earonly )
do
 echo $celltype
 for chr in   {1..22}
 do
  echo $chr
  outname=$(echo $celltype | sed 's/_genes.txt//g' | sed 's/lui_/lui./g')
  /mnt/ukbb/software/ldsc/make_annot.py --gene-set-file "$celltype" --gene-coord-file entrezid.coords --windowsize 100 --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC."$chr".bim --annot-file  lui_ldsc/"$outname"."$chr".gz
  /mnt/ukbb/software/ldsc/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC."$chr" --ld-wind-kb 100 --annot lui_ldsc/"$outname"."$chr".gz --thin-annot --out lui_ldsc/"$outname"."$chr" 
 done
done



#Get positions of all genes

##LDSC analysis

ls | grep lui | grep genes | grep -v avg | sed 's/_genes.txt//g' | sed s'/lui_/lui./g' | awk '{print $1, "lui_ldsc/"$1".,lui_ldsc/lui.genome_avg."}' > lui.ldcts

ls | grep lui | grep earonly  | sed 's/_genes.txt//g' | sed s'/lui_/lui./g' | awk '{print $1, "lui_ldsc/"$1".,lui_ldsc/lui.genome_avg."}' > lui_earonly.ldcts


 /mnt/ukbb/software/ldsc/ldsc.py \
    --h2-cts /mnt/ukbb/adam/tinnitus_gwas/paper2/metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz  \
    --ref-ld-chr /mnt/ukbb/software/ldsc/baseline_v1.2/baseline. \
    --out ldsc_cts/$celltype \
    --ref-ld-chr-cts lui.ldcts \
    --w-ld-chr /mnt/ukbb/software/ldsc/weights_hm3_no_hla/weights.
	
   /mnt/ukbb/software/ldsc/ldsc.py \
    --h2-cts /mnt/ukbb/adam/tinnitus_gwas/paper2/metal_results/ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.munge.gz.sumstats.gz  \
    --ref-ld-chr /mnt/ukbb/software/ldsc/baseline_v1.2/baseline. \
    --out ldsc_cts/$celltype \
    --ref-ld-chr-cts lui_earonly.ldcts \
    --w-ld-chr /mnt/ukbb/software/ldsc/weights_hm3_no_hla/weights.
	
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 ##MAGMA analysis

#Setp one: do not condition on covariates
/mnt/ukbb/software/magma --gene-results ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based.genes.raw \
--gene-covar mouseexpressionlui_human_mapped.txt  \
--out lui_unconditioned

cat ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based.genes.raw.gsa1.gsa.out

#Condition on mean of other cell types
/mnt/ukbb/software/magma --gene-results ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.gene_based.genes.raw \
--gene-covar mouseexpressionlui_human_mapped.txt  \
--model condition-hide=genome_avg  \
--out lui_genome_avg
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#Identify unexpressed genes (who to cite?)
 unexpressed <- which(d1$pillars_mean == 0 & d1$deiters_mean == 0 & d1$ohc_mean == 0 & d1$ihc_mean == 0 & d1$melanocytes_mean ==0)
 head(d1[unexpressed,])
 
 
 d2 <- merge(d1[-unexpressed,],genesV2,by="mousegene") #notice that we REMOVE the unexpressed genes
 d3 <- d2[,c("humangene","ihc_mean","ohc_mean","deiters_mean","pillars_mean","melanocytes_mean")]
 
 #For duplicates, I exclude both
 library(plyr)
 dupgenes <- table(d3$humangene)>=2
 duplicated_genes <- which(dupgenes == TRUE)
 duplicated_genes_ids <- names(dupgenes)[duplicated_genes]
 
 d3 <- d3[-which(d3$humangene %in% duplicated_genes_ids),]
 
 d3$pillars_specificity = d3$pillars_mean / (d3$ohc_mean + d3$ihc_mean + d3$pillars_mean + d3$deiters_mean+ d3$melanocytes_mean)
 d3$deiters_specificity = d3$deiters_mean / (d3$ohc_mean + d3$ihc_mean + d3$pillars_mean + d3$deiters_mean+ d3$melanocytes_mean)
 d3$ohcs_specificity = d3$ohc_mean / (d3$ohc_mean + d3$ihc_mean + d3$pillars_mean + d3$deiters_mean+ d3$melanocytes_mean)
 d3$ihcs_specificity = d3$ihc_mean / (d3$ohc_mean + d3$ihc_mean + d3$pillars_mean + d3$deiters_mean+ d3$melanocytes_mean)
 d3$melanocytes_specificity = d3$ihc_mean / (d3$ohc_mean + d3$ihc_mean + d3$pillars_mean + d3$deiters_mean + d3$melanocytes_mean)
  
 d3$pillars_specificity_top10 <- as.numeric(d3$pillars_specificity >= quantile(d3$pillars_specificity,0.9,na.rm=T))
 d3$deiters_specificity_top10 <- as.numeric(d3$deiters_specificity >= quantile(d3$deiters_specificity,0.9,na.rm=T))
 d3$ohcs_specificity_top10 <- as.numeric(d3$ohcs_specificity >= quantile(d3$ohcs_specificity,0.9,na.rm=T))
 d3$ihcs_specificity_top10 <- as.numeric(d3$ihcs_specificity >= quantile(d3$ihcs_specificity,0.9,na.rm=T))
 d3$melanocytes_specificity_top10 <- as.numeric(d3$melanocytes_specificity >= quantile(d3$melanocytes_specificity,0.9,na.rm=T))
 
 pillars_sub <- t(subset(d3,pillars_specificity_top10==TRUE, select=c(humangene)))
 pillars_sub <- c("pillars",pillars_sub)
 write.table(t(pillars_sub),"pillars_top10_julien.txt",row.names=F,quote=F,col.names=F)
 
 
 #Cederroth method:

We processed and calculated the expression specificity as previously described.36 
Briefly, in each tissue expression dataset (i.e., organ of Corti,38 stria vascularis/SGNs,37 
and neural tissue39), we first aggregated the count per gene per cell type and excluded
 genes that are (1) not expressed in any cell type, (2) with duplicated identifier, or 
 (3) not 1:1 orthologous between mouse and human. We then normalized the expression to 1
 TPM (transcripts per million) per cell type. Next, gene expression specificity was 
 calculated per gene per cell type as:
