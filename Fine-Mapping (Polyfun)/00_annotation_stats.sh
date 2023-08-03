library(data.table)

d1 <- fread('finemapping_annotated.csv',data.table=F)

library(plyr)

#Exclude loci that do not finemap

d1 <- subset(d1$finemapped==1)

d1$rdb_lt5 <- d1$RDB
d1[which(d1$RDB %in% c("1a",   "1b" ,  "1c" ,  "1d"  , "1f" , "2a" ,  "2b" ,  "2c" ,  "3a" ,  "3b" ,   "4")),]$rdb_lt5 <- 1
d1[which(d1$RDB %in% c("5","6","7")),]$rdb_lt5 <- 0
d1$rdb_lt5 <- as.numeric(d1$rdb_lt5)


#
t.test(d1$CADD ~ d1$credible_set)

table(d1$credible_set,d1$func)
table(d1$credible_set,d1$rdb_lt5)

#R table behavior is completely fucked up now, it has a null column..

 chisq.test(table(d1$credible_set,d1$rdb_lt5)[1:2,2:3])
 
 
#Compared to SNPs outside of the credible set, SNPs within credible sets were more likely to be 
n_outside <- table(d1$credible_set)[1]
n_credible <-  table(d1$credible_set)[2]
2058 SNPs in creidlge set
11281 outside
1125

 chisq.test(table(d1$credible_set,d1$func == "exonic" ))
 chisq.test(table(d1$credible_set,d1$func == "intergenic" ))
 chisq.test(table(d1$credible_set,d1$func == "UTR3" ))
 chisq.test(table(d1$credible_set,d1$func == "upstream" ))
 chisq.test(table(d1$credible_set,d1$func == "downstream" ))
 chisq.test(table(d1$credible_set,d1$func == "ncRNA_exonic" ))
 chisq.test(table(d1$credible_set,d1$func == "ncRNA_intronic" ))
  chisq.test(table(d1$credible_set,d1$func == "intronic" ))
    
 
#Get N variants in credible set of each locus
ddply(d1, ~ GenomicLocus , colwise(.fun=sum, .cols= "credible_set"))

table( ddply(d1, ~ GenomicLocus , colwise(.fun=sum, .cols= "credible_set"))[,2])

#locus 19 - 1 variant
#locus 45, 49, 53 -2 variants
53 - 2 variants

#Get N variants total (counts total variants, but does not count ones tested
ddply(d1, ~ GenomicLocus , colwise(.fun=length, .cols= "credible_set"))


 UTR3 UTR5
 
 
 p_exonic_credible <- 39/n_credible 
 p_exonic_outside  <- 48/n_outside
 

 #Take median cadd score versus take median within each locus then compute..
 
 
 #eQTL
 chisq.test(table(d1$eqtlMapFilt,d1$credible_set))
 chisq.test(table(d1$ciMapFilt,d1$credible_set))
 #39 exonic
 
 t.test(d1$dist ~ d1$cr
 
 #you have to correct for gene length when you do this...
 
#is adjusting for the locus the same as adjusting for the length? I think it may be better, the length adjustment is because they CANT adjust for things as factors
 summary(lm(CADD ~ as.factor(GenomicLocus) + credible_set,data=d1))
 summary(lm(commonChrState ~ as.factor(GenomicLocus) + credible_set,data=d1))
 summary(lm(minChrState ~ as.factor(GenomicLocus) + credible_set,data=d1))
 
 by(d1$commonChrState,d1$credible_set,median,na.rm=T)
  by(d1$minChrState,d1$credible_set,median,na.rm=T)
  
table(d1$RDB,d1$credible_set)

#perhaps load lengths as a file
  summary(glm(rdb_lt5 ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="exonic") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="UTR3") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="UTR5") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="intron") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))

  summary(glm(I(func=="ncRNA_exonic") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="ncRNA_exonic:splicing") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="ncRNA_intronic") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="splicing") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="upstream") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="downstream") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
  summary(glm(I(func=="upstream:downstream") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
    
 
 
  summary(glm(I(func=="intergenic") ~ as.factor(GenomicLocus) + credible_set,data=d1,family='binomial'))
 