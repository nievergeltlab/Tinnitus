#Tinnitus studymax phenotype construction
#last rev april 2, 2019
#Phenotype construction notes:

#From ukb24405.tab:
#4803 is self reported tinnitus
#21003-0.0 age when at center
#31-0.0 is sex (self report)
#22009.0.x are PCs
#22027 are samples who failed QC
#22000 - genotyping center

head -n1 /run/media/royce/ukb/pheno/ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.4803.0.0")
      printf("Column %d is f.4803.0.0\n", i-1)

  }
  exit 0
}
' 

head -n1 /run/media/royce/ukb/pheno/ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.21003.0.0")
      printf("Column %d is f.21003.0.0\n", i-1)

  }
  exit 0
}
' 


head -n1 ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.22009.0.1")
      printf("Column %d is f.22009.0.1\n", i-1)

  }
  exit 0
}
' 

head -n1 ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.22027.0.0")
      printf("Column %d is f.22027.0.1\n", i-1)

  }
  exit 0
}
' 

head -n1 ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.22027.0.0")
      printf("Column %d is f.22000.0.1\n", i-1)

  }
  exit 0
}
' 


head -n1 ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.22000.0.0")
      printf("Column %d is f.22000.0.1\n", i-1)

  }
  exit 0
}
' 

head -n1 ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.54.0.0")
      printf("Column %d is f.54.0.1\n", i-1)

  }
  exit 0
}
' 


head -n1 ukb24405.tab | awk  '
{
  for(i=1;i<=NF;i++) {
    if($i == "f.22001.0.0")
      printf("Column %d is f.22001.0.1\n", i-1)

  }
  exit 0
}
' 

#f.4803.0.0      f.4803.1.0      f.4803.2.0 : tinnitus in one or both ears
#f.21003.0.0     f.21003.1.0     f.21003.2.0 age at assessment center
#31 sex (column 2)
#2247 is hearing problems. Note there are deaf people! the notes say 4803 was not collected in teh deaf thou
#22009.0.1 - are PCs
#22027 - failed sample. value 1 means exclude!
#22000.0.0 - gentoyping batch
#54 - assessment center
#22001.0.0 - genetic sex


#Take the important columns
awk '{print $1,$2,$15,$16,$17, $945,$946,$947,$1070,$1071,$1072,$1082,$1083,$1084,$1085,$1086,$1087,$1088,$1089,$1090,$1091,$1130,$1073,$1074,$8}'  ukb24405.tab > UKBB_tinnitus_alldata_apr2_2019.txt

#Construct phenotype in R
R
library(data.table)
d1 <- fread('UKBB_tinnitus_alldata_apr2_2019.txt',data.table=F)
dim(d1)
d1$FID <- d1$f.eid

genodat <- fread('/mnt/ukbb/adam/preimp/ukb_chr12_v2.fam',data.table=F)
genodat$FID <- genodat$V1
#Immediately remove samples NOT in genotype data.. 
d1 <- subset(d1,FID %in% genodat$FID)
dim(d1)

d1 <- subset(d1,is.na(f.22027.0.0))
dim(d1)
#Sex mismatch check
table(d1$f.22001.0.0,d1$f.31.0.0)
d1 <- subset(d1,f.22001.0.0 == f.31.0.0)

table(d1$f.4803.0.0 )
table(d1$f.4803.1.0 )
table(d1$f.4803.2.0 )
#11 = yes now or most all of the time
#12 = yes now a lot of the time
#13 = some of the time
#14 = yes but not now, have in the past
#-3 is prefer to not answer
#-1 dont know
#0 no, never

#Recode this to studymax
 d1$f.4803.0.0_recode <- NA
 d1[which(d1$f.4803.0.0 == 11),]$f.4803.0.0_recode <- 4
 d1[which(d1$f.4803.0.0 == 12),]$f.4803.0.0_recode <- 3
 d1[which(d1$f.4803.0.0 == 13),]$f.4803.0.0_recode <- 2
 d1[which(d1$f.4803.0.0 == 14),]$f.4803.0.0_recode <- 1
 d1[which(d1$f.4803.0.0 == 0),]$f.4803.0.0_recode <- 0
 d1[which(d1$f.4803.0.0 == -1),]$f.4803.0.0_recode <- -1
 d1[which(d1$f.4803.0.0 == -3),]$f.4803.0.0_recode <- -3

 d1$f.4803.1.0_recode <- NA
 d1[which(d1$f.4803.1.0 == 11),]$f.4803.1.0_recode <- 4
 d1[which(d1$f.4803.1.0 == 12),]$f.4803.1.0_recode <- 3
 d1[which(d1$f.4803.1.0 == 13),]$f.4803.1.0_recode <- 2
 d1[which(d1$f.4803.1.0 == 14),]$f.4803.1.0_recode <- 1
 d1[which(d1$f.4803.1.0 == 0),]$f.4803.1.0_recode <- 0
 d1[which(d1$f.4803.1.0 == -1),]$f.4803.1.0_recode <- -1
 d1[which(d1$f.4803.1.0 == -3),]$f.4803.1.0_recode <- -3

 d1$f.4803.2.0_recode <- NA
 d1[which(d1$f.4803.2.0 == 11),]$f.4803.2.0_recode <- 4
 d1[which(d1$f.4803.2.0 == 12),]$f.4803.2.0_recode <- 3
 d1[which(d1$f.4803.2.0 == 13),]$f.4803.2.0_recode <- 2
 d1[which(d1$f.4803.2.0 == 14),]$f.4803.2.0_recode <- 1
 d1[which(d1$f.4803.2.0 == 0),]$f.4803.2.0_recode <- 0
 d1[which(d1$f.4803.2.0 == -1),]$f.4803.2.0_recode <- -1
 d1[which(d1$f.4803.2.0 == -3),]$f.4803.2.0_recode <- -3

 d1$f.4803.whichmax <- apply(d1[,c("f.4803.0.0","f.4803.1.0","f.4803.2.0")],1,which.max)
 
 takecol <- function(x)
 {
  column_to_pick <- x[4]
  return(x[column_to_pick])
 }
 
 d1$f.4803.max <- apply(d1[,c("f.4803.0.0_recode","f.4803.1.0_recode","f.4803.2.0_recode","f.4803.whichmax")],1,takecol)
 
 #Assign the age at the maximum tinnitus severity
 d1$f.21003.max <- apply(d1[,c("f.21003.0.0","f.21003.1.0","f.21003.2.0","f.4803.whichmax")],1,takecol)
 #And hearing loss status at this timepoint
 d1$f.2247.max <- apply(d1[,c("f.2247.0.0","f.2247.1.0","f.2247.2.0","f.4803.whichmax")],1,takecol)

 table(d1$f.2247.0.0)
 #1 = yes
 #0 = no
 #-1 = dont know 
 #-3 = dont want to answer
 #Dont need to recode, fine as is

table(d1$f.2247.max,d1$f.4803.max)

table(d1$f.2247.max,d1$f.4803.max,useNA="always")

#Ok, it's looking great.
#I'll do the highest level vs. lowest
#And any v no tinnitus codings
 d1$f.4803.max_coding1 <- NA
 d1[which(d1$f.4803.max == 4),]$f.4803.max_coding1 <- 1
 d1[which(d1$f.4803.max == 0),]$f.4803.max_coding1 <- 0

 d1$f.4803.max_coding2 <- NA
 d1[which(d1$f.4803.max %in% c(1,2,3,4)),]$f.4803.max_coding2 <- 1
 d1[which(d1$f.4803.max == 0),]$f.4803.max_coding2 <- 0

 #third coding just getting rid of the -1 and -3 people
 d1$f.4803.max_coding3 <- d1$f.4803.max
 d1[which(d1$f.4803.max %in% c(-1,-3)),]$f.4803.max_coding3 <- NA
 
 #Fourth coding where its additive but you combine the top 2 categories
 d1$f.4803.max_coding4 <- d1$f.4803.max
 d1[which(d1$f.4803.max == 4),]$f.4803.max_coding4 <- 3
 d1[which(d1$f.4803.max %in% c(-1,-3)),]$f.4803.max_coding4 <- NA
 
 table(d1$f.4803.max_coding3,d1$f.4803.max_coding4 )
 #Hearing loss coding
 d1$f.2247.max_coding1 <- NA
 d1[which(d1$f.2247.max == 1),]$f.2247.max_coding1 <- 1
 d1[which(d1$f.2247.max == 0),]$f.2247.max_coding1 <- 0

 table(d1$f.4803.max_coding1,d1$f.2247.max_coding1 )
 table(d1$f.4803.max_coding2,d1$f.2247.max_coding1 )

 #Missing age data??? 
 table(is.na(d1$f.21003.max),is.na(d1$f.4803.max_coding1)) #68 people, i dont give a shit.
 
 #load in ancestry information
 ancestries <- fread('ukbb_ancestry_allchr.predpc_oneweek.header_adam',data.table=F)

 d1a <- merge(ancestries,d1,by="FID") #Some subjects are getting lost due to not having genetic data??

 d1_eur <- subset(d1a,bestpop_oneweek_reclass2=="eur")

 table(d1_eur$f.4803.max_coding1) #So we have at least 10,000 cases

 #Load in relatedness (should come after ancestry, as there may be cross-ancestry relateds)
 #Note that only related pairs are in here.
 #According to the manual:
 #Close relatives can be inferred fairly reliably based on the estimated kinship coefficients as shown in the following simple algorithm: 
 #an estimated kinship coefficient range >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] 
 #corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively.
 
 relatedness0 <-  fread('ukb40951_rel_chr1_s488292.dat',data.table=F)
 
 #Filter to only european subjects (both subjects must be european for me to care about this pruning, otherwise they wont be analyzed!!!
 #And on degree of relationship
 relatedness <- subset(relatedness0, ID1 %in% d1_eur$FID & ID2 %in% d1_eur$FID & Kinship >= 0.0442)
 
 #Filter down relatedness data to only care about > N degree related
 #Discussion with caroline says we filter out third or more degree related.
 
  
 #Question: I have two tinnitus phenotypes. It seems like I should join on the less restrictive tinnitus phenotype to find more cases
 #Then use tinnitus status 1 as a tiebreaker. However status 1 does not even consider intermediates..
 
 #Ill just go by max symptoms
 
 
 #Note: this has to be iterative, I have to merge the relatedness again with the pruned list until all relateds are gone!!
 #Hence I make a df called d1_eur2, which keeps getting updated in the loop 
 
 d1_eur2 <- d1_eur
 
#I use a while loop to see if there is any remaining overlap between my dataset and the relatedness dataset
#By seeing if the dimension of the dataset decreases
#if the data doesnt prune any more, stop

 d1preprunedim <- dim(d1_eur2)[1] #Dimension of data before pruning
 
 d1_eur_s1 <- subset(d1_eur2,select=c(FID,f.4803.max))
 d1_eur_s2 <- d1_eur_s1
 
 names(d1_eur_s1)[1] <- "ID1"
 names(d1_eur_s2)[1] <- "ID2"
 
 #
 r1 <- merge(relatedness,d1_eur_s1,by="ID1",suffixes=c("","_id1"))
 r2 <- merge(r1,d1_eur_s2,by="ID2",suffixes=c("","_id2")) 
 
 #Note who will be pruned
 r2$prune <- NA
  for (i in 1:dim(r2)[1])
 {
  #Look at each pair.
  #Note: if there is a tie, will just take first one
  #Note I am taking the minimum as this is a list of who to toss
  #If one of these is NA, just toss the NA
  #if both are NA, we dont have to care
  if(is.na(r2[i,]$f.4803.max))
  {
   r2[i,]$prune <- r2[i,]$ID1
  } else if (is.na(r2[i,]$f.4803.max_id2))
  {
   r2[i,]$prune <- r2[i,]$ID2
  } else {
   minsub <- which.min(c(r2[i,]$f.4803.max,r2[i,]$f.4803.max_id2))
   r2[i,]$prune <- r2[i,minsub] #Since the first column is hte first subject, second column is second subject, this is adequate for pruning
   rm(minsub)} #Since this is a loop I dont want weird shit happening because this doesnt get declared
  }
  
  d1_eur2 <- subset(d1_eur2, !(FID %in% r2$prune))
  
  dim(d1_eur)[1]  #Dimension of data before pruning
  dim(d1_eur2)[1] #Dimension of data after pruning

 #For each related pair, prefer to retain a tinnitus case.
 #~10k cases after adjustment for 

 #Table of assesement center, batch
 table(d1_eur$f.22000.0.0,d1_eur$f.54.0.0)
 
 #Table 
 
#Write phenotype files (unrelated and related included)
write.table(subset(d1_eur2,select=c(FID,IID,f.4803.max_coding1,f.4803.max_coding2,f.4803.max_coding3,f.4803.max_coding4,f.2247.max_coding1,f.21003.max,f.31.0.0,f.22009.0.1,f.22009.0.2,f.22009.0.3,f.22009.0.4, f.22009.0.5, f.22009.0.6,f.22000.0.0,f.54.0.0)),"UKB_tinnitus_eur_unrelated_may2_2019.pheno",row.names=F,quote=F)
write.table(subset(d1_eur,select=c(FID,IID,f.4803.max_coding1,f.4803.max_coding2,f.4803.max_coding3,f.4803.max_coding4,f.2247.max_coding1,f.21003.max,f.31.0.0,f.22009.0.1,f.22009.0.2,f.22009.0.3,f.22009.0.4, f.22009.0.5,f.22009.0.6,f.22000.0.0,f.54.0.0)),"UKB_tinnitus_eur_related_may2_2019.pheno",row.names=F,quote=F)

# DIagnostic code:: Compare versions of phenotype file..
library(data.table)
d1 <- fread('UKB_tinnitus_eur_unrelated_april2_2019.pheno')
d2 <- fread('UKB_tinnitus_eur_unrelated_feb27_2019.pheno')
dm <- merge(d1,d2,by="FID")
dim(d1)
dim(d2) #should be slightly different due to no filter on quality
table(dm$f.2247.max_coding1.y , dm$f.2247.max_coding1.x) #this info should be hte same
table(dm$f.4803.max_coding3.y , dm$f.4803.max_coding3.y) 
#OK looks fine..

summary(lm(f.4803.max_coding3 ~ as.factor(f.22000.0.0) + as.factor(f.54.0.0),data=d1))

 anova(lm(f.4803.max_coding3 ~f.22009.0.1 + f.22009.0.2 + f.22009.0.3+ f.22009.0.4+  f.22009.0.5 + f.21003.max + as.factor(f.22000.0.0) + as.factor(f.54.0.0)  ,data=d1))
 