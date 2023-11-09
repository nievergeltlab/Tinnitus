# Running July 21, 2021 
# Author: Marianna Gasperi, Ph.D. 
# Purpose: Evaluating the genetic relationship between Tinnitus 
# (updated META 07/21) and HD META from MVP and
# UKBB and 11 psychiatric conditions from PGC

# This script: 1) munges the data 2) runs multivariate LD score regression for
# all, and odd, and even chromosomes 3) saves the output for LD score regression
# for all three ways I prefer to do this part on a lab machine due to file size
# and run time.
#
# EFA and CFA are conducted using another script that I typically run locally.
# The only output that is needed from this script are the three LDSC files and
# the log files saved out at the end. The munged files can be kept and used for
# other analyses. Note, GSEM does not like to overwrite log files, so if you
# have to rerun something, delete appropriate log before rerunning.

--------------------------------------
   #If needed
   # install.packages("devtools")
   # library(devtools)
   # install_github("MichelNivard/GenomicSEM")
--------------------------------------

#can ignore the 24 warnings




# Fix Summary Stat files -------------------------------------------------------

#This is example code from other changes made with column names/rs.

#rename columns in AN  - pgcAN2.2019-07.vcf.tsv.gz from REF to A2 and ALT to A1 $$$THESE ARE REVERSED BECASUE THEY ARE INCORRECT IN THE ORIGINAL FILE
# pgcAN2.2019.07.vcf.tsv <- fread("pgcAN2.2019-07.vcf.tsv.gz", data.table=F, quote="\"")
# names(pgcAN2.2019.07.vcf.tsv)[names(pgcAN2.2019.07.vcf.tsv) == "REF"] <- "A2"
# names(pgcAN2.2019.07.vcf.tsv)[names(pgcAN2.2019.07.vcf.tsv) == "ALT"] <- "A1"
# write.table(pgcAN2.2019.07.vcf.tsv, file = "pgc_sumstats/pgcAN2.2019-07.vcf.tsv",
            # sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE, append = FALSE)

# #rename effect column in ALC Dep - pgc_alcdep.eur_discovery.aug2018_release.txt from Z to Zscore
# pgc_alcdep.eur_discovery.aug2018_release.txt <- fread("zcat pgc_alcdep.eur_discovery.aug2018_release.txt.gz", data.table=F)
# names(pgc_alcdep.eur_discovery.aug2018_release.txt)[names(pgc_alcdep.eur_discovery.aug2018_release.txt) == "Z"] <- "Zscore"
# pgc_alcdep.eur_discovery.aug2018_release.txt$SNP<-sapply(strsplit(pgc_alcdep.eur_discovery.aug2018_release.txt$SNP, ":"), `[`, 1)
# write.table(pgc_alcdep.eur_discovery.aug2018_release.txt, file = "pgc_sumstats/pgc_alcdep.eur_discovery.aug2018_release.txt-wide.txt",
            # sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE, append = FALSE)

# MUNGE -------------------------------------------------------------------

# All  files being combined and munged - make sure that everything is in the right order!!!
# If you get "Error: cannot allocate vector....." the job is too big, easiest solution: split up into chunks (next section)

# MUNGE -------------------------------------------------------------------
# munging all files at once, might run out of memory, so can break this step up into two runs (see below)
# # # All  files being combined and munged - make sure that everything is in the right order!!!

# echo "SNP CHR BP A1 A2 BETA SE P NSTUD  HETP" | cat  - pgc_sumstats/GPC-2.NEUROTICISM.full.txt > pgc_sumstats/GPC-2.NEUROTICISM.full.txt2

# #I have to unfuck the stupid summary format of these UKBB data
# zcat /mnt/ukbb/adam/tinnitus_gwas/mendelian_randomization_inputs/variants.tsv.bgz | awk '{if (NR==1 ||( $13 >= 0.01 && $13 <= 0.99) )print $1,$4,$5,$6,$10,$13}' |LC_ALL=C sort -k1b,1 > variants.tsv.bgz.txt 

# 6159_3.gwas.imputed_v3.both_sexes.tsv.bgz 
# for traits in 6154_100.gwas.imputed_v3.both_sexes.tsv.bgz 135.gwas.imputed_v3.both_sexes.tsv.bgz 137.gwas.imputed_v3.both_sexes.tsv.bgz 20002_1265.gwas.imputed_v3.both_sexes.tsv.bgz 6154_3.gwas.imputed_v3.both_sexes.tsv.bgz 6145_100.gwas.imputed_v3.both_sexes.tsv.bgz 6159_100.gwas.imputed_v3.both_sexes.tsv.bgz 2080.gwas.imputed_v3.both_sexes.tsv.bgz 2188.gwas.imputed_v3.both_sexes.tsv.bgz
# do
 # LC_ALL=C join variants.tsv.bgz.txt  <(zcat $traits | sed 's/pval/P/g' | LC_ALL=C sort -k1b,1 ) | sort -g -k6 | sed 's/minor_AF/MAF/' | sed 's/ > "$traits".fixed
# done

#ANXIETY AND SUICIDALITY DROPPED, HEARING EXCLUDED (Per CAROLINE, dec 17,2021)
require(data.table)
library(GenomicSEM)

files=scan(what=character())
ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz
SWB_Full.txt.gz
GPC-2.NEUROTICISM.full.txt2
6159_3.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
6154_100.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
135.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
137.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
20002_1265.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
6154_3.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
6145_100.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
6159_100.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
2080.gwas.imputed_v3.both_sexes.tsv.bgz.fixed
2188.gwas.imputed_v3.both_sexes.tsv.bgz.fixed


trait.names <- scan(what=character())
TIN
SWB
NEUR
PAIN1
MED1
ILLN
TXN
MIG
MED2
STRESS
PAIN2
TIRED
LILL
SLEEP


Ns <- scan(what=numeric())
467369
298420
160671
360391
357084
361141
361141
361141
357084
358836
360391
350580
352798
359020

sample.prev <- scan(what=numeric())
NA
NA
NA
0.228
0.555
NA
NA
0.029
0.218
0.563
0.404
NA
0.325
NA

population.prev <- scan(what=numeric())
NA
NA
NA
0.05
0.55
NA
NA
0.03
0.22
0.06
0.05
NA
0.33
NA

# munge(files=files[4],
      # hm3 = "/mnt/ukbb/adam/tinnitus_gwas/w_hm3.noMHC.snplist",
      # trait.names=traits[4],
      # N=Ns[4],
      # info.filter = 0.9, maf.filter = 0.01)

# Run multivariable LDSC --------------------------------------------------

# sample.prev: A vector of sample prevalences of length equal to the number of
# traits. Sample prevalence is calculated as the number of cases over the total
# number of participants (i.e., cases + controls). Possible range = 0-1. If the
# trait is continuous, the values should equal NA. population.prev: A vector of
# population prevalences. These estimates can be obtained from a number of
# sources, such as large scale epidemiological studies. Possible range = 0-1.
# Again, if the trait is continuous the values should equal NA.

#HEARING EXCLUDED!
traits <- paste(trait.names,".sumstats.gz",sep="")


ld <- "/mnt/ukbb/royce/eur_w_ld_chr/"
wld <- "/mnt/ukbb/royce/eur_w_ld_chr/"


#Run the LDSC three ways. All for reporting sample genetic correlations and models not needing EFA; odd for EFA; even for CFA (if needed)
LDSCoutput.all <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, ldsc.log="LOG_ldsc_TINADJ_ldhub_all")
LDSCoutput.odd <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, select = "ODD", ldsc.log="LOG_ldsc_TINADJ_ldhub_odd")
LDSCoutput.evn <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, select = "EVEN", ldsc.log="LOG_ldsc_TINADJ_ldhub_even")

# LDSC Output -------------------------------------------------------------

# The output (named LDSCoutput here) is a list with X named variables in it:
# LDSCoutput$S is the covariance matrix (on the liability scale for case/control
# designs). LDSCoutput$V which is the sampling covariance matrix in the format
# expected by lavaan. LDSCoutput$I is the matrix of LDSC intercepts and
# cross-trait (i.e., bivariate) intercepts. LDSCoutput$N contains the sample
# sizes (N) for the heritabilities and sqrt(N1N2) for the co-heritabilities.
# These are the sample sizes provided in the munging process. LDSCoutput$m is
# the number of SNPs used to construct the LD score.

##optional command to save the ldsc output in case you want to use it in a later R session.
#These are the files used for EFA/CFA
save(LDSCoutput.all, file="ldsc_ldhub_jan3_2022_all.RData")
save(LDSCoutput.odd, file="ldsc_ldhub_jan3_2022_odd.RData")
save(LDSCoutput.evn, file="ldsc_ldhub_jan3_2022_even.RData")


load(file="ldsc_ldhub_jan3_2022_all.RData")
load(file="ldsc_ldhub_jan3_2022_odd.RData")
load(file="ldsc_ldhub_jan3_2022_even.RData")


#PLOTS
require(corrplot)

rownames(LDSCoutput.all$S_Stand) <- trait.names

pdf('tinnitus_ldhub_jan3_2022.pdf',12,12)
corrplot(LDSCoutput.all$S_Stand,  
         method="color",
         #title = "Genetic Correlations",
         order ="hclust",
         #addrect = 2, 
         tl.col = "black",
         tl.srt = 0,
         tl.offset = 1.,
         addgrid.col = "black",
         addCoef.col = "black")
         dev.off()

#############################################################
#USING ODD AND EVEN
#############################################################

#smooth the S matrix for EFA using the nearPD function in the Matrix package. 
require(Matrix)
Ssmooth.odd<-as.matrix((nearPD(LDSCoutput.odd$S, corr = FALSE))$mat)
Ssmooth.even<-as.matrix((nearPD(LDSCoutput.evn$S, corr = FALSE))$mat)
Ssmooth.all<-as.matrix((nearPD(LDSCoutput.all$S, corr = FALSE))$mat)


#############################################################
# Common Factor Model
#############################################################


# chisq: The model chi-square, reflecting index of exact fit to observed data, with lower values indicating better fit. 
# df and p_chisq: The degrees of freedom and p-value for the model chi-square. 
# AIC: Akaike Information Criterion. Can be used to compare models regardless of whether they are nested.
# CFI: Comparative Fit Index. Higher = better. > .90 = acceptable fit; > .95 = good model fit 
# SRMR: Standardized Room Mean Square Residual. Lower = better. < .10 = acceptable fit; < .05 = good fit	

# COMMON Factor Model using built in function -----------------------------

#Using DWLS estimation
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput.evn, estimation="DWLS")
CommonFactor_DWLS

# $modelfit  I LIKE TO PASTE FIT STATS HERE SO I CAN QUICKLY SCAN THESE AS I BUILD MODELS
# chisq df      p_chisq      AIC       CFI      SRMR
# df 423.4066 20 2.811742e-77 455.4066 0.7969974 0.1237625

# COMMON Factor model using usermodel Lavaan code-------------------------- 
#This should yield the exact results as above - run both and check that you 
#get the coding and output. 

#Define model using Lavaan 
F1.a.model<-'
F1 =~ TINADJ+SWB+NEUR+PAIN1+MED1+ILLN+TXN+MIG+MED2+STRESS+PAIN2+TIRED+LILL +SLEEP
'
F1.a.fit.evn<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F1.a.model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F1.a.fit.evn




# $modelfit
# chisq df      p_chisq      AIC       CFI      SRMR
# df 423.4066 20 2.811742e-77 455.4066 0.7969974 0.1237625

# Common factor models don't require EFA since there is just one factor. For other models, we do EFA, then CFA

#############################################################
# TWO Factors 
#############################################################

### EFA --------------------------------- 
EFA2.p.odd<-factanal(covmat = Ssmooth.odd,factors = 2, rotation = "promax")  #can use varimax
EFA2.p.odd

# Loadings:
# TINADJ  0.202   0.252
# SWB     0.371  -1.160
# NEUR            0.858
# PAIN1   0.719
# MED1   -1.075   0.214
# ILLN    0.747
# TXN     0.898
# MIG     0.500  -0.131
# MED2    1.001  -0.199
# STRESS -0.583
# PAIN2  -0.792
# TIRED   0.448   0.477
# LILL    0.646   0.141
# SLEEP  -0.212  -0.131

###CFA ---------------------------------   

F2.a.model.even<-'
F1 =~ SWB + PAIN1 + MED1 + ILLN+TXN+MIG+MED2+STRESS+PAIN2+TIRED+LILL
F2 =~ TINADJ + SWB + NEUR +TIRED
SWB ~~ a*SWB;
a > 0.001
'

F2.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F2.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F2.a.fit.even


F3.a.model.even<-'
F1 =~ PAIN1 + ILLN +TXN + STRESS+PAIN2+TIRED+LILL+SLEEP
F2 =~ SWB+MED1+TXN+MIG+MED2+PAIN2
F3 =~ TINADJ + SWB + NEUR
MED1 ~~ a*MED1
a > 0.001
'

F3.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F3.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F3.a.fit.even


F4.a.model.even<-'
F1 =~ TINADJ + SWB + NEUR 
F2 =~ PAIN1 + STRESS + PAIN2 + TIRED + LILL + SLEEP
F3 =~ ILLN + TXN + LILL
F4 =~ MIG + MED1 + MED2
MED1 ~~ a*MED1
a > 0.001
TXN ~~ a*TXN
a > 0.001


'

F4.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F4.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F4.a.fit.even


F5.a.model.even<-'
F1 =~ SWB + NEUR 
F2 =~ PAIN1 + STRESS + PAIN2 + TIRED + LILL + SLEEP
F3 =~ ILLN + TXN + LILL
F4 =~ MIG + MED1 + MED2
F5 =~ TINADJ 
MED1 ~~ a*MED1
a > 0.001
TXN ~~ a*TXN
a > 0.001
NEUR ~~ a*NEUR
a > 0.001
'


F5.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F5.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F5.a.fit.even


F4.a.model.alt <-'
F1 =~  SWB + NEUR 
F2 =~ PAIN1 + STRESS + PAIN2 + TIRED + LILL + SLEEP
F3 =~ ILLN + TXN + LILL
F4 =~ MIG + MED1 + MED2
FTin =~ TINADJ
TINADJ ~~0*TINADJ
MED1 ~~ a*MED1
a > 0.001
TXN ~~ a*TXN
a > 0.001
NEUR ~~ a*NEUR
a > 0.001
'

F4.a.fit.even.alt<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F4.a.model.alt, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F4.a.fit.even.alt



# $modelfit
# chisq df      p_chisq      AIC       CFI       SRMR
# df 200.6988 18 7.262435e-33 236.6988 0.9080621 0.08739148

# This is an iterative process.  Have to inspect CFA output for loadings, make sure there are no Heywood cases.
# For Heywood, add the following two lines for each Heywood case:
# ANX ~~ a*ANX;
# a>.001
#have to do a separate constant letter for each variable, otherwise you will fix their error variance to be the same. 
#inspect for Heywood again, sometimes it is a bit of Whac-A-Mole

EFA2.p.odd<-factanal(covmat = Ssmooth.odd,factors = 2, rotation = "promax")  #can use varimax
EFA2.p.odd

EFA3.p.odd<-factanal(covmat = Ssmooth.odd,factors = 3, rotation = "promax")  #can use varimax
EFA3.p.odd

EFA4.p.odd<-factanal(covmat = Ssmooth.odd,factors = 4, rotation = "promax")  #can use varimax
EFA4.p.odd

EFA5.p.odd<-factanal(covmat = Ssmooth.odd,factors = 5, rotation = "promax")  #can use varimax
EFA5.p.odd


EFA2v.p.odd<-factanal(covmat = Ssmooth.odd,factors = 2, rotation = "varimax")  #can use varimax
EFA2v.p.odd

EFA3v.p.odd<-factanal(covmat = Ssmooth.odd,factors = 3, rotation = "varimax")  #can use varimax
EFA3v.p.odd

EFA4v.p.odd<-factanal(covmat = Ssmooth.odd,factors = 4, rotation = "varimax")  #can use varimax
EFA4v.p.odd

EFA5v.p.odd<-factanal(covmat = Ssmooth.odd,factors = 5, rotation = "varimax")  #can use varimax
EFA5v.p.odd

write.table(EFA2.p.odd$loadings,file="ldhub_efa2.txt",row.names=T)
write.table(EFA3.p.odd$loadings,file="ldhub_efa3.txt",row.names=T)
write.table(EFA4.p.odd$loadings,file="ldhub_efa4.txt",row.names=T)
write.table(EFA5.p.odd$loadings,file="ldhub_efa5.txt",row.names=T)
write.table(EFA2v.p.odd$loadings,file="ldhub_efa2v.txt",row.names=T)
write.table(EFA3v.p.odd$loadings,file="ldhub_efa3v.txt",row.names=T)
write.table(EFA4v.p.odd$loadings,file="ldhub_efa4v.txt",row.names=T)
write.table(EFA5v.p.odd$loadings,file="ldhub_efa5v.txt",row.names=T)

