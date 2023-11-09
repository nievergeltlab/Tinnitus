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


#ANXIETY AND SUICIDALITY DROPPED, HEARING EXCLUDED (Per CAROLINE, dec 17,2021)
require(data.table)
library(GenomicSEM)



trait.names<-c("TINADJ","ADHD",
               "ALCH", "AN","AUT",
               "BP", "MDD", "OCD",
               "PTSD","SCZ","TS")


#HEARING EXCLUDED!
traits <- c("TINADJ.sumstats.gz", "ADHD.sumstats.gz",
            "ALCH.sumstats.gz","AN.sumstats.gz", "AUT.sumstats.gz",
            "BP.sumstats.gz", "MDD.sumstats.gz",  "OCD.sumstats.gz",
            "PTSD.sumstats.gz", "SCZ.sumstats.gz",  "TS.sumstats.gz")

sample.prev <- c(NA, 0.358,
                 0.248, 0.234,0.397,
                 0.394, 0.341,0.276,
                 NA, 0.418,0.337) #NA for continuous  

population.prev <- c(NA,  0.05,
                     0.159, 0.009, 0.019,
                     0.044, 0.353, 0.025,
                     NA, 0.01, 0.008) #NA for continuous  

ld <- "/mnt/ukbb/royce/eur_w_ld_chr/"
wld <- "/mnt/ukbb/royce/eur_w_ld_chr/"



#Run the LDSC three ways. All for reporting sample genetic correlations and models not needing EFA; odd for EFA; even for CFA (if needed)
LDSCoutput.all <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, ldsc.log="LOG_ldsc_TINHD_PSYCHPGC_all")
LDSCoutput.odd <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, select = "ODD", ldsc.log="LOG_ldsc_TINHD_PSYCHPGC_odd")
LDSCoutput.evn <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, select = "EVEN", ldsc.log="LOG_ldsc_TINHD_PSYCHPGC_even")



##optional command to save the ldsc output in case you want to use it in a later R session.
#These are the files used for EFA/CFA
save(LDSCoutput.all, file="ldsc_psych_jan3_2022_all.RData")
save(LDSCoutput.odd, file="ldsc_psych_jan3_2022_odd.RData")
save(LDSCoutput.evn, file="ldsc_psych_jan3_2022_even.RData")


load(file="ldsc_psych_jan3_2022_all.RData")
load(file="ldsc_psych_jan3_2022_odd.RData")
load(file="ldsc_psych_jan3_2022_even.RData")


#PLOTS
require(corrplot)

rownames(LDSCoutput.all$S_Stand) <- trait.names

pdf('tinnitus_psych_jan3_2022.pdf',12,12)
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

write.table(EFA2.p.odd$loadings,file="psych_efa2.txt",row.names=T)
write.table(EFA3.p.odd$loadings,file="psych_efa3.txt",row.names=T)
write.table(EFA4.p.odd$loadings,file="psych_efa4.txt",row.names=T)
write.table(EFA5.p.odd$loadings,file="psych_efa5.txt",row.names=T)
write.table(EFA2v.p.odd$loadings,file="psych_efa2v.txt",row.names=T)
write.table(EFA3v.p.odd$loadings,file="psych_efa3v.txt",row.names=T)
write.table(EFA4v.p.odd$loadings,file="psych_efa4v.txt",row.names=T)
write.table(EFA5v.p.odd$loadings,file="psych_efa5v.txt",row.names=T)



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
F1 =~ TINadj+ADHD+ALCH+AN+AUT+BP+MDD+OCD+PTSDc+SCZ+TS
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
F1 =~ TINADJ + ADHD + ALCH + AUT + BP +MDD + OCD+ PTSD+SCZ
F2 =~ ALCH + AN + MDD + OCD + TS
OCD ~~ a*OCD
a > 0.001
'

F2.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F2.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F2.a.fit.even


F3.a.model.even<-'
F1 =~ TINADJ + ADHD + ALCH +AUT + MDD+ PTSD
F2 =~ ALCH + AN + MDD + OCD + TS
F3 =~ BP + SCZ

'

F3.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F3.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F3.a.fit.even


F4.a.model.even<-'
F1 =~ ADHD + ALCH + MDD+ PTSD
F2 =~ AN + MDD + OCD + TS
F3 =~ TINADJ + ALCH +AUT
F4 =~ BP + SCZ

ALCH~~ a*ALCH
a > 0.001
'

F4.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F4.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F4.a.fit.even


F5.a.model.even<-'
F1 =~ ADHD + ALCH + MDD+  PTSD
F2 =~ AN + MDD + OCD + TS
F3 =~ BP + SCZ
F4 =~ AUT
F5 =~ TINADJ 


'


F5.a.fit.even<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F5.a.model.even, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F5.a.fit.even


F4.a.model.alt <-'
F1 =~  ADHD + ALCH + MDD+ PTSD
F2 =~ AN+MDD + OCD + TS
F3 =~ ALCH + AUT
F4 =~ BP + SCZ
FTin =~ TINADJ
TINADJ ~~0*TINADJ

'

F4.a.fit.even.alt<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F4.a.model.alt, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F4.a.fit.even.alt #FAILS TO CONVERGE!!



F3.a.model.alt <-'
F1 =~ ADHD + ALCH +AUT + MDD+ PTSD
F2 =~ ALCH + AN + MDD + OCD + TS
F3 =~ BP + SCZ
FTin =~ TINADJ
TINADJ ~~0*TINADJ
'

F3.a.fit.even.alt<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F3.a.model.alt, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F3.a.fit.even.alt


F3.a.model.alt2 <-'
F1 =~ ADHD + ALCH + MDD+ PTSD
F2 =~ ALCH + AN + MDD + OCD + TS
F3 =~ BP + SCZ
FTin =~ TINADJ
TINADJ ~~0*TINADJ
FAUT =~ AUT
AUT ~~0*AUT

'

F3.a.fit.even.alt2<-usermodel(LDSCoutput.evn, estimation = "DWLS", model = F3.a.model.alt2, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
F3.a.fit.even.alt2





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

write.table(EFA2.p.odd$loadings,file="psych_efa2.txt",row.names=T)
write.table(EFA3.p.odd$loadings,file="psych_efa3.txt",row.names=T)
write.table(EFA4.p.odd$loadings,file="psych_efa4.txt",row.names=T)
write.table(EFA5.p.odd$loadings,file="psych_efa5.txt",row.names=T)
write.table(EFA2v.p.odd$loadings,file="psych_efa2v.txt",row.names=T)
write.table(EFA3v.p.odd$loadings,file="psych_efa3v.txt",row.names=T)
write.table(EFA4v.p.odd$loadings,file="psych_efa4v.txt",row.names=T)
write.table(EFA5v.p.odd$loadings,file="psych_efa5v.txt",row.names=T)

