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


# Summary Stat Data File Description --------------------------------------

# 1) ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz            #updated 07/21
# PGC META Tinnitus
# meta
# Sample prev=continuous
# Pop prev=0.143
# N=482031
# 
# 2) ukbb_hd_mvp_hearing_dec17_20201.tbl.fuma.gz
# HDmeta
# METAMETA HD
# binary
# Sample prev=0.359
# Pop prev=0.2
# N=490462
#
# 3) adhd_eur_jun2017.gz
# PGC ADHD
# binary
# Sample prev=0.358
# Pop prev=0.05
# N=53293
# 
# 4)pgc_alcdep.eur_discovery.aug2018_release.txt.gz
# PGC Alcohol Dependence
# binary
# Sample prev=0.248
# Pop prev=0.159
# N=46568
# 
# 5) iPSYCH-PGC_ASD_Nov2017.gz
# PGC Autism
# binary
# Sample prev=0.397
# Pop prev=0.012
# N=46351
# 
# 6)PGC_UKB_depression_genome-wide.txt
# PGC MDD
# binary
# Sample prev=0.341
# Pop prev=0.353
# N=500199
# 
# 7)eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz
# PGC PTSD continuous
# continuous
# Sample prev=
# Pop prev=
# N=182,199
# 
# 8)anxiety.meta.full.cc.tbl.gz
# PGC Anxiety
# binary
# Sample prev=0.317
# Pop prev=0.1
# N=17310

# 9) pgcAN2.2019-07.vcf.tsv.gz
# AN
# PGCAnorexia Nervosa
# binary
# Sample prev=0.234
# Pop prev=0.009
# N=72517
#
# 10) daner_PGC_BIP32b_mds7a_0416a.gz
# BP
# PGC Bipolar
# binary
# Sample prev=0.394
# Pop prev=0.044
# N=51710
# 
# 11)ocd_aug2017.gz
# OCD
# PGCOCD
# binary
# Sample prev=0.276
# Pop prev=0.025
# N=9725
# 
# 12)PGC3_SCZ_wave3_public.v2.tsv.gz
# SCZ
# PGCSchizophrenia
# binary
# Sample prev=0.418
# Pop prev=0.01
# N=161405
# 
# 13)TS_Oct2018.gz
# TS
# PGCTourettes
# binary
# Sample prev=0.337
# Pop prev=0.008
# N=14307
# 


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
# # All  files being combined and munged - make sure that everything is in the right order!!!


#ANXIETY AND SUICIDALITY DROPPED, HEARING EXCLUDED (Per CAROLINE, dec 17,2021)
require(data.table)
library(GenomicSEM)

munge(files=c("ukbbcoding3relatedhearcov_mvpanytinnitusanyhear_ssw1.tbl.fuma.gz", 
              "pgc_sumstats/adhd_eur_jun2017.gz",
              "pgc_sumstats/pgc_alcdep.eur_discovery.aug2018_release.txt-wide.txt.gz",
              "pgc_sumstats/pgcAN2.2019-07.vcf.tsv.gz",
              "pgc_sumstats/iPSYCH-PGC_ASD_Nov2017.gz",
              "pgc_sumstats/daner_PGC_BIP32b_mds7a_0416a.gz",
              "pgc_sumstats/PGC_UKB_depression_genome-wide.txt",
              "pgc_sumstats/ocd_aug2017.gz",
              "pgc_sumstats/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz",
              "pgc_sumstats/PGC3_SCZ_wave3_public.v2.tsv.gz",
              "pgc_sumstats/TS_Oct2018.gz"),
      hm3 = "/mnt/ukbb/adam/tinnitus_gwas/w_hm3.noMHC.snplist",
      trait.names=c("TINadj", "ADHD",
                    "ALCH","AN", "AUT",
                    "BP","MDD","OCD",
                    "PTSDc", "SCZ","TS"
                    ),
      N=c(467369,53293,
          46568,72517,46351,
          51710,500199,9725,
          182199,161405,14307),
      info.filter = 0.9, maf.filter = 0.01)

# Run multivariable LDSC --------------------------------------------------

# sample.prev: A vector of sample prevalences of length equal to the number of
# traits. Sample prevalence is calculated as the number of cases over the total
# number of participants (i.e., cases + controls). Possible range = 0-1. If the
# trait is continuous, the values should equal NA. population.prev: A vector of
# population prevalences. These estimates can be obtained from a number of
# sources, such as large scale epidemiological studies. Possible range = 0-1.
# Again, if the trait is continuous the values should equal NA.

#HEARING EXCLUDED!
traits <- c("TINadj.sumstats.gz", "ADHD.sumstats.gz",
            "ALCH.sumstats.gz","AN.sumstats.gz", "AUT.sumstats.gz",
            "BP.sumstats.gz", "MDD.sumstats.gz",  "OCD.sumstats.gz",
            "PTSDc.sumstats.gz", "SCZ.sumstats.gz",  "TS.sumstats.gz")

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

trait.names<-c("TINadj","ADHD",
               "ALCH", "AN","AUT",
               "BP", "MDD", "OCD",
               "PTSDc","SCZ","TS")

#Run the LDSC three ways. All for reporting sample genetic correlations and models not needing EFA; odd for EFA; even for CFA (if needed)
LDSCoutput.all <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, ldsc.log="LOG_ldsc_TINHD_PSYCHPGC_all")
LDSCoutput.odd <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, select = "ODD", ldsc.log="LOG_ldsc_TINHD_PSYCHPGC_odd")
LDSCoutput.evn <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand=TRUE, select = "EVEN", ldsc.log="LOG_ldsc_TINHD_PSYCHPGC_even")

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
save(LDSCoutput.all, file="ldsc_TBI_2021120_all.RData")
save(LDSCoutput.odd, file="ldsc_TBI_2021120_odd.RData")
save(LDSCoutput.evn, file="ldsc_TBI_2021120_even.RData")


