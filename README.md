# Tinnitus  
Codebase for 2023 Tinnitus GWAS manuscript

Description of contents

CC GWAS
Code for CC-GWAS analysis comparing tinnitus to hearing difficulty

CCGWAS_strand.R: CC GWAS normally throws out variants with ambiguous alleles. Here we have hacked CC GWAS function to skip this step (our compared data are on the same strand, we prefer to retain these variants).
00_cc_gwas.sh: Analysis code for CC GWAS (data reformat and command)
