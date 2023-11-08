# Tinnitus  
Codebase for 2023 Tinnitus GWAS manuscript (Clifford, Maihofer,..., Nievergelt)

Most analyses follow a main script that has been commented. Read for details.

### Meta analysis
Perform meta analyses

01_mvp_ukbb_meta.sh: Calls to METAL. Details which meta analysis script corresponds to which analysis.
metal_input_scripts: Scripts used in METAL


### CC GWAS
Code for CC-GWAS analysis comparing tinnitus to hearing difficulty

00_cc_gwas.sh: Analysis code for CC GWAS (data reformat and command)
CCGWAS_strand.R: CC GWAS normally throws out variants with ambiguous alleles. Here we have hacked CC GWAS function to skip this step (our compared data are on the same strand, we prefer to retain these variants).

### Fine-mapping
Polyfun based finemapping for summary statistics

00_polyfun.sh: Perform polyfun analysis for tinnitus
00_hearing_adjusted_polyfun.sh: Polyfun analysis for tinnitus adjusted (follows methods in first script)
01_credible_sets.sh: Calculate size of credible sets from polyfun outputs
02_wide_format_genes.txt: Manipulate genes.txt from FUMA to obtain a wide format list of genes within the risk locus
ldfiles_polyfun_tinnitus.xlsx: Polyfun requires huge LD reference files. The presupplied files cover specific windows. This file is used to determine which LD reference dataset file should be downloaded



### inversion
Inversion calling for chromosome 8

### MiXeR
Genetic architecture (univaraite and bivariate)
00_prepare_files.sh: Format files into MiXeR format
01_mixer_univariate.sh: Do univariate tests
02_mixer_univariate_plots.sh: Make univariate plots
03_mixer_bivariate.sh: Do bivariate tests
04_mixer_bivariate_plots.sh: Make bivariate plots


00_chr8_inversion_calling.sh: Calling in the UKBB data
01_gwas_stratified_by_inversion.sh: Perform tinnitus GWAS stratified by inversion status
inversion.mi - script to perform meta-analysis of inversion stratified GWAS
01_inversion_status_z.txt: meta-analysis of inversion results to determine whether the inversion itself is associated with tinnitus

### PRS
Polygenic risk score calculation


### rg
Genetic correlations of tinnitus to other phenotypes

00_rg_to_psychiatric_traits.sh: LDSC rg calculations with PGC psychiatric disorders
