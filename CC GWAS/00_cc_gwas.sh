#reformat sumstats
zcat metal_results/ukbb_linprobit_mvp_tinnitus_ivw1.tbl.fuma.gz                 | awk 'OFS="\t" {beta="OR";se="SE"; if(NR>1) {beta=exp($7); se=$8; Neff="511443"} ; if (NR==1) {$1="CHR";$2="BP";$3="SNP"; $4="EA";$5="NEA";$6="FRQ";Neff="Neff";$9="P"; beta="OR"; se="SE"}; print $3,$1,$2,toupper($4),toupper($5),$6,beta,se,$9,Neff} ' | gzip > metal_results/ukbb_linprobit_mvp_tinnitus_ivw1.tbl.fuma.gz.ccgwas.gz
zcat metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_2020_vw1.tbl.fuma.gz    | awk 'OFS="\t" {beta="OR";se="SE"; if(NR>1) {beta=exp($7); se=$8; Neff="450950"} ; if (NR==1) {$1="CHR";$2="BP";$3="SNP"; $4="EA";$5="NEA";$6="FRQ";Neff="Neff";$9="P"; beta="OR"; se="SE"}; print $3,$1,$2,toupper($4),toupper($5),$6,beta,se,$9,Neff} ' | gzip > metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_2020_vw1.tbl.fuma.gz.ccgwas.gz


R
library(data.table)
library(R.utils)
library(devtools)
#install_github("wouterpeyrot/CCGWAS")
library(CCGWAS)

#to not remove strand ambiguous vars
 source('CCGWAS_strand.R')
  
CCGWAS_strand( outcome_file = "ccgwas_tin_v_hear_metas_more_overlap_ambig" , A_name = "tinnitusmeta" , B_name = "hearingmeta" , 
        sumstats_fileA1A0 = "metal_results/ukbb_linprobit_mvp_tinnitus_ivw1.tbl.fuma.gz.ccgwas.gz" ,
        sumstats_fileB1B0 = "metal_results/hearing_meta/ukbb_hd_mvp_hearing_dec17_2020_vw1.tbl.fuma.gz.ccgwas.gz" ,
        K_A1A0 = 0.35 , K_A1A0_high = 0.50 , K_A1A0_low = 0.2 ,  
        K_B1B0 = 0.31 , K_B1B0_high = 0.50  , K_B1B0_low = 0.2, 
        h2l_A1A0 = 0.07 , h2l_B1B0 = 0.07 , rg_A1A0_B1B0 = 0.5701 , intercept_A1A0_B1B0 =0.2202 , m = 1e4 ,  
        N_A1 = 173703  , N_B1 = 175838 , N_A0 = 308171 , N_B0 = 391364  , N_overlap_A0B0 = 169495 ,save.all=TRUE ) #for ukb, the a0b0 overlap is 90k..55%. I'll keep this overlap for MVP.
        
