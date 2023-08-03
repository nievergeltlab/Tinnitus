
#!/bin/bash

 
 module load 2019
 module load Boost.Python/1.67.0-intel-2019b-Python-3.6.6
 module load Tk/8.6.8-GCCcore-8.3.0 
 
 #When you start LISA:
 #LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 
 
 cd  /home/maihofer/ehr_mixer/results_filtered
 
##Convert sumstats

#Studies:

 study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
 study2=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz
 
 study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
 study2=TotalPCL_MVP_eur.gz

 study1=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz
 study2=TotalPCL_MVP_eur.gz

# study1=eur_ptsd_pcs_v4_aug3_2021.fuma.gz 
# study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
# study1=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz 
# study1=TotalPCL_MVP_eur.gz
# study1=eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz


# study2=pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz
  
#study1=


#

study2=pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz

eur_ptsd_pcs_v4_aug3_2021.fuma.gz 
eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz
 TotalPCL_MVP_eur.gz
 
 #Singular matrix for EHR, case/control
 
for study1 in eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz   eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz
do
echo $study1
studycomb="$study1"_"$study2"

#Bivarate plot based on fit
  python /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$studycomb".rep@.json --out      "$studycomb".rep.fit         
  #python /home/maihofer/mixer/precimed/mixer_figures.py two --json "$studycomb".rep.fit.json  --out "$studycomb".rep.combine  --statistic mean std
 
# Bivariate plot based on total    
  python /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$studycomb".test.rep@.json --out "$studycomb".test.rep.fit
  #python /home/maihofer/mixer/precimed/mixer_figures.py two --json "$studycomb".test.rep.fit.json  --out "$studycomb".test.rep.combine  --statistic mean std
     

  python /home/maihofer/mixer/precimed/mixer_figures.py two --json-fit "$studycomb".rep.fit.json --json-test "$studycomb".test.rep.fit.json  --out "$studycomb".test.rep.combine  --statistic mean std
     
   done


# python precimed/mixer_figures.py combine --json <prefix>.fit.rep@.json --out <prefix>.fit
# python precimed/mixer_figures.py two --json <prefix>.fit.json --out <out_file> --statistic mean std