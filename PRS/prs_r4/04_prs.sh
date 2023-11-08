#Load PRScs and PLINK
module load PRScs/0878d30
module load plink

#awk '{print $3, toupper($4), toupper($5), $8, $9}' tinprs.txt > tinprs.txt.prscs

#Make an output directory
scratchpath=/scratch/scratch13/mvp039/vhasdcmaihoa/tinnitus_paper2_prs
mkdir "$scratchpath"

#Filter genotypes to the list of SNPs in the PRS file -- you need to link to the ducking IDs.
mvpprspath=/home/home1/vhasdcmaihoa/mvp039/AM/tinnitus_paper2_prs

#need to link our variants to the database,rename our variants to theirs.

zcat /data/data1/Public_Common_data/Annotation/dbSNP/db149/GRCh37p13/All_20161121.vcf.gz | grep -v "#" | awk '{print $3, $1":"$2":"$4":"$5}' | LC_ALL=C sort -k1b,1 > "$scratchpath"/All_20161121.vcf.gz.rs.pos.txt

cat "$mvpprspath"//tinprs.txt.prscs | LC_ALL=C sort -k1b,1 > "$mvpprspath"//tinprs.txt.prscs.sorted

LC_ALL=C join "$scratchpath"/All_20161121.vcf.gz.rs.pos.txt "$mvpprspath"//tinprs.txt.prscs.sorted > "$mvpprspath"//tinprs.txt.prscs.sorted.linked


echo "SNP A1 A2 BETA P" > "$mvpprspath"//PRScs.header
awk '{print $2,$3,$4,$5,$6}' "$mvpprspath"//tinprs.txt.prscs.sorted.linked | cat "$mvpprspath"//PRScs.header - > "$mvpprspath"//tinprs.txt.prscs.sorted.linked.header

awk '{print $1,$3,$4,$5,$6}' "$mvpprspath"//tinprs.txt.prscs.sorted.linked | cat "$mvpprspath"//PRScs.header - > "$mvpprspath"//tinprs.txt.prscs.sorted.linked.header1


echo "SNP MVP A1 A2 BETA P" > "$mvpprspath"//PRScs.header2
awk '{print $1,$2,$3,$4,$5,$6}' "$mvpprspath"//tinprs.txt.prscs.sorted.linked | cat "$mvpprspath"//PRScs.header2 - > "$mvpprspath"//tinprs.txt.prscs.sorted.linked.header.bothids


cp "$mvpprspath"//tinprs.txt.prscs.sorted.linked.header "$scratchpath"/.
cp "$mvpprspath"//tinprs.txt.prscs.sorted.linked.header1 "$scratchpath"/.

#Now must subset genomics data list to SNPlist... It will be almost all
#Need to keep in the original IDs for this step!!!!

cat /data/data1/mvp039/mvp_imputed/Release4_PGEN/chr*/*.pvar | grep -v "#" | awk '{print $1,$3,$2,$2,$4,$5}' >  "$scratchpath"/mvp_all.bim

cat "$scratchpath"/mvp_all.bim | LC_ALL=C sort -k2b,2  > "$scratchpath"/mvp_all.bim.sorted

LC_ALL=C join -1 2 -2 2 "$scratchpath"/mvp_all.bim.sorted <(cat "$mvpprspath"//tinprs.txt.prscs.sorted.linked.header.bothids | awk '{print $1,$2}' | LC_ALL=C sort -k2b,2)  > "$scratchpath"//mvp_all.bim.sorted.prscsin.bim.almost
awk '{print $2,$7,$3,$4,$5,$6}' "$scratchpath"//mvp_all.bim.sorted.prscsin.bim.almost > "$scratchpath"//mvp_all.bim.sorted.prscsin


head "$scratchpath"//mvp_all.bim.sorted.prscsin
wc -l "$scratchpath"//mvp_all.bim.sorted.prscsin #should be around 6.6 million

