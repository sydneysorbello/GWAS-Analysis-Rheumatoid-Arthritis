## FINAL PROJECT
# Sydney Sorbello
# Due: 4/7


# let set the path to the NARAC data
BASE=/projectnb/bs859/data/RheumatoidArthritis/final_project/narac_hg19


# lets load some key packages
module load R
module load plink/1.90b6.27
module load eigensoft
module load metal


# 1: Preprocessing and PCA analysis
# First, let's filter the dataset based on the following requirements
# MAF 0.01
# Fewer than 5% missing genotypes
# hwe p-value greater than 0.000001
# individuals with less than 5% genotype misssingness
plink --bfile $BASE --maf 0.01 --geno 0.05 --hwe 1e-6 --mind 0.05 --make-bed --out narac_filtered

# here we can check the sex assignmnet aligns with the genotypes
plink --bfile narac_filtered --check-sex --out narac_sexcheck
awk '$5 == "PROBLEM"' narac_sexcheck.sexcheck > failed_sexcheck.txt

# the sex check revealed 7 individuals with problems during the sex check, lets remove them from the dataset
cut -f1,2 failed_sexcheck.txt > remove_sexcheck.txt
plink --bfile narac_filtered --remove remove_sexcheck.txt --make-bed --out narac_sexcheck_filtered

# now we check for relatedness. the data description noted that the participants should not be related, but we check
regardless
plink --bfile narac_sexcheck_filtered --genome --out narac_ibd
awk '$10 > 0.185' narac_ibd.genome > related_pairs.txt
wc -l related_pairs.txt

# now we prune LD from the cleaned data
plink --bfile narac_sexcheck_filtered --indep-pairwise 10000kb 1 0.2 --out narac_pruned

# finally we extract the pruned snps
plink --bfile narac_sexcheck_filtered --extract narac_pruned.prune.in --make-bed --out narac_pca

# now we perform PCA
smartpca -p narac_pca.par > narac_pca.log

# lets plot the first two principal components
Rscript --vanilla plotPCs.R narac_pca.evec 1 2 10
Rscript --vanilla plotPCs.R narac_pca.evec 1 4 10
Rscript --vanilla plotPCs.R narac_pca.evec 2 4 10

# lets find the outliers in the PC1 cluster < -0.067
awk '$2 < -0.065 { split($1, id, ":"); print id[1], id[2] }' narac_pca.evec > pca_outliers.txt

# and remove them from the dataset
plink --bfile narac_sexcheck_filtered --remove pca_outliers.txt --make-bed --out narac_pcaout_filtered
plink --bfile narac_pca --remove pca_outliers.txt --make-bed --out narac_pcaout_filtered


# Part 2A: Sex Stratified Analysis
# the covariate file for the dataset has already been supplied
# however, we need to align this file with the individuals we are still working with in the dataset
# this awk command with create a new covariate file consistent with the current individuals
(head -n 1 narac.cov && awk 'NR==FNR {keep[$1,$2]; next} ($1,$2) in keep' narac_pcaout_filtered.fam narac.cov) >
narac_cov_synced.cov

# we pick out PC1, PC2, and PC4 so that we can add them to the covariate file
awk '{ split($1, id, ":"); print id[1], id[2], $2, $3, $5 }' narac_pca.evec > pca_selected.txt

# save and remove header
## ONLY WORKS directly in terminal
head -n 1 narac_cov_synced.cov > cov_header.txt
tail -n +2 narac_cov_synced.cov > narac_cov_noheader.cov

# make sure the files are formatted similarly
sed -i 's/\r$//' pca_selected.txt
sed -i 's/\r$//' narac_cov_noheader.cov

# and merge the files
{ echo -e "FID IID RA sex SEN antiCCP IgM PC1 PC2 PC4"; awk 'NR==FNR && $1 !~ /^#/ {a[$1]=$3"\t"$4"\t"$5; next} ($1 in
a) {print $0 "\t" a[$1]}' pca_selected.txt narac_cov_noheader.cov; } > merged_output.txt

# in order to continue with the sex-stratified analysis, we split the covariate file by the sex
# Female-only covariates
awk '$4 == 2' merged_output.txt > narac_females_cov.cov

# Male-only covariates
awk '$4 == 1' merged_output.txt > narac_males_cov.cov

# set the new header for the file
header="FID IID RA sex SEN antiCCP IgM PC1 PC2 PC4"

# Add to female covariates
{ echo -e "$header"; cat narac_females_cov.cov; } > narac_females_cov_with_header.cov

# Add to male covariates
{ echo -e "$header"; cat narac_males_cov.cov; } > narac_males_cov_with_header.cov

# now lets perform GWAS using the selected PCs for females
plink --bfile narac_pcaout_filtered --covar narac_females_cov_with_header.cov --covar-name PC1,PC2,PC4 --logistic beta
hide-covar --ci .95 --out narac_female_gwas

# and the same for the male population
plink --bfile narac_pcaout_filtered --covar narac_males_cov_with_header.cov --covar-name PC1,PC2,PC4 --logistic beta
hide-covar --ci .95 --out narac_male_gwas

# Lets show the qqplot for bth analyses
Rscript --vanilla qqplot.R narac_female_gwas.assoc.logistic PC124adj_Female ADD
Rscript --vanilla qqplot.R narac_male_gwas.assoc.logistic PC124adj_Male ADD

# lets take a look at significant SNPs
awk '$9 < 5e-8' narac_female_gwas.assoc.logistic > female_sig_snps.txt
awk '$9 < 5e-8' narac_male_gwas.assoc.logistic > male_sig_snps.txt

# and check the size of these files
wc female_sig_snps.txt
wc male_sig_snps.txt

# lets make a manhattan plot
Rscript --vanilla gwaplot_pgen.R narac_female_gwas.assoc.logistic "Manhattan Plot for PC1,2,4 Female"
manhattan_pc124_female2
Rscript --vanilla gwaplot_pgen.R narac_male_gwas.assoc.logistic "Manhattan Plot for PC1,2,4 Male" manhattan_pc124_male


# PART 2B: Meta analysis
# first we have to reformat the plink files for METAL
awk 'BEGIN {print "A2"} {print $6}' narac_pcaout_filtered.bim > A2.txt

# lets add A2 to the female data
awk 'NR==FNR {a[NR]=$1; next} {print $1, $2, $3, $4, a[FNR], $5, $6, $7, $8, $9, $10, $11, $12}' A2.txt
narac_female_gwas.assoc.logistic > narac_female_final.assoc.logistic

# and to the male dataset
awk 'NR==FNR {a[NR]=$1; next} {print $1, $2, $3, $4, a[FNR], $5, $6, $7, $8, $9, $10, $11, $12}' A2.txt
narac_male_gwas.assoc.logistic > narac_male_final.assoc.logistic

# and perform the meta analysis
metal metal.txt > metal.log

# lets take a look at the most significant gene and compare between the GWAS analyses
awk 'NR==1||$1=="rs660895"{print $0}' METAANALYSIS1.TBL
awk 'NR==1||$2=="rs660895"{print $0}' narac_female_gwas.assoc.logistic
awk 'NR==1||$2=="rs660895"{print $0}' narac_male_gwas.assoc.logistic

# lets take a look at the 9 most significant SNPs for all analyses
(head -n 1 METAANALYSIS1.TBL && tail -n +2 METAANALYSIS1.TBL | sort -k6,6g) | head
(head -n 1 narac_female_gwas.assoc.logistic && tail -n +2 narac_female_gwas.assoc.logistic | awk '$12 != "NA"' | sort
-k12,12g) | head
(head -n 1 narac_male_gwas.assoc.logistic && tail -n +2 narac_male_gwas.assoc.logistic | awk '$12 != "NA"' | sort
-k12,12g) | head

# Lets create a file for the meta data that can be used to visualize the data
awk 'NR==1 {print "SNP P"} {print $1, $6}' METAANALYSIS1.TBL > meta_p.txt
awk 'NR != 2' meta_p.txt > tmp && mv tmp meta_p.txt

# remove the P value column so we can insert the P from the meta analysis
awk '{
for (i = 1; i <= NF; i++) {
if (i != 12) {
printf "%s%s", $i, (i == NF || i == 11 ? ORS : OFS)
}
}
}' narac_female_gwas.assoc.logistic > meta_without_p.txt

# and add the p value fromt he meta analysis
awk '
NR==FNR && FNR > 1 { pval[$1] = $2; next }
FNR==1 { print $0, "P"; next }
{ print $0, (pval[$2] ? pval[$2] : "NA") }
' meta_p.txt meta_without_p.txt > meta_joined.txt

# remove the first row
awk 'NR > 1' meta_joined.txt > metaanalysis.txt

# finaize the file
sed -i 's/\r$//' metaanalysis.txt

# make a qqplot
Rscript --vanilla qqplot.R metaanalysis.txt Metaanalysis_QQPlot ADD

# and a manhattan plot
Rscript --vanilla gwaplot_pgen.R metaanalysis.txt "Manhattan Plot for Metaanalysis" manhattan_metaanalysis


### Task 3C: LD Score Regression
# first lets initilialize the path to the final folder
OKADA=/projectnb/bs859/data/RheumatoidArthritis/final_project

# and import the necessary packages
module load python2
module load ldsc

### ALL
# lets take a look at the data from both populations
zcat $OKADA/RA_GWASmeta_TransEthnic_v2.txt.gz|head
zcat $OKADA/RA_GWASmeta_TransEthnic_v2.txt.gz|wc

# format OKADA summary statistics for ldsc
# the reformatting process is repeated for the three analyses
zcat $OKADA/RA_GWASmeta_TransEthnic_v2.txt.gz | awk 'BEGIN {OFS="\t"} {
if ($1 ~ /:/) {
split($1, a, ":");
$1 = "rs" a[2]
}
print
}' > RA_ALL_RS.txt

# add a column for BETA
gawk 'BEGIN {OFS="\t"}
NR==1 {print $0, "Beta"; next}
{print $0, log($6)}' RA_ALL_RS.txt > RA_ALL_BETA.txt

# and the sample size
awk 'BEGIN {OFS="\t"}
NR==1 {print $0, "N"; next}
{print $0, 80799}' RA_ALL_BETA.txt > RA_ready.txt

# rename the column OR_A1 as just OR
awk 'BEGIN {OFS="\t"} NR==1 { $6 = "OR" } { print }' RA_ready.txt > RA_ALL.txt

# intialize the path to the LD files
export LDDIR='/projectnb/bs859/data/ldscore_files'

# SNPID Chr Position_hg19 A1 A2 OR_A1 OR_95CIlow OR_95CIup P-val
# reformat the file for ldsc with munge sumstats
munge_sumstats.py \
--sumstats RA_ALL.txt \
--snp SNPID \
--a1 A1 \
--a2 A2 \
--p P-val \
--signed-sumstats Beta,0 \
--N-col N \
--merge-alleles $LDDIR/w_hm3.snplist \
--out RA_ALL_DONE

# now we can run ldsc
ldsc.py \
--h2 RA_ALL_DONE.sumstats.gz \
--ref-ld $LDDIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDDIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out RA_ALL_LDSC

#### EUR
# lets take a look at the data from both populations
head $OKADA/RA_GWASmeta_European_v2.txt
wc $OKADA/RA_GWASmeta_European_v2.txt

# format OKADA summary statistics for ldsc
awk 'BEGIN {OFS="\t"} {
if ($1 ~ /:/) {
split($1, a, ":");
$1 = "rs" a[2]
}
print
}' $OKADA/RA_GWASmeta_European_v2.txt > RA_EUR_RS.txt
gawk 'BEGIN {OFS="\t"}
NR==1 {print $0, "Beta"; next}
{print $0, log($6)}' RA_EUR_RS.txt > RA_EUR_BETA.txt
awk 'BEGIN {OFS="\t"}
NR==1 {print $0, "N"; next}
{print $0, 80799}' RA_EUR_BETA.txt > RA_EUR_ready.txt
awk 'BEGIN {OFS="\t"} NR==1 { $6 = "OR" } { print }' RA_EUR_ready.txt > RA_EUR.txt
export LDDIR='/projectnb/bs859/data/ldscore_files'

# SNPID Chr Position_hg19 A1 A2 OR_A1 OR_95CIlow OR_95CIup P-val
munge_sumstats.py \
--sumstats RA_EUR.txt \
--snp SNPID \
--a1 A1 \
--a2 A2 \
--p P-val \
--signed-sumstats Beta,0 \
--N-col N \
--merge-alleles $LDDIR/w_hm3.snplist \
--out RA_EUR_DONE

# now we can run ldsc
ldsc.py \
--h2 RA_EUR_DONE.sumstats.gz \
--ref-ld $LDDIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--w-ld $LDDIR/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out RA_EUR_LDSC

#### EAS
# lets take a look at the data from both populations
zcat $OKADA/RA_GWASmeta_Asian_v2.txt.gz|head
zcat $OKADA/RA_GWASmeta_Asian_v2.txt.gz|wc

# format OKADA summary statistics for ldsc
zcat $OKADA/RA_GWASmeta_Asian_v2.txt.gz | awk 'BEGIN {OFS="\t"} {
if ($1 ~ /:/) {
split($1, a, ":");
$1 = "rs" a[2]
}
print
}' > RA_EAS_RS.txt
gawk 'BEGIN {OFS="\t"}
NR==1 {print $0, "Beta"; next}
{print $0, log($6)}' RA_EAS_RS.txt > RA_EAS_BETA.txt
awk 'BEGIN {OFS="\t"}
NR==1 {print $0, "N"; next}
{print $0, 80799}' RA_EAS_BETA.txt > RA_ready_EAS.txt
awk 'BEGIN {OFS="\t"} NR==1 { $6 = "OR" } { print }' RA_ready_EAS.txt > RA_EAS.txt
export LDDIR='/projectnb/bs859/data/ldscore_files'
munge_sumstats.py \
--sumstats RA_EAS.txt \
--snp SNPID \
--a1 A1 \
--a2 A2 \
--p P-val \
--signed-sumstats Beta,0 \
--N-col N \
--merge-alleles $LDDIR/w_hm3.snplist \
--out RA_EAS_DONE

# now we can run ldsc
ldsc.py \
--h2 RA_EAS_DONE.sumstats.gz \
--ref-ld $LDDIR/UKBB.ALL.ldscore/UKBB.EAS.rsid \
--w-ld $LDDIR/UKBB.ALL.ldscore/UKBB.EAS.rsid \
--out RA_EAS_LDSC