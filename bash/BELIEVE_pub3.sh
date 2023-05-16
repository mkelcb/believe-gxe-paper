# Main UKBB vs BELIEVE GxE analysis, v2, +whole genome GWAS

#################################

###################
# 0 preliminary QC: create reduced input for R-scripts:
###################


# Keep lists: Relationships_TWIST-GxS_README_v2.0.txt (these are all have the IIDs from the .fams)
# CAMBRIDGE-BELIEVE_Freeze_One.GxS.commonsnps_samples_ancestries.txt: hapmap populations(if $4 != "SAS", then exclude)
# CAMBRIDGE-BELIEVE_Freeze_One.GxS.genome.FILTERED.genome.3rd_degree_unrelated: list of unrelated indis (~33K)
# CAMBRIDGE-BELIEVE_Freeze_One.GxS.commonsnps.eigenvec: 2 cols for IDs, rest is for the 20PCs
# CAMBRIDGE-BELIEVE_Freeze_One.GxS.commonsnps.SNPs_used_for_PCA.bim: SNPs used for PCA (only 5K)

# remove ancestry mismatched indis from the 3rd degree unrelateds

awk 'FNR == NR { 
if($4 != "SAS") { file1[$1] = $1; }
next; 
} FNR <= NR { 
if ( $1 in file1 == 0) {  print $0"\t"$0 }
 }' $relatedBaseLoc$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.commonsnps_samples_ancestries.txt'  $relatedBaseLoc$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.genome.FILTERED.genome.3rd_degree_unrelated' > $believeScratchLoc$'indi_keeplist'
head $believeScratchLoc$'indi_keeplist'
wc -l $relatedBaseLoc$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.genome.FILTERED.genome.3rd_degree_unrelated'
wc -l $believeScratchLoc$'indi_keeplist' # 33712

head '/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/genotype/genomewide/plink/may_2022/CAMBRIDGE-BELIEVE_Freeze_One.GxS.bim'
wc -l '/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/genotype/genomewide/plink/may_2022/CAMBRIDGE-BELIEVE_Freeze_One.GxS.bim' # 1,414,695 variants
wc -l '/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/genotype/genomewide/plink/may_2022/CAMBRIDGE-BELIEVE_Freeze_One.GxS.fam' # 5,4620 participants

###########################################################################
# I. Liftover BELIEVE from hg38 back to hg19 to make BELIEVE compatible with UKBB
###########################################################################

mkdir -p $believeRawLoc$'liftover/'
cd $believeRawLoc$'liftover/'

# liftover
# https://www.biostars.org/p/252938/
# download database and script
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/GscPythonUtility/master/liftOverPlink.py
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/ibdqc.pl

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz


head -n 2 $believeScratchLoc$'indi_keeplist' > $believeScratchLoc$'indi_keeplist_2'

# create a dummy plink file with 2 indis, as we will be transcoding to text plink files
# rebuild plink file to avoid chromsome-miss-order problem
arguments=' --memory 43000 --bfile '$believeGxS$' --keep '$believeScratchLoc$'indi_keeplist_2 --make-bed --out '$believeScratchLoc$'liftOverDummy.sort --allow-extra-chr --allow-no-sex'
$plink $arguments

# space to tab to generate bed files for liftOver from hg18 to hg19
$plink --bfile $believeScratchLoc$'liftOverDummy.sort' --recode tab --out $believeScratchLoc$'liftOverDummy.sort.tab'

# apply liftOverPlink.py to update hg18 to hg19 or hg38
# Only works in BIRC10, Not HPC, Caused by Python version
mkdir liftOver
chmod 0777 liftOver
python2 liftOverPlink.py -m $believeScratchLoc$'liftOverDummy.sort.tab.map' -p $believeScratchLoc$'liftOverDummy.sort.tab.ped' -o $believeScratchLoc$'liftOverDummy.hg19' -c hg38ToHg19.over.chain.gz -e ./liftOver

# convert back to binary bim
arguments=' --memory 43000 --file '$believeScratchLoc$'liftOverDummy.hg19 --make-bed --out '$believeScratchLoc$'liftedOver --allow-extra-chr --allow-no-sex'
$plink $arguments
head $believeScratchLoc$'liftedOver.bim'
wc -l $believeScratchLoc$'liftedOver.bim'

# the above may not keep the A1/A2, as we only kept 2 indis...
# we find the original A1/A2 and match them via the SNPid
awk 'FNR == NR { 
file1[ $2 ] = $5"\t"$6;  next; }
FNR <= NR {  
if( $2 in file1 ) {print $1"\t"$2"\t"$3"\t"$4"\t"file1[$2]} }
' $believeGxS$'.bim' $believeScratchLoc$'liftedOver.bim' > $believeScratchLoc$'liftedOver_origA1A2.bim'
head $believeScratchLoc$'liftedOver_origA1A2.bim'


# now get the matching RSIds
awk 'FNR == NR { 
file1[ "chr"$1"_"$4"_"$5"_"$6 ] = $2; file2[ "chr"$1"_"$4"_"$6"_"$5 ] = $2; next; }
FNR <= NR {  
A1=$5
A2=$6
if(A1 == "G") {A1_rev="C"}
else if(A1 == "C") {A1_rev="G"}
else if(A1 == "A") {A1_rev="T"}
else if(A1 == "T") {A1_rev="A"}

if(A2 == "G") {A2_rev="C"}
else if(A2 == "C") {A2_rev="G"}
else if(A2 == "A") {A2_rev="T"}
else if(A2 == "T") {A2_rev="A"}

lookup1="chr"$1"_"$4"_"A1"_"A2
lookup2="chr"$1"_"$4"_"A1_rev"_"A2_rev

if( lookup1 in file1 ) {print $2"\t"file1[lookup1] } 
else if(lookup1 in file2 ) {print $2"\t"file2[lookup1] } 
else if(lookup2 in file1 ) {print $2"\t"file1[lookup2] } 
else if(lookup2 in file2) {print $2"\t"file2[lookup2] } }
' $hapmap3_b37bim $believeScratchLoc$'liftedOver_origA1A2.bim' > $believeRawLoc$'GxS.rsids'
head $believeRawLoc$'GxS.rsids'
wc -l $believeRawLoc$'GxS.rsids' # 195, 216


awk '{print $2}' $believeRawLoc$'GxS.rsids' > $believeRawLoc$'GxS.rsidsjustRSid' 
head $believeRawLoc$'GxS.rsidsjustRSid'

awk '{print $1}' $believeRawLoc$'GxS.rsids' > $believeRawLoc$'GxS.rsidsjustGxSIds' 
head $believeRawLoc$'GxS.rsidsjustGxSIds' 


# extract GxS PCA subset
arguments=' --memory 43000 --bfile '$believeGxS$' --keep '$believeScratchLoc$'indi_keeplist --extract '$believeRawLoc$'GxS.rsidsjustGxSIds --make-bed --out '$believeScratchLoc$'PCA_subset --allow-extra-chr --allow-no-sex'
$plink $arguments


# apply the liftover
awk 'FNR == NR { 
file1[ $2 ] = $1"\t"$2"\t"$3"\t"$4;next; } # only take the first 4 cols, as the A1/A2 were dummy
FNR <= NR {  
if( $2 in file1) { print file1[$2]"\t"$5"\t"$6 }
else {print $0} }
' $believeScratchLoc$'liftedOver_origA1A2.bim' $believeScratchLoc$'PCA_subset.bim' > $believeScratchLoc$'PCA_subset_liftover'
head $believeScratchLoc$'PCA_subset_liftover'

# overwrite the GxS IDs with RSids
awk 'FNR == NR { 
file1[ $1 ] = $2;next; } 
FNR <= NR {  
if( $2 in file1) { print $1"\t"file1[$2]"\t"$3"\t"$4"\t"$5"\t"$6 }
else {print $0} }
' $believeRawLoc$'GxS.rsids' $believeScratchLoc$'PCA_subset_liftover' > $believeScratchLoc$'PCA_subset_liftover_rsid.bim'
head $believeScratchLoc$'PCA_subset_liftover_rsid.bim'

# overwrite .bim
mv $believeScratchLoc$'PCA_subset.bim' $believeScratchLoc$'PCA_subset.old'
cp $believeScratchLoc$'PCA_subset_liftover_rsid.bim' $believeScratchLoc$'PCA_subset.bim' # 195216


###########################################################################
# II. extract UKBB PCA subset
###########################################################################

# extract UKB PC data
awk '{print $1"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36"\t"$37"\t"$38"\t"$39"\t"$40}' $QCfile > $believeRawLoc$'UKBB_PCA16'
head $believeRawLoc$'UKBB_PCA16'

# sex discordant or low quality genotype ($9: Submitted.Gender        $10: Inferred.Gender, or $18:het.missing.outliers)
awk '{ if ($9 != $10 || $18 == 1) {print $1} }' $QCfile > $SEX_DISC_or_LQ
head $SEX_DISC_or_LQ
wc -l $SEX_DISC_or_LQ # 374

cat $WITHDRAWALS $SEX_DISC_or_LQ  > $believeRawLoc$'excluded_keeprelateds' ## 550, instead of 40227


# get the keeplists with inferred ancestry by using the pop centroids by Prive
mkdir -p $believeRawLoc$'PCA_related/'
arguments='/home/mk907/scripts/R/BelievePCA.R '$believeRawLoc$'UKBB_PCA16 '$believeRawLoc$'excluded_keeprelateds '$believeRawLoc$'PCA_related/PCA_related 33712'
Rscript $arguments
wc -l $believeRawLoc$'PCA_related/PCA_related_PCAIndis' # 45,777


for ((i=1; i<=$numChroms; i++)); do
pgen='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr'$i
arguments=' --pfile '$pgen$' --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --make-bed  --keep '$believeRawLoc$'PCA_related/PCA_related_PCAIndis --extract '$believeRawLoc$'GxS.rsidsjustRSid --out '$believeScratchLoc$'UKBB_PCA_subset_'$i$' --allow-extra-chr --snps-only --bp-space 1'
$plink2 $arguments
#sbatch --mem ${PLINKRAM} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account SHARRISON-SL3-CPU --wrap "$plink2 $arguments"
done

# merge UKBB
plinkFileList=$believeScratchLoc$'plinkFileList'
rm -rf $plinkFileList
for ((i=1; i<=$numChroms; i++)); do
echo $believeScratchLoc$'UKBB_PCA_subset_'$i >> ${plinkFileList}
done # end of chrom loop

arguments=' --memory '$PLINKRAM$' --merge-list '$plinkFileList$' --make-bed --out '$believeRawLoc$'UKBB_PCA_subset_all --allow-extra-chr --allow-no-sex'
$plink $arguments

# try merging UKBB and BELIEVE
rm -rf $plinkFileList
echo $believeRawLoc$'UKBB_PCA_subset_all' >> ${plinkFileList}
echo $believeScratchLoc$'PCA_subset' >> ${plinkFileList}
arguments=' --memory '$PLINKRAM$' --merge-list '$plinkFileList$' --make-bed --out '$believeRawLoc$'UKBB_BELIEVE --allow-extra-chr --allow-no-sex'
$plink $arguments

# this fails because of multi allelic variants, and variants with mismatched positions
# create list of variants that have mismatched positions between UKBB and BELIEVE
awk 'FNR == NR { 
file1[ $2 ] = $4; }
FNR <= NR {  
if(  file1[$2] != $4 ) { print $2"\t"$4"\t"file1[$2] }
 }
' $believeRawLoc$'UKBB_PCA_subset_all.bim' $believeScratchLoc$'PCA_subset.bim' > $believeScratchLoc$'mismatchedSNPs'
head $believeScratchLoc$'mismatchedSNPs'
wc -l $believeScratchLoc$'mismatchedSNPs'  # 384

# add it to list of SNPs with multi allelic vars:
cat /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/raw/UKBB_BELIEVE-merge.missnp  $believeScratchLoc$'mismatchedSNPs' > $believeScratchLoc$'allbadSNPs'

# remove these from both
arguments=' --bfile '$believeRawLoc$'UKBB_PCA_subset_all --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --make-bed  --exclude '$believeScratchLoc$'allbadSNPs --out '$believeScratchLoc$'UKBB_PCA_subset_all_QC --allow-extra-chr --snps-only --bp-space 1'
$plink2 $arguments

arguments=' --bfile '$believeScratchLoc$'PCA_subset --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --make-bed  --exclude '$believeScratchLoc$'allbadSNPs --out '$believeScratchLoc$'PCA_subset_QC --allow-extra-chr --snps-only --bp-space 1'
$plink2 $arguments

# try merge again UKBB and BELIEVE
rm -rf $plinkFileList
echo $believeScratchLoc$'UKBB_PCA_subset_all_QC' >> ${plinkFileList}
echo $believeScratchLoc$'PCA_subset_QC' >> ${plinkFileList}
arguments=' --memory '$PLINKRAM$' --merge-list '$plinkFileList$' --make-bed --out '$believeRawLoc$'UKBB_BELIEVE --allow-extra-chr --allow-no-sex'
$plink $arguments

##################
# Perform QC: 

# create backups of the original files
mv $believeRawLoc$'UKBB_BELIEVE.bim' $believeRawLoc$'UKBB_BELIEVE_noQC.bim'
mv $believeRawLoc$'UKBB_BELIEVE.fam' $believeRawLoc$'UKBB_BELIEVE_noQC.fam'
mv $believeRawLoc$'UKBB_BELIEVE.bed' $believeRawLoc$'UKBB_BELIEVE_noQC.bed'


# 1. remove ambiguous SNPs
remove_ambiguous_alleles_bim $believeRawLoc$'UKBB_BELIEVE_noQC.bim'

arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE_noQC --make-bed --out '$believeScratchLoc$'UKBB_BELIEVE_noambig --extract '$believeRawLoc$'UKBB_BELIEVE_noQC.bim_noAmbiguousAlleles_SNPlist --allow-extra-chr --allow-no-sex'
$plink $arguments
# 179675 variants and 79488 people pass filters and QC.

# 2. exclude QC failed SNPs in the UKBB
arguments=' --memory '$PLINKRAM$' --bfile '$believeScratchLoc$'UKBB_BELIEVE_noambig --make-bed --out '$believeScratchLoc$'UKBB_BELIEVE_UKBQC --exclude '$scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC --allow-extra-chr --allow-no-sex'
$plink $arguments
# 171287 variants and 79488 people pass filters and QC. # if we exclude INFOS



# 3. apply MAF /HWE filters in each cohort separately 

# do NOT do  --hwe 1e-10 midp, as the UKB is a non-homogenous dataset. They did HWE already in the $GENOTYPE_QC file, rationale: (https://biobank.ctsu.ox.ac.uk/crystal/crystal/docs/genotyping_qc.pdf)
# UKB --maf 0.001 --geno 0.05
arguments=' --memory '$PLINKRAM$' --bfile '$believeScratchLoc$'UKBB_BELIEVE_UKBQC --keep '$believeRawLoc$'UKBB_PCA_subset_all.fam --write-snplist --out '$believeScratchLoc$'UKBB_BELIEVE_QC_UKB --maf 0.001 --geno 0.05 --allow-extra-chr --allow-no-sex'
$plink $arguments


arguments=' --memory '$PLINKRAM$' --bfile '$believeScratchLoc$'UKBB_BELIEVE_UKBQC --remove '$believeRawLoc$'UKBB_PCA_subset_all.fam --write-snplist --out '$believeScratchLoc$'UKBB_BELIEVE_QC_BELIEVE --maf 0.001 --geno 0.05 --hwe 1e-10 midp --allow-extra-chr --allow-no-sex'
$plink $arguments
# 156726 variants and 33712 people pass filters and QC.

# intersect 2 lsits
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeScratchLoc$'UKBB_BELIEVE_QC_UKB.snplist' $believeScratchLoc$'UKBB_BELIEVE_QC_BELIEVE.snplist' > $believeScratchLoc$'UKBB_BELIEVE_QCCOMMON.snplist'
head $believeScratchLoc$'UKBB_BELIEVE_QCCOMMON.snplist'
wc -l $believeScratchLoc$'UKBB_BELIEVE_QCCOMMON.snplist' # 127,315


arguments=' --memory '$PLINKRAM$' --bfile '$believeScratchLoc$'UKBB_BELIEVE_UKBQC  --make-bed --out '$believeRawLoc$'UKBB_BELIEVE --extract '$believeScratchLoc$'UKBB_BELIEVE_QCCOMMON.snplist --allow-extra-chr --allow-no-sex'
$plink $arguments




wc -l $believeRawLoc$'UKBB_BELIEVE.bim'
###########################################################################
# III. PCA:
###########################################################################

resLoc=$believeResultsLoc$'UKBB_BELIEVE/'
mkdir -p $resLoc
rm -rf $resLoc$'res.err'
rm -rf $resLoc$'res.out'

arguments='/home/mk907/scripts/R/PCA_ncMCE.R '$believeRawLoc$'UKBB_BELIEVE.bed '$believeResultsLoc$'UKBB_BELIEVE/ '$NCORES_LDPred2$' 20'
sbatch --mem 160000 --cpus-per-task $NCORES_LDPred2 --time 12:0:0 --partition cclake --account danesh-sl3-cpu -e $resLoc$'res.err' -o $resLoc$'res.out' --wrap "Rscript $arguments"

 wc -l $believeRawLoc$'UKBB_BELIEVE.bim' # 127,315
  wc -l $believeRawLoc$'UKBB_BELIEVE.fam' # 79,488
  
# visualise w ethnic group
AddEthnicLabels $believeRawLoc$'PCA_related/PCA_related_PCAIndis' $believeResultsLoc$'UKBB_BELIEVE/_PC_PCA' $believeResultsLoc$'ukbb_believe_ethnic'

# visualise first 4 PCs
echo -e "1\tV2\n2\tV3\n3\tV4\n4\tV5" > $believeResultsLoc$'ethnic_PCs'

arguments='/home/mk907/scripts/R/PCA_visualiser.R '$believeResultsLoc$'UKBB_BELIEVE/_PC_PCA '$believeResultsLoc$'ethnic_PCs '$believeResultsLoc$'ukbb_believe_ethnic '$believeResultsLoc$'UKBB_BELIEVE/ethnic_PCA_ PC'
Rscript $arguments

# this generates the coloured PCA plot that shows where the BELIEVE indies are located relative to the UKB known populations (slide 8)

###########################################################################
# IV. Find subset of indis in UKBB that are genetically matched to BELIEVE:
###########################################################################

mkdir -p $believeResultsLoc$'ancestry_tests/'
# Project excluded indis into the existing PCA space to increase numbers
arguments='/home/mk907/scripts/R/ProjectPCA.R '$believeResultsLoc$'UKBB_BELIEVE/ '$believeRawLoc$'UKBB_BELIEVE.bed '$believeResultsLoc$'ancestry_tests/PCAProj'
Rscript $arguments
# "written  PCs for  8436  new indis


# identify as 'BELIEVE-MATCHED' in the UKBB, as those which are at most the furthest distance away on each PC
AddEthnicLabels $believeRawLoc$'PCA_related/PCA_related_PCAIndis' $believeResultsLoc$'ancestry_tests/PCAProj' $believeResultsLoc$'ancestry_tests/Related_ukbb_believe_ethnic'

arguments='/home/mk907/scripts/R/BelieveMatcher.R '$believeResultsLoc$'ancestry_tests/PCAProj '$believeResultsLoc$'ancestry_tests/Related_ukbb_believe_ethnic '$believeResultsLoc$'ancestry_tests/UKBB_related_matched'
Rscript $arguments
#India Other
#   24    69
#India Other
#   60   109

arguments='/home/mk907/scripts/R/BelieveMatcher.R '$believeResultsLoc$'UKBB_BELIEVE/_PC_PCA '$believeResultsLoc$'ukbb_believe_ethnic '$believeResultsLoc$'ancestry_tests/UKBB_matched'
Rscript $arguments
#India Other
# 1626  1099
#India Other
# 5122  2431


# these produce 2 lists of indis: (UKBB_matched_matched and UKBB_matched_matchedPC4), which have people in the UKBB that are matched either precisely (up to 20PCs) or more relaxed (up to 4PCs) to BELIEVE

# produce combined lists that have both the related and unrelated indis:
cat $believeResultsLoc$'ancestry_tests/UKBB_matched_matched' $believeResultsLoc$'ancestry_tests/UKBB_related_matched_matched' > $believeResultsLoc$'ancestry_tests/UKBB_all_matched'
cat $believeResultsLoc$'ancestry_tests/UKBB_matched_matchedPC4' $believeResultsLoc$'ancestry_tests/UKBB_related_matched_matchedPC4' > $believeResultsLoc$'ancestry_tests/UKBB_all_matchedPC4'
cat $believeResultsLoc$'UKBB_BELIEVE/_PC_PCA' $believeResultsLoc$'ancestry_tests/PCAProj' > $believeResultsLoc$'ancestry_tests/ALL_PCA'

wc -l $believeResultsLoc$'ancestry_tests/UKBB_all_matched' # 2,818 
wc -l $believeResultsLoc$'ancestry_tests/UKBB_all_matchedPC4' # 7,722 

###############

# Precise match (up to PC20)
awk '{print $1"\t"$2}' $believeScratchLoc$'PCA_subset_QC.fam' > $believeScratchLoc$'BELIEVE_ONLY'
awk '{print $1"\t"$1}' $believeResultsLoc$'ancestry_tests/UKBB_all_matched' > $believeScratchLoc$'MATCHED_ALL_ONLY' # 2,818 
cat $believeScratchLoc$'BELIEVE_ONLY' $believeScratchLoc$'MATCHED_ALL_ONLY' >  $believeScratchLoc$'BELIEVE_MATCHED_ALL' # this is all the BELIEVE indis plus  the UKBB that are ancestrally matched

# get master keeplist: a fam file
awk 'FNR == NR { file1[ $2 ] = $2;next; } 
FNR <= NR {  if( $2 in file1 ) {  print $0 } }
' $believeScratchLoc$'BELIEVE_MATCHED_ALL' $believeRawLoc$'UKBB_BELIEVE.fam' > $believeRawLoc$'BELIEVE_MATCHED_ALL.fam'
head $believeRawLoc$'BELIEVE_MATCHED_ALL.fam'
wc -l $believeRawLoc$'BELIEVE_MATCHED_ALL.fam'

arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE --keep '$believeRawLoc$'BELIEVE_MATCHED_ALL.fam --make-bed --out '$believeRawLoc$'BELIEVE_MATCHED_ALL --allow-extra-chr --allow-no-sex'
$plink $arguments

# perform new pca to see how this worked
resLoc=$believeResultsLoc$'UKBB_BELIEVE_MATCHED/'
mkdir -p $resLoc
rm -rf $resLoc$'res.err'
rm -rf $resLoc$'res.out'

arguments='/home/mk907/scripts/R/PCA_ncMCE.R '$believeRawLoc$'BELIEVE_MATCHED_ALL.bed '$resLoc$' '$NCORES_LDPred2$' 20'
sbatch --mem 160000 --cpus-per-task $NCORES_LDPred2 --time 12:0:0 --partition cclake --account danesh-sl3-cpu -e $resLoc$'res.err' -o $resLoc$'res.out' --wrap "Rscript $arguments"

###############
# Relaxed match (up to PC4)
awk '{print $1"\t"$1}' $believeResultsLoc$'ancestry_tests/UKBB_all_matchedPC4' > $believeScratchLoc$'MATCHED_ALL_ONLYPC4'
cat $believeScratchLoc$'BELIEVE_ONLY' $believeScratchLoc$'MATCHED_ALL_ONLYPC4' >  $believeScratchLoc$'BELIEVE_MATCHED_ALLPC4'

# get master keeplist: a fam file
awk 'FNR == NR { file1[ $2 ] = $2;next; } 
FNR <= NR {  if( $2 in file1 ) {  print $0 } }
' $believeScratchLoc$'BELIEVE_MATCHED_ALLPC4' $believeRawLoc$'UKBB_BELIEVE.fam' > $believeRawLoc$'BELIEVE_MATCHED_ALLPC4.fam'
head $believeRawLoc$'BELIEVE_MATCHED_ALLPC4.fam'
wc -l $believeRawLoc$'BELIEVE_MATCHED_ALLPC4.fam'

arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE --keep '$believeRawLoc$'BELIEVE_MATCHED_ALLPC4.fam --make-bed --out '$believeRawLoc$'BELIEVE_MATCHED_ALLPC4 --allow-extra-chr --allow-no-sex'
$plink $arguments

# perform new pca
resLoc=$believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/'
mkdir -p $resLoc
rm -rf $resLoc$'res.err'
rm -rf $resLoc$'res.out'

arguments='/home/mk907/scripts/R/PCA_ncMCE.R '$believeRawLoc$'BELIEVE_MATCHED_ALLPC4.bed '$resLoc$' '$NCORES_LDPred2$' 20'
sbatch --mem 160000 --cpus-per-task $NCORES_LDPred2 --time 12:0:0 --partition cclake --account danesh-sl3-cpu -e $resLoc$'res.err' -o $resLoc$'res.out' --wrap "Rscript $arguments"



# visualise results of above

###############

# Precise match
# prepare groups: black BELIEVE, otherwise blue
awk 'FNR == NR { file1[ $1 ] = $1;next; } 
FNR <= NR {  
if(FNR == 1) {print "IID\tgroup\tcolor"}
if( $1 in file1 ) { print $1"\tUKBB_MATCH\tblue"} 
else {print $1"\tBELIEVE\tblack"} }
' $believeResultsLoc$'ancestry_tests/UKBB_all_matched' $believeResultsLoc$'UKBB_BELIEVE_MATCHED/_PC_PCA' > $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels2'
head $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels2'  
wc -l $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels2'

arguments='/home/mk907/scripts/R/PCA_visualiser.R '$believeResultsLoc$'UKBB_BELIEVE_MATCHED/_PC_PCA '$believeResultsLoc$'ethnic_PCs '$believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels2 '$believeResultsLoc$'UKBB_BELIEVE_MATCHED/MATCHED PC'
Rscript $arguments
# produces the UKB-matched + BEIEVE indis (in black) (Slide 9)

###############
# Relaxed up to PC4 matches
# prepare groups: black BELIEVE, otherwise blue
awk 'FNR == NR { file1[ $1 ] = $1;next; } 
FNR <= NR {  
if(FNR == 1) {print "IID\tgroup\tcolor"}
if( $1 in file1 ) { print $1"\tUKBB_MATCH\tblue"} 
else {print $1"\tBELIEVE\tblack"} }
' $believeResultsLoc$'ancestry_tests/UKBB_all_matchedPC4' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/_PC_PCA' > $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labels'
head $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labels'
wc -l $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labels'

arguments='/home/mk907/scripts/R/PCA_visualiser.R '$believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/_PC_PCA '$believeResultsLoc$'ethnic_PCs '$believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labels '$believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/MATCHED PC'
Rscript $arguments
# produces the UKB-matched + BEIEVE indis (in black) (Slide 9)
############################


# prepare genotype files: subset to GBP, UKBB_Matched, and a quarter of BELIEVE, so that total is around ~30K

# get white british
awk '{if ($2 == "UK") {print $1"\t"$1}}' $believeResultsLoc$'ukbb_believe_ethnic' > $believeResultsLoc$'ukbb_british_only'

# get all BELIEVE_Matched, and related matched and EUR
cat $believeScratchLoc$'BELIEVE_MATCHED_ALL' $believeResultsLoc$'ukbb_british_only' > $believeScratchLoc$'BELIEVE_MATCHED_EUR'
cat $believeScratchLoc$'BELIEVE_MATCHED_ALLPC4' $believeResultsLoc$'ukbb_british_only' > $believeScratchLoc$'BELIEVE_MATCHED_EURPC4'

# create subset genotype files
arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE --keep '$believeScratchLoc$'BELIEVE_MATCHED_EUR --make-bed --out '$believeRawLoc$'BELIEVE_MATCHED_EUR --allow-extra-chr --allow-no-sex'
$plink $arguments

arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE --keep '$believeScratchLoc$'BELIEVE_MATCHED_EURPC4 --make-bed --out '$believeRawLoc$'BELIEVE_MATCHED_EURPC4 --allow-extra-chr --allow-no-sex'
$plink $arguments

###########################################################################
# V. Prepare keeplists/labels:
###########################################################################
awk -v RS='\r?\n' 'FNR == NR { 
stem="CAMBRIDGE-BELIEVE_"
id=substr($3,2); # remove the first character, as in the ID mapping file it has a 0 prefix, where as in the fam files it does not
famID=stem id"_"id
file1[$2] =famID;
next; 
} FNR <= NR { 
if ( FNR > 1 && $2 in file1 && $46 != "") {  print file1[$2]"\t"file1[$2]"\t"$46  }
 }' OFS="\t" FS="," $phenoMap  FS="," $phenoRaw > $believeDataLoc$'ethnic_all.phe'
head $believeDataLoc$'ethnic_all.phe'
wc -l $believeDataLoc$'ethnic_all.phe' # 52028
# check how many cases/controls
awk '{count[$3]++} END {for (word in count) print word, count[word]}' $believeDataLoc$'ethnic_all.phe'

# ethnic groups,  1 = Bangali; 2 = Bihari; 3 = Other
awk '{
color = "gray"
group = "Other"
if (FNR == 1) print "IID\tgroup\tcolor"
if ($3 == "1") {color = "green"; group = "Bangali"; }
else if ($3 == "2") {color = "blue"; group = "Bihari"; }
print $2"\t"group"\t"color
}' $believeDataLoc$'ethnic_all.phe' > $believeResultsLoc$'ethnic_groups'
head $believeResultsLoc$'ethnic_groups'


# prepare pop files with header:
#IID     pop     superpop
awk '{print $1"\tGBP\tEUR"}' $believeResultsLoc$'ukbb_british_only' > $believeResultsLoc$'ukbb_british_only_list'

# precise
awk '{print $1"\tUK_SEA\tSEA"}' $believeResultsLoc$'ancestry_tests/UKBB_matched_matched' > $believeResultsLoc$'ancestry_tests/UKBB_matched_matched_list'
head $believeResultsLoc$'ancestry_tests/UKBB_matched_matched_list'
cat $believeResultsLoc$'ukbb_british_only_list' $believeResultsLoc$'ancestry_tests/UKBB_matched_matched_list' > $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_list'
head $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_list'
wc -l $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_list'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_list'
# GBP 6297
# UK_SEA 2725

# relaxed
awk '{print $1"\tUK_SEA\tSEA"}' $believeResultsLoc$'ancestry_tests/UKBB_matched_matchedPC4' > $believeResultsLoc$'ancestry_tests/UKBB_matched_matched_listPC4'
head $believeResultsLoc$'ancestry_tests/UKBB_matched_matched_list'
cat $believeResultsLoc$'ukbb_british_only_list' $believeResultsLoc$'ancestry_tests/UKBB_matched_matched_listPC4' > $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_listPC4'
head $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_listPC4'
wc -l $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_listPC4'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_listPC4'
#GBP 6297
#UK_SEA 7553



# Problem:
# the issue is that there were a few people that got excluded back when I've done QC from the "non BELIEVE" list, which list we then used to assign "BELIEVE" labels, which then would mean
# that a few people who were UKBB got incorrectly assigned to be BELIEVE
# so we prepare the final label lists by using the "$believeScratchLoc$'BELIEVE_ONLY'" to mmake sure all BELIEVE indis are correctly labelled

# precise
awk 'FNR == NR { file1[ $1 ] = "BELIEVE\tSEA";next; } 
FNR <= NR {  
if(FNR == 1) {print "IID\tpop\tsuperpop"}
if( $1 in file1 ) { print $1"\t"file1[$1]} 
else {print $1"\tQCFAIL\tQCFAIL"} }
' $believeScratchLoc$'BELIEVE_ONLY' $believeRawLoc$'BELIEVE_MATCHED_EUR.fam' > $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels_v1'
head $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels_v1'
wc -l $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels_v1' # 41611
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels_v1'
#QCFAIL 9115
#pop 1
#BELIEVE 33712

awk 'FNR == NR { file1[ $1 ] = $2"\t"$3;next; } 
FNR <= NR {  
if(FNR == 1) {print $0}
else if( $1 in file1 ) { print $1"\t"file1[$1]} 
else {print $0} }
' $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_list' $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels_v1' > $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels'
head $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels'
wc -l $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels' # 41611
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels'
# QCFAIL 93
# pop 1
# GBP 6297
# BELIEVE 33712
# UK_SEA 2725


# relaxed
awk 'FNR == NR { file1[ $1 ] = "BELIEVE\tSEA";next; } 
FNR <= NR {  
if(FNR == 1) {print "IID\tpop\tsuperpop"}
if( $1 in file1 ) { print $1"\t"file1[$1]} 
else {print $1"\tQCFAIL\tQCFAIL"} }
' $believeScratchLoc$'BELIEVE_ONLY' $believeRawLoc$'BELIEVE_MATCHED_EURPC4.fam' > $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_v1'

awk 'FNR == NR { file1[ $1 ] = $2"\t"$3;next; } 
FNR <= NR {  
if(FNR == 1) {print $0}
else if( $1 in file1 ) { print $1"\t"file1[$1]} 
else {print $0} }
' $believeResultsLoc$'ancestry_tests/UKBB_matched_GBP_listPC4' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_v1' > $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4'
head $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4'
wc -l $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4' # 41611
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4'
#QCFAIL 169
#pop 1
#GBP 6297
#BELIEVE 33712
#UK_SEA 7553


###########################################################################
# VI. Prepare phenotypes:
###########################################################################

######################
# Prepare covariates
# the mapping file, the 3rd col is the one we want
#HseIdentifier,Identifier,GeneticID
#1000007,100000701,0334272530

# new variable "Age" col77, Gender is now col 78
# create overall covariate matrix of: age age^2 and sex
# the "HseIdentifier" and	"Identifier" are different between releases, so we need to match via the GeneticID,
# which is the $519 col in the new phenotype files

# the pheno files are dirty, they have "," within fields, so awk cannot process them. Have to use R:
R
# without specifying the colClasses, it would read GeneticID as number, and remove the leading zeroes, which would mismatch the .fam files # , colClasses=c("GeneticID"="character")
phenoMap=read.csv("/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/BELIEVEdata_P5031_20220310.csv" )
phenoMap$id = paste("CAMBRIDGE-BELIEVE_" , as.character(phenoMap$GeneticID),"_" ,as.character(phenoMap$GeneticID), sep="" )
phenoMap_exist = phenoMap[which(phenoMap$Age != "NA" & phenoMap$Age != ""& is.na(phenoMap$Age) == F & phenoMap$Gender != "NA" & phenoMap$Gender != ""& is.na(phenoMap$Gender) == F),]
write.table(cbind(phenoMap_exist$id,phenoMap_exist$id, phenoMap$Age, phenoMap$Age^2,  phenoMap$Gender), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/covars_new.phe", quote = F, col.names = F, row.names = F);
quit()

# UKBB: get sex + age data for ukbb
R
library(readstata13)
library('data.table')
ceu <- read.dta13("/rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/20210901/STATA/analysis.dta")
setDT(ceu)
vars = attributes(ceu)$var.labels
names(vars) = names(ceu)
myData=ceu[,.(eid=idno,ages, ages^2, sex)]
write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/scratch/ukbb_covs_raw2", sep="\t", row.names = F, col.names = T, quote = F)
quit()


#######################
# BELIEVE: extract height, BMI, Edu Years and high cholesterol
# use R, as awk has hard time dealing with dirty csv files that come with the BELIEVE data
R
phenoMap=read.csv("/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/archive/BELIEVEdata_GeneticIDMapping_P5031_20210826.csv")
# as the above reads the IDs as numbers, it removes the leading zeroes
phenoMap$id = paste("CAMBRIDGE-BELIEVE_" , as.character(phenoMap$GeneticID),"_" ,as.character(phenoMap$GeneticID), sep="" )
pheno_all=read.csv("/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/archive/BELIEVEdata_P5031_20210826.csv")

# extract height
pheno=pheno_all[,c(2,82)]
pheno =  na.omit(pheno)
pheno_ID = merge(pheno, phenoMap, by="Identifier")
pheno_ID = pheno_ID[,c(5,5,2)]


# exclude those below 100, and those over 200 as outliers (the max height was 192cm)
min(pheno_ID$Height) # 0
max(pheno_ID$Height) # 192
outliers = which( pheno_ID$Height < 100 | pheno_ID$Height > 200) # there were 11 of these, with several people claiming 15 cm heights
length(outliers)
pheno_ID$Height[outliers]
pheno_ID = pheno_ID[-outliers,]
nrow(pheno_ID)
write.table(pheno_ID,"/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/height_all.phe", sep="\t",col.names = F, row.names = F, quote = F)

# extract BMI
pheno=pheno_all[,c(2,443)]
pheno =  na.omit(pheno)
pheno_ID = merge(pheno, phenoMap, by="Identifier")
pheno_ID = pheno_ID[,c(5,5,2)]
min(pheno_ID$BMI) #  12.1 # this is very low but still realistic
max(pheno_ID$BMI) # 66.84
# no exclusions

write.table(pheno_ID,"/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/bmi_all.phe", sep="\t",col.names = F, row.names = F, quote = F)

# extract Edu years (FullTimeEdYears)
pheno=pheno_all[,c(2,194)]
pheno =  na.omit(pheno)
pheno_ID = merge(pheno, phenoMap, by="Identifier")
pheno_ID = pheno_ID[,c(5,5,2)]

min(pheno_ID$FullTimeEdYears) # 1 # in Bangladesh some people may never go to school
max(pheno_ID$FullTimeEdYears) # 60

outliers = which(pheno_ID$FullTimeEdYears > 30) # 6670 indis with 1 year of  education, is this realistic?, there is one individual with 60 years of education ( pheno_ID$FullTimeEdYears < 2 | )
#outliers = which( pheno_ID$FullTimeEdYears < 2 )
#outliers = which( pheno_ID$FullTimeEdYears > 20 )
#pheno_ID$FullTimeEdYears[outliers]
#length(outliers)
pheno_ID = pheno_ID[-outliers,]
nrow(pheno_ID)
write.table(pheno_ID,"/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/edu_all.phe", sep="\t",col.names = F, row.names = F, quote = F)
quit()


# High cholesterol was extracted differently (due to having done this work earlier)


# fix DOS line endings: https://stackoverflow.com/questions/29612237/awk-script-doesnt-print-last-column-of-csv-file
awk -v RS='\r?\n' 'FNR == NR { 
stem="CAMBRIDGE-BELIEVE_"
id=substr($3,2); # remove the first character, as in the ID mapping file it has a 0 prefix, where as in the fam files it does not
famID=stem""id"_"id
file1[$2] =famID;
next; 
} FNR <= NR { 
if ( FNR > 1 && $2 in file1 && $94 != "") {  print file1[$2]"\t"file1[$2]"\t"$94  }
 }' OFS="\t" FS="," $phenoMap  FS="," $phenoRaw > $believeDataLoc$'highCholesterol_all.phe'
head $believeDataLoc$'highCholesterol_all.phe'
wc -l $believeDataLoc$'highCholesterol_all.phe'

# check how many cases/controls
awk '{count[$3]++} END {for (word in count) print word, count[word]}' $believeDataLoc$'highCholesterol_all.phe'
#1 6701
#2 45326
# 6701/ (6701 + 45326) =  0.1287985, so 13% are cases

# need to reverse cholesterol to be 1=control, 2= to be case to be consistent with UKBB
cp $believeDataLoc$'highCholesterol_all.phe' $believeDataLoc$'highCholesterol_all.phe.old'
awk '{if ($3 == "1") {print $1"\t"$2"\t2"} else {print $1"\t"$2"\t1"}}' $believeDataLoc$'highCholesterol_all.phe.old' > $believeDataLoc$'highCholesterol_all.phe'
awk '{count[$3]++} END {for (word in count) print word, count[word]}' $believeDataLoc$'highCholesterol_all.phe'



###############
# UKBB: extract same phenos
R
library(readstata13)
library('data.table')

ceu <- read.dta13("/rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/20210901/STATA/analysis.dta")


setDT(ceu)
vars = attributes(ceu)$var.labels
names(vars) = names(ceu)

# height
myData=ceu[,.(eid=idno,eid=idno, Height=ht)]
myData =  na.omit(myData)

nrow(myData) # 499921
outliers = which( myData$Height < 100 | myData$Height > 200) # there were 11 of these, with several people claiming 15 cm heights
length(outliers)
myData$Height[outliers]
myData = myData[-outliers,]
nrow(myData) # 499848


write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_height.phe", row.names = F, col.names = F, quote = F, sep="\t")

# BMI
myData=ceu[,.(eid=idno,eid=idno, BMI=bmi)]
myData =  na.omit(myData)
min(myData$BMI) #  12.12
max(myData$BMI) #  74.68
write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_bmi.phe", row.names = F, col.names = F, quote = F, sep="\t")

# Edu
myData=ceu[,.(eid=idno,eid=idno, FullTimeEdYears=edyrs)]
myData =  na.omit(myData)
min(myData$FullTimeEdYears) # 0
max(myData$FullTimeEdYears) # 30
write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_edu.phe", row.names = F, col.names = F, quote = F, sep="\t")

# high cholesterol binary, https://www.mayoclinic.org/diseases-conditions/high-blood-cholesterol/diagnosis-treatment/drc-20350806
# tchol Total cholesterol (mmol/l)
myData=ceu[,.(eid=idno,eid=idno, Totalcholesterol=tchol)]
myData =  na.omit(myData)

myData$HighCholesterol = 1 # we flipped it to be the usual 1=control, 2= case for BELIEVE now too
highCholesterolIndices = which(myData$Totalcholesterol > 6.2) 
length(highCholesterolIndices) # 147630 out of 469548
#myData$HighCholesterol[highCholesterolIndices] = 1
myData$HighCholesterol[highCholesterolIndices] = 2
myData = myData[,-3]
str(myData)
write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_cholesterol.phe", row.names = F, col.names = F, quote = F, sep="\t")

quit()

################



# generate phenotypes for the 1500 subset. This is NEVER directly used, as we only ever use the 1,500 subset to run a simulation
# to show how much more accurate it is, relative the the imputed version
# BLV 
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'highCholesterol_all.phe' $believeScratchLoc$'cholest_cov2' 1
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'height_all.phe' $believeScratchLoc$'height_cov2' 0
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'bmi_all.phe' $believeScratchLoc$'bmi_cov2' 0
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'edu_all.phe' $believeScratchLoc$'edu_cov2' 0

# UKBB
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_cholesterol.phe' $believeScratchLoc$'ukbb_cholesterol_cov2' 1
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_height.phe' $believeScratchLoc$'ukbb_height_cov2' 0
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_bmi.phe' $believeScratchLoc$'ukbb_bmi_cov2' 0
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_edu.phe' $believeScratchLoc$'ukbb_edu_cov2' 0

# cat ukbb and believe
cat $believeScratchLoc$'cholest_cov2_res_id' $believeScratchLoc$'ukbb_cholesterol_cov2_res_id' > $believeDataLoc$'highCholesterol_all2'
cat $believeScratchLoc$'height_cov2_res_id' $believeScratchLoc$'ukbb_height_cov2_res_id' > $believeDataLoc$'height_all2'
cat $believeScratchLoc$'bmi_cov2_res_id' $believeScratchLoc$'ukbb_bmi_cov2_res_id' > $believeDataLoc$'bmi_all2'
cat $believeScratchLoc$'edu_cov2_res_id' $believeScratchLoc$'ukbb_edu_cov2_res_id' > $believeDataLoc$'edu_all2'



# Subset plink to only have the precise indis of each group:
# This is only used for the simulations involving  the precise-set of indis
subsetPhen $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels' $believeDataLoc$'bmi_all2' $believeScratchLoc$'1'
subsetPhen $believeScratchLoc$'1' $believeDataLoc$'height_all2' $believeScratchLoc$'2'
subsetPhen $believeScratchLoc$'2' $believeDataLoc$'highCholesterol_all2'  $believeScratchLoc$'3'

awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'1'

# we loose half of our data with _edu...
# just work with bmi/height/cholesterol, and edu will take a hit at the time when I fit the model

subsetPhen $believeScratchLoc$'3' $believeDataLoc$'edu_all' $believeScratchLoc$'has_edu'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'has_edu'
#QCFAIL 37
#pop 1
#GBP 3959
#BELIEVE 21392
#UK_SEA 1245



awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'3'
#QCFAIL 85
#pop 1
#GBP 5966
#BELIEVE 25112
#UK_SEA 2564


# this is for bmi,height, cholesterol
cp $believeScratchLoc$'3' $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels_hasphe'


awk '{if ($2 == "UK_SEA") print $1"\t"$1"\t"$2 }' $believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels_hasphe' > $believeScratchLoc$'UK_SEA_keep'
head $believeScratchLoc$'UK_SEA_keep'
wc -l $believeScratchLoc$'UK_SEA_keep'
awk '{count[$3]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'UK_SEA_keep'
# UK_SEA 2564

# problem is that we want to preferentially keep those that have eduY, even among the GBP and the BELIEVE indis
# so while we keep everyone from UK_SEA, we only work off the smaller $believeScratchLoc$'has_edu' list, that has eduY for the others

# pick numHasEdu indis from BELIEVE  
awk '{if ($2 == "BELIEVE") print $1"\t"$1"\t"$2 }' $believeScratchLoc$'has_edu' > $believeScratchLoc$'BELIEVE_all'
wc -l $believeScratchLoc$'BELIEVE_all' # 21392
numToP=$(wc -l < $believeScratchLoc$'BELIEVE_all')
echo $numToP

# reproducible random number gen
arguments='/home/mk907/scripts/R/randomNums.R '$numToP$' 42 '$believeScratchLoc$'randomNums'
Rscript $arguments


# pick the first 2564
awk '{if (FNR <= 2564) print $0}' $believeScratchLoc$'randomNums' > $believeScratchLoc$'randomNums_indis'
head $believeScratchLoc$'randomNums_indis'

# pick the actual sample ids 
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $0 } }
' $believeScratchLoc$'randomNums_indis' $believeScratchLoc$'BELIEVE_all' > $believeScratchLoc$'BELIEVE_keep' 
head $believeScratchLoc$'BELIEVE_keep' 
wc -l $believeScratchLoc$'BELIEVE_keep' 


# pick  white british too
awk '{if ($2 == "GBP") print $1"\t"$1"\t"$2 }' $believeScratchLoc$'has_edu' > $believeScratchLoc$'GBP_all'
wc -l $believeScratchLoc$'GBP_all'  # 3959
numToP=$(wc -l < $believeScratchLoc$'GBP_all')
echo $numToP

# reproducible random number gen
arguments='/home/mk907/scripts/R/randomNums.R '$numToP$' 42 '$believeScratchLoc$'randomNums'
Rscript $arguments

# pick the first 2564
awk '{if (FNR <= 2564) print $0}' $believeScratchLoc$'randomNums' > $believeScratchLoc$'randomNums_indis'
head $believeScratchLoc$'randomNums_indis'

# pick the actual sample ids 
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $0 } }
' $believeScratchLoc$'randomNums_indis' $believeScratchLoc$'GBP_all' > $believeScratchLoc$'GBP_keep' 
head $believeScratchLoc$'GBP_keep' 
wc -l $believeScratchLoc$'GBP_keep' 

cat  $believeScratchLoc$'BELIEVE_keep' $believeScratchLoc$'UK_SEA_keep' $believeScratchLoc$'GBP_keep' > $believeScratchLoc$'2564_keep'
wc -l $believeScratchLoc$'2564_keep' # 7692

# extract the genotypes
arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'BELIEVE_MATCHED_EUR --keep '$believeScratchLoc$'2564_keep --make-bed --out '$believeRawLoc$'2564_ALL --allow-extra-chr --allow-no-sex'
$plink $arguments


# create a labels file, with matching order to the .fam
awk 'FNR == NR { test[ $1 ] = $3; next; } FNR <= NR { if( $2 in test ) {print $2"\t"test[$2] } }
' $believeScratchLoc$'2564_keep' $believeRawLoc$'2564_ALL.fam' > $believeRawLoc$'2564_ALL_labels'
head $believeRawLoc$'2564_ALL_labels'
wc -l $believeRawLoc$'2564_ALL_labels'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeRawLoc$'2564_ALL_labels'


awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if( $2 in test ) {print $2"\t"$2"\t"test[$2] } }
' $believeDataLoc$'height_all2' $believeRawLoc$'2564_ALL.fam' > $believeDataLoc$'height_2564'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $2 in test ) {print $0"\t"test[$2] } }
' $believeRawLoc$'2564_ALL_labels'  $believeDataLoc$'height_2564' > $believeDataLoc$'height_2564_pop'

# generate the final real phenotype files that match the above (no z-scoring)
wc -l $believeRawLoc$'2564_ALL.fam' # 7692

###########################################################################
# VII. Simulation of precisely-matched indis 
###########################################################################

resultLoc=$believeResultsLoc$'UKBB_BELIEVE_MATCHED/BLV/'
mkdir -p $resultLoc

# simulation: 2 pops 
# run it manually (this is faster for some reason than submitting to the cluster...
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'2564_ALL '$believeDataLoc$'height_2564_pop '$resultLoc$'Sims_res_v2 /home/mk907/scripts/R/pheDecomp_functs_lowmem.R 0 1 20 0 BELIEVE UK_SEA'
Rscript $arguments


# move results to somewhere accessible 
mkdir -p $resultLoc$'finals/'
cp $resultLoc$'Sims_res_v2_VG.txt' $resultLoc$'finals/Sims_res_VG.txt'
cp $resultLoc$'Sims_res_v2_VG_Env.txt' $resultLoc$'finals/Sims_res_v2_VG_Env.txt'
cp $resultLoc$'Sims_res_v2_VE_ENV_GXE.txt' $resultLoc$'finals/Sims_res_VE_ENV_GXE.txt'
cp $resultLoc$'Sims_res_v2_VG_sub.txt' $resultLoc$'finals/Sims_res_VG_sub.txt'
cp $resultLoc$'Sims_res_v2_VG_Env_sub.txt' $resultLoc$'finals/Sims_res_VG_Env_sub.txt'
cp $resultLoc$'Sims_res_v2_VE_ENV_GXE_sub.txt' $resultLoc$'finals/Sims_res_VE_ENV_GXE_sub.txt'
cp $resultLoc$'Sims_res_v2_GV.txt' $resultLoc$'finals/Sims_res_GV.txt'

# run plotting 
cp $resultLoc$'Sims_res_v2_VG.txt' $resultLoc$'finals/Sims_res_VG.txt'
#cp $resultLoc$'Sims_res_v2_VG_Env.txt' $resultLoc$'finals/Sims_res_VG_Env.txt'
#cp $resultLoc$'Sims_res_v2_VE_ENV_GXE.txt' $resultLoc$'finals/Sims_res_VE_ENV_GXE.txt'

cp $resultLoc$'Sims_res_v2_VG_sub.txt' $resultLoc$'finals/Sims_res_VG_sub.txt'
cp $resultLoc$'Sims_res_v2_VG_Env_sub.txt' $resultLoc$'finals/Sims_res_VG_Env_sub.txt'
cp $resultLoc$'Sims_res_v2_VE_ENV_GXE_sub.txt' $resultLoc$'finals/Sims_res_VG_ENV_GXE_sub.txt'

cp $resultLoc$'Sims_res_v2_GV.txt' $resultLoc$'finals/Sims_res_GV.txt'

# I mixed up the VG_ENV and VG_ENV_GXE models
cp $resultLoc$'Sims_res_v2_VG_Env.txt' $resultLoc$'finals/Sims_res_VG_ENV_GXE.txt'
cp $resultLoc$'Sims_res_v2_VE_ENV_GXE.txt' $resultLoc$'finals/Sims_res_VG_Env.txt'

# truth scenario VG Env + GxE
arguments='/home/mk907/scripts/R/pheDecomp_plotter.R '$resultLoc$'finals/Sims_res_VG_ENV_GXE.txt '$resultLoc$'finals/Sims_res_VG_ENV_GXE_sub.txt '$resultLoc$'finals/Sims_res_GV.txt '$resultLoc$'finals/Sims_VG_ENV_GXE 0.33 1 0.33 0.0016'
Rscript $arguments


# truth scenario VG Env + GxE
arguments='/home/mk907/scripts/R/pheDecomp_plotter.R '$resultLoc$'finals/Sims_res_VG_Env.txt '$resultLoc$'finals/Sims_res_VG_Env_sub.txt '$resultLoc$'finals/Sims_res_GV.txt '$resultLoc$'finals/Sims_VG_ENV 0.33 1 0 0.0016'
Rscript $arguments


# truth scenario VG
arguments='/home/mk907/scripts/R/pheDecomp_plotter.R '$resultLoc$'finals/Sims_res_VG.txt '$resultLoc$'finals/Sims_res_VG_sub.txt '$resultLoc$'finals/Sims_res_GV.txt '$resultLoc$'finals/Sims_VG 0.33 0 0 0.0016'
Rscript $arguments


###########################################################################
# VIII. TOPMed imputation: expand the 200K directly genotyped to include 800K hapmap3
###########################################################################
# imputation description: /rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/genotype/imputed/imputation.v0.1.1.README.txt
genotypeFolder="/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/genotype/imputed/"
wc -l $genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.GT.pvar'  # 307,949,575
wc -l $genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.HDS.pvar' # 307,949,575


# convert hapmap3 coords into Build38, this will be easier to subset than the other way around
head $UKBB_PLINK1$'45K_IMPUTED.bim'
wc -l $UKBB_PLINK1$'45K_IMPUTED.bim' # 13,048,567 

# create a dummy plink file with 2 indis, as we will be transcoding to text plink files
head -n 2 $UKBB_PLINK1$'45K_IMPUTED.fam' > $believeScratchLoc$'UKBB_45K_indi_keeplist_2'
# rebuild plink file to avoid chromsome-miss-order problem
arguments=' --memory 43000 --bfile '$UKBB_PLINK1$'45K_IMPUTED --keep '$believeScratchLoc$'UKBB_45K_indi_keeplist_2 --make-bed --out '$believeScratchLoc$'UKBB_45K_liftOverDummy.sort --allow-extra-chr --allow-no-sex'
$plink $arguments


# space to tab to generate bed files for liftOver from hg19 to hg38
$plink --bfile $believeScratchLoc$'UKBB_45K_liftOverDummy.sort' --recode tab --out $believeScratchLoc$'UKBB_45K_liftOverDummy.sort.tab'

cd $believeRawLoc$'liftover/'
python2 liftOverPlink.py -m $believeScratchLoc$'UKBB_45K_liftOverDummy.sort.tab.map' -p $believeScratchLoc$'UKBB_45K_liftOverDummy.sort.tab.ped' -o $believeScratchLoc$'UKBB_45K_liftOverDummy.hg38' -c hg19ToHg38.over.chain.gz -e ./liftOver

# convert back to binary bim
arguments=' --memory 43000 --file '$believeScratchLoc$'UKBB_45K_liftOverDummy.hg38 --make-bed --out '$believeScratchLoc$'UKBB_45K_liftedOver --allow-extra-chr --allow-no-sex'
$plink $arguments


# the above may not keep the A1/A2, as we only kept 2 indis...
# we find the original A1/A2 and match them via the SNPid
awk 'FNR == NR { 
file1[ $2 ] = $5"\t"$6;  next; }
FNR <= NR {  
if( $2 in file1 ) {print $1"\t"$2"\t"$3"\t"$4"\t"file1[$2]} }
' $UKBB_PLINK1$'45K_IMPUTED.bim' $believeScratchLoc$'UKBB_45K_liftedOver.bim' > $believeScratchLoc$'UKBB_45K_liftedOver_origA1A2.bim'
head $believeScratchLoc$'UKBB_45K_liftedOver_origA1A2.bim'
wc -l $believeScratchLoc$'UKBB_45K_liftedOver_origA1A2.bim' # 13,044,615

# create keeplist for the pvar files
awk 'FNR == NR { 
file1[ $1":"$4 ] = $1":"$4; A1[ $1":"$4 ] = $5; A2[ $1":"$4 ] = $6;  next; }
FNR <= NR {  
lookup=$1":"$2
if( lookup in file1 && ($4 == A1[lookup] && $5 == A2[lookup] ||  $4 == A2[lookup] && $5 == A1[lookup]) ) {print $0 } }
' $believeScratchLoc$'UKBB_45K_liftedOver_origA1A2.bim' $genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.GT.pvar' > $believeScratchLoc$'topmed_'
head $believeScratchLoc$'topmed_'
wc -l $believeScratchLoc$'topmed_' # 12,712,741

# filter remaining to only INFO > 0.9 MAF > 0.001
zless $genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.snpstats.txt.gz'
#1	        2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24
#alternate_ids	rsid	chromosome	position	alleleA	alleleB	comment	alleleA_count	alleleB_count	alleleA_frequency	alleleB_frequency	minor_allele_frequency	minor_allele	major_allele	info	impute_info	missing_proportion	A	B	AA	AB	BB	NULL	total
#1:10481:A:G	1:10481:A:G	1	10481	A	G	NA	109232	0.089	0.999999	8.1478E-07	8.1478E-07	G	A	0.00428009	0.00428009	0	0	0	54615.9	0.089	0	0	54616
zcat $genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.snpstats.txt.gz' | awk '{if ($16 > 0.9 && $16 != "NA" && $16 != "" && $12 > 0.001) {print $3"\t"$4"\tchr"$2"\t"$5"\t"$6"\t"$16} }' > $believeScratchLoc$'topmed_09_MAF001'
head $believeScratchLoc$'topmed_09_MAF001'
wc -l $believeScratchLoc$'topmed_09_MAF001' # 7,714,085

zcat $genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.snpstats.txt.gz' | awk '{{print $3"\t"$4"\tchr"$2"\t"$5"\t"$6"\t"$16} }' > $believeScratchLoc$'topmed_all'

zcat $genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.snpstats.txt.gz' | awk '{{print $3":"$4} }' > $believeScratchLoc$'topmed_all_pos'



# subset to UKB
awk 'FNR == NR { file1[ $3 ] = $3;  next; }
FNR <= NR {  if( $3 in file1 ) {print $0 } }
' $believeScratchLoc$'topmed_' $believeScratchLoc$'topmed_09_MAF001' > $believeScratchLoc$'topmed_UKB'
head $believeScratchLoc$'topmed_UKB'
wc -l $believeScratchLoc$'topmed_UKB' # 5,825,182

# extract the 0.9 ones 
arguments=' --pgen '$genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.GT.pgen  --pvar '$genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.GT.pvar --psam '$genotypeFolder$'CAMBRIDGE-BELIEVE_Freeze_One.GxS.TOPMED_dosages.GT.COLLAB.psam --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeScratchLoc$'PCA_subset.fam --extract '$believeScratchLoc$'topmed_UKB  --maf 0.001 --geno 0.05 --hwe 1e-10 midp --make-bed --out '$believeDataLoc$'believe_imputed --allow-extra-chr --snps-only --bp-space 1'
$plink2 $arguments
# 5818177 variants remaining after main filters.
# CAMBRIDGE-BELIEVE_337860107_337860107  -> this guy is missing from the imputed dataset...


# overwrite the SNP ids in BELIEVE with the proper rsIDs.
cp $believeDataLoc$'believe_imputed.bim' $believeScratchLoc$'believe_imputed_oldbim'

# find the rsIds, for build 38
awk 'FNR == NR { 
file1[ $1":"$4 ] = $2;next; } 
FNR <= NR {  
if( $1":"$4 in file1) { print $1"\t"file1[$1":"$4]"\t"$3"\t"$4"\t"$5"\t"$6 }
else {print $0} }
' $believeScratchLoc$'UKBB_45K_liftedOver.bim' $believeScratchLoc$'believe_imputed_oldbim' > $believeScratchLoc$'believe_imputed_B38.bim'
head $believeScratchLoc$'believe_imputed_B38.bim'
wc -l $believeScratchLoc$'believe_imputed_B38.bim' # 5,818,177

# now also map the coordinates back to B37, as UKBB is in B37, so the positions would not match
awk 'FNR == NR { 
file1[ $2 ] = $4;next; } 
FNR <= NR {  
if( $2 in file1) { print $1"\t"$2"\t"$3"\t"file1[$2]"\t"$5"\t"$6 }
else {print $0} }
' $UKBB_PLINK1$'45K_IMPUTED.bim' $believeScratchLoc$'believe_imputed_B38.bim' > $believeDataLoc$'believe_imputed.bim'
head $believeDataLoc$'believe_imputed.bim'
wc -l $believeDataLoc$'believe_imputed.bim' # 5818177
# $UKBB_PLINK1$'45K_IMPUTED.bim' # B37
# $believeScratchLoc$'UKBB_45K_liftedOver.bim' # B38

# sanity check that the surviving SNPs are truly common with no missing ones between UKBB and BELIEVE
# awk 'FNR == NR { 
# file1[ $2 ] = $2;next; } 
# FNR <= NR {  
# if( $2 in file1) { print $0 }
 # }
# ' $UKBB_PLINK1$'45K_IMPUTED.bim' $believeDataLoc$'believe_imputed.bim' > $believeScratchLoc$'believe_ukbb_commonRSid.bim'
# wc -l $believeDataLoc$'believe_imputed.bim' # 5818177
# wc -l $believeScratchLoc$'believe_ukbb_commonRSid.bim' # 5818177

#######################################

# Extra QC step: Exclude excess hets from BLV (for the UKB the sample QC that came with the dataset has already done this
arguments=' --bfile '$believeDataLoc$'believe_imputed --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$'  --extract '$w_hm3_snplist$' --het --out '$believeDataLoc$'blv_hets --allow-extra-chr --allow-no-sex'
$plink $arguments

mkdir -p $believeDataLoc$'hets_PCA/'

# need a PCA for the heterozygosity step
# create a lower number SNP subset just for the PCA, using the previously selected variants suitable for PCA
arguments=' --bfile '$believeDataLoc$'believe_imputed --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$'  --extract '$believeRawLoc$'UKBB_BELIEVE.bim --make-bed --out '$believeDataLoc$'hets_PCA/blv_plink_pca --allow-extra-chr --allow-no-sex'
$plink $arguments

# perform PCA on the above 
arguments='/home/mk907/scripts/R/PCA_ncMCE.R '$believeDataLoc$'hets_PCA/blv_plink_pca.bed '$believeDataLoc$'hets_PCA/ '$NCORES_LDPred2$' 10'
Rscript $arguments

# Project excluded indis into the existing PCA space to increase numbers
arguments='/home/mk907/scripts/R/ProjectPCA.R '$believeDataLoc$'hets_PCA/ '$believeDataLoc$'hets_PCA/blv_plink_pca.bed '$believeDataLoc$'hets_PCA/PCAProj'
Rscript $arguments


cat $believeDataLoc$'hets_PCA/_PC_PCA' $believeDataLoc$'hets_PCA/PCAProj' > $believeDataLoc$'hets_PCA/PCAProj_ALL'
$believeDataLoc$'blv_hets'
######################################

# https://choishingwan.github.io/PRS-Tutorial/target/
R
baseloc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/"
filename="blv_hets.het"
dat <- read.table(paste0(baseloc,filename), header=T) # Read in the .het file, specify it has header

# correct for PCs as per the UKB documentation: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0579-z/MediaObjects/41586_2018_579_MOESM1_ESM.pdf
PCs <- read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/hets_PCA/PCAProj_ALL", header=F)

PCs = PCs[,c(1:7)] # select the first 6 PCs only
colnames(PCs) = c("IID", paste0("PC", 1:6))
str(PCs)

# merge PCs with the hets rates
dat_PCs = merge(dat,PCs, by = "IID")

# add quadratic terms into the DF
dat_PCs$PC1_sq = dat_PCs$PC1^2
dat_PCs$PC2_sq = dat_PCs$PC2^2
dat_PCs$PC3_sq = dat_PCs$PC3^2
dat_PCs$PC4_sq = dat_PCs$PC4^2
dat_PCs$PC5_sq = dat_PCs$PC5^2
dat_PCs$PC6_sq = dat_PCs$PC6^2
# fit all terms, quadratic terms, and all 2-way interactions
model = lm (F ~ (PC1 + PC2 + PC3 + PC4 + PC5 + PC6)^2 + PC1_sq + PC2_sq + PC3_sq + PC4_sq + PC5_sq + PC6_sq, data = dat_PCs )
summary(model)
# dat_PCs$F_adjusted = model$residuals # actually we want the fitted values, 
dat_PCs$F_adjusted =fitted(model)

# plot residuals against original
plotTitle=paste0(filename, " Heterozygosity adjusted")
png(filename=paste0(baseloc, filename, "_PCadjusted.png") , width=1280, height=720, res=128)
plot(dat_PCs$F, dat_PCs$F_adjusted, xlab = "raw heterozygosity", ylab = "PC adjusted heterozygosity", main = plotTitle )
dev.off()


# repeat QC procedure...
m <- mean(dat_PCs$F_adjusted) # Calculate the mean  
s <- sd(dat_PCs$F_adjusted) # Calculate the SD
UB= m+3*s
LB=m-3*s
print(paste0(nrow(subset(dat_PCs, F_adjusted >= UB )), " / ", nrow(dat_PCs), " are > 3 SD hets (DNA sample contamination?)"))
print(paste0(nrow(subset(dat_PCs, F_adjusted <= LB )), " / ", nrow(dat_PCs), " are < 3 SD hets (inbreeding?)"))
# "423 / 33710 are > 3 SD hets (DNA sample contamination?)"
#  "0 / 33710 are < 3 SD hets (inbreeding?)"

validIndices = which(dat_PCs$F_adjusted <= UB & dat_PCs$F_adjusted >= LB )
valid <- dat_PCs[validIndices,] # Get any samples with F coefficient within 3 SD of the population mean
print(paste0(nrow(valid), " / ", nrow(dat_PCs), " are within 3 SD hets"))


plotTitle=paste0(filename, "Adjusted Heterozygosity")
png(filename=paste0(baseloc, filename, "_adjusted.png") , width=1280, height=720, res=128)
hist(dat_PCs$F_adjusted, xlab = "heterozygosity", breaks = "Scott", main = plotTitle )
abline(v = UB, col="red", lwd=3, lty=2)
abline(v = LB, col="red", lwd=3, lty=2)
dev.off()

# too high hets: DNA sample contamination
# too low hets: inbreeding
write.table(valid[,c(1,2)], paste0(baseloc, filename, "_adjusted.valid"), quote=F, row.names=F, col.names=F) # print FID and IID for valid samples
write.table(dat_PCs[-validIndices,c(1,2)], paste0(baseloc, filename, "_adjusted.invalid"), quote=F, row.names=F, col.names=F) # print FID and IID for INVvalid samples



q() # exit R

wc -l $believeDataLoc$'blv_hets.het_adjusted.valid' # 33287

######################################
mv $believeDataLoc$'believe_imputed.fam' $believeDataLoc$'believe_imputed_noqc.fam'
mv $believeDataLoc$'believe_imputed.bim' $believeDataLoc$'believe_imputed_noqc.bim'
mv $believeDataLoc$'believe_imputed.bed' $believeDataLoc$'believe_imputed_noqc.bed'


arguments=' --bfile '$believeDataLoc$'believe_imputed_noqc --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeDataLoc$'blv_hets.het_adjusted.valid --make-bed --out '$believeDataLoc$'believe_imputed --allow-extra-chr --allow-no-sex'
$plink $arguments



####################
# export relevant subset of SNPs and indis for the UKB


# first, find the high quality SNPs
GENOTYPE_QC='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/pre_qc_data/affy_ukbiobank_array/raw_data/showcase_release_19Jul/snpqc/ukb_snp_qc.txt'
scratchLoc='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/scratch/'

head -n 1 $GENOTYPE_QC
#rs_id affymetrix_snp_id affymetrix_probeset_id chromosome position allele1_ref allele2_alt strand array Batch_b001_qc Batch_b002_qc Batch_b003_qc Batch_b004_qc Batch_b005_qc Batch_b006_qc Batch_b007_qc Batch_b008_qc Batch_b009_qc Batch_b010_qc Batch_b011_qc Batch_b012_qc Batch_b013_qc Batch_b014_qc Batch_b015_qc Batch_b016_qc Batch_b017_qc Batch_b018_qc Batch_b019_qc Batch_b020_qc Batch_b021_qc Batch_b022_qc Batch_b023_qc Batch_b024_qc Batch_b025_qc Batch_b026_qc Batch_b027_qc Batch_b028_qc Batch_b029_qc Batch_b030_qc Batch_b031_qc Batch_b032_qc Batch_b033_qc Batch_b034_qc Batch_b035_qc Batch_b036_qc Batch_b037_qc Batch_b038_qc Batch_b039_qc Batch_b040_qc Batch_b041_qc Batch_b042_qc Batch_b043_qc Batch_b044_qc Batch_b045_qc Batch_b046_qc Batch_b047_qc Batch_b048_qc Batch_b049_qc Batch_b050_qc Batch_b051_qc Batch_b052_qc Batch_b053_qc Batch_b054_qc Batch_b055_qc Batch_b056_qc Batch_b057_qc Batch_b058_qc Batch_b059_qc Batch_b060_qc Batch_b061_qc Batch_b062_qc Batch_b063_qc Batch_b064_qc Batch_b065_qc Batch_b066_qc Batch_b067_qc Batch_b068_qc Batch_b069_qc Batch_b070_qc Batch_b071_qc Batch_b072_qc Batch_b073_qc Batch_b074_qc Batch_b075_qc Batch_b076_qc Batch_b077_qc Batch_b078_qc Batch_b079_qc Batch_b080_qc Batch_b081_qc Batch_b082_qc Batch_b083_qc Batch_b084_qc Batch_b085_qc Batch_b086_qc Batch_b087_qc Batch_b088_qc Batch_b089_qc Batch_b090_qc Batch_b091_qc Batch_b092_qc Batch_b093_qc Batch_b094_qc Batch_b095_qc UKBiLEVEAX_b1_qc UKBiLEVEAX_b2_qc UKBiLEVEAX_b3_qc UKBiLEVEAX_b4_qc UKBiLEVEAX_b5_qc UKBiLEVEAX_b6_qc UKBiLEVEAX_b7_qc UKBiLEVEAX_b8_qc UKBiLEVEAX_b9_qc UKBiLEVEAX_b10_qc UKBiLEVEAX_b11_qc in_HetMiss in_Relatedness in_PCA PC1_loading PC2_loading PC3_loading PC4_loading PC5_loading PC6_loading PC7_loading PC8_loading PC9_loading PC10_loading PC11_loading PC12_loading PC13_loading PC14_loading PC15_loading PC16_loading PC17_loading PC18_loading PC9_loading PC20_loading PC21_loading PC22_loading PC23_loading PC24_loading PC25_loading PC26_loading PC27_loading PC28_loading PC9_loading PC30_loading PC31_loading PC32_loading PC33_loading PC34_loading PC35_loading PC36_loading PC37_loading PC38_loading PC9_loading PC40_loading in_Phasing_Input
#rs28659788 Affx-13546538 AX-32115783 1 723307 C G + 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA 0

# find SNPs that failed in any of the batches
awk '{
for (i =10; i<=115 ; i++)  {
 if( $i == "0") 
 { print $1; next; }
 } }' $GENOTYPE_QC > $scratchLoc$'GENOTYPE_FAILED_QC'
 
head $scratchLoc$'GENOTYPE_FAILED_QC'
wc -l $scratchLoc$'GENOTYPE_FAILED_QC' # 51,707 , 51K SNPs failed QC in at least 1 of the batches


# UKBB SNP QC file format: $14 = MAF, $18 = INFO
#     1               2              3             4          5       6      7              8               9               10             11                 12                  13                    14                   15           16             17          18
# alternate_ids     rsid          chromosome    position   alleleA alleleB comment   HW_exact_p_value HW_lrt_p_value   alleleA_count   alleleB_count   alleleA_frequency   alleleB_frequency     minor_allele_frequency minor_allele   major_allele     info     impute_info   missing_proportion A B AA AB BB NULL total
# 21:9411239_G_A    rs559462325      21        9411239     G          A       NA            1           0.97576         974788         29.9137             0.999969              3.06865e-05          3.06865e-05          A             G            0.260815     0.260815       5.97113e-16 0 0 487379 29.8824 0.0156863 2.91038e-10 487409

rm -rf $scratchLoc$'ALL_FAIL_INFO_09'
for ((i=1; i<=$numChroms; i++)); do

chromQCFile='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/reference_files/ukb_impv3_chr'$i$'_snpstats.txt'
# b) exclude SNPs that have low imputation quality INFO < 0.9

# need to skip lines with comments, otherwise awk will default to non-numeric comparisons
awk '/^[^#]/ { print $0 }' $chromQCFile | awk '{ if ($18 < 0.9 ) {print $2} }'  > $scratchLoc$i$'_FAIL_INFO_09'

cat $scratchLoc$i$'_FAIL_INFO_09' >>  $scratchLoc$'ALL_FAIL_INFO_09'

done

wc -l $scratchLoc$'ALL_FAIL_INFO_09' # 75,887,604

# also add in the genotype failed SNPs
cat $scratchLoc$'GENOTYPE_FAILED_QC' $scratchLoc$'ALL_FAIL_INFO_09' >  $scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC'

# find out how many we are going to exclude
wc -l $scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC' # 75939311 


########

# Extract the final subset of imputed SNPs from the UKBB
# generate PRS per each chrom
mkdir -p $UKBB_PLINK1

for ((i=1; i<=22; i++)); do

# only submit if it doesn't exist yet
outfile=$UKBB_PLINK1$i$'_45K.bim'
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
pgen='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr'$i
arguments=' --pfile '$pgen$' --memory 43000  --threads '$SLURM_CPUS_ON_NODE$' --maf 0.001 --geno 0.05 --hard-call-threshold 0.1 --keep '$believeRawLoc$'PCA_related/PCA_related_PCAIndis --exclude '$scratchLoc$'ALL_FAIL_INFO_09_GENOTYPEQC --make-bed --allow-extra-chr --allow-no-sex --snps-only --bp-space 1 --out '$UKBB_PLINK1$i$'_45K'
#sbatch --mem 43000 --cpus-per-task 16 --time 4:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$plink2 $arguments"
$plink2 $arguments

fi
done

# merge chroms into one
plinkFileList=$believeScratchLoc'/mergelist.txt'
rm -f $plinkFileList

for ((j=1; j<=22; j++)); do
echo $UKBB_PLINK1$j$'_45K' >> ${plinkFileList}

done

arguments=' --memory  43000 --merge-list '$plinkFileList$' --make-bed --out '$UKBB_PLINK1$'45K_IMPUTED --allow-extra-chr --allow-no-sex'
$plink $arguments


############
# Perform QC on UKB
# remove ambiguous SNPs
remove_ambiguous_alleles_bim $UKBB_PLINK1$'45K_IMPUTED.bim'

# do NOT do  --hwe 1e-10 midp, as the UKB is a non-homogenous dataset. They did HWE already in the $GENOTYPE_QC file, rationale: (https://biobank.ctsu.ox.ac.uk/crystal/crystal/docs/genotyping_qc.pdf)
# UKB --maf 0.001 --geno 0.05
arguments=' --memory '$PLINKRAM$' --bfile '$UKBB_PLINK1$'45K_IMPUTED --keep '$believeRawLoc$'UKBB_PCA_subset_all.fam --write-snplist --out '$believeScratchLoc$'45K_IMPUTED_QC --extract '$UKBB_PLINK1$'45K_IMPUTED.bim_noAmbiguousAlleles_SNPlist --maf 0.001 --geno 0.05 --allow-extra-chr --allow-no-sex'
$plink $arguments

# the above is NOT enough, as the maf etc thresholds are for the 45K total, which includes GBP. so there would be SNPs, like rs1589393, that will have > 1% total, but 0% MAF when just looking at the UK_SEA. 
# so we need to perform a QC just on the UK SEA pop
awk '{if ($2 == "UK_SEA") print $1"\t"$1}' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4' > $believeRawLoc$'UK_SEA_intermediatelist'
head $believeRawLoc$'UK_SEA_intermediatelist'
wc -l $believeRawLoc$'UK_SEA_intermediatelist'

arguments=' --memory '$PLINKRAM$' --bfile '$UKBB_PLINK1$'45K_IMPUTED --keep '$believeRawLoc$'UK_SEA_intermediatelist --write-snplist --extract '$UKBB_PLINK1$'45K_IMPUTED.bim_noAmbiguousAlleles_SNPlist --out '$believeScratchLoc$'45K_IMPUTED_QC_2 --maf 0.001 --geno 0.05 --allow-extra-chr --allow-no-sex'
$plink $arguments
grep "rs201676649" $believeScratchLoc$'45K_IMPUTED_QC_2.snplist'
grep "rs1589393" $believeScratchLoc$'45K_IMPUTED_QC_2.snplist'

# also subset to GBP and AFR too
awk '{if ($2 == "GBP") print $1"\t"$1}' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4' > $believeRawLoc$'GBP_intermediatelist'
arguments=' --memory '$PLINKRAM$' --bfile '$UKBB_PLINK1$'45K_IMPUTED --keep '$believeRawLoc$'GBP_intermediatelist --write-snplist --extract '$UKBB_PLINK1$'45K_IMPUTED.bim_noAmbiguousAlleles_SNPlist --out '$believeScratchLoc$'45K_IMPUTED_QC_3 --maf 0.001 --geno 0.05 --allow-extra-chr --allow-no-sex'
$plink $arguments


head $ancestryL$'ukbb_believe_ethnic'
wc -l $ancestryL$'ukbb_believe_ethnic'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $ancestryL$'ukbb_believe_ethnic'
# Nigeria 3735
# Caribbean 2261
# take an afro caribbian keeplist as an outgroup
awk '{if( $2 == "Nigeria" || $2 == "Caribbean") print $1"\t"$1}' $ancestryL$'ukbb_believe_ethnic' > $believeDataLoc$'AFR_keeplist'
head $believeDataLoc$'AFR_keeplist'
wc -l $believeDataLoc$'AFR_keeplist' # 5996. this is just enough to be the same as the ones used for the other simulations
arguments=' --memory '$PLINKRAM$' --bfile '$UKBB_PLINK1$'45K_IMPUTED --keep '$believeDataLoc$'AFR_keeplist --write-snplist --extract '$UKBB_PLINK1$'45K_IMPUTED.bim_noAmbiguousAlleles_SNPlist --out '$believeScratchLoc$'45K_IMPUTED_QC_4 --maf 0.001 --geno 0.05 --allow-extra-chr --allow-no-sex'
$plink $arguments

$believeScratchLoc$'45K_IMPUTED_QC.snplist' # original list: with everyne at once
$believeScratchLoc$'45K_IMPUTED_QC_2.snplist' # UK_SEA split out
$believeScratchLoc$'45K_IMPUTED_QC_3.snplist' # GBP split out
$believeScratchLoc$'45K_IMPUTED_QC_4.snplist' # AFR split out

# intersect the 4 keeplists
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeScratchLoc$'45K_IMPUTED_QC.snplist' $believeScratchLoc$'45K_IMPUTED_QC_2.snplist' > $believeScratchLoc$'45K_IMPUTED_QC_12.snplist'
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeScratchLoc$'45K_IMPUTED_QC_12.snplist' $believeScratchLoc$'45K_IMPUTED_QC_3.snplist' > $believeScratchLoc$'45K_IMPUTED_QC_123.snplist'
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeScratchLoc$'45K_IMPUTED_QC_123.snplist' $believeScratchLoc$'45K_IMPUTED_QC_4.snplist' > $believeScratchLoc$'45K_IMPUTED_QC_1234.snplist'
head $believeScratchLoc$'45K_IMPUTED_QC_1234.snplist'
wc -l $believeScratchLoc$'45K_IMPUTED_QC_1234.snplist'  # 4340304, so we loose a lot of SNPs if we intersect all 4...

####################
# intersect the UKB QC-d snp list, with what we ended up with in BELIEVE, 
# if we intersected then with the BELIEVE SNP list, then we would only have 3.2 mil SNPs
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeDataLoc$'believe_imputed.snplist' $believeScratchLoc$'45K_IMPUTED_QC_1234.snplist' > $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON2.snplist'
head $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON2.snplist'
wc -l $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON2.snplist' # 3,262,091

#awk '{print $2}' $believeDataLoc$'believe_imputed.bim' > $believeDataLoc$'believe_imputed.snplist'
#awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeDataLoc$'believe_imputed.snplist' $believeScratchLoc$'45K_IMPUTED_QC.snplist' > $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON.snplist'
#head $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON.snplist'
#wc -l $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON.snplist' # 4,910,709 # 4.9 mil SNPs in total


# lets see if we intersect with just the UK_SEA
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeDataLoc$'believe_imputed.snplist' $believeScratchLoc$'45K_IMPUTED_QC_12.snplist' > $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON_UKSEA.snplist'
head $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON_UKSEA.snplist'
wc -l $believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON_UKSEA.snplist'  # 4,610,701, ie 4.6 mil SNPs
# if we intersected without the AFR and GBP ones... we keep almost all SNPs. and since the GWAS will NOT use anyone except UK_SEA and BELIEVE, we should go with this

# merge UKBB and BELIEVE
plinkFileList=$believeScratchLoc$'plinkFileList'
rm -rf $plinkFileList
echo $UKBB_PLINK1$'45K_IMPUTED' >> ${plinkFileList}
echo $believeDataLoc$'believe_imputed' >> ${plinkFileList}


# then extract these SNPs that are common, otherwise we would have the full union set of SNPs, with the missing SNPs set to 0 for the other cohort
arguments=' --memory '$PLINKRAM$' --merge-list '$plinkFileList$' --make-bed --out '$believeRawLoc$'UKBB_BELIEVE_FULL --extract '$believeScratchLoc$'UKB_BELIEVE_IMPUTED_QC_QCCOMMON_UKSEA.snplist --allow-extra-chr --allow-no-sex'
sbatch --mem ${PLINKRAM} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$plink $arguments"

wc -l $believeRawLoc$'UKBB_BELIEVE_FULL.bim' # 4,610,701 SNPs in total
wc -l $believeRawLoc$'UKBB_BELIEVE_FULL.fam' # 79,063 indis in total


###########################################################################
# IX. Expand phenotypes to include more cardio metabolic traits
###########################################################################

# find phenos in BELIEVE questionnaire that can be matched to the UKBB curated list (get infos about how to QC them)
# extract all UKBB phenos
#BELIEVE       UKBB
#Mi            hxmi
#Stroke        hxstroke
#Type2Diab     hxdiabbin
#Angina        hxang
#SystolicBp1   sbp
#DiastolicBp1  dbp

# correct for the effect of anti hypertensive medication, I would add +10 to DBP, and +15 SBP to those who are taking meds. 
#In the UKBB I would use "hypdbin", and  for BELIEVE it would be"AntiHyperT

# for BELIEVE SBP/DBP:
# make sure none are missing
# check for  high/low values, for DBP>SBP, for first and second measurements that are highly discrepant
# then take average of surviving (SystolicBp1+SystolicBp2)/2 and (DiastolicBp1+DiastolicBp2)/2

# what cutoff to use for sbp/dbp
# Not sure physiologically (there may be some cutoffs used on other CEU studies or UK Biobank that we could look for?) but for example 
# CEU Data Managers often use 777 or 999 as missing codes so look out for those, plus for other variables in BELIEVE we found an over-representation of the values 1, 2 and 10, which we think are where the RA has answered the wrong question in the form.

# UKBB: 
R
library(readstata13)
library('data.table')
ceu <- read.dta13("/rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/20210901/STATA/analysis.dta")
setDT(ceu)
vars = attributes(ceu)$var.labels
names(vars) = names(ceu)
myData=ceu[,.(eid=idno,MI=hxmi, stroke=hxstroke, T2D=hxdiabbin, angina=hxang, sbp,dbp, hypdbin)]

# set indis with no clear yes/no category to missing (these will be excluded later, we don't want to exclude people who are missing anything as that would result in too few people due to lowest common denominator effects for many phenos and we only have like 6K SAS indis)
levels(myData$MI)
myData$MI= as.character(myData$MI) # need this otherwise I wouldn't be able to transcode the characters into numeric values
myData$MI[which(myData$MI=="Poss/susp")] = NA
myData$MI[which(myData$MI=="Yes")] = 2
myData$MI[which(myData$MI=="No")] = 1

levels(myData$stroke)
myData$stroke= as.character(myData$stroke) # need this otherwise I wouldn't be able to transcode the characters into numeric values
myData$stroke[which(myData$stroke=="Poss/susp")] = NA
myData$stroke[which(myData$stroke=="Yes")] = 2
myData$stroke[which(myData$stroke=="No")] = 1

levels(myData$T2D)
myData$T2D= as.character(myData$T2D) # need this otherwise I wouldn't be able to transcode the characters into numeric values
myData$T2D[which(myData$T2D=="Definite diabetic")] = 2
myData$T2D[which(myData$T2D=="Other")] = 1
table(myData$T2D)

levels(myData$angina)
myData$angina= as.character(myData$angina) # need this otherwise I wouldn't be able to transcode the characters into numeric values
myData$angina[which(myData$angina=="Poss/susp")] = NA
myData$angina[which(myData$angina=="Yes")] = 2
myData$angina[which(myData$angina=="No")] = 1
table(myData$angina)

# correct SBP/DBP for being on antihypertensive meds
levels(myData$hypdbin)
onMeds =which(myData$hypdbin=="Current")
myData$sbp[onMeds] = myData$sbp[onMeds] +15
myData$dbp[onMeds] = myData$dbp[onMeds] +10

myData <- subset(myData, select = -c(hypdbin))

write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/scratch/ukbb_outcomes_raw", sep="\t", row.names = F, col.names = T, quote = F)

subs= myData[which(is.na(myData$MI) ==F  ),]
write.table(cbind(subs$eid,subs$eid,subs$MI), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_MI.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= myData[which(is.na(myData$stroke) ==F  ),]
write.table(cbind(subs$eid,subs$eid,subs$stroke), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_stroke.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= myData[which(is.na(myData$T2D) ==F  ),]
write.table(cbind(subs$eid,subs$eid,subs$T2D), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_t2d.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= myData[which(is.na(myData$angina) ==F  ),]
write.table(cbind(subs$eid,subs$eid,subs$angina), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_angina.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= myData[which(is.na(myData$sbp) ==F  ),]
write.table(cbind(subs$eid,subs$eid,subs$sbp), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_sbp.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= myData[which(is.na(myData$dbp) ==F  ),]
write.table(cbind(subs$eid,subs$eid,subs$dbp), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/ukbb_dbp.phe", sep="\t", row.names = F, col.names = F, quote = F)

quit()


##############
# BELIEVE
# the pheno files are dirty, they have "," within fields, so awk cannot process them. Have to use R:
R
# without specifying the colClasses, it would read GeneticID as number, and remove the leading zeroes, which would mismatch the .fam files # , colClasses=c("GeneticID"="character")
phenoMap=read.csv("/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/BELIEVEdata_P5031_20220310.csv" )
phenoMap$id = paste("CAMBRIDGE-BELIEVE_" , as.character(phenoMap$GeneticID),"_" ,as.character(phenoMap$GeneticID), sep="" )
phenoMap_exist = phenoMap[which(phenoMap$Age != "NA" & phenoMap$Age != ""& is.na(phenoMap$Age) == F & phenoMap$Gender != "NA" & phenoMap$Gender != ""& is.na(phenoMap$Gender) == F),]


#  SBP=50-250, DBP=30-150
validPhenoIndices = which(phenoMap_exist$SystolicBp1 != "" & is.na(phenoMap_exist$SystolicBp1) == F & phenoMap_exist$SystolicBp1 >= 50  & phenoMap_exist$SystolicBp1 <= 250 & 

phenoMap_exist$SystolicBp2 != "" & is.na(phenoMap_exist$SystolicBp2) == F & phenoMap_exist$SystolicBp2 >= 50  & phenoMap_exist$SystolicBp2 <= 250 & 

phenoMap_exist$DiastolicBp1 != "" & is.na(phenoMap_exist$DiastolicBp1) == F & phenoMap_exist$DiastolicBp1 >= 30  & phenoMap_exist$DiastolicBp1 <= 150 & 

phenoMap_exist$DiastolicBp2 != "" & is.na(phenoMap_exist$DiastolicBp2) == F & phenoMap_exist$DiastolicBp2 >= 30  & phenoMap_exist$DiastolicBp2 <= 150 & 

phenoMap_exist$DiastolicBp1 < phenoMap_exist$SystolicBp1 & phenoMap_exist$DiastolicBp2 < phenoMap_exist$SystolicBp2 & 

phenoMap_exist$Mi != "" & is.na(phenoMap_exist$Mi) == F & phenoMap_exist$Stroke != "" & is.na(phenoMap_exist$Stroke) == F  & phenoMap_exist$Type2Diab != "" & is.na(phenoMap_exist$Type2Diab) == F & phenoMap_exist$Angina != "" & is.na(phenoMap_exist$Angina)  == F
)

length(validPhenoIndices)
nrow(phenoMap_exist)

phenos = phenoMap_exist[validPhenoIndices,]

# BELIEVE has case/control status in reverse, 1=yes, 2=no???
# need to check cholesterol: it was OK, there I made the opposite decision, made UKBB 2=control and 1=case too.
# but for the rest let's flip it to 1=control, and 2=case, as that is more common (eg PLINK)

caseIndices= which(phenos$Mi == 1)
controlIndices= which(phenos$Mi == 2)
phenos$Mi[caseIndices] = 2
phenos$Mi[controlIndices] = 1


caseIndices= which(phenos$Stroke == 1)
controlIndices= which(phenos$Stroke == 2)
phenos$Stroke[caseIndices] = 2
phenos$Stroke[controlIndices] = 1

caseIndices= which(phenos$Type2Diab == 1)
controlIndices= which(phenos$Type2Diab == 2)
phenos$Type2Diab[caseIndices] = 2
phenos$Type2Diab[controlIndices] = 1

caseIndices= which(phenos$Angina == 1)
controlIndices= which(phenos$Angina == 2)
phenos$Angina[caseIndices] = 2
phenos$Angina[controlIndices] = 1


phenos$SystolicBp = (phenos$SystolicBp1 + phenos$SystolicBp2)/2
phenos$DiastolicBp = (phenos$DiastolicBp1 + phenos$DiastolicBp2)/2


subs= phenos[which(is.na(phenos$Mi) ==F  ),]
write.table(cbind(subs$id,subs$id,subs$Mi), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/MI.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= phenos[which(is.na(phenos$Stroke) ==F  ),]
write.table(cbind(subs$id,subs$id,subs$Stroke), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/stroke.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= phenos[which(is.na(phenos$Type2Diab) ==F  ),]
write.table(cbind(subs$id,subs$id,subs$Type2Diab), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/t2d.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= phenos[which(is.na(phenos$Angina) ==F  ),]
write.table(cbind(subs$id,subs$id,subs$Angina), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/angina.phe", sep="\t", row.names = F, col.names = F, quote = F)


subs= phenos[which(is.na(phenos$SystolicBp) ==F  ),]
write.table(cbind(subs$id,subs$id,subs$SystolicBp), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/sbp.phe", sep="\t", row.names = F, col.names = F, quote = F)

subs= phenos[which(is.na(phenos$DiastolicBp) ==F  ),]
write.table(cbind(subs$id,subs$id,subs$DiastolicBp), "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/dbp.phe", sep="\t", row.names = F, col.names = F, quote = F)


####
# get prevelances:
table(phenos$HighBloodChol)[1] / (table(phenos$HighBloodChol)[1] + table(phenos$HighBloodChol)[2])
# cholesterol=0.128921

table(phenos$Mi)[2] / (table(phenos$Mi)[1] + table(phenos$Mi)[2])
# MI=0.02070714

table(phenos$Stroke)[2] / (table(phenos$Stroke)[1] + table(phenos$Stroke)[2])
# Stroke=0.02253369

table(phenos$Angina)[2] / (table(phenos$Angina)[1] + table(phenos$Angina)[2])
# Angina=0.04143291


table(phenos$Type2Diab)[2] / (table(phenos$Type2Diab)[1] + table(phenos$Type2Diab)[2])
# Type2Diab=0.1703539
quit()


# BLV 
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'MI.phe' $believeScratchLoc$'MI_cov2' 1
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'stroke.phe' $believeScratchLoc$'stroke_cov2' 1
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'t2d.phe' $believeScratchLoc$'t2d_cov2' 1
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'angina.phe' $believeScratchLoc$'angina_cov2' 1
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'sbp.phe' $believeScratchLoc$'sbp_cov2' 0
GenBLV2 $believeDataLoc$'covars_new.phe' $believeDataLoc$'dbp.phe' $believeScratchLoc$'dbp_cov2' 0

# UKBB
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_MI.phe' $believeScratchLoc$'ukbb_MI_cov2' 1
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_stroke.phe' $believeScratchLoc$'ukbb_stroke_cov2' 1
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_t2d.phe' $believeScratchLoc$'ukbb_t2d_cov2' 1
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_angina.phe' $believeScratchLoc$'ukbb_angina_cov2' 1
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_sbp.phe' $believeScratchLoc$'ukbb_sbp_cov2' 0
GenUKBB2 $believeScratchLoc$'ukbb_covs_raw2' $believeDataLoc$'ukbb_dbp.phe' $believeScratchLoc$'ukbb_dbp_cov2' 0

# cat ukbb and believe
cat $believeScratchLoc$'MI_cov2_res_id' $believeScratchLoc$'ukbb_MI_cov2_res_id' > $believeDataLoc$'MI_all2'
cat $believeScratchLoc$'stroke_cov2_res_id' $believeScratchLoc$'ukbb_stroke_cov2_res_id' > $believeDataLoc$'stroke_all2'
cat $believeScratchLoc$'t2d_cov2_res_id' $believeScratchLoc$'ukbb_t2d_cov2_res_id' > $believeDataLoc$'t2d_all2'
cat $believeScratchLoc$'angina_cov2_res_id' $believeScratchLoc$'ukbb_angina_cov2_res_id' > $believeDataLoc$'angina_all2'

cat $believeScratchLoc$'sbp_cov2_res_id' $believeScratchLoc$'ukbb_sbp_cov2_res_id' > $believeDataLoc$'sbp_all2'
cat $believeScratchLoc$'dbp_cov2_res_id' $believeScratchLoc$'ukbb_dbp_cov2_res_id' > $believeDataLoc$'dbp_all2'




############## 
# phenotypes (post adjustment for covariates)
#cp $believeScratchLoc$'labelsPC4_height' $believeScratchLoc$'labelsPC4_height_old'

# we need to remove any BELIEVE individual who is NOT in the imputed dataset to prevent picking it (eg: "CAMBRIDGE-BELIEVE_337860107_337860107")
awk 'FNR == NR { test[ $2 ] = $2; next; } FNR <= NR { if( $2 != "BELIEVE" || $1 in test) {print $0 } }
' $believeDataLoc$'believe_imputed.fam' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4' > $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_IMPUTED'
wc -l $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4' # 47732
wc -l $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_IMPUTED' # 47307 , so we lost some more due to the hets QC

#cp $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_IMPUTED' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_IMPUTED_old' # the 47730 list

# create subsets of indis that have phenos for all 4 phenos (first check edu attainment)
awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'height_all2' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_IMPUTED' > $believeScratchLoc$'labelsPC4_height'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_height' # check how many of each group left
#QCFAIL 169
#pop 1
#GBP 6280
#BELIEVE 31763
#UK_SEA 7386

awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'bmi_all2' $believeScratchLoc$'labelsPC4_height' > $believeScratchLoc$'labelsPC4_heightbmi'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_heightbmi' # check how many of each group left
# QCFAIL 169
# pop 1
# GBP 6275
# BELIEVE 24785 # we loose here 7K BELIEVE
# UK_SEA 7377

awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'highCholesterol_all2' $believeScratchLoc$'labelsPC4_heightbmi' > $believeScratchLoc$'labelsPC4_heightbmichol'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_heightbmichol' # check how many of each group left


awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'MI_all2' $believeScratchLoc$'labelsPC4_heightbmichol' > $believeScratchLoc$'labelsPC4_p2_mi'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_p2_mi' # check how many of each group left

awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'stroke_all2' $believeScratchLoc$'labelsPC4_p2_mi' > $believeScratchLoc$'labelsPC4_p2_mistr'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_p2_mistr' # check how many of each group left

awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'t2d_all2' $believeScratchLoc$'labelsPC4_p2_mistr' > $believeScratchLoc$'labelsPC4_p2_mistrt2d'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_p2_mistrt2d' # check how many of each group left

awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'angina_all2' $believeScratchLoc$'labelsPC4_p2_mistrt2d' > $believeScratchLoc$'labelsPC4_p2_mistrt2dang'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_p2_mistrt2dang' # check how many of each group left

awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'sbp_all2' $believeScratchLoc$'labelsPC4_p2_mistrt2dang' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbp'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbp' # check how many of each group left

awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'dbp_all2' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbp' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp' # check how many of each group left

# QCFAIL 155
# pop 1
# GBP 5921
# BELIEVE 24760
# UK_SEA 6664
# so we will have 5921 indis in total?? but why not 6664? Why use the GBP as the metric??


awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $1 in test ) {print $0 } }
' $believeDataLoc$'edu_all2' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu' # check how many of each group left
# QCFAIL 78
# pop 1
# GBP 3925
# BELIEVE 21148
# UK_SEA 3453

# UK_SEA 3453 # so edu is the lowest common denom again

# first randomize the order for all, so that we can just use a master list 
# find out how many has edu
numHasEdu=$(awk '{if($2 == "UK_SEA")hasedu++} END {print hasedu}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu')
echo $numHasEdu

numHasEduGBP=$(awk '{if($2 == "GBP")hasedu++} END {print hasedu}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu')
echo $numHasEduGBP # 3925

# randomize the 6664 list
awk '{if (FNR > 1) print $0 }'  $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2' # 37500
numToP=$(wc -l < $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2')
echo $numToP


arguments='/home/mk907/scripts/R/randomNums.R '$numToP$' 42 '$believeScratchLoc$'randomNumspc4'
Rscript $arguments


awk 'FNR == NR { test[ NR ] = $0; next; } FNR <= NR { print test[$1]  }
' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2' $believeScratchLoc$'randomNumspc4' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand'
head $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand' # 37500

# randomize the has_edu list
awk '{if (FNR > 1) print $0 }'  $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu2'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu2' # 28829
numToP=$(wc -l < $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu2')
echo $numToP # 28604

arguments='/home/mk907/scripts/R/randomNums.R '$numToP$' 42 '$believeScratchLoc$'randomNumspc4_hasEdu'
Rscript $arguments

awk 'FNR == NR { test[ NR ] = $0; next; } FNR <= NR { print test[$1]  }
' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu2' $believeScratchLoc$'randomNumspc4_hasEdu' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_rand'
head $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_rand'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_rand' # 28604


# first select from the edu subset hasEdu list
awk -v numHasEduGBP=$numHasEduGBP '{if ($2 == "GBP" && count < numHasEduGBP) {print $0; count++;} }' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_rand' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_GBP_hasEdu'
head $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_GBP_hasEdu'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_GBP_hasEdu' # 3925
# DONT CARE about GBP now

awk -v numHasEdu=$numHasEdu '{if ($2 == "BELIEVE" && count < numHasEdu) {print $0; count++;} }' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_rand' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_BELIEVE_hasEdu'
head $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_BELIEVE_hasEdu'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_BELIEVE_hasEdu' # 3453

awk -v numHasEdu=$numHasEdu '{if ($2 == "UK_SEA" && count < numHasEdu) {print $0; count++;} }' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_rand' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_UK_SEA_hasEdu'
head $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_UK_SEA_hasEdu'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_UK_SEA_hasEdu' # 3453

# concat them to get an exclude list # 
cat $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_GBP_hasEdu' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_BELIEVE_hasEdu' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_UK_SEA_hasEdu' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_All3_hasEdu'
#cat $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_BELIEVE_hasEdu' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_UK_SEA_hasEdu' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_All3_hasEdu'


# exclude the above from the rest of the phenos list
awk 'FNR == NR { test[ $1 ] = $1; next; } FNR <= NR { if ($1 in test == 0) print $0 }
' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_All3_hasEdu' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand' > $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand3'
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand'  # 37500
wc -l $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand3' # 26669         30594


# now from the remaining list, pick 6664 - 3453 = 3211 indis from each category
numTotal=$(awk '{if($2 == "UK_SEA")numTotal++} END {print numTotal}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp')
remainder=$(awk -v numHasEdu=$numHasEdu -v numTotal=$numTotal 'BEGIN {print numTotal-numHasEdu}')
echo $remainder # 3211

remainderGBP=$(awk -v numHasEduGBP=$numHasEduGBP -v numTotal=$numTotal 'BEGIN {print numTotal-numHasEduGBP}')
echo $remainderGBP # 2739


# (this accomplishes: we first ranomdly picked the hasEdu indis that have all 4 phenos, then we pick an additional remainder indis who have phenos just for the others)
awk -v remainderGBP=$remainderGBP '{if ($2 == "GBP" && count < remainderGBP) {print $0; count++;} }' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand3' > $believeScratchLoc$'labelsPC4_otherphenos_GBP_remainder'
head $believeScratchLoc$'labelsPC4_otherphenos_GBP_remainder'
wc -l $believeScratchLoc$'labelsPC4_otherphenos_GBP_remainder' # 1996

awk -v remainder=$remainder '{if ($2 == "BELIEVE" && count < remainder) {print $0; count++;} }' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand3' > $believeScratchLoc$'labelsPC4_otherphenos_BELIEVE_remainder'
wc -l $believeScratchLoc$'labelsPC4_otherphenos_BELIEVE_remainder'  # 3211

awk -v remainder=$remainder '{if ($2 == "UK_SEA" && count < remainder) {print $0; count++;} }' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp2_rand3' > $believeScratchLoc$'labelsPC4_otherphenos_UK_SEA_remainder'
wc -l $believeScratchLoc$'labelsPC4_otherphenos_UK_SEA_remainder'  # 3211

# now concat the above 3, together with the all3 hasEdu list
cat $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu_All3_hasEdu' $believeScratchLoc$'labelsPC4_otherphenos_GBP_remainder' $believeScratchLoc$'labelsPC4_otherphenos_BELIEVE_remainder' $believeScratchLoc$'labelsPC4_otherphenos_UK_SEA_remainder' > $believeScratchLoc$'labelsPC4_ALL'
head $believeScratchLoc$'labelsPC4_ALL'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_ALL'
#GBP 5921 # so we have less for GBP, but that doesnt matter, as we won't use that for anything
#UK_SEA 6664
#BELIEVE 6664

##########
# export the imputed genotype data for the above indis
# exclude GBP here, as we are not using them for anything from here onwards
awk '{if ($2 != "GBP") { print $1"\t"$1 } }' $believeScratchLoc$'labelsPC4_ALL' > $believeScratchLoc$'labelsPC4_ALL_keep'
head $believeScratchLoc$'labelsPC4_ALL_keep'
wc -l $believeScratchLoc$'labelsPC4_ALL_keep' # 13328

arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL --keep '$believeScratchLoc$'labelsPC4_ALL_keep --make-bed --out '$believeRawLoc$'UKBB_BELIEVE_FULL_6664 --allow-extra-chr --allow-no-sex'
$plink $arguments


###############

# create a labels file, with matching order to the .fam
awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if( $2 in test ) {print $2"\t"test[$2] } }
' $believeScratchLoc$'labelsPC4_ALL' $believeRawLoc$'UKBB_BELIEVE_FULL_6664.fam' > $believeRawLoc$'6664_ALL_labels'
head $believeRawLoc$'6664_ALL_labels'
wc -l $believeRawLoc$'6664_ALL_labels'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeRawLoc$'6664_ALL_labels'


# merge the phenotypes too with the pop informations, to make them aligned exactly to the genotype order in the fams
awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'bmi_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'bmi_6664_pop'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'height_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'height_6664_pop'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'highCholesterol_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'highCholesterol_6664_pop'

# for edu, we want to make sure all 3 pops are  the same size, the below would result in UK_SAS being too small
awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'edu_all2' $believeRawLoc$'6664_ALL_labels'  > $believeDataLoc$'edu_6664_pop2'

# so we subset it so that GBP/BELIEVE all have a max of only what we have UK_SEA
numHasEdu=$(awk '{if($2 == "UK_SEA")hasedu++} END {print hasedu}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbpedu')
echo $numHasEdu
awk -v numHasEdu=$numHasEdu '{if (FNR == 1 || $4 == "UK_SEA" || $4 == "BELIEVE" && BELIEVEcount < numHasEdu) {print $0; if ($4 == "BELIEVE") {BELIEVEcount++}  } }' $believeDataLoc$'edu_6664_pop2' > $believeDataLoc$'edu_6664_pop'
head $believeDataLoc$'edu_6664_pop'
wc -l $believeDataLoc$'edu_6664_pop'
awk '{count[$4]++} END {for (word in count) print word, count[word]}' $believeDataLoc$'edu_6664_pop'
# BELIEVE 3453
# UK_SEA 3453


awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'MI_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'MI_6664_pop'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'stroke_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'stroke_6664_pop'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'t2d_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'t2d_6664_pop'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'angina_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'angina_6664_pop'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'sbp_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'sbp_6664_pop'

awk 'FNR == NR { test[ $1 ] = $2; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } if( $1 in test ) {print $1"\t"$1"\t"test[$1]"\t"$2 } }
' $believeDataLoc$'dbp_all2' $believeRawLoc$'6664_ALL_labels' > $believeDataLoc$'dbp_6664_pop'


# prepare phenotypes/covariates for the binary traits, so that we could have a model that uses covariates (quantitative traits are pre-adjusted)
# harmonize covars between UKBB and BLV and concat them
# BLV, 1=male, 2=Female
awk '{if ($5==2) {sex="Female"} else {sex="Male"} print $1"\t"$3"\t"$4"\t"sex  }' $believeDataLoc$'covars_new.phe' > $believeScratchLoc$'covars_harmonized.phe'
head $believeScratchLoc$'covars_harmonized.phe'
wc -l $believeScratchLoc$'covars_harmonized.phe'

awk '{if (FNR > 1) {print $1"\t"$2"\t"$3"\t"$4 } }' $believeScratchLoc$'ukbb_covs_raw2' > $believeScratchLoc$'ukbb_covs_harmonized.phe'
head $believeScratchLoc$'ukbb_covs_harmonized.phe'
wc -l $believeScratchLoc$'ukbb_covs_harmonized.phe'

cat $believeScratchLoc$'ukbb_covs_harmonized.phe' $believeScratchLoc$'covars_harmonized.phe' > $believeScratchLoc$'ukbb_blv_covars_harmonized.phe'
wc -l $believeScratchLoc$'ukbb_blv_covars_harmonized.phe'

cat $believeDataLoc$'MI.phe' $believeDataLoc$'ukbb_MI.phe' > $believeScratchLoc$'ukbb_blv_MI_binary.phe'
cat $believeDataLoc$'stroke.phe' $believeDataLoc$'ukbb_stroke.phe' > $believeScratchLoc$'ukbb_blv_stroke_binary.phe'
cat $believeDataLoc$'t2d.phe' $believeDataLoc$'ukbb_t2d.phe' > $believeScratchLoc$'ukbb_blv_t2d_binary.phe'
cat $believeDataLoc$'angina.phe' $believeDataLoc$'ukbb_angina.phe' > $believeScratchLoc$'ukbb_blv_angina_binary.phe'
cat $believeDataLoc$'highCholesterol_all.phe' $believeDataLoc$'ukbb_cholesterol.phe' > $believeScratchLoc$'ukbb_blv_cholesterol_binary.phe'

binarPheno $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_MI_binary.phe' $believeDataLoc$'MI_6664_pop'
binarPheno $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_stroke_binary.phe' $believeDataLoc$'stroke_6664_pop'
binarPheno $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_t2d_binary.phe' $believeDataLoc$'t2d_6664_pop'
binarPheno $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_angina_binary.phe' $believeDataLoc$'angina_6664_pop'
binarPheno $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_cholesterol_binary.phe' $believeDataLoc$'highCholesterol_6664_pop'


#############
# get disease prevelances for binary traits for all pops
prevMI_GBP=$(findprevalence "GBP" $believeDataLoc$'ukbb_MI.phe') 
echo $prevMI_GBP # 0.0243202
prevStroke_GBP=$(findprevalence "GBP" $believeDataLoc$'ukbb_stroke.phe')
echo $prevStroke_GBP # 0.0182402
prevT2D_GBP=$(findprevalence "GBP" $believeDataLoc$'ukbb_t2d.phe')
echo $prevT2D_GBP # 0.0439115
prevAngina_GBP=$(findprevalence "GBP" $believeDataLoc$'ukbb_angina.phe')
echo $prevAngina_GBP # 0.0344536
prevChol_GBP=$(findprevalence "GBP" $believeDataLoc$'ukbb_cholesterol.phe')
echo $prevChol_GBP # 0.320554

prevMI_UKSAS=$(findprevalence "UK_SEA" $believeDataLoc$'ukbb_MI.phe')
echo $prevMI_UKSAS # 0.0381152
prevStroke_UKSAS=$(findprevalence "UK_SEA" $believeDataLoc$'ukbb_stroke.phe')
echo $prevStroke_UKSAS # 0.0168067
prevT2D_UKSAS=$(findprevalence "UK_SEA" $believeDataLoc$'ukbb_t2d.phe')
echo $prevT2D_UKSAS # 0.184574
prevAngina_UKSAS=$(findprevalence "UK_SEA" $believeDataLoc$'ukbb_angina.phe')
echo $prevAngina_UKSAS # 0.0660264
prevChol_UKSAS=$(findprevalence "UK_SEA" $believeDataLoc$'ukbb_cholesterol.phe')
echo $prevChol_UKSAS # 0.203782

# BLV 
prevMI_BLV=$(findprevalence "BELIEVE" $believeDataLoc$'MI.phe')
echo $prevMI_BLV #   0.0217587
prevStroke_BLV=$(findprevalence "BELIEVE" $believeDataLoc$'stroke.phe')
echo $prevStroke_BLV # 0.0204082
prevT2D_BLV=$(findprevalence "BELIEVE" $believeDataLoc$'t2d.phe')
echo $prevT2D_BLV # 0.176921
prevAngina_BLV=$(findprevalence "BELIEVE" $believeDataLoc$'angina.phe')
echo $prevAngina_BLV # 0.0484694
prevChol_BLV=$(findprevalence "BELIEVE" $believeDataLoc$'highCholesterol_all.phe')
echo $prevChol_BLV # 0.138205


###########################################################################
# X. 6664 Scenarios
###########################################################################

# need to create a HM3 version of the datasets, as the full datasets take 325GB of RAM, and that is not possible:
arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL_6664 --make-bed --out '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 --extract '$w_hm3_snplist$' --allow-extra-chr --allow-no-sex'
$plink $arguments
# --extract: 710505 variants remaining.



# repeat everything 2x, one for GCTA, and a second time via PopKin
outputArray=( '_popKin' '' )
popKinArray=( '1' '0' )
arraylength=${#popKinArray[@]}
c=1

for (( c=1; c<${arraylength}+1; c++ )); do
isPopKin=${popKinArray[$c-1]} 
out=${outputArray[$c-1]} 

resultLoc=$believeResultsLoc$'UKBB_BELIEVE_PC46664_bin'$out$'/'
mkdir -p $resultLoc
mkdir -p $resultLoc$'finals/'

arguments='/home/mk907/scripts/R/phenoDecompositionReal_bin.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'highCholesterol_6664_pop_bin '$resultLoc$'highCholesterol_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA '$prevChol_BLV$' '$prevChol_UKSAS
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal_bin.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'MI_6664_pop_bin '$resultLoc$'MI_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA '$prevMI_BLV$' '$prevMI_UKSAS
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal_bin.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'stroke_6664_pop_bin '$resultLoc$'stroke_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA '$prevStroke_BLV$' '$prevStroke_UKSAS
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal_bin.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'t2d_6664_pop_bin '$resultLoc$'t2d_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA '$prevT2D_BLV$' '$prevT2D_UKSAS
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal_bin.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'angina_6664_pop_bin '$resultLoc$'angina_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA '$prevAngina_BLV$' '$prevAngina_UKSAS
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"


# continuous phenos (as I am testing the LRT for the GxE component explicitly (edu Y wasn't significant as it was the first variance component that wasn't significant which is expected if it was 0)
arguments='/home/mk907/scripts/R/phenoDecompositionReal2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'height_6664_pop '$resultLoc$'height_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'bmi_6664_pop '$resultLoc$'bmi_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'edu_6664_pop '$resultLoc$'edu_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'sbp_6664_pop '$resultLoc$'sbp_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"

arguments='/home/mk907/scripts/R/phenoDecompositionReal2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'dbp_6664_pop '$resultLoc$'dbp_res /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"


#done # end of loop for isPopKin
# wait until above finishes on the cluster, then execute below:


myarray=( "height" "bmi" "highCholesterol" "edu" "angina" "t2d" "stroke" "MI" "sbp" "dbp" ) # 
myarray_bin=( "0" "0" "1" "0" "1" "1" "1" "1" "0" "0" ) # 
arraylength=${#myarray[@]}
j=8

for (( j=1; j<${arraylength}+1; j++ )); do
trait=${myarray[$j-1]}
isBin=${myarray_bin[$j-1]}
echo $trait
resultLoc=$believeResultsLoc$'UKBB_BELIEVE_PC46664_bin'$out$'/'
# get the p-value for VGxE ,h2, and subpop h2s, and add them into the final results file
VGxE_Pval=$(awk -v binary="$isBin" '{ if(binary == "1") {line=20} else {line =14}  if (FNR == line) {print $2}}' $resultLoc$trait$'_res.hsq')
h2_Pval=$(awk -v binary="$isBin" '{ if(binary == "1") {line=13} else {line ="10"}  if (FNR == line) {print $2}}' $resultLoc$trait$'_res_env.hsq')
BELIEVE_Pval=$(awk -v binary="$isBin" '{ if(binary == "1") {line=13} else {line =10}  if (FNR == line) {print $2}}' $resultLoc$trait$'_res_sub1.hsq')
UKBB_Pval=$(awk -v binary="$isBin" '{ if(binary == "1") {line=13} else {line =10}  if (FNR == line) {print $2}}' $resultLoc$trait$'_res_sub2.hsq')
echo $VGxE_Pval
echo $h2_Pval
echo $BELIEVE_Pval
echo $UKBB_Pval
# for PopKin MI, the basic h2 did not converge as it was 0, so we use a p of 1
#VGxE_Pval=1

awk -v VGxE_Pval="$VGxE_Pval" -v h2_Pval="$h2_Pval" -v BELIEVE_Pval="$BELIEVE_Pval" -v UKBB_Pval="$UKBB_Pval" '{print $0} END { print "VGxE_Pval: "VGxE_Pval"\nh2_Pval: "h2_Pval"\nBELIEVE_Pval: "BELIEVE_Pval"\nUKBB_Pval: "UKBB_Pval}' $resultLoc$trait$'_res_results.txt' >  $resultLoc$'finals/'$trait$'_res_results.txt'


done # end of loop for results aggregation

# generate copy & paste table
arguments='/home/mk907/scripts/R/realResTableGenerator.R '$resultLoc$'finals/ '$resultLoc$'finals/table'$out
Rscript $arguments


done # end of loop for isPopKin

###########################################################################
# XI. Simulations: Basic Sims and also AFR sims check what happens to simulations if we include genetically divergent group, 
###########################################################################

########################
# Prepare data for Simulations:  

# 1. check for more outgroups with larger Fst, eg AFR: just for the simulations, add the AFR indis, to see how bad it gets

numTotal=$(awk '{if($2 == "UK_SEA")numTotal++} END {print numTotal}' $believeScratchLoc$'labelsPC4_p2_mistrt2dangsbpdbp')
numAFR=$(wc -l < $believeDataLoc$'AFR_keeplist')
echo $numAFR
remainderToPick=$(awk -v numAFR=$numAFR -v numTotal=$numTotal 'BEGIN {print numAFR-numTotal}')
echo $remainderToPick # -668


awk -v remainderToPick=$remainderToPick '{if (FNR > remainderToPick) {print $0}}' $believeDataLoc$'AFR_keeplist' > $believeDataLoc$'AFR_keeplist_6664'
wc -l $believeDataLoc$'AFR_keeplist_6664' # 5996

# we actually have fewer africans than UK_SEA, so to keep pop sizes the same, we need to downsize the latter
numAFR_existing=$(wc -l < $believeDataLoc$'AFR_keeplist_6664')

# extract these indis into a separate dataset that includes just the BELIEVE and the AFR ones
# re-run simulations of BELIEVE vs AFR to see if the larger Fst difference makes any difference
arguments=' --memory 43000 --bfile '$UKBB_PLINK1$'45K_IMPUTED --keep '$believeDataLoc$'AFR_keeplist_6664 --extract '$believeRawLoc$'UKBB_BELIEVE_FULL_6664.bim --make-bed --out '$believeDataLoc$'AFR_FULL --allow-extra-chr --allow-no-sex'
$plink $arguments

awk -v numAFR_existing=$numAFR_existing '{ if($4 == "BELIEVE" && count < numAFR_existing ) {print $1"\t"$1; count++}}' $believeDataLoc$'height_6664_pop' > $believeDataLoc$'BELIEVE_keeplist_FULL'
head $believeDataLoc$'BELIEVE_keeplist_FULL'
wc -l $believeDataLoc$'BELIEVE_keeplist_FULL'

arguments=' --memory 43000 --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL_6664 --keep '$believeDataLoc$'BELIEVE_keeplist_FULL --make-bed --out '$believeDataLoc$'BELIEVE_only_FULL --allow-extra-chr --allow-no-sex'
$plink $arguments

wc -l $believeDataLoc$'BELIEVE_only_FULL.fam' # 5997


# now merge AFR with the BELIEVE 6664
plinkFileList=$believeScratchLoc$'plinkFileList'
rm -rf $plinkFileList
echo $believeDataLoc$'BELIEVE_only_FULL' >> ${plinkFileList}
echo $believeDataLoc$'AFR_FULL' >> ${plinkFileList}

arguments=' --memory '$PLINKRAM$' --merge-list '$plinkFileList$' --make-bed --out '$believeRawLoc$'AFR_BELIEVE_FULL --allow-extra-chr --allow-no-sex'
$plink $arguments

# create poplabels for AFR
awk '{ print $1"\t"$1"\t-1\tAFR"}' $believeDataLoc$'AFR_FULL.fam' > $believeScratchLoc$'AFR_lab_FULL'
awk '{ if ($4 =="BELIEVE") print $0}' $believeDataLoc$'height_6664_pop' > $believeScratchLoc$'BELIEVE_lab_FULL'
head $believeScratchLoc$'BELIEVE_lab_FULL'
cat $believeScratchLoc$'AFR_lab_FULL' $believeScratchLoc$'BELIEVE_lab_FULL' > $believeScratchLoc$'AFR_BELIEVE_raw_lab'

awk 'FNR == NR { test[ $1 ] = $0; next; } FNR <= NR { if(FNR == 1) {print "FID\tIID\tpheno\tpop" } print test[$2] }
' $believeScratchLoc$'AFR_BELIEVE_raw_lab' $believeRawLoc$'AFR_BELIEVE_FULL.fam' > $believeDataLoc$'AFR_BELIEVE_6664_pop'
head $believeDataLoc$'AFR_BELIEVE_6664_pop'
wc -l $believeDataLoc$'AFR_BELIEVE_6664_pop'

# also use AFR QC-d list of SNPs
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $w_hm3_snplist $believeScratchLoc$'45K_IMPUTED_QC_4.snplist' > $believeScratchLoc$'HM3_AFR_QC_SNPs'

# need to create a HM3 version of the datasets, as the full datasets take 325GB of RAM, and that is not possible:
arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'AFR_BELIEVE_FULL --make-bed --out '$believeRawLoc$'AFR_BELIEVE_HM3_6664 --extract '$w_hm3_snplist$' --allow-extra-chr --allow-no-sex'
$plink $arguments
# 710505 variants and 11992 people pass filters and QC. 


########################
# Run Simulations:  (do it broken into multiple jobs, to avoid overruning limit)


# repeat everything 2x, one for GCTA, and a second time via PopKin
outputArray=( '_popKin' '' )
popKinArray=( '1' '0' )
numSumSplits=( '21' '4' ) # how many ways we split the file (+1, as we loop <)
arraylength=${#popKinArray[@]}
c=1

for (( c=1; c<${arraylength}+1; c++ )); do
isPopKin=${popKinArray[$c-1]} 
out=${outputArray[$c-1]} 
looplength=${numSumSplits[$c-1]} 



# submit simulations for UK AFR vs BELIEVE  
resultLoc=$believeResultsLoc$'UKBB_BELIEVE_PC46664_bin'$out$'/'
resultOutLoc=$resultLoc$'Sims/out/'
mkdir -p $resultLoc$'finals/'

outfile=$resultLoc$'Sims_BELIEVE_AFR_1_VE_ENV_GXE.txt'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'AFR_BELIEVE_HM3_6664 '$believeDataLoc$'AFR_BELIEVE_6664_pop '$resultLoc$'Sims_BELIEVE_AFR_1 /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' 1 7 0 BELIEVE AFR'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"
fi

outfile=$resultLoc$'Sims_BELIEVE_AFR_2_VE_ENV_GXE.txt'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'AFR_BELIEVE_HM3_6664 '$believeDataLoc$'AFR_BELIEVE_6664_pop '$resultLoc$'Sims_BELIEVE_AFR_2 /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' 8 15 0 BELIEVE AFR'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"
fi

outfile=$resultLoc$'Sims_BELIEVE_AFR_3_VE_ENV_GXE.txt'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'AFR_BELIEVE_HM3_6664 '$believeDataLoc$'AFR_BELIEVE_6664_pop '$resultLoc$'Sims_BELIEVE_AFR_3 /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' 16 20 0 BELIEVE AFR'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"
fi



# for some reason, the "BELIEVE"    "UK_SEA"  PopKin runs were very slow, only completing 1 simulation iteration before timing out
# to solve this, resubmit these as 20 jobs, one for each 
#(technically, I've had to run this 1x and fail, as we are using the precalculated kinship matrix for the 20 parts)
if [ $isPopKin == "1" ] ; then 
echo "popkin, we submit in 20 parts" 
for (( k=1; k<$looplength; k++ )); do
echo $k

outfile=$resultLoc$'Sims_HM3_'$k$'_VE_ENV_GXE.txt'
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile


if [ $k == "1" ] ; then 
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'height_6664_pop '$resultLoc$'Sims_HM3_'$k$' /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' '$k$' '$k$' 0 BELIEVE UK_SEA'
#echo "FIRST ONE"
dumm=1
else
# WAIT FOR THE FIRST ONE TO FINISH, so that kinship matrix will be done 
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'height_6664_pop '$resultLoc$'Sims_HM3_'$k$' /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' '$k$' '$k$' '$resultLoc$'Sims_HM3_1 BELIEVE UK_SEA'
#echo "NOT FIRST ONE"
fi # end of if first one

sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"
fi # end of if file exists

done # end of looplength


else
echo "not popkin"
# submit simulations for UK SAS and BELIEVE
outfile=$resultLoc$'Sims_HM3_1_VE_ENV_GXE.txt'
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'height_6664_pop '$resultLoc$'Sims_HM3_1 /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' 1 7 0  BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"
fi


outfile=$resultLoc$'Sims_HM3_2_VE_ENV_GXE.txt'
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'height_6664_pop '$resultLoc$'Sims_HM3_2 /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' 8 15 0 BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"
fi

outfile=$resultLoc$'Sims_HM3_3_VE_ENV_GXE.txt'
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments='/home/mk907/scripts/R/phenoDecompositionSim2.R '$believeRawLoc$'UKBB_BELIEVE_HM3_6664 '$believeDataLoc$'height_6664_pop '$resultLoc$'Sims_HM3_3 /home/mk907/scripts/R/pheDecomp_functs_lowmem.R '$isPopKin$' 16 20 0 BELIEVE UK_SEA'
sbatch --mem 160000 --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "Rscript $arguments"
fi

fi # end of if not popkin


#done # end of loop for isPopKin

# wait until above finishes on the cluster, then execute below:


# plot all 
resultOutLoc=$resultLoc$'out/'
mkdir -p $resultOutLoc

simsresarray=( 'BELIEVE_AFR' 'HM3' ) # 
arraylength2=${#simsresarray[@]}
j=1


# concat the results of the above
for (( j=1; j<${arraylength2}+1; j++ )); do
sim=${simsresarray[$j-1]}

if [ $isPopKin == "1" ] && [ $sim == "HM3" ] ; then 
echo "popkin BELIEVE UK SEA"

# generate input for the function out of all the 20 splits
args1=()
args2=()
args3=()
args4=()
args5=()
args6=()
args7=()
for (( k=1; k<$looplength; k++ )); do
args1+=($resultLoc$'Sims_'$sim$'_'$k$'_VG.txt')
args2+=($resultLoc$'Sims_'$sim$'_'$k$'_VG_Env.txt')
args3+=($resultLoc$'Sims_'$sim$'_'$k$'_VE_ENV_GXE.txt')
args4+=($resultLoc$'Sims_'$sim$'_'$k$'_VG_sub.txt')
args5+=($resultLoc$'Sims_'$sim$'_'$k$'_VG_Env_sub.txt')
args6+=($resultLoc$'Sims_'$sim$'_'$k$'_VE_ENV_GXE_sub.txt')
args7+=($resultLoc$'Sims_'$sim$'_'$k$'_GV.txt')
done 

concatResLoop $resultLoc$'Sims_'$sim$'_VG.txt' "${args1[@]}"
concatResLoop $resultLoc$'Sims_'$sim$'_VG_Env.txt' "${args2[@]}"
concatResLoop $resultLoc$'Sims_'$sim$'_VE_ENV_GXE.txt' "${args3[@]}"
concatResLoop $resultLoc$'Sims_'$sim$'_VG_sub.txt' "${args4[@]}"
concatResLoop $resultLoc$'Sims_'$sim$'_VG_Env_sub.txt' "${args5[@]}"
concatResLoop $resultLoc$'Sims_'$sim$'_VE_ENV_GXE_sub.txt' "${args6[@]}"
concatResLoop $resultLoc$'Sims_'$sim$'_GV.txt' "${args7[@]}"



else 
concatRes $resultLoc$'Sims_'$sim$'_1_VG.txt' $resultLoc$'Sims_'$sim$'_2_VG.txt' $resultLoc$'Sims_'$sim$'_3_VG.txt' $resultLoc$'Sims_'$sim$'_VG.txt'
concatRes $resultLoc$'Sims_'$sim$'_1_VG_Env.txt' $resultLoc$'Sims_'$sim$'_2_VG_Env.txt' $resultLoc$'Sims_'$sim$'_3_VG_Env.txt' $resultLoc$'Sims_'$sim$'_VG_Env.txt'
concatRes $resultLoc$'Sims_'$sim$'_1_VE_ENV_GXE.txt' $resultLoc$'Sims_'$sim$'_2_VE_ENV_GXE.txt' $resultLoc$'Sims_'$sim$'_3_VE_ENV_GXE.txt' $resultLoc$'Sims_'$sim$'_VE_ENV_GXE.txt'
concatRes $resultLoc$'Sims_'$sim$'_1_VG_sub.txt' $resultLoc$'Sims_'$sim$'_2_VG_sub.txt' $resultLoc$'Sims_'$sim$'_3_VG_sub.txt' $resultLoc$'Sims_'$sim$'_VG_sub.txt'
concatRes $resultLoc$'Sims_'$sim$'_1_VG_Env_sub.txt' $resultLoc$'Sims_'$sim$'_2_VG_Env_sub.txt' $resultLoc$'Sims_'$sim$'_3_VG_Env_sub.txt' $resultLoc$'Sims_'$sim$'_VG_Env_sub.txt'
concatRes $resultLoc$'Sims_'$sim$'_1_VE_ENV_GXE_sub.txt' $resultLoc$'Sims_'$sim$'_2_VE_ENV_GXE_sub.txt' $resultLoc$'Sims_'$sim$'_3_VE_ENV_GXE_sub.txt' $resultLoc$'Sims_'$sim$'_VE_ENV_GXE_sub.txt'
concatRes $resultLoc$'Sims_'$sim$'_1_GV.txt' $resultLoc$'Sims_'$sim$'_2_GV.txt' $resultLoc$'Sims_'$sim$'_3_GV.txt' $resultLoc$'Sims_'$sim$'_GV.txt'

fi # end of if not popkin





### PLot Sims
cp $resultLoc$'Sims_'$sim$'_VG.txt' $resultLoc$'finals/Sims_'$sim$'res_VG.txt'
# I mixed up the VG_ENV and VG_ENV_GXE models in the R script, as I wrote out the GxE results with filename "_VG_Env, so swap it back"
cp $resultLoc$'Sims_'$sim$'_VG_Env.txt' $resultLoc$'finals/Sims_'$sim$'res_VG_ENV_GXE.txt'
cp $resultLoc$'Sims_'$sim$'_VE_ENV_GXE.txt' $resultLoc$'finals/Sims_'$sim$'res_VG_Env.txt'

cp $resultLoc$'Sims_'$sim$'_VG_sub.txt' $resultLoc$'finals/Sims_'$sim$'res_VG_sub.txt'
cp $resultLoc$'Sims_'$sim$'_VG_Env_sub.txt' $resultLoc$'finals/Sims_'$sim$'res_VG_Env_sub.txt'
cp $resultLoc$'Sims_'$sim$'_VE_ENV_GXE_sub.txt' $resultLoc$'finals/Sims_'$sim$'res_VG_ENV_GXE_sub.txt'

cp $resultLoc$'Sims_'$sim$'_GV.txt' $resultLoc$'finals/Sims_'$sim$'res_GV.txt'

# find Fst, col 2 of the last line
fst=$(awk '{if (FNR == 3) {print $2}}' $resultLoc$'Sims_'$sim$'_1_popFSt_prive.txt')

# truth scenario VG Env + GxE
arguments='/home/mk907/scripts/R/pheDecomp_plotter.R '$resultLoc$'finals/Sims_'$sim$'res_VG_ENV_GXE.txt '$resultLoc$'finals/Sims_'$sim$'res_VG_ENV_GXE_sub.txt '$resultLoc$'finals/Sims_'$sim$'res_GV.txt '$resultOutLoc$'/Sims_'$sim$'VG_ENV_GXE 0.33 1 0.33 '$fst
Rscript $arguments

# truth scenario VG + Env
arguments='/home/mk907/scripts/R/pheDecomp_plotter.R '$resultLoc$'finals/Sims_'$sim$'res_VG_Env.txt '$resultLoc$'finals/Sims_'$sim$'res_VG_Env_sub.txt '$resultLoc$'finals/Sims_'$sim$'res_GV.txt '$resultOutLoc$'/Sims_'$sim$'VG_ENV 0.33 1 0 '$fst
Rscript $arguments

# truth scenario VG
arguments='/home/mk907/scripts/R/pheDecomp_plotter.R '$resultLoc$'finals/Sims_'$sim$'res_VG.txt '$resultLoc$'finals/Sims_'$sim$'res_VG_sub.txt '$resultLoc$'finals/Sims_'$sim$'res_GV.txt '$resultOutLoc$'/Sims_'$sim$'VG 0.33 0 0 '$fst
Rscript $arguments

done # end of loop for results aggregation

done # end of loop for isPopKin



###########################################################################
# XII. Fst Analyses
###########################################################################


######################
# re-run Fst analysis using the relaxed PC4 subsets
awk 'FNR == NR { 
file1[ $1 ] = $2; next; } 
FNR <= NR {  
if(FNR == 1) {print $0}
else if( $1 in file1 ) { print $1"\tBLV_"file1[$1]"\t"$3} 
else if( $2 == "BELIEVE") { print $1"\tBLV_Other\t"$3} 
else {print $0} }
' $believeResultsLoc$'ethnic_groups' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4' > $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_bang'
head $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_bang'
wc -l $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_bang'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_bang'
# QCFAIL 169
# pop 1
# BLV_Other 1533
# BLV_Bihari 1309
# BLV_Bangali 30870
# GBP 6297
# UK_SEA 7553

mkdir -p $believeResultsLoc$'UKBB_BELIEVE_MATCHED_BANGPC4/'



#########################
# Recalcualte Fsts on console (used for the Table of how well we have done at matching the pops
R
library(Matrix)
library(genio)
#library(GxEMM)
library(popkin)
library(bigstatsr)
library(bigsnpr)
##################

# Relaxed (Top 4PCs)
plinkGenotypesLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/raw/BELIEVE_MATCHED_EURPC4"
phenoLabelsLoc = "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/BELIEVE_MATCHED_ALLPC4/labelsPC4_bang"
outLoc = "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED_BANGPC4/pop_"

# also do the precisely-matched TOP 20 PC ones
plinkGenotypesLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/raw/BELIEVE_MATCHED_EUR"
phenoLabelsLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED/labels_bang"
outLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED_BANGPC4/pop_1500_"



# Relaxed (Top 4PCs) and # also do the precisely matched TOP 20 PC ones
plinks=c("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/raw/BELIEVE_MATCHED_EURPC4", "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/raw/BELIEVE_MATCHED_EUR")
phens=c("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/BELIEVE_MATCHED_ALLPC4/labelsPC4_bang", "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED/labels_bang")
outs=c("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED_BANGPC4/pop_", "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED_BANGPC4/pop_precise_")


pops = c("BLV_Other","BLV_Bihari","BLV_Bangali","UK_SEA","GBP", "QCFAIL")
simFunctLoc= "/home/mk907/scripts/R/pheDecomp_functs_lowmem.R"
source(simFunctLoc)

# write loop to do both relaxed and precise, and also WC and PopKin
i=2
for (i in 1:2) {
plinkGenotypesLoc=plinks[i]
phenoLabelsLoc = phens[i]
outLoc = outs[i]

print(plinkGenotypesLoc)
print(phenoLabelsLoc)
print(outLoc)


 
phenoLabels=read.table(phenoLabelsLoc, header=T)
#phenoLabels$indi_index=1:nrow(phenoLabels)
tempDirLoc=paste0(outLoc,"temp/")
unlink(tempDirLoc, recursive = TRUE)
  
(obj.bed <- bed(paste0(plinkGenotypesLoc,".bed") )) # this is only needed for the Fst calculation 
  

snpData = snp_readBed2(paste0(plinkGenotypesLoc,".bed"), # ind.row = ind.row,
                         backingfile = tempDirLoc, # sub_bed(bedfile)
                         ncores = 1) # using 16 cores causes segmentation error for some reason
						 
snpAttached <- snp_attach(snpData)	
fam2=snpAttached$fam
fam2$indi_index=1:nrow(fam2)
fam_pheno = merge (fam2, phenoLabels,by.x="sample.ID", by.y="IID" )

calcFst_WC(obj.bed, fam_pheno, outLoc, pops )

unlink(tempDirLoc, recursive = TRUE)
######################
# run PopKin version too by loading in the Kinship matrix? (may have to create the kinship matrix!

# as popkin needs to load in the whole genotype matrix and calculate an n x n kinship matrix from that we don't want to load in too many indis

maxNum = 1538 # only load in a certain max number, this is arbitrary just to stop OOM
pops = names(table(fam_pheno$pop))
popsNumTotals = table(fam_pheno$pop)
fam_pheno$index = 1:nrow(fam_pheno)
numtoPick=c()
set.seed(42)
all_picked_indices = c()
  for (j in 1:length(pops)) {
  if(popsNumTotals[j] <= maxNum) { numtoPick= c(numtoPick, popsNumTotals[j]) 
  } else {numtoPick= c(numtoPick,maxNum)}

  popNumber = as.numeric(numtoPick[j]) 
  print(popNumber)
  
  all_indices = fam_pheno$index[which(fam_pheno$pop ==pops[j])]
  picked_indices = sample(all_indices,popNumber)
  # want to sort these back into original order

   picked_indices = picked_indices[order(picked_indices)]
   all_picked_indices = c(all_picked_indices,picked_indices )
   }
   

fam_pheno_subset = fam_pheno[all_picked_indices,]

snpData_subset = snp_readBed2(paste0(plinkGenotypesLoc,".bed"), # ind.row = ind.row,
                         backingfile = tempDirLoc, ind.row = all_picked_indices,
                         ncores = 1) # if 16, then this causes segmentation errors
						 
snpAttached <- snp_attach(snpData_subset)	
X_all = snpAttached$genotypes
X_all = X_all[]

K_standard =  calcKinshipViaPopKin(X_all,fam_pheno_subset, paste0(outLoc, "_popkin") )
K_standard = K_standard /2 # for the Fst calculations we need to undo the 2 * correction...
calcFst_popKin(K_standard, fam_pheno_subset, paste0(outLoc, "_popkinFst"))

##########################
unlink(tempDirLoc, recursive = TRUE)


}

quit()

# the finished fst files (for popkin)
/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED_BANGPC4/pop_precise__popkinFst_popFSt_prive.txt
/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED_BANGPC4/pop_precise__popFSt_prive.txt

###########################################################################
# XIII. GWAS, Additive and GxE (REGENIE)
###########################################################################


# !!! HERE
# for the GWAS, use all indis, even those with missing phenos, and code them up as NA, and let REGENIE handle that
# - for REGENIE, need to work off from this list:
$believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4_IMPUTED'
$believeRawLoc$'UKBB_BELIEVE_FULL.bim'
$believeRawLoc$'UKBB_BELIEVE_FULL.fam'
$believeDataLoc$'MI_all2'
$believeDataLoc$'stroke_all2'
$believeDataLoc$'t2d_all2'
$believeDataLoc$'angina_all2'
$believeDataLoc$'sbp_all2'
$believeDataLoc$'dbp_all2'
$believeDataLoc$'highCholesterol_all2'
$believeDataLoc$'height_all2'
$believeDataLoc$'bmi_all2'
$believeDataLoc$'edu_all2'

# split plink datas into UK_SAS and BELIEVE
awk 'FNR == NR { 
file1[ $2 ] = $2; next; } 
FNR <= NR {  
if(FNR == 1) {print $0}
else if( $1 in file1 ) { print $0} }
' $believeRawLoc$'UKBB_BELIEVE_FULL.fam' $believeResultsLoc$'BELIEVE_MATCHED_ALLPC4/labelsPC4' > $believeScratchLoc$'labelsPC4_re'
head $believeScratchLoc$'labelsPC4_re'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_re'
# QCFAIL 169
# pop 1
# GBP 6297
# BELIEVE 33287
# UK_SEA 7553


awk 'FNR == NR { 
file1[ $1 ] = $1; next; } 
FNR <= NR {  
if(FNR == 1) {print $0}
else if( $1 in file1 ) { print $0} }
' $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'labelsPC4_re'  > $believeScratchLoc$'labelsPC4_reg'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $believeScratchLoc$'labelsPC4_reg'
# QCFAIL 169
# pop 1
# GBP 6296
# BELIEVE 31776
# UK_SEA 7552


wc -l $believeScratchLoc$'labelsPC4_reg'




#



#cp $believeScratchLoc$'blv_keeplist' $believeScratchLoc$'blv_keeplist_old'
#cp $believeScratchLoc$'ukbb_keeplist' $believeScratchLoc$'ukbb_keeplist_old'
#cp $believeScratchLoc$'ukbb_blv_keeplist' $believeScratchLoc$'ukbb_blv_keeplist_old'

awk '{if ($2 == "BELIEVE") {print $1"\t"$1}}' $believeScratchLoc$'labelsPC4_reg' > $believeScratchLoc$'blv_keeplist'
head $believeScratchLoc$'blv_keeplist'
wc -l $believeScratchLoc$'blv_keeplist' # 33287

awk '{if ($2 == "UK_SEA") {print $1"\t"$1}}' $believeScratchLoc$'labelsPC4_reg' > $believeScratchLoc$'ukbb_keeplist'
head $believeScratchLoc$'ukbb_keeplist'
wc -l $believeScratchLoc$'ukbb_keeplist' # 7553



cat $believeScratchLoc$'blv_keeplist' $believeScratchLoc$'ukbb_keeplist' > $believeScratchLoc$'ukbb_blv_keeplist'

# prepare covars in plink format for binary phenos
awk '{if (FNR == 1) {print "FID\tIID\tage\tage2\tsex"} print $1"\t"$0}' $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' > $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS'
head $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS'
wc -l $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS'

#################################

# prepare  phenos in regenie format ( IE all phenos in a single file: FID IID Y1 Y2 )
binarPheno_all $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_MI_binary.phe' $believeDataLoc$'MI_all2_bin'
binarPheno_all $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_stroke_binary.phe' $believeDataLoc$'stroke_all2_bin'
binarPheno_all $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_t2d_binary.phe' $believeDataLoc$'t2d_all2_bin'
binarPheno_all $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_angina_binary.phe' $believeDataLoc$'angina_all2_bin'
binarPheno_all $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' $believeScratchLoc$'ukbb_blv_cholesterol_binary.phe' $believeDataLoc$'highCholesterol_all2_bin'


myarray=( "height" "bmi" "highCholesterol" "edu" "angina" "t2d" "stroke" "MI" "sbp" "dbp" ) # 
myarray_bin=( "" "" "_bin" "" "_bin" "_bin" "_bin" "_bin" "" "" ) # 
arraylength=${#myarray[@]}

j=3

awk '{if(FNR == 1) {print "FID\tIID"} print $0}' $believeScratchLoc$'ukbb_blv_keeplist' > $believeDataLoc$'binary_traits.txt'
awk '{if(FNR == 1) {print "FID\tIID"} print $0}' $believeScratchLoc$'ukbb_blv_keeplist' > $believeDataLoc$'quantitative_traits.txt'

for (( j=1; j<${arraylength}+1; j++ )); do

trait=${myarray[$j-1]}
isBin=${myarray_bin[$j-1]}
echo $trait

if [[ "$isBin" == "_bin" ]]; then
echo "binary pheno"
traitFile=$believeDataLoc$'binary_traits.txt'
else
echo "continuoius pheno"
traitFile=$believeDataLoc$'quantitative_traits.txt'
fi

cp $traitFile $traitFile$'_temp'
awk -v trait=$trait 'FNR == NR { 
file1[ $1 ] = $2; next; } 
FNR <= NR {  
if(FNR == 1) {print $0"\t"trait}
else if( $2 in file1 ) { print $0"\t"file1[$2]}  else { print $0"\tNA"}}
' $believeDataLoc$trait$'_all2'$isBin $traitFile$'_temp' > $traitFile
# Missing values must be coded as NA: https://rgcgithub.github.io/regenie/options/#phenotype-file-format

done

head $believeDataLoc$'binary_traits.txt'
wc -l $believeDataLoc$'binary_traits.txt' # 39329

head $believeDataLoc$'quantitative_traits.txt'
wc -l $believeDataLoc$'quantitative_traits.txt' # 39329



# prepare covars for everyone
awk 'FNR == NR { 
file1[ $2 ] = $2; next; } 
FNR <= NR {  
if(FNR == 1) {print "FID\tIID\tage\tage2\tsex"}
if( $1 in file1 ) { print $1"\t"$0  } }
' $believeScratchLoc$'ukbb_blv_keeplist' $believeScratchLoc$'ukbb_blv_covars_harmonized.phe' > $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS'
head $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS'
wc -l $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS'


# add the Country as a categorical covariate for the binary covariates file
awk 'FNR == NR { 
file1[ $2 ] = $2; next; } 
FNR <= NR {  
if(FNR == 1) {print $0"\tcountry"}
else if( $2 in file1 ) { print $0"\tBangladesh"}  else { print $0"\tUK"}}
' $believeScratchLoc$'blv_keeplist' $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS' > $believeScratchLoc$'ukbb_blv_covars_binary'
head $believeScratchLoc$'ukbb_blv_covars_binary'
wc -l $believeScratchLoc$'ukbb_blv_covars_binary'

# create a covariates file for the continuos trait GxE
awk '{ print $1"\t"$2"\t"$6}' $believeScratchLoc$'ukbb_blv_covars_binary' > $believeScratchLoc$'ukbb_blv_covars_quant'
head $believeScratchLoc$'ukbb_blv_covars_quant'
wc -l $believeScratchLoc$'ukbb_blv_covars_quant' # 39329


head $believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS'

# covariates file to use for the separate cohort binary GWAS, (as this will NOT have the country covariate)
$believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS' # --catCovarList sex

# covariates file to use for the combined cohort binary GWAS (both Additive and GxE)
$believeScratchLoc$'ukbb_blv_covars_binary' # --catCovarList sex,country

# covariates file to use for the separate cohort quantitative GWAS
# NONE: as the phenos were pre-adjusted, and when fit separately, we won't have a country covariate

# covariates file to use for the combined cohort quant GWAS (both Additive and GxE)
$believeScratchLoc$'ukbb_blv_covars_quant' # --catCovarList country



############################
# A) Step I of REGENIE 
# doc says we should use 500K SNPs that are directly genotyed: https://rgcgithub.github.io/regenie/faq/#:~:text=We%20recommend%20to%20use%20a,in%20the%20level%201%20models.
# we only have around 120K of those, so we will use 380K random HM3 SNps with MAF > 5%

# make sure that these criteria work for BOTH cohorts separately, so we fit it twice
arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL --keep '$believeScratchLoc$'blv_keeplist --maf 0.05 --mac 100 --write-snplist --out '$believeScratchLoc$'BLV_HM3_GWAS --extract '$w_hm3_snplist$' --allow-extra-chr --allow-no-sex'
$plink $arguments
arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL --keep '$believeScratchLoc$'ukbb_keeplist --maf 0.05 --mac 100 --write-snplist --out '$believeScratchLoc$'UKB_HM3_GWAS --extract '$w_hm3_snplist$' --allow-extra-chr --allow-no-sex'
$plink $arguments


arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL --keep '$believeScratchLoc$'blv_keeplist --maf 0.05 --mac 100 --write-snplist --out '$believeScratchLoc$'BLV_PCA_GWAS --extract '$believeRawLoc$'2564_ALL.bim --allow-extra-chr --allow-no-sex'
$plink $arguments
arguments=' --memory '$PLINKRAM$' --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL --keep '$believeScratchLoc$'ukbb_keeplist --maf 0.05 --mac 100 --write-snplist --out '$believeScratchLoc$'UKB_PCA_GWAS --extract '$believeRawLoc$'2564_ALL.bim --allow-extra-chr --allow-no-sex'
$plink $arguments



# and then intersect the snp lists
awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeScratchLoc$'BLV_HM3_GWAS.snplist' $believeScratchLoc$'UKB_HM3_GWAS.snplist' > $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON.snplist'
head $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON.snplist'
wc -l $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON.snplist' # 624421

awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $believeScratchLoc$'BLV_PCA_GWAS.snplist' $believeScratchLoc$'UKB_PCA_GWAS.snplist' > $believeScratchLoc$'UKB_BELIEVE_PCA_GWAS_QCCOMMON.snplist'
head $believeScratchLoc$'UKB_BELIEVE_PCA_GWAS_QCCOMMON.snplist'
wc -l $believeScratchLoc$'UKB_BELIEVE_PCA_GWAS_QCCOMMON.snplist' # 92067


# remove all PCA SNPs from HM3 (as those will be selected regardless)
awk  'FNR == NR { file1[ $1 ] = $1; next; } 
FNR <= NR {if( $1 in file1 == 0) { print $0}  }
' $believeScratchLoc$'UKB_BELIEVE_PCA_GWAS_QCCOMMON.snplist' $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON.snplist' >  $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON_NOPCA.snplist'
head $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON_NOPCA.snplist'
wc -l $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON_NOPCA.snplist' # 535716

# now randomly select the right number from the above list to make up a total of 500K SNPs
numPCASNPs=$(wc -l < $believeScratchLoc$'UKB_BELIEVE_PCA_GWAS_QCCOMMON.snplist')
echo $numPCASNPs
numHM3SNPs=$(wc -l < $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON_NOPCA.snplist')
echo $numHM3SNPs

# num to pick: 500K - numPCASNPs
numNeeded=500000
numHM3ToPick=$(awk -v numPCASNPs=$numPCASNPs -v numNeeded=$numNeeded 'BEGIN {print numNeeded-numPCASNPs}')
echo $numHM3ToPick # 407933

# now pick this many random SNPs from the UKB_BELIEVE_HM3_GWAS_QCCOMMON_NOPCA
arguments='/home/mk907/scripts/R/randomNums.R '$numHM3SNPs$' 42 '$believeScratchLoc$'numHM3ToPick'
Rscript $arguments

# pick the first numHM3ToPick
awk  -v numHM3ToPick=$numHM3ToPick  '{if (FNR <= numHM3ToPick) print $0}' $believeScratchLoc$'numHM3ToPick' > $believeScratchLoc$'numHM3ToPick_indices'
wc -l $believeScratchLoc$'numHM3ToPick_indices'


awk 'FNR == NR {  test[ NR ] = $0; next; } FNR <= NR { print test[$1]  }
' $believeScratchLoc$'UKB_BELIEVE_HM3_GWAS_QCCOMMON_NOPCA.snplist' $believeScratchLoc$'numHM3ToPick_indices' > $believeScratchLoc$'numHM3ToPick_SNPs'
head $believeScratchLoc$'numHM3ToPick_SNPs'
wc -l $believeScratchLoc$'numHM3ToPick_SNPs' 

# combine these lists to get the final 500K
cat $believeScratchLoc$'UKB_BELIEVE_PCA_GWAS_QCCOMMON.snplist' $believeScratchLoc$'numHM3ToPick_SNPs' > $believeScratchLoc$'UKB_BELIEVE_STEP1.snplist'
head $believeScratchLoc$'UKB_BELIEVE_STEP1.snplist'
wc -l $believeScratchLoc$'UKB_BELIEVE_STEP1.snplist' # 500000

################
# REGENIE:  
resultLoc=$believeResultsLoc$'GWAS/'
mkdir -p $resultLoc$'bin/'
mkdir -p $resultLoc$'quant/'
mkdir -p $resultLoc$'temp/'

###############
# BELIEVE GWAS

# binary
arguments=' --step 1 --bt --cc12 --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'binary_traits.txt --covarFile '$believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS --catCovarList sex --lowmem-prefix '$resultLoc$'temp/ --extract '$believeScratchLoc$'UKB_BELIEVE_STEP1.snplist --keep '$believeScratchLoc$'blv_keeplist --out '$resultLoc$'bin/BLV_BT_step1 --bsize 1000 --threads '$SLURM_CPUS_ON_NODE
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"


   # - 'highCholesterol': 4457 cases and 27318 controls
   # - 'angina': 1452 cases and 30289 controls
   # - 't2d': 5844 cases and 25897 controls
   # - 'stroke': 767 cases and 30974 controls
   # - 'MI': 723 cases and 31018 controls

# quant #  --lowmem
arguments=' --step 1 --qt --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'quantitative_traits.txt --lowmem-prefix '$resultLoc$'temp/ --extract '$believeScratchLoc$'UKB_BELIEVE_STEP1.snplist --keep '$believeScratchLoc$'blv_keeplist --out '$resultLoc$'quant/BLV_QT_step1 --bsize 1000 --threads '$SLURM_CPUS_ON_NODE
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"



###############
# UKB GWAS  --lowmem

# binary
arguments=' --step 1 --bt --cc12 --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'binary_traits.txt --covarFile '$believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS --catCovarList sex --lowmem-prefix '$resultLoc$'temp/ --extract '$believeScratchLoc$'UKB_BELIEVE_STEP1.snplist --keep '$believeScratchLoc$'ukbb_keeplist --out '$resultLoc$'bin/UKB_BT_step1 --bsize 1000 --threads '$SLURM_CPUS_ON_NODE
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"

# quant
arguments=' --step 1 --qt --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'quantitative_traits.txt --lowmem-prefix '$resultLoc$'temp/ --extract '$believeScratchLoc$'UKB_BELIEVE_STEP1.snplist --keep '$believeScratchLoc$'ukbb_keeplist --out '$resultLoc$'quant/UKB_QT_step1 --bsize 1000 --threads '$SLURM_CPUS_ON_NODE
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"


###############
# COMBINED (UKB+BELIEVE) + GXE
# Additive SNP effect, when country is in the model "ADD-CONDTL"
#--interaction country[Bangladesh] (e.g. --interaction VAR[BASELEVEL]).

# binary
arguments=' --step 1 --bt --cc12 --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'binary_traits.txt --covarFile '$believeScratchLoc$'ukbb_blv_covars_binary --catCovarList sex,country --lowmem-prefix '$resultLoc$'temp/ --extract '$believeScratchLoc$'UKB_BELIEVE_STEP1.snplist --keep '$believeScratchLoc$'ukbb_blv_keeplist --out '$resultLoc$'bin/UKB_BLV_BT_step1 --bsize 1000 --threads '$SLURM_CPUS_ON_NODE
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"


# quant
arguments=' --step 1 --qt --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'quantitative_traits.txt --covarFile '$believeScratchLoc$'ukbb_blv_covars_quant --catCovarList country --lowmem-prefix '$resultLoc$'temp/ --extract '$believeScratchLoc$'UKB_BELIEVE_STEP1.snplist --keep '$believeScratchLoc$'ukbb_blv_keeplist --out '$resultLoc$'quant/UKB_BLV_QT_step1 --bsize 1000 --threads '$SLURM_CPUS_ON_NODE
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"


# !!! HERE

############################
# B) Step II of REGENIE (using the full set of SNPs)
# bsize in step2 isn't that important, using 400 is fine: https://github.com/rgcgithub/regenie/issues/239

###############
resultLoc=$believeResultsLoc$'GWAS/'
chr=21
for ((chr=1; chr<=$numChroms; chr++)); do

###############
# BELIEVE GWAS


# binary
outfile=$resultLoc$'bin/BLV_BT_step2_'$chr$'_highCholesterol.regenie'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments=' --chr '$chr$' --step 2 --bt --cc12 --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'binary_traits.txt --covarFile '$believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS --catCovarList sex --lowmem-prefix '$resultLoc$'temp/ --minMAC 100 --firth --approx --pThresh 0.01 --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeScratchLoc$'blv_keeplist --pred '$resultLoc$'bin/BLV_BT_step1_pred.list --out '$resultLoc$'bin/BLV_BT_step2_'$chr$' --bsize 400 --lowmem'
#$regenie $arguments

sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"

fi


# quant
outfile=$resultLoc$'quant/BLV_QT_step2_'$chr$'_height.regenie'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments=' --chr '$chr$' --step 2 --qt --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'quantitative_traits.txt --lowmem-prefix '$resultLoc$'temp/ --minMAC 100 --firth --approx --pThresh 0.01 --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeScratchLoc$'blv_keeplist --pred '$resultLoc$'quant/BLV_QT_step1_pred.list --out '$resultLoc$'quant/BLV_QT_step2_'$chr$' --bsize 400 --lowmem'
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"
fi


###############
# UKB GWAS

# binary
outfile=$resultLoc$'bin/UKB_BT_step2_'$chr$'_highCholesterol.regenie'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments=' --chr '$chr$' --step 2 --bt --cc12 --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'binary_traits.txt --covarFile '$believeScratchLoc$'ukbb_blv_covars_harmonized_GWAS --catCovarList sex --lowmem-prefix '$resultLoc$'temp/ --minMAC 100 --firth --approx --pThresh 0.01 --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeScratchLoc$'ukbb_keeplist --pred '$resultLoc$'bin/UKB_BT_step1_pred.list --out '$resultLoc$'bin/UKB_BT_step2_'$chr$' --bsize 400 --lowmem'
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"
fi



# quant
outfile=$resultLoc$'quant/UKB_QT_step2_'$chr$'_height.regenie'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments=' --chr '$chr$' --step 2 --qt --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'quantitative_traits.txt --lowmem-prefix '$resultLoc$'temp/ --minMAC 100 --firth --approx --pThresh 0.01 --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeScratchLoc$'ukbb_keeplist --pred '$resultLoc$'quant/UKB_QT_step1_pred.list --out '$resultLoc$'quant/UKB_QT_step2_'$chr$' --bsize 400 --lowmem'
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"
fi

###############
# COMBINED (UKB+BELIEVE) + GXE

# binary
outfile=$resultLoc$'bin/UKB_BLV_BT_step2_'$chr$'_highCholesterol.regenie'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments=' --chr '$chr$' --step 2 --bt --cc12 --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'binary_traits.txt  --interaction country[Bangladesh] --covarFile '$believeScratchLoc$'ukbb_blv_covars_binary --catCovarList sex,country --lowmem-prefix '$resultLoc$'temp/ --minMAC 100 --firth --approx --pThresh 0.01 --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeScratchLoc$'ukbb_blv_keeplist --pred '$resultLoc$'bin/UKB_BLV_BT_step1_pred.list --out '$resultLoc$'bin/UKB_BLV_BT_step2_'$chr$' --bsize 400 --lowmem'
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"
fi

# quant
outfile=$resultLoc$'quant/UKB_BLV_QT_step2_'$chr$'_height.regenie'
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
echo "submitting: "$outfile
arguments=' --chr '$chr$' --step 2 --qt --bed '$believeRawLoc$'UKBB_BELIEVE_FULL --phenoFile '$believeDataLoc$'quantitative_traits.txt --interaction country[Bangladesh] --covarFile '$believeScratchLoc$'ukbb_blv_covars_quant --catCovarList country --lowmem-prefix '$resultLoc$'temp/ --minMAC 100 --firth --approx --pThresh 0.01 --threads '$SLURM_CPUS_ON_NODE$' --keep '$believeScratchLoc$'ukbb_blv_keeplist --pred '$resultLoc$'quant/UKB_BLV_QT_step1_pred.list --out '$resultLoc$'quant/UKB_BLV_QT_step2_'$chr$' --bsize 400 --lowmem'
#$regenie $arguments
sbatch --mem ${PLINKRAM_SMALL} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$regenie $arguments"
fi

done # end of looping chroms

##################


# consider --apply-rint for quant phenos, maybe? 

################################
# concat all results across all chroms
# BLV
 # * case-control counts for each trait:
   # - 'highCholesterol': 4457 cases and 27318 controls
   # - 'angina': 1452 cases and 30289 controls
   # - 't2d': 5844 cases and 25897 controls
   # - 'stroke': 767 cases and 30974 controls
   # - 'MI': 723 cases and 31018 controls
# * number of observations for each trait:
   # - 'height': 31763 observations
   # - 'bmi': 24791 observations
   # - 'edu': 27566 observations
   # - 'sbp': 31741 observations
   # - 'dbp': 31741 observations
   
# # UKB
# * case-control counts for each trait:
   # - 'highCholesterol': 1446 cases and 5718 controls
   # - 'angina': 513 cases and 6861 controls
   # - 't2d': 1389 cases and 5946 controls
   # - 'stroke': 125 cases and 7159 controls
   # - 'MI': 305 cases and 6989 controls
 # * number of observations for each trait:
   # - 'height': 7386 observations
   # - 'bmi': 7377 observations
   # - 'edu': 3873 observations
   # - 'sbp': 7525 observations
   # - 'dbp': 7525 observations


#awk: fatal: cannot open file `/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/GWAS/bin/UKB_BT_step2_23.log' for reading (No such file or directory)



resultLoc=$believeResultsLoc$'GWAS/'
myarray=( "height" "bmi" "highCholesterol" "edu" "angina" "t2d" "stroke" "MI" "sbp" "dbp" ) # 
myarray_bin=( "quant" "quant" "bin" "quant" "bin" "bin" "bin" "bin" "quant" "quant" ) # 
filenam_bin=( "QT" "QT" "BT" "QT" "BT" "BT" "BT" "BT" "QT" "QT" )
arraylength=${#myarray[@]}

GWAS_type=( "BLV" "UKB" "UKB_BLV" )
GWAS_arraylength=${#GWAS_type[@]}

k=2
j=3
for (( k=1; k<${GWAS_arraylength}+1; k++ )); do
GWAS=${GWAS_type[$k-1]}
echo $GWAS


for (( j=1; j<${arraylength}+1; j++ )); do
trait=${myarray[$j-1]}
fileBin=${filenam_bin[$j-1]}
isBin=${myarray_bin[$j-1]}

outfile=$resultLoc$isBin$'/'$GWAS$'_'$fileBin$'_step2_'
echo $trait " / "$outfile


# extract the sample sizes from the logfile
awk 'BEGIN {start=0; end =0} {
if( index($0, "* number of observations for each trait:") != 0  || index($0, "* case-control counts for each trait:") != 0 ) {start=1}
else if(index($0, "* LOCO predictions") != 0) {end=1}
if(start == 1 && end != 1) {print $0}
}' $outfile$'1.log' > $outfile$'_sampleSizes.txt'
#head $outfile$'_sampleSizes.txt'


chr=21
# concat all chrom results without the header
# create header
head -n 1 $outfile$'1_'$trait$'.regenie' > $resultLoc$GWAS$'_'$trait$'.regenie'
for ((chr=1; chr<=$numChroms; chr++)); do
chrFile=$outfile$chr$'_'$trait$'.regenie'

if [ -s "$chrFile" ] ; then
dumm=1
awk '{if(FNR > 1) print $0}' $chrFile >> $resultLoc$GWAS$'_'$trait$'.regenie'
else
echo "DID NOT EXIST: "$chrFile
fi




done # end of chrom loop

# split the additive and GxE results for the UKB_BLVresults
if [[ "$GWAS" == "UKB_BLV" ]]; then
echo "plotting GxE"
awk '{if(FNR == 1 || $8 == "ADD-INT_SNPxcountry=UK") {print $0}}' $resultLoc$GWAS$'_'$trait$'.regenie' > $resultLoc$GWAS$'_'$trait$'_GXE.regenie'
awk '{if(FNR == 1 || $8 == "ADD-CONDTL") {print $0}}' $resultLoc$GWAS$'_'$trait$'.regenie' > $resultLoc$GWAS$'_'$trait$'_Add.regenie'

produce_manhattan_REGENIE $resultLoc$GWAS$'_'$trait$'_GXE.regenie' $resultLoc$GWAS$'_'$trait$'_GXE' $trait$'_(GxE)' 7.30103
produce_manhattan_REGENIE $resultLoc$GWAS$'_'$trait$'_Add.regenie' $resultLoc$GWAS$'_'$trait $trait 7.30103

# find GxE sig ones for eah trait
awk '{if (FNR == 1) {print "CHR\tSNP\tBP\tP"}  else if( $12 >= 7.30103) {print $1"\t"$3"\t"$2"\t"10^(-$12)}}' $resultLoc$GWAS$'_'$trait$'_GXE.regenie' > $resultLoc$GWAS$'_'$trait$'_GXE_sig'

# Ld clump those, but only if there are any (we have at least 1 line, due to the header)
numSig=$(wc -l < $resultLoc$GWAS$'_'$trait$'_GXE_sig')

if [ "$numSig" -gt 1 ]; then
rm -rf $resultLoc$GWAS$'_'$trait$'_GXE_Clump.clumped'
arguments=' --bfile '$believeRawLoc$'UKBB_BELIEVE_FULL --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$'  --keep '$believeScratchLoc$'ukbb_blv_keeplist --clump '$resultLoc$GWAS$'_'$trait$'_GXE_sig --clump-kb 500 --clump-r2 0.10  --out '$resultLoc$GWAS$'_'$trait$'_GXE_Clump  --allow-extra-chr --snps-only --bp-space 1 --allow-no-sex'
$plink $arguments
rm -rf $resultLoc$GWAS$'_'$trait$'_GXE_Clump.log'
rm -rf $resultLoc$GWAS$'_'$trait$'_GXE_Clump.nosex'
dumm=3
fi



else 
echo "plotting Additive"
produce_manhattan_REGENIE $resultLoc$GWAS$'_'$trait$'.regenie' $resultLoc$GWAS$'_'$trait $trait 7.30103
fi

done # end of trait loop
done # end of GWAS loop




#####################################

# Find out genetic correlation between UK_SAS and BELIEVE via LDSC
resultLoc=$believeResultsLoc$'GWAS/'
echo -e "\nrG" > $resultLoc$'_RGALL.csv'
myarray=( "height" "bmi" "highCholesterol" "edu" "angina" "t2d" "stroke" "MI" "sbp" "dbp" ) # 
myarray_bin=( "" "" "_bin" "" "_bin" "_bin" "_bin" "_bin" "" "" ) # 
arraylength=${#myarray[@]}
j=6
for (( j=1; j<${arraylength}+1; j++ )); do
trait=${myarray[$j-1]}
isBin=${myarray_bin[$j-1]}
echo $trait

# grab the association files
ukbbres=$resultLoc$'UKB_'$trait$'.regenie'
blvres=$resultLoc$'BLV_'$trait$'.regenie'

# need to convert log10p to just P for munge...

# we also need to remove N from the data columns, otherwise Munge will use that from the sumstats, instead of the one I supplied
# https://github.com/bulik/ldsc/blob/master/munge_sumstats.py  , line 347: "if 'N' not in dat.columns:"
# this is not right, as for case/control studies, we should use effective sampel size, and not the N reported by REGENIE
awk '{line=$1" "$2" "$3" "$4" "$5" "$6" "$8" "$9" "$10" "$11" "$12" "$13; if (FNR == 1) {print line" P"} else {print line" "10^(-$12)}}' $ukbbres > $ukbbres$'_P'
awk '{line=$1" "$2" "$3" "$4" "$5" "$6" "$8" "$9" "$10" "$11" "$12" "$13; if (FNR == 1) {print line" P"} else {print line" "10^(-$12)}}' $blvres > $blvres$'_P'

# awk '{if (FNR == 1) {print $0" P"} else {print $0" "10^(-$12)}}' $blvres > $blvres$'_P'

# find sample sizes used for the GWAS
if [[ "$isBin" == "_bin" ]]; then
echo "binary pheno"

N_eff_UKBB=$(GetN_eff_regenie $trait $believeDataLoc$'binary_traits.txt' $believeScratchLoc$'ukbb_keeplist')
N_eff_BLV=$(GetN_eff_regenie $trait $believeDataLoc$'binary_traits.txt' $believeScratchLoc$'blv_keeplist')


else
echo "continuoius pheno"
# read sample size from the assoc file and just average it, as there it is simply the nunber used
N_eff_UKBB=$(awk '{ if(FNR > 1) {sum += $7; n++ } } END { printf("%0.f", sum / n) ; }' $ukbbres)
N_eff_BLV=$(awk '{ if(FNR > 1) {sum += $7; n++ } } END { printf("%0.f", sum / n) ; }' $blvres)
fi

# REGENIE FORMAT: 
#   1    2     3    4       5       6   7   8   9   10  11      12   13
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
#  reference allele (allele 0), alternative allele (allele 1),

# for munge:
#A1 : reference allele
#A2 : alternative allele
# regenie records the total sample size, not the effective sample size, which is different for case/control phenos

# generate LDSC format sumstats
$munge_sumstats --sumstats $ukbbres$'_P' --snp ID --a1 ALLELE0 --a2 ALLELE1 --frq A1FREQ --N $N_eff_UKBB --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $resultLoc$trait$'_UKBB'
$munge_sumstats --sumstats $blvres$'_P' --snp ID --a1 ALLELE0 --a2 ALLELE1 --frq A1FREQ --N $N_eff_BLV --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $resultLoc$trait$'_BLV'

# get rG between BLV and UKBB
$lsdc --rg $resultLoc$trait$'_UKBB.sumstats.gz',$resultLoc$trait$'_BLV.sumstats.gz' --ref-ld $csa_ref --w-ld $csa_ref --out $resultLoc$trait$'_rG'
# needed to rename UKBB.CSA.rsid.l2.ldscore.gz
# to needed to rename UKBB.CSA.l2.ldscore.gz
# as the old "UKBB.CSA.l2.ldscore.gz" did NOT use the rsids, so it could not be merged with the summary stats

# extract rG, and add asterisk if the p value is significant (0.05)
rGtext=$(awk 'BEGIN {finaltext =""} {
if( index($0, "Genetic Correlation:") != 0 ) {finaltext=$3 }
else if( index($0, "P:") != 0 && $2 < 0.05) { finaltext = finaltext"*" }
} END {print finaltext }' $resultLoc$trait$'_rG.log')

cp $resultLoc$'_RGALL.csv' $resultLoc$'_RGALL.csv_temp'
awk -v trait=$trait -v rGtext=$rGtext '{if (FNR == 1) {print $0"\t"trait} else {print $0"\t"rGtext}}' $resultLoc$'_RGALL.txt_temp' > $resultLoc$'_RGALL.csv'


done

# MI, T2D and stroke have nans as sample sizes were too low 
#nano $resultLoc$'MI_rG.log'
# nano $resultLoc$'edu_rG.log'  -0.1271, and p val is not significant 
#####################################
# Aggregate the GWAS results into a table 

# create header for results file with signature:
#        GxE                     BELIEVE+UKB                       BELIEVE                   UKB
#SNP     p       BETA            p       BETA       MAF            p       BETA       MAF    p       BETA       MAF
echo -e "\t\tGxE\t\tBELIEVE+UKB\t\t\tBELIEVE\t\t\tUKB\t\t\t" > $resultLoc$GWAS$'_results.csv'
echo -e "trait\tSNP\tp\tbeta\tp\tbeta\tmaf\tp\tbeta\tmaf\tp\tbeta\tmaf" >> $resultLoc$GWAS$'_results.csv'
head $resultLoc$GWAS$'_results.csv'

resultLoc=$believeResultsLoc$'GWAS/'
j=6
for (( j=1; j<${arraylength}+1; j++ )); do
trait=${myarray[$j-1]}
fileBin=${filenam_bin[$j-1]}
isBin=${myarray_bin[$j-1]}
GWAS="UKB_BLV"
outfile=$resultLoc$isBin$'/'$GWAS$'_'$fileBin$'_step2_'
echo $trait " / "$outfile

# find out if there is at least 1 significant result
numSig=$(wc -l < $resultLoc$GWAS$'_'$trait$'_GXE_sig')
if [ "$numSig" -gt 1 ]; then


#head $resultLoc$GWAS$'_'$trait$'_results.txt'

# loop the list of significant clumped SNPs

counter=0
while read p; do # go through each SNP

counter=$((counter+1))

if [ $counter != "1" ]; then  # skip header
p=$(echo $p | sed -e 's/\r//g')  # the last column contains a newline character, which would fuck up awk string comparison and create unreadable files etc
IFS=' ' read -r -a arrIN <<< "$p"  # https://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
targetVar=${arrIN[2]} 
targetVar_P=${arrIN[4]}

#targetVar="rs236554"        
#targetVar_P="2.01e-08"

# we also get blank lines in the clumped file for some reason, so we want to skip those too
if [[ $targetVar != "" ]]; then 
# add this to a running total 
echo -e $targetVar$"\t"$targetVar_P

# get same SNPs from the additive associations for Combined, and just within BELIEVE and UKB cohorts too

# REGENIE FORMAT: 
#   1    2     3    4       5       6   7   8   9   10  11      12   13
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
# write found values into an intermediate file
awk -v targetVar=$targetVar '{if($3 == targetVar) {print "targetVar_GXEBeta="$9; }}' $resultLoc$GWAS$'_'$trait$'_GXE.regenie' > $resultLoc$GWAS$'_'$trait$'_GXE.regenie_intermediate'
head $resultLoc$GWAS$'_'$trait$'_GXE.regenie_intermediate'


awk -v targetVar=$targetVar '{if($3 == targetVar) {print "targetVar_MAF="$6; print "targetVar_AddP="10^(-$12); print "targetVar_AddBeta="$9; }}' $resultLoc$GWAS$'_'$trait$'_Add.regenie' > $resultLoc$GWAS$'_'$trait$'_ADD_Clump.clumped_intermediate'
awk -v targetVar=$targetVar '{if($3 == targetVar) {print "targetVarUKB_MAF="$6; print "targetVarUKB_AddP="10^(-$12); print "targetVarUKB_AddBeta="$9; }}' $resultLoc$'UKB_'$trait$'.regenie' > $resultLoc$'UKB_'$trait$'.regenie_intermediate'
awk -v targetVar=$targetVar '{if($3 == targetVar) {print "targetVarBLV_MAF="$6; print "targetVarBLV_AddP="10^(-$12); print "targetVarBLV_AddBeta="$9; }}' $resultLoc$'BLV_'$trait$'.regenie' > $resultLoc$'BLV_'$trait$'.regenie_intermediate'

# load it back into bash
source $resultLoc$GWAS$'_'$trait$'_GXE.regenie_intermediate'
source $resultLoc$GWAS$'_'$trait$'_ADD_Clump.clumped_intermediate'
source $resultLoc$'UKB_'$trait$'.regenie_intermediate'
source $resultLoc$'BLV_'$trait$'.regenie_intermediate'

echo $targetVar_GXEBeta

echo $targetVar_MAF
echo $targetVar_AddP
echo $targetVar_AddBeta


echo $targetVarUKB_MAF
echo $targetVarUKB_AddP
echo $targetVarUKB_AddBeta

echo $targetVarBLV_MAF
echo $targetVarBLV_AddP
echo $targetVarBLV_AddBeta

# add to final table
traitText=$trait
if [ $counter != "2" ]; then  # only the first SNP found gets the "trait" mentioned, the rest will be just a tab
traitText=""
fi # end of if first lead SNP
echo $traitText

echo -e $traitText$"\t"$targetVar$"\t"$targetVar_P$"\t"$targetVar_GXEBeta$"\t"$targetVar_AddP$"\t"$targetVar_AddBeta$"\t"$targetVar_MAF$"\t"$targetVarBLV_AddP$"\t"$targetVarBLV_AddBeta$"\t"$targetVarBLV_MAF$"\t"$targetVarUKB_AddP$"\t"$targetVarUKB_AddBeta$"\t"$targetVarUKB_MAF >> $resultLoc$GWAS$'_results.csv'

fi # end of if not blank line
fi # end of if header

done <$resultLoc$GWAS$'_'$trait$'_GXE_Clump.clumped'
fi # end of if there were results




done # end of trait loop



############################
resultLoc=$believeResultsLoc$'GWAS/'
head $resultLoc$GWAS$'_t2d_GXE_Clump.clumped'
#CHR    F             SNP         BP        P    TOTAL   NSIG    S05    S01   S001  S0001    SP2
#  17    1        rs236554   68248158   2.01e-08        1      0      0      0      0      1 rs236555(1)

# REGENIE FORMAT: 
#   1    2     3    4       5       6   7   8   9   10  11      12   13
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
  
# locus zoom: https://statgen.github.io/localzoom/
# extract area, and also add the P column
gxeres=$resultLoc$GWAS$'_t2d_GXE.regenie'
mkdir -p $believeResultsLoc$'GWAS/locuszoom/'

# 100K +- around target
#68248158-100000
#68248158+100000
awk '{line=$0; if (FNR == 1) {print line" P"} else if ($1 == 17 && $2 >= 68148158 && $2 <= 68348158) {print line" "10^(-$12)}}' $gxeres > $gxeres$'_P'
head $gxeres$'_P'
wc -l $gxeres$'_P' # 312

grep "rs236554"  $gxeres$'_P'


# need to also liftover SNPs from hg19 to hg38, as SAS is only available in hg38
awk '{if(FNR > 1) print $1"\t"$3"\t0\t"$2}' $gxeres$'_P' > $gxeres$'_P_dummy.map'
head $gxeres$'_P_dummy.map'
awk 'BEGIN {line="dumm1\tdumm1\tdumm3\tdumm4\t1\t0"} {if(FNR > 1) {line = line"\t"$5"\t"$4}} END {print line}' $gxeres$'_P' > $gxeres$'_P_dummy.ped'
head $gxeres$'_P_dummy.ped'
wc -l $gxeres$'_P_dummy.ped'
awk '{print NF}' $gxeres$'_P_dummy.ped' | sort -nu | tail -n 1  # 628 cols, so, (628 -6)/2, so 311 SNPs, which is correct


cd $believeRawLoc$'liftover/'
python2 liftOverPlink.py -m $gxeres$'_P_dummy.map' -p $gxeres$'_P_dummy.ped' -o $believeResultsLoc$'GWAS/locuszoom/hg38' -c hg19ToHg38.over.chain.gz -e ./liftOver

# now map the new coords onto the sumstats file
awk 'FNR == NR {  test[ $2 ] = $4; next; } FNR <= NR { if(FNR == 1) {print } else if ($3 in test) {  $2=test[$3]; print $0 } }
' OFS=" " FS="\t" $believeResultsLoc$'GWAS/locuszoom/hg38.map' FS=" " $gxeres$'_P' > $gxeres$'_P_38'
head $gxeres$'_P_38'
wc -l $gxeres$'_P_38'

## add plink to the path temporarily: https://codereview.stackexchange.com/questions/88236/unix-shell-function-for-adding-directories-to-path
# add2path() {
  # if ! echo $PATH | egrep "(^|:)$1(:|\$)" > /dev/null ; then
    # if [[ $2 = "front" ]]; then
      # PATH="$1:$PATH"
    # else
      # PATH="$PATH:$1"
    # fi
    # export PATH
  # fi
# }

# add2path "/home/mk907/software/plink/" "front"



locuszoom="/rds/project/jmmh2/rds-jmmh2-public_databases/locuszoom/bin/locuszoom"
$locuszoom --source 1000G_Nov2014 --build hg38 --pop SAS --refsnp rs236554 --metal $gxeres$'_P_38' --delim space title="UKB-BELIEVE-GXE" --markercol ID --pvalcol P --chr 17 --start 70152209 --end 70349437 --no-date --plotonly --prefix=BELIEVE --rundir $believeResultsLoc$'GWAS/locuszoom/'


grep "rs236554" $gxeres$'_P_38'


