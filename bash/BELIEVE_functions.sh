
# Reusable command line functions for the BELIEVE project
# screen -r -D 245278.pts-76.login-e-10 # 2X sims

# screen -r -D 33240.pts-9.login-e-15

believeBaseLoc='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/'
believeScratchLoc=$believeBaseLoc$'scratch/'
believeDataLoc=$believeBaseLoc$'data/'
believeRawLoc=$believeDataLoc$'raw/'
believeResultsLoc=$believeBaseLoc$'results/'

believeGxS='/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/genotype/genomewide/CAMBRIDGE-BELIEVE_Freeze_One.GxS'
believeWES='/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/genotype/wes/CAMBRIDGE-BELIEVE_Freeze_One.norm'

phenoBaseLoc='/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/'
plink="/home/mk907/software/plink/plink"
plink2_alpha='/home/mk907/software/plink2_alpha/plink2'
relatedBaseLoc='/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/reference_files/genetic/relationships/'

phenoMap="/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/BELIEVEdata_GeneticIDMapping_P5031_20210826.csv"
phenoRaw="/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/BELIEVEdata_P5031_20210826.csv"
phenoRaw2="/rds/project/jmmh2/rds-jmmh2-post_qc_data/believe/phenotype/BELIEVEdata_P5031_20220310.csv"
PLINKRAM=97000
PLINKRAM_SMALL=12000
SLURM_CPUS_ON_NODE=32
NCORES_LDPred2=16

hapmap3_b37bim='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim'
plink2='/home/mk907/software/plink2/plink2'
ancestryL="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/"
QCfile='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/reference_files/ukb_sqc_v2.txt' # eid matched to our app


SNPQCfile='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/reference_files/ukb_snp_qc.txt' # 
#   1            2                 3                 4          5         6         7         8       9       10             11            12            13            14          15              16           17             18              19          20             21          22            23             24            25            26            27           28            29             30            31            32            33          34             35            36           37              38           39             40
# rs_id affymetrix_snp_id affymetrix_probeset_id chromosome position allele1_ref allele2_alt strand array Batch_b001_qc Batch_b002_qc Batch_b003_qc Batch_b004_qc Batch_b005_qc Batch_b006_qc Batch_b007_qc Batch_b008_qc Batch_b009_qc Batch_b010_qc Batch_b011_qc Batch_b012_qc Batch_b013_qc Batch_b014_qc Batch_b015_qc Batch_b016_qc Batch_b017_qc Batch_b018_qc Batch_b019_qc Batch_b020_qc Batch_b021_qc Batch_b022_qc Batch_b023_qc Batch_b024_qc Batch_b025_qc Batch_b026_qc Batch_b027_qc Batch_b028_qc Batch_b029_qc Batch_b030_qc Batch_b031_qc Batch_b032_qc Batch_b033_qc Batch_b034_qc Batch_b035_qc Batch_b036_qc Batch_b037_qc Batch_b038_qc Batch_b039_qc Batch_b040_qc Batch_b041_qc Batch_b042_qc Batch_b043_qc Batch_b044_qc Batch_b045_qc Batch_b046_qc Batch_b047_qc Batch_b048_qc Batch_b049_qc Batch_b050_qc Batch_b051_qc Batch_b052_qc Batch_b053_qc Batch_b054_qc Batch_b055_qc Batch_b056_qc Batch_b057_qc Batch_b058_qc Batch_b059_qc Batch_b060_qc Batch_b061_qc Batch_b062_qc Batch_b063_qc Batch_b064_qc Batch_b065_qc Batch_b066_qc Batch_b067_qc Batch_b068_qc Batch_b069_qc Batch_b070_qc Batch_b071_qc Batch_b072_qc Batch_b073_qc Batch_b074_qc Batch_b075_qc Batch_b076_qc Batch_b077_qc Batch_b078_qc Batch_b079_qc Batch_b080_qc Batch_b081_qc Batch_b082_qc Batch_b083_qc Batch_b084_qc Batch_b085_qc Batch_b086_qc Batch_b087_qc Batch_b088_qc Batch_b089_qc Batch_b090_qc Batch_b091_qc Batch_b092_qc Batch_b093_qc Batch_b094_qc Batch_b095_qc UKBiLEVEAX_b1_qc UKBiLEVEAX_b2_qc UKBiLEVEAX_b3_qc UKBiLEVEAX_b4_qc UKBiLEVEAX_b5_qc UKBiLEVEAX_b6_qc UKBiLEVEAX_b7_qc UKBiLEVEAX_b8_qc UKBiLEVEAX_b9_qc UKBiLEVEAX_b10_qc UKBiLEVEAX_b11_qc in_HetMiss in_Relatedness in_PCA PC1_loading PC2_loading PC3_loading PC4_loading PC5_loading PC6_loading PC7_loading PC8_loading PC9_loading PC10_loading PC11_loading PC12_loading PC13_loading PC14_loading PC15_loading PC16_loading PC17_loading PC18_loading PC9_loading PC20_loading PC21_loading PC22_loading PC23_loading PC24_loading PC25_loading PC26_loading PC27_loading PC28_loading PC9_loading PC30_loading PC31_loading PC32_loading PC33_loading PC34_loading PC35_loading PC36_loading PC37_loading PC38_loading PC9_loading PC40_loading in_Phasing_Input
#head -n 1 $SNPQCfile
#awk '{if(FNR ==1) print $116 }' $SNPQCfile # in_HetMiss

rawLocV3='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/'
WITHDRAWALS=$rawLocV3$'WITHDRAWALS' # people no longer in the UKBB
SEX_DISC_or_LQ=$rawLocV3$'SEX_DISC_or_LQ'  # Sex discordant or those with low quality genotypes ( too many missing or excess heterozygosity)
RELATEDS=$rawLocV3$'RELATEDS' # pairs of people who are too closely related 


# signature:
#  1              2                 3          4            5         6            7               8                   9                        10                  11              12                 13                     14                 15                  16                    17                       18                              19                              20                                           21            22                              23                                24                 25      26      27      28      29      30      31      32      33       34      35      36   37       38      39      40
# eid     genotyping.array        Batch   Plate.Name      Well    Cluster.CR      dQC     Internal.Pico..ng.uL.   Submitted.Gender        Inferred.Gender       X.intensity     Y.intensity     Submitted.Plate.Name    Submitted.Well  sample.qc.missing.rate  heterozygosity  heterozygosity.pc.corrected  het.missing.outliers     putative.sex.chromosome.aneuploidy      in.kinship.table        excluded.from.kinship.inference    excess.relatives        in.white.British.ancestry.subset      used.in.pca.calculation     PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10    PC11    PC12 PC13     PC14    PC15    PC16    PC17    PC18    PC19    PC20    PC21    PC22    PC23    PC24    PC25    PC26    PC27    PC28    PC29    PC30    PC31 PC32     PC33    PC34    PC35    PC36    PC37    PC38    PC39    PC40    in.Phasing.Input.chr1_22        in.Phasing.Input.chrX   in.Phasing.Input.chrXY
baseLoc='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/'
dataLoc=$baseLoc$'data/'
UKBB_PLINK1=$dataLoc$'UKBB_PLINK1/'

numChroms=22
module load R/4.0.3


module load miniconda/2
source activate ldsc
# module load ceuadmin/regenie/3.2.5
module load ceuadmin/regenie/3.2.5.3
#regenie="/usr/local/Cluster-Apps/ceuadmin/regenie/3.2.5/regenie"
regenie="/usr/local/Cluster-Apps/ceuadmin/regenie/3.2.5.3/regenie"

csa_ref="/home/mk907/software/ldsc/panUKBB/UKBB.ALL.ldscore/UKBB.CSA"
ldscDir='/home/mk907/software/ldsc/'
w_hm3_snplist='/home/mk907/software/ldsc/w_hm3.snplist'
eur_w_ld_chr='/home/mk907/software/ldsc/eur_w_ld_chr/'
lsdc=$ldscDir$'ldsc/ldsc.py'
munge_sumstats=$ldscDir$'ldsc/munge_sumstats.py'



# get the UKBB PCs
rawLocV3='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/'
WITHDRAWALS=$rawLocV3$'WITHDRAWALS' # people no longer in the UKBB
SEX_DISC_or_LQ=$rawLocV3$'SEX_DISC_or_LQ'  # Sex discordant or those with low quality genotypes ( too many missing or excess heterozygosity)
RELATEDS=$rawLocV3$'RELATEDS' # pairs of people who are too closely related 


mkdir -p $believeScratchLoc
mkdir -p $believeRawLoc
mkdir -p $believeResultsLoc
###########################################################################

#################################



# 3) exclude ambiguous alleles:  A/T or C/G from a .bim formatted file
function remove_ambiguous_alleles_bim { 
sumstatsLoc=$1
#  1     2     3    4     5
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
awk ' 
FNR <= NR {
if (FNR == 1) {print $0 }
else { 
if( toupper($6) == "A" && toupper($5) == "T" || toupper($5) == "A" && toupper($6) == "T" || toupper($6) == "G" && toupper($5) == "C" || toupper($5) == "G" && toupper($6) == "C") {}
else {print $0}
}
}
' $sumstatsLoc > $sumstatsLoc$'_noAmbiguousAlleles'

head $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc
}
export -f remove_ambiguous_alleles_bim # this makes local functions executable when bsubbed


function AddEthnicLabels { 
labels=$1
PCAres=$2
outLocc=$3

awk 'FNR == NR { 
file1[ $2 ] = $3;next; } 
FNR <= NR {  
if(FNR == 1) {print "IID\tgroup\tcolor"}

if( $1 in file1) { 
part1=$1"\t"file1[$1]
if ( file1[$1] == "Ashkenazi") { color = "#88CCEE"}
else if ( file1[$1] == "Caribbean") { color = "#CC6677"}
else if ( file1[$1] == "China"    ) { color = "#DDCC77"} 
else if ( file1[$1] == "India"  ) { color = "#117733"}   
else if ( file1[$1] == "Iran"  ) { color = "#332288"}    
else if ( file1[$1] == "Italy"  ) { color = "#AA4499"}   
else if ( file1[$1] == "Nigeria" ) { color = "#44AA99"}  
else if ( file1[$1] == "Other"   ) { color = "#888888"} 
else if ( file1[$1] == "Poland"   ) { color = "#882255"} 
else if ( file1[$1] == "UK") { color = "#661100"}
print part1"\t"color
} 
else {print $1"\tBELIEVE\tblack"} }
' $labels $PCAres > $outLocc
head $outLocc
wc -l $outLocc
}
export -f AddEthnicLabels # this makes local functions executable when bsubbed



# labelsfile=$believeResultsLoc$'UKBB_BELIEVE_MATCHED/labels'
# phefile=$believeDataLoc$'bmi_all'
function subsetPhen { # removes those from the labels that do not have a pheno
labelsfile=$1
phefile=$2
outLo=$3
awk  'FNR == NR {  file1[ $1 ] = $1; next;  } 
FNR <= NR {  if(FNR == 1) {print $0}
else if( $1 in file1 ) 
{ print $0 } }' $phefile $labelsfile > $outLo
head $outLo
wc -l $outLo
wc -l $labelsfile
}




function GenUKBB { # the UKBB already has "Male" and "Female" so we dont have to convert those
covFile=$1
phefile=$2
outputN=$3
binary=$4

awk -v binary="$binary" 'FNR == NR { if($2 == "Male" || $2 == "Female") {file1[ $1 ] = $2;} next;  } 
FNR <= NR {  if(FNR == 1) {print "IID\tpheno\tSEX"}
if( $1 in file1 ) 
{ 
if(binary == 0 || binary ==1 && $3 !=0) {print $1"\t"$3"\t"file1[$1]}
} }' $covFile $phefile > $outputN
head $outputN

GenRes $outputN $binary
}


# covFile=$believeDataLoc$'covars_all.phe'
# phefile=$believeDataLoc$'highCholesterol_all.phe'
# outputN=$believeScratchLoc$'cholest_cov'
# outputName=$outputN
# binary=1
function GenBLV { 
covFile=$1
phefile=$2
outputN=$3
binary=$4

awk -v binary="$binary" 'FNR == NR { 
if($3 == "1") {file1[ $1 ] = "Male";}
else if($3 == "2") {file1[ $1 ] = "Female";}
next;  } 
FNR <= NR {  if(FNR == 1) {print "IID\tpheno\tSEX"}
if( $1 in file1 ) 
{ 
if(binary == 0 || binary ==1 && $3 !=0) {print $1"\t"$3"\t"file1[$1]}
} }' $covFile $phefile > $outputN
head $outputN

GenRes $outputN $binary
}

function GenRes { 
outputName=$1
bin=$2

awk -v binary="$bin" '{ if(FNR > 1) { if(binary == 1) {print $2-1} else {print $2} } }' $outputName > $outputName$'_phe'
awk '{ print $3 }' $outputName > $outputName$'_cov'
awk '{ if(FNR > 1) {  print $1}  }' $outputName > $outputName$'_id'

arguments='/home/mk907/scripts/R/phenoRegress_backward.R '$outputName$'_phe '$outputName$' _res '$outputName$'_cov 0 0 '$bin$' 0'
Rscript $arguments 

# then re-merge the pheno residuals with their pheno pheno IDs
paste $outputName$'_id' $outputName$'_res_z' > $outputName$'_res_z_id'
head $outputName$'_res_z_id'
wc -l $outputName$'_res_z_id'

paste $outputName$'_id' $outputName$'_res' > $outputName$'_res_id'
head $outputName$'_res_id'
wc -l $outputName$'_res_id'
}


#### Pheno file prep for the 1,500 indis runs (only kept for backwards compatibiltiy, these are not used in the final analysis)

# covFile=$believeDataLoc$'covars_all.phe'
# phefile=$believeDataLoc$'highCholesterol_all.phe'
# outputN=$believeScratchLoc$'cholest_cov'
# outputName=$outputN
# binary=1
function GenUKBB2 { # the UKBB already has "Male" and "Female" so we dont have to convert those
covFile=$1
phefile=$2
outputN=$3
binary=$4

awk -v binary="$binary" 'FNR == NR { if($4 == "Male" || $4 == "Female") {file1[ $1 ] = $2"\t"$3"\t"$4;} next;  } 
FNR <= NR {  if(FNR == 1) {print "IID\tpheno\tAGE\tAGE2\tSEX"}
if( $1 in file1 ) 
{ 
if(binary == 0 || binary ==1 && $3 !=0) {print $1"\t"$3"\t"file1[$1]}
} }' $covFile $phefile > $outputN
head $outputN

GenRes2 $outputN $binary
}


# covFile=$believeDataLoc$'covars_new.phe'
# phefile=$believeDataLoc$'highCholesterol_all.phe'
# outputN= $believeScratchLoc$'cholest_cov2'
# binary=1
covFile=$believeDataLoc$'covars_new.phe'
phefile=$believeDataLoc$'edu_all.phe'
outputN=$believeScratchLoc$'edu_cov2'
binary=0
function GenBLV2 { 
covFile=$1
phefile=$2
outputN=$3
binary=$4

awk -v binary="$binary" 'FNR == NR { 
if($5 == "1") {file1[ $1 ] = $3"\t"$4"\tMale";}
else if($5 == "2") {file1[ $1 ] = $3"\t"$4"\tFemale";}
next;  } 
FNR <= NR {  if(FNR == 1) {print "IID\tpheno\tAGE\tAGE2\tSEX"}
if( $1 in file1 ) 
{ 
if(binary == 0 || binary ==1 && $3 !=0) {print $1"\t"$3"\t"file1[$1]}
} }' $covFile $phefile > $outputN
head $outputN

GenRes2 $outputN $binary
}


outputName=$outputN
bin=$binary
function GenRes2 { 
outputName=$1
bin=$2

awk -v binary="$bin" '{ if(FNR > 1) { if(binary == 1) {print $2-1} else {print $2} } }' $outputName > $outputName$'_phe'
awk '{ print $3"\t"$4"\t"$5 }' $outputName > $outputName$'_cov'
awk '{ if(FNR > 1) {  print $1}  }' $outputName > $outputName$'_id'

arguments='/home/mk907/scripts/R/phenoRegress_backward.R '$outputName$'_phe '$outputName$' _res '$outputName$'_cov 0 0 '$bin$' 0'
Rscript $arguments 


# phenoLoc =  "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/scratch/edu_cov2_phe"
# outputLoc= "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/scratch/edu_cov2"
# outputName = "_res"
# covLoc = "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/scratch/edu_cov2_cov"
# onlySigPredictors = 0
# bonf = 0
# binary = 0
# numFactors = 0

#### Pheno file prep for the 1,500 indis runs (only kept for backwards compatibiltiy, these are not used in the final analysis)




# then re-merge the pheno residuals with their pheno pheno IDs
paste $outputName$'_id' $outputName$'_res' > $outputName$'_res_id'
head $outputName$'_res_id'
wc -l $outputName$'_res_id'
}




# UKBB calculate the prevalence
function findprevalence { 
pop=$1
currentphe=$2
awk -v pop=$pop 'FNR == NR { if($2 == pop) {file1[ $1 ] = $2;} next;  } 
FNR <= NR {  
if( $1 in file1 ) { if ($3 == 1) {controls++} else {cases++}   } } END { print cases / (controls + cases)}
' $believeScratchLoc$'labelsPC4_ALL' $currentphe 
}


#outFile=$resultLoc$'Sims_'$sim$'_VG.txt'
#files=( "${args1[@]}" )
function concatResLoop { # concatones the results of any number of simulations
outFile=$1
shift            # Shift all arguments to the left (original $1 gets lost)
files=("$@") # Rebuild the array with rest of arguments # https://askubuntu.com/questions/674333/how-to-pass-an-array-as-function-argument   and https://stackoverflow.com/questions/12303974/copying-a-bash-array-fails

arraylength3=${#files[@]}
echo $arraylength3
echo $outFile
# get header
head -n 1 ${files[0]} > $outFile
for (( z=1; z<${arraylength3}+1; z++ )); do
echo ${files[$z-1]}

awk '{if (FNR > 1 && $1 != "NA") {print $0}}' ${files[$z-1]} > ${files[$z-1]}'_clean'
cat ${files[$z-1]}'_clean' >>  $outFile
rm -rf ${files[$z-1]}'_clean'
done
}



function concatRes { # concatones the results of 3 simulations
file1=$1
file2=$2
file3=$3
outFile=$4

head -n 1 $file1 > $outFile$'_header'
awk '{if (FNR > 1 && $1 != "NA") {print $0}}' $file1 > $file1$'_clean'
awk '{if (FNR > 1 && $1 != "NA") {print $0}}' $file2 > $file2$'_clean'
awk '{if (FNR > 1 && $1 != "NA") {print $0}}' $file3 > $file3$'_clean'

cat $outFile$'_header' $file1$'_clean' $file2$'_clean' $file3$'_clean' > $outFile

rm -rf $outFile$'_header'
rm -rf $file1$'_clean'
rm -rf $file2$'_clean'
rm -rf $file3$'_clean'
}



function binarPheno { 
covs=$1
binPhe=$2
currentph=$3

awk 'FNR == NR { {file1[ $2 ] = $3;} next;  } 
FNR <= NR {  
if (FNR == 1) {print $0}
else { print $1"\t"$2"\t"file1[$1]"\t"$4} }
' $binPhe $currentph > $currentph$'_bin2'


awk 'FNR == NR { {file1[ $1 ] = $2"\t"$3"\t"$4;} next;  } 
FNR <= NR {  
if (FNR == 1) {print $0"\tage\tage2\tsex"}
else { print $0"\t"file1[$1]} }
' $covs $currentph$'_bin2' > $currentph$'_bin'

rm -rf $currentph$'_bin2'

head $currentph$'_bin'
}

# simple reformats phenos into IID Pheno format, subsetting them to what we have in the covars
function binarPheno_all { 
covs=$1
binPhe=$2
currentph=$3


awk 'FNR == NR { {file1[ $1 ] = $2"\t"$3"\t"$4;} next;  } 
FNR <= NR {  
if($2 in file1) { print $2"\t"$3} }
' $covs $binPhe > $currentph

head $currentph
}

trait="MI"
pheFile=$believeDataLoc$'binary_traits.txt' 
keeplist=$believeScratchLoc$'ukbb_keeplist' 




# calculates N_eff based on regenie multiple pheno file, and a keeplist for a given trait
function GetN_eff_regenie {
trait=$1
pheFile=$2 # stores the cases/controls counts
keeplist=$3 # stores the actualy people used for training (this may be more than the phefile)

# 3) filter for QC failed in UKBB
awk -v trait=$trait 'FNR == NR { file1[$2] = $2; next; }
FNR <= NR { 
# find column number based on specified trait
if(FNR == 1) {
for(i=1; i<=NF; i++) {if($i==trait) traitCol=i;}
}
# then look at the rest of the file, 
# we only care about individuals on the keeplist
if($2 in file1) {

if($traitCol == 1) controls++
else if ($traitCol == 2) cases++
else if ($traitCol == "NA") NAs++
 }
 } END { 
# print "numCases "cases" / numcontrols: "controls" / NAs: "NAs
 # calculate N_eff
N_eff = 4 / (1 / cases + 1 / controls) 
print int(N_eff) 
} ' $keeplist $pheFile

}



# calculates N_eff based on a .phe file and a list of individuals used for assoc
function GetN_eff {
pheFile=$1 # stores the cases/controls counts
totalsFile=$2 # stores the actualy people used for training (this may be more than the phefile)

# 3) filter for QC failed in UKBB
awk 'FNR == NR { file1[$2] = $2; next; }
FNR <= NR { 

if($2 in file1) {

if($3 == 1) controls++
else if ($3 == 2) cases++
 }
 } END { 
# print "numCases "cases" / numcontrols: "controls
 # calculate N_eff
N_eff = 4 / (1 / cases + 1 / controls) 
print int(N_eff) 
} ' $totalsFile $pheFile
}



# generates manhattan plot from a REGENIE GWAS association file
function produce_manhattan_REGENIE { 
plinkAssocFile=$1 # this is the standard plink .phe.glm.logistic.hybrid file
output_loc=$2
plottitle=$3
sigTreshold=$4

# map the REGENIE output format 

#  1    2     3    4       5       6   7   8   9   10  11      12   13
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
# TO
# SNP CHR BP P
awk '{ 
if ( FNR == 1) {print "SNP\tCHR\tBP\tP"}
else  { print $3"\t"$1"\t"$2"\t"$12 }
 }' $plinkAssocFile > $output_loc$'_gwasResults'


# call Rscript to produce plot
arguments='/home/mk907/scripts/R/make_manhattan_fixed.R '$output_loc$'_gwasResults '$output_loc$' '$plottitle$' '$sigTreshold$' 0 0'
Rscript $arguments # has to be run on farm5, as dgx-server has the "X11 is not available" error

rm -rf $output_loc$'_gwasResults'
}
export -f produce_manhattan_REGENIE # this makes local functions executable when bsubbed





# generates manhattan plot from a PLINK2 GWAS association file
function produce_manhattan_PLINK2 { 
plinkAssocFile=$1 # this is the standard plink .phe.glm.logistic.hybrid file
output_loc=$2
plottitle=$3
sigTreshold=$4
linear=$5


# map the PLINK2 output format depending on if results file is binary or linear
if [[ "$linear" == '1' ]]; then
echo 'linear'
#   1 		2     3        4       5       6           7       8       9        10    11       12       13      14       15    16  
#CHROM	   POS	  ID	  REF	   ALT	   A1   	A1_FREQ	 TEST	 OBS_CT	  BETA	  SE	  L95	    U95	 T_STAT	  P	    ERRCODE
# TO
# extract summary stats into format into the following signature: # SNP CHR BP         P    zscore
awk '{ 
if ( FNR == 1) {print "SNP\tCHR\tBP\tP\tzscore"}
else  { print $3"\t"$1"\t"$2"\t"$15"\t"$14 }
 }' $plinkAssocFile > $output_loc$'_gwasResults'
else
echo 'binary'
# 1      2         3         4   5  6      7         8        9      10       11            12        13          14          15         16             17
#CHROM	POS	      ID	    REF	ALT	A1	A1_FREQ	   FIRTH?	TEST	OBS_CT	OR	        LOG(OR)_SE	L95       	U95	        Z_STAT	     P	        ERRCODE 
#9	  22115286	rs944797	C	T	T	0.493031	N	     ADD	190724	0.875327	0.0254157	0.832792	0.920035	-5.2392 	1.61274e-07	. 
# TO
# extract summary stats into format into the following signature: # SNP CHR BP         P    zscore
awk '{ 
if ( FNR == 1) {print "SNP\tCHR\tBP\tP\tzscore"}
else  { print $3"\t"$1"\t"$2"\t"$16"\t"$15 }
 }' $plinkAssocFile > $output_loc$'_gwasResults'
fi

# call Rscript to produce plot
arguments='/home/mk907/scripts/R/make_manhattan.R '$output_loc$'_gwasResults '$output_loc$' '$plottitle$' '$sigTreshold
Rscript $arguments # has to be run on farm5, as dgx-server has the "X11 is not available" error

#rm -rf $output_loc$'_gwasResults'
}
export -f produce_manhattan_PLINK2 # this makes local functions executable when bsubbed


# generates manhattan plot from a PLINK1 formatted file
function produce_manhattan_PLINK1 { 
ssFile=$1 # this is the standard plink .phe.glm.logistic.hybrid file
output_loc=$2
plottitle=$3
sigTreshold=$4


# map the PLINK1 output format depending on if results file is binary or linear
#   1 		2     3        4       5       6           7       8       9        10    11       12
# CHR      SNP     BP     A1       TEST    NMISS       BETA    SE      L95      U95   STAT      P
# TO
# extract summary stats into format into the following signature: # SNP CHR BP         P    zscore
awk '{ 
if ( FNR == 1) {print "SNP\tCHR\tBP\tP\tzscore"}
else  { print $2"\t"$1"\t"$3"\t"$12"\t"$11 }
 }' $ssFile > $output_loc$'_gwasResults'


# call Rscript to produce plot
arguments='/home/mk907/scripts/R/make_manhattan.R '$output_loc$'_gwasResults '$output_loc$' '$plottitle$' '$sigTreshold
Rscript $arguments 

#rm -rf $output_loc$'_gwasResults'
}
export -f produce_manhattan_PLINK1 # this makes local functions executable when bsubbed



# 3) exclude ambiguous alleles:  A/T or C/G from a .bim formatted file
function remove_ambiguous_alleles_bim { 
sumstatsLoc=$1
#  1     2     3    4     5
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
awk ' 
FNR <= NR {
if (FNR == 1) {print $0 }
else { 
if( toupper($6) == "A" && toupper($5) == "T" || toupper($5) == "A" && toupper($6) == "T" || toupper($6) == "G" && toupper($5) == "C" || toupper($5) == "G" && toupper($6) == "C") {}
else {print $0}
}
}
' $sumstatsLoc > $sumstatsLoc$'_noAmbiguousAlleles'

awk '{print $2}' $sumstatsLoc$'_noAmbiguousAlleles' > $sumstatsLoc$'_noAmbiguousAlleles_SNPlist'

head $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc
}
export -f remove_ambiguous_alleles_bim # this makes local functions executable when bsubbed

