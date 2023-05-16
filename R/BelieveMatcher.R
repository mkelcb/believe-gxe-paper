#############################
# Matches UKBB individuals to the BELIEVE data
#############################



#options(error=traceback)
# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")}


# arguments
inputPCsLoc=args[1]
groupLabelsLoc=args[2] # list of predictors who match the colnames of the first arg
outLoc = args[3]


# arguments
#inputPCsLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE/_PC_PCA"
#PCsToVisualiseLoc=args[2] # list of predictors who match the colnames of the first arg
#groupLabelsLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/ukbb_believe_ethnic"



# load PCs (or MCEs)
PCs = read.table(inputPCsLoc, header=F) 
groupLabels = read.table(groupLabelsLoc, sep="\t", header=T, comment.char = "") # need this otherwise the colors with #123456 wouldn't be able to be read


PCs_needed = PCs


renamedCols = colnames(PCs_needed)
renamedCols = as.numeric( gsub("V", "", renamedCols) ) -1
#renamedCols
colnames(PCs_needed) = c("IID", renamedCols[2:length(renamedCols)])
#colnames(PCs_needed) 
# load group labels

# merge group labels with the PCs
labels_PCs = merge(groupLabels, PCs_needed, by = "IID")

legendNames = names( table(labels_PCs$group) )
#legendNames

BELIEVE_minPC = c()
BELIEVE_maxPC = c()
existingIndis = c()
existingUptoPC4 = c()
# go through all PCs, subset to BELIEVE indis, and get their min/max values for their PCs
i=4
for (i in 4:ncol(labels_PCs) ) {
  
  # get the BELIEVE indis
  BELIEVE_PC = labels_PCs[labels_PCs$group == "BELIEVE",i]
  
  minPC = min(BELIEVE_PC)
  maxPC = max(BELIEVE_PC)
  BELIEVE_minPC = c(BELIEVE_minPC, minPC )
  BELIEVE_maxPC = c(BELIEVE_maxPC, maxPC )
  
  # get all indis who are within the BELIEVE PC range from those indis that are NOT from BELIEVE (ie from the UKBB)
  NON_BELIEVE_PC = labels_PCs[labels_PCs$group != "BELIEVE",c(1,i)]
  indis = which(NON_BELIEVE_PC[,2] <= maxPC & NON_BELIEVE_PC[,2] >= minPC )
  print(paste0("found ",length(indis), " within PC", (i-3) ) )
  
  if(i == 4) {
    existingIndis = indis
  } else {
    existingIndis = intersect(existingIndis, indis)
    if(i-3 == 4) {existingUptoPC4 = existingIndis} 
    print(paste0("found ",length(existingIndis), " so far up to PC", (i-3) ) )
  }
  
}

# get all non BELIEVE Indis
NON_BELIEVE = labels_PCs[labels_PCs$group != "BELIEVE",c(1,2)]

# new
#  India Other
#  614   924
# why is India lower? probably due to being excluded being related to the "Other"s

UKBB_MATCHED = NON_BELIEVE[existingIndis,]
table(UKBB_MATCHED$group)

outFileLoc=paste0(outLoc,"_matched")
write.table(UKBB_MATCHED, outFileLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("written ", nrow(UKBB_MATCHED), " matched indis to: ", outFileLoc))



# get more people who are less well matched
UKBB_MATCHED = NON_BELIEVE[existingUptoPC4,]
table(UKBB_MATCHED$group)

outFileLoc=paste0(outLoc,"_matchedPC4")
write.table(UKBB_MATCHED, outFileLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("written ", nrow(UKBB_MATCHED), " matched (up to PC4) indis to: ", outFileLoc))



