# from UKBB, this extracts a subset of indis for PCA, as well as a larger subset of whites

##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")}

PC_UKBB_IdLoc = args[1]
all_excludedLoc = args[2]
outLoc = args[3]
numBELIEVE = as.numeric(args[4])

# Debug vars
# PC_UKBB_IdLoc <- read.table("C:/0LocalHPC/believe/UKBB_PCA16", header=T)
# all_excludedLoc <- read.table("C:/0LocalHPC/believe/all_excluded", header=F)
# outLoc = "C:/0LocalHPC/believe/PCA"
# numBELIEVE=33712

# get the to be excluded people


PC_UKBB_Id <- read.table(PC_UKBB_IdLoc, header=T)
all_excluded <- read.table(all_excludedLoc, header=F)
to_be_excluded_indis = match(all_excluded$V1 ,PC_UKBB_Id$eid)
#to_be_excluded_indis = to_be_excluded_indis[is.na(to_be_excluded_indis) == F]

PC_UKBB = PC_UKBB_Id[-1]
all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)

# compute square distances from the pop center for each PC
all_sq_dist <- apply(all_centers[-1], 1, function(one_center) {
  rowSums(sweep(PC_UKBB, 2, one_center, '-')^2)
})


thr_sq_dist <- max(dist(all_centers[-1])^2) * 0.002 / 0.16
group <- apply(all_sq_dist, 1, function(x) {
  grp <- NA
  ind <- which.min(x)
  if (isTRUE(x[ind] < thr_sq_dist)) {
    grp <- all_centers$Ancestry[ind]
    # We used a more stringent cutoff for the Ashkenazi group
    if (grp == "Ashkenazi" && x[ind] > 12.5^2) grp <- NA
  }
  grp
})


#head(group)
table(group, exclude = NULL)
#Ashkenazi      Caribbean          China          India           Iran          Italy        Nigeria 
#2493           2657                1852           6676           1224           6809           4087 
#Poland United Kingdom           <NA> 
#  4308         445944          11355 


# get table that has people and their pops
PC_UKBB_Id_group = cbind.data.frame(PC_UKBB_Id$eid,group)
colnames(PC_UKBB_Id_group) = c("IID","pop")

# from these exclude the ones that need to be excluded
PC_UKBB_Id_excluded = PC_UKBB_Id_group[-to_be_excluded_indis,]
PC_UKBB_Id_excluded$pop[is.na(PC_UKBB_Id_excluded$pop)] = "Other"
PC_UKBB_Id_excluded$pop[PC_UKBB_Id_excluded$pop == "United Kingdom"] = "UK"
allGroups = table(PC_UKBB_Id_excluded$pop, exclude = NULL)
PC_UKBB_Id_excluded = cbind(PC_UKBB_Id_excluded$IID,PC_UKBB_Id_excluded) # duplicate the IID col, so that it can be used as a keeplist in plink without any changes


#Ashkenazi      Caribbean          China          India           Iran          Italy        Nigeria 
#2354           2470           1813                6327           1202           6492           3944 
#Poland United Kingdom           <NA> 
#  4148         411614          10865

# pick a random (max) India number people from all groups
countryNames = names(allGroups)
india_index = which(countryNames == "India")
numIndia = allGroups[india_index] # 6327 

set.seed(42)
i=1
allPCA = NULL
for (i in 1:length(countryNames)) {
  numAvailable  = allGroups[i]
  num = numIndia
  if(countryNames[i] == "Other"){
    num = nrow(PC_UKBB_Id_excluded)
    print(paste0("For Other, we pick ", num))
  }
  numToPick = min(numAvailable,num)
  
  allIndisOfGroup = PC_UKBB_Id_excluded[PC_UKBB_Id_excluded$pop == countryNames[i],]
  index = 1:nrow(allIndisOfGroup)
  indisToPick = sample(index,numToPick)
  indis = allIndisOfGroup[indisToPick,]
  
  if(is.null(allPCA) ) {  allPCA = indis  } else {
    allPCA = rbind(allPCA,indis)
  }
}

outFileLoc=paste0(outLoc,"_PCAIndis")
write.table(allPCA, outFileLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("written ", nrow(allPCA), " PCA indis to: ", outFileLoc))

# now write just the indians
justIndians = allPCA[which(allPCA$pop == "India"),]
outFileLoc=paste0(outLoc,"_India")
write.table(justIndians, outFileLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("written ", nrow(justIndians), " Indian indis to: ", outFileLoc))


# write white british,
justBritish = allPCA[which(allPCA$pop == "UK"),]
outFileLoc=paste0(outLoc,"_British")
write.table(justBritish, outFileLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("written ", nrow(justBritish), " British indis to: ", outFileLoc))

justOther = allPCA[which(allPCA$pop == "Other"),]
outFileLoc=paste0(outLoc,"_Other")
write.table(justOther, outFileLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("written ", nrow(justOther), " Other indis to: ", outFileLoc))



# write white british of the same sample size BELIEVE (33712)
allIndisOfGroup = PC_UKBB_Id_excluded[PC_UKBB_Id_excluded$pop == "UK",]
index = 1:nrow(allIndisOfGroup)
indisToPick = sample(index,numBELIEVE)
indis = allIndisOfGroup[indisToPick,]

outFileLoc=paste0(outLoc,"_British_BELIEVE")
write.table(indis, outFileLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("written ", nrow(indis), " British indis to: ", outFileLoc))


# cols = rainbow(length(countryNames))
# pie(rep(1, 8), col = cols)
# safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
#                              "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
# 
# pie(rep(1, 8), col = safe_colorblind_palette)
# 

