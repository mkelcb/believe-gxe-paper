# Projects new indis into an existing PC space

################################

library(bigstatsr)
library(bigsnpr)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 


dataLoc=args[1]
inputBedLoc=args[2]
outLoc=args[3]

#usedIndisLoc = paste0(dataLoc,"usedIndis.txt")
usedIndis_idxLoc = paste0(dataLoc,"usedIndis_idx.txt")
#usedSNPs_idxLoc =  paste0(dataLoc,"usedSNPs_idx.txt")

#DEBUG VARS:
#dataLoc = '/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/GxS_common/_svd'
obj.svd2 = readRDS(file = paste0(dataLoc,"_svd") ) 
(obj.bed <- bed(inputBedLoc)) # this needs to be loaded for both PCa and ncMCE

#usedIndis=unlist(read.table(usedIndisLoc,  header = FALSE))
ind.row = unlist(read.table(usedIndis_idxLoc,  header = FALSE) )
#ind.col2 = unlist(read.table(usedSNPs_idxLoc,  header = FALSE) )
newIndis=obj.bed$.fam$sample.ID[-ind.row]

proj <- bed_projectSelfPCA(obj.svd2, obj.bed,
                           ind.row = rows_along(obj.bed)[-ind.row],
                           ncores = 1) # useless -> too few individuals

PCmat_indis = cbind(newIndis,proj$OADP_proj) # concat it with the indi IDs

write.table(PCmat_indis,outLoc , sep = "\t", row.names = F, col.names = F, quote = FALSE)
print(paste("written  PCs for ",nrow(PCmat_indis)," new indis to", outLoc))

# also generate sanity check plot



PCs <- matrix(NA, nrow(obj.bed), ncol(obj.svd2$u))
PCs[ind.row, ] <- predict(obj.svd2)
PCs[-ind.row, ] <- proj$OADP_proj


plotTitle=paste0("PCA: New Indis projection")
png(filename=paste0(outLoc,"newIndis.png")  , width=1280, height=720, res=128)
plot(PCs[ind.row, 1:2], pch = 20, xlab = "PC1", ylab = "PC2")
points(PCs[-ind.row, 1:2], pch = 20, col = "blue")
dev.off()



