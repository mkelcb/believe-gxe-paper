# produces Manhattan plot from GWAS/TWAS/PWAS association statistics
if (!require("qqman")) install.packages("qqman")
library("qqman")

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")} 
 
gwasLoc = args[1]
outputLoc= args[2]
plotName= args[3]
sigtrshold= as.numeric(args[4]) # 7.30103

#-log10(5e-8) # GWAS # 7.30103
# -log10(5.76*10e-6) #  4.239578
# -log10(5*10e-6)# 4.30103

#
# TWAS
#p=0.0000000000001
#p_log= -log10(p)
#p_log

# Debug wars
# gwasLoc = "C:/softwares/Cluster/GIANT/miniPRS/manhattan/_gwasResults"
# outputLoc = "C:/softwares/Cluster/GIANT/miniPRS/manhattan/"
# plotName = "PWAS results"

# 1. load phenos
gwasRes= read.table(gwasLoc, header = T) 
gwasRes$CHR = as.numeric(gwasRes$CHR)
minAcceptablePval = 5*10e-100
#gwasRes$P[1] = NA
#gwasRes$P
# clamp the lowest p val to something sensible
print( paste("BEFORE minium p val is: ", min(gwasRes$P) ) )
gwasRes$P[which(gwasRes$P == 0) ] = minAcceptablePval

gwasRes <- na.omit(gwasRes) 

print( paste("AFTER minium p val is: ", min(gwasRes$P) ) )

# max(gwasRes$CHR)
# min(gwasRes$P)
# 
# CHR_23_hits = gwasRes[ which(gwasRes$CHR==23),]
# CHR_24_hits = gwasRes[ which(gwasRes$CHR==24),]
# 
# CHR_15_hits = gwasRes[ which(gwasRes$CHR==15),]


filename = paste(outputLoc,plotName, "_manhattan.pdf" , sep="")
pdf(filename, width=12.8 , height=6.4);
manhattan(gwasRes, main = plotName, col = c("blue4", "orange3"), suggestiveline=FALSE, genomewideline=sigtrshold)
dev.off()


#options(bitmapType='cairo')
filename = paste(outputLoc,plotName, "_manhattan.png" , sep="")
png(filename, width=1280 , height=640, res =128);
manhattan(gwasRes, main = plotName, col = c("blue4", "orange3"), suggestiveline=FALSE, genomewideline=sigtrshold)
dev.off()




print(paste("saved manhattan plot to: ", filename ))



