# produces Manhattan plot from GWAS/TWAS/PWAS association statistics
if (!require("qqman")) install.packages("qqman")
library("qqman")

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 6) {stop("not enough Arguments received")} 
 
gwasLoc = args[1]
outputLoc= args[2]
plotName= args[3]
sigtrshold= as.numeric(args[4]) # 7.30103
setLog10P=F # if false, then we assume that P vals are already -log10
if (args[5] == "1") { isLogP = T}
makePDF = F
if (args[6] == "1") { makePDF = T}

#-log10(5e-8) # GWAS # 7.30103
# -log10(5.76*10e-6) #  4.239578
# -log10(5*10e-6)# 4.30103

#
# TWAS
#p=0.0000000000001
#p_log= -log10(p)
#p_log

# Debug wars
# gwasLoc = "C:/0LocalHPC/results/AAD_any/GWAS_HM3_gwasResults"
# outputLoc = "C:/0LocalHPC/results/AAD_any/GWAS_HM3"
# plotName = "AAD_any_GWAS"
# sigtrshold= 7.30103


# 1. load phenos
gwasRes= read.table(gwasLoc, header = T) 
gwasRes$CHR = as.numeric(gwasRes$CHR)
if(setLog10P) {
  

minAcceptablePval = 5*10e-100
#gwasRes$P[1] = NA
#gwasRes$P
# clamp the lowest p val to something sensible
print( paste("BEFORE minium p val is: ", min(gwasRes$P) ) )
gwasRes$P[which(gwasRes$P == 0) ] = minAcceptablePval
gwasRes <- na.omit(gwasRes) 
print( paste("AFTER minium p val is: ", min(gwasRes$P) ) )
} else {
  gwasRes <- na.omit(gwasRes) 
  print( paste("AFTER minium p val is: ", max(gwasRes$P) ) )
}
#gwasRes[which(gwasRes$P == min(gwasRes$P)), ]





# max(gwasRes$CHR)
# min(gwasRes$P)
# 
# CHR_23_hits = gwasRes[ which(gwasRes$CHR==23),]
# CHR_24_hits = gwasRes[ which(gwasRes$CHR==24),]
# 
# CHR_15_hits = gwasRes[ which(gwasRes$CHR==15),]


mycolors =c("blue4", "orange3")
if(makePDF) {
filename = paste(outputLoc, "_manhattan.pdf" , sep="")
pdf(filename, width=12.8 , height=6.4);
manhattan(gwasRes, main = gsub("_", " ", plotName), col = mycolors, suggestiveline=FALSE, genomewideline=sigtrshold, logp = setLog10P) # , cex.axis = 0.5
dev.off()
}

#options(bitmapType='cairo')
factor=1
filename = paste(outputLoc, "_manhattan.png" , sep="")
png(filename, width=1280 *factor , height=800*factor, res =128);
manhattan(gwasRes, main = gsub("_", " ", plotName), col = mycolors, suggestiveline=FALSE, genomewideline=sigtrshold, logp = setLog10P) # , cex.axis = 0.5
dev.off()

# how to force display tick labels: use lower cex.axis ... or custom place them
# https://stackoverflow.com/questions/8688854/r-draw-all-axis-labels-prevent-some-from-being-skipped
# ?axis

print(paste("saved manhattan plot to: ", filename ))



