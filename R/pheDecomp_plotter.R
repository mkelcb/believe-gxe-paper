# Plots loaded data from pheno decomp simulations

library(Matrix)


##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################


set.seed(42)
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 8) {stop("not enough Arguments received")}

scenarioLoc=args[1]
subpopsLoc=args[2]
GVsLoc=args[3]
outLoc=args[4]
h2=as.numeric(args[5])
EnvFixedBeta=as.numeric(args[6])
VGxE_true=as.numeric(args[7])
fst_loaded=as.numeric(args[8])



# 
# # DEBUG DATA
# h2 = 0.25
# EnvFixedBeta = 1
# VGxE_true = 0.25
# 
# # create dummy data
# numEstimates=4
# numModels=(2+3+4)
# inputDatas_mat = matrix(rnorm(numTests*numModels,0.5,0.2),ncol=numModels, nrow=numTests)
# 
# subpop = as.data.frame( matrix(rnorm(numTests*length(pops),0.5,0.2), ncol=length(pops), nrow=numTests) )
# colnames(subpop) = pops
# 
# GVs = as.data.frame( matrix(rnorm(numTests*length(pops),0.5,0.2), ncol=length(pops), nrow=numTests) )
# colnames(GVs) = pops
# fst_loaded = 0.1
# pops = c("BELIEVE", "UK_SEA")

# load data
inputDatas_mat =read.table(scenarioLoc, header=T)
subpop=read.table(subpopsLoc, header=T)
GVs=read.table(GVsLoc, header=T)
pops = colnames(subpop)

####################

# plot main results
colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
allColours=c(colorBlind7[1:3],colorBlind7[6]) # pick 4 nice colors for the 4 estimates (Vg, Ve, VGxE, Fenv)
allColours_rep = c(allColours[1:2], allColours[1:3], allColours[1:4]) # the first model only has Vg and Ve, the second Vg,Ve,Env and the third has Vg,Ve,Env,VGxE

plotName= paste0("truth: Vg: ",h2," / Env: ",EnvFixedBeta,", GxE: ",VGxE_true) #, in different colors (green, blue, orange) #"" # plotname disabled for publication
minVal = 0
#maxVal = max((h2 + VGxE_true), EnvFixedBeta)  * 1.1 # max y is the larger of either h2 or the env, whichever is greater plus some headroom
maxVal = 1  * 1.1  # always use the max scale of 1, to ensure same scale across the plots
boxNames = c("VG only", "VG+Env", "VG+Env+GxE") # each box will have the name of the model (ie what we are testing, NOT the truth)
pixelScale = 100
yaxislab = "estimate"


allMeans = list()
for(i in 1:ncol(inputDatas_mat)) { allMeans[[i]] =mean(inputDatas_mat[,i]) }


# create true value list
truevals= list(h2,1-h2,  h2,1-h2,EnvFixedBeta,  h2,1-h2-VGxE_true,EnvFixedBeta,VGxE_true)

inputDatas_mat = as.data.frame(inputDatas_mat)
# insert 2 blank columns in the middle to separate out the 3 models: https://stackoverflow.com/questions/58708529/how-to-create-spaces-between-groups-and-control-size-of-axis-labels-in-boxplot
#n_ = ncol(inputDatas_mat)
#VEC =  seq(1, n_/4, length.out=n_)*4 - c(0, .2)
VEC = 1:ncol(inputDatas_mat)
VEC[3:length(VEC)] = VEC[3:length(VEC)] +1
VEC[6:length(VEC)] = VEC[6:length(VEC)] +1


filen =paste(outLoc,".png", sep="" )
png(filen, width=3.5* length(boxNames) * pixelScale , height=6.4 * pixelScale);

boxplot(allMeans , outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0,xaxt="n", border =allColours_rep, xlab="", ylab=parse(text=yaxislab), main=plotName, ylim=c(minVal, maxVal), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8, at=VEC)
stripchart( inputDatas_mat, vertical = TRUE, pch = 16, col = allColours_rep, add = TRUE, method = "jitter", cex = 2,at=VEC) #
# add true vals
boxplot(truevals ,notch =T, outline = FALSE, boxlty =0, whisklty = 0, add = TRUE, staplelty = 0,xaxt="n", xlab="", ylab="", main="", ylim=c(minVal, maxVal), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8, at=VEC)


unit = ( ncol(inputDatas_mat) / 3) /2 + 0.5
axis(1,at=c(1.5, 5,5.5 +4),labels=boxNames, cex.axis = 2,mgp=c(3,1.5,0))  # , tck=0
legend("topright", legend=c("Vg","Ve","Env","VGxE","truth"),col=c(allColours,"black"), text.col= c(allColours,"black"), lty=1, lwd=3)
title(xlab = "models", line = 3.7, cex.lab = 1.5)

dev.off()


#######

# plot subpops
allMeans = list()
truevals= list() # create true value list
for(i in 1:ncol(subpop)) { 
  allMeans[[i]] =mean(subpop[,i]) 
  truevals[[i]] = h2
}

allColours=colorBlind7[1] # pick 4 nice colors for the 4 estimates (Vg, Ve, VGxE, Fenv)
allColours_rep = rep(allColours, length(pops)) # the first model only has Vg and Ve, the second Vg,Ve,Env and the third has Vg,Ve,Env,VGxE
maxVal=1

plotName = "VG if estimated in subpops"
boxNames = pops

filen =paste(outLoc,"_sub.png", sep="" )
png(filen, width=3.5* length(boxNames) * pixelScale , height=6.4 * pixelScale);

boxplot(allMeans , outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0,xaxt="n", border =allColours_rep, xlab="", ylab=parse(text=yaxislab), main=plotName, ylim=c(minVal, maxVal), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8)

stripchart( subpop, vertical = TRUE, pch = 16, col = allColours_rep, add = TRUE, method = "jitter", cex = 2) #
# add true vals
boxplot(truevals ,notch =T, outline = FALSE, boxlty =0, whisklty = 0, add = TRUE, staplelty = 0,xaxt="n", xlab="", ylab="", main="", ylim=c(minVal, maxVal), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8)


axis(1, at =1:length(pops),labels=boxNames, cex.axis = 2,mgp=c(3,1.5,0))  # , tck=0

legend("topright", legend=c("Vg","truth"),col=c(allColours,"black"), text.col= c(allColours,"black"), lty=1, lwd=3)
title(xlab = "populations", line = 3.7, cex.lab = 1.5)

dev.off()


# plot Genetic Values
plotName=paste0("True Genetic Values for populations (Fst:",round(fst_loaded,3),")")

allMeans = list()
for(i in 1:ncol(GVs)) { allMeans[[i]] =mean(GVs[,i]) }


filen =paste(outLoc,"_GV.png", sep="" )
png(filen, width=3.5* length(boxNames) * pixelScale , height=6.4 * pixelScale);

boxplot(allMeans , outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0,xaxt="n", border =allColours_rep, xlab="", ylab=parse(text=yaxislab), main=plotName, ylim=c(minVal, maxVal), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8)

stripchart( GVs, vertical = TRUE, pch = 16, col = allColours_rep, add = TRUE, method = "jitter", cex = 2) #


axis(1, at =1:length(pops),labels=boxNames, cex.axis = 2,mgp=c(3,1.5,0))  # , tck=0

legend("topright", legend=c("GV"),col=allColours, text.col= allColours, lty=1, lwd=3)
title(xlab = "populations", line = 3.7, cex.lab = 1.5)

dev.off()



  
