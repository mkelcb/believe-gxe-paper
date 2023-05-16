# takes phenotypes ( binary or continuous)
# and regresses out a set of covariates, and write the residuals to disk
# can enforce categorical/factor predictors by specifying their indices (number of these is the 7th argument, and the 8th argument onwards this is )
# uses AOV instead of lm/glm as we want to regress out all significant covars and only aov can do that for many leveled factors (?)

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 5) {stop("not enough Arguments received")} 
 
phenoLoc = args[1]
outputLoc= args[2]
outputName = args[3]
covLoc = args[4]
onlySigPredictors = as.numeric(args[5]) # 1= for only significant predictors to be regressed out, 0 for not
bonf = as.numeric(args[6]) # 1= for bonferroni correct, 0 for not
binary = as.numeric(args[7])

numFactors = as.numeric(args[8]) # an array length for all the factors
factors= vector(length=0 ) # these are the (1-based) indices of the factors

print(paste("numFactors is", numFactors))
if ( numFactors > 0 ) {
factors = vector(length = (length(args) -8) ) # 7 is the number of arguments so far
counter = 1
for (i in 9:length(args)) { # loop through the rest of the arguments, where each argument is an index of a factor
  factors[counter] = as.numeric(args[i])
  counter = counter +1
  print(as.numeric(args[i]))
}
}

# TODO:
# by checking the entries of the first row of the covariates, determine if they need to be treated as factors or numerics
# cant do this as the 'centre' variable is numeric like '11011' so it would not look like a factor even though it is
# decide based on number of values for the pheno if it is an lm or an AOV (only 2 possible values mean case control)

# WHy am I fitting an AOV in the first place??: because when obtaining the residuals, an lm with factors would be difficult


# Debug wars
# phenoLoc = "C:/softwares/Cluster/GIANT/miniPRS/BackwardSelection/pheno_short"
# outputLoc = "C:/softwares/Cluster/GIANT/miniPRS/BackwardSelection/"
# outputName = "residualsbackwards"
# covLoc = "C:/softwares/Cluster/GIANT/miniPRS/BackwardSelection/covariates_short"
# onlySigPredictors = 1
# bonf = 1
# numFactors = 1
# factors = c(4) # the 4th index is a factor
# binary =0

# 
# phenoLoc = "C:/softwares/Cluster/UKBB2/height/pheno_short"
# outputLoc = "C:/softwares/Cluster/UKBB2/"
# outputName = "height_all"
# covLoc = "C:/softwares/Cluster/UKBB2/height/covariates_short"
# onlySigPredictors = 1
# bonf = 1
# numFactors = 1
# factors = c(4) # the 4th index is a factor
# binary =0
# 

# 1. load phenos
pheFile = as.matrix( read.table(phenoLoc) )
pheno = as.matrix( unlist(as.numeric(pheFile[,1]) ))

# load the covariates
covs =  read.table(covLoc, header=TRUE) 

# force cast things as factors that may be interpreted as numeric
names = colnames(covs)
#print(names)
if ( numFactors > 0 ) {
for (i in 1:numFactors) { 
  print(paste( names[ factors[i] ], " is a factor") )
  covs[,factors[i] ] = as.character(covs[,factors[i] ])
  
}
}

#covs$CENTRE = as.character(allvars_df$CENTRE) 
allvars_df = cbind(pheno,covs)
#allvars_df = allvars_df[1:5000,]
#myModel_glm = lm(pheno ~.,allvars_df )


# Do backward selection
if (binary == 1) {
  print("binary phenotype, fitting GLM")
  # want to know more details about the direction of effect
  inner_glm = glm(pheno ~.,allvars_df,family=binomial(logit) )
  #summary(inner_glm)
  
  myModel_glm = aov( inner_glm )
  
} else {
  inner_glm = lm(pheno ~.,allvars_df )
  myModel_glm = aov( inner_glm  )
}

origColnames = colnames(covs) # save the original col names and their indices



outFileLoc = paste(outputLoc,outputName, "_modelSummary_inner", sep="")
capture.output(summary(inner_glm), file = outFileLoc)
print(paste("written inner model summary to: ", outFileLoc))

outFileLoc = paste(outputLoc,outputName, "_modelSummary", sep="")
capture.output(summary(myModel_glm), file = outFileLoc)
print(paste("written model summary to: ", outFileLoc))




############################

allColsToRemove <- list() # will hold all the cols to be removed (indices mapped to their original positions)

# first round
AIC_model =  AIC(myModel_glm)

#step(myModel_glm) # this errors out sometimes, and does nothing actually?
#step(inner_glm)


summary(myModel_glm)
# get all p values
p_values_all = as.numeric ( summary(myModel_glm)[[1]][["Pr(>F)"]] )

# exclude the intercept
p_values_touse = p_values_all[1:(length(p_values_all) -1) ]

# find the largest p value
largest_p_index =which(p_values_touse == max(p_values_touse)) # the index of this is not mapped to the original list

largestPval =p_values_touse[largest_p_index]

# establish Sig threshold on the first run
alpha=0.05 # apply bonferroni correction if needed
if ( bonf == 1 ) { alpha = alpha / length(p_values_touse)}
print( paste("alpha is", alpha, "after number of corrections:", length(p_values_touse)) )
if ( onlySigPredictors == 0 ) { 
  alpha = 1
  print("setting alpha to 1, as we want to keep all covariates" )
  }

while (largestPval > alpha) {

# map this index to the column names IN the model ( this may be different than the total, as we may already have removed some)
currentColnames= colnames(myModel_glm$model) # this will be different as we keep removing terms
currentColnames = currentColnames[2:length(currentColnames)]
toRemoveCol = currentColnames[largest_p_index]
print(paste("removing term",toRemoveCol ," with p-val of:",largestPval, "model AIC:", round(AIC_model, 3)  ))
toRemoveCol_origIndex = which(origColnames == toRemoveCol) # map the index of the surviving term to the original index


# add the column to be removed to the total
allColsToRemove <- append(allColsToRemove, toRemoveCol_origIndex)


# remove this from the data and refit
allvars_df = cbind(pheno,covs[,-unlist(allColsToRemove)])
#colnames(allvars_df)

if (binary == 1) {
  print("binary phenotype, fitting GLM")
  # want to know more details about the direction of effect
  inner_glm = glm(pheno ~.,allvars_df,family=binomial(logit) )
  #summary(inner_glm)
  
  myModel_glm = aov( inner_glm )
  
} else {
  inner_glm = lm(pheno ~.,allvars_df )
  myModel_glm = aov( inner_glm  )
}
#summary(myModel_glm)
AIC_model =  AIC(myModel_glm)

# get all p values
p_values_all = as.numeric ( summary(myModel_glm)[[1]][["Pr(>F)"]] )

# exclude the intercept
p_values_touse = p_values_all[1:(length(p_values_all) -1) ]
#print( length(p_values_all) )
# find the largest p value
largest_p_index =which(p_values_touse == max(p_values_touse)) # the index of this is not mapped to the original list

largestPval =p_values_touse[largest_p_index]

}


# Covariates kept:
if( length(unlist(allColsToRemove) ) > 0 ) { # crappy R will freak out with 0 length vectors
sigCols = colnames(covs[,-unlist(allColsToRemove)])
} else { sigCols = colnames(covs)}

outFileLoc = paste(outputLoc,outputName, "_significantCovariates.txt", sep="")
write.table(sigCols, outFileLoc, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
print(paste("written ", length(sigCols), " significant covariates to: ", outFileLoc))
print( sigCols ) 


print( paste("out of the ", (length(origColnames)), " predictors, ",(length(sigCols)), " were significant:" , sep="") )



outFileLoc = paste(outputLoc,outputName, "_modelSummary_inner_REDUCED", sep="")
capture.output(summary(inner_glm), file = outFileLoc)
print(paste("written _REDUCED inner model summary to: ", outFileLoc))

outFileLoc = paste(outputLoc,outputName, "_modelSummary_REDUCED", sep="")
capture.output(summary(myModel_glm), file = outFileLoc)
print(paste("written _REDUCED model summary to: ", outFileLoc))






phenoResiduals = myModel_glm$residuals

#View(as.matrix(phenoResiduals))




# export out the residuals WITH the intercept,
# as that may be different between populations
filename = paste(outputLoc, outputName , sep="")
write.table(myModel_glm$coefficients[1] +phenoResiduals, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE) 
print(paste("written Phenotype residuals WITH intercept to: ", filename ))

# and without the intercept
filename = paste(outputLoc, outputName,"_noInt" , sep="")
write.table(phenoResiduals, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE) 
print(paste("written Phenotype residuals to: ", filename ))


filename = paste(outputLoc, outputName,"_z" , sep="")
write.table(scale(phenoResiduals), filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE) 

print(paste("written z-scored Phenotype residuals to: ", filename ))


# , breaks="Scott"
pdf(paste(outputLoc,"adjusted",".pdf", sep="" ), width=12.8 , height=6.4);
hist(phenoResiduals, col="red", probability = T,main ="Distribution of adjusted phenotypes", xlab = "phenotype", breaks="Scott")
dev.off()


pdf(paste(outputLoc,"original",".pdf", sep="" ), width=12.8 , height=6.4);
hist(pheno, col="red", probability = T,main ="Distribution of original phenotypes", xlab = "phenotype", breaks="Scott")
dev.off()
