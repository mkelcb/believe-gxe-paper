# Performs pheno decomposition analysis on real data 
# 
# load these BEFORE, we load the custom ones
library(Matrix)
library(genio)
#library(GxEMM)
library(popkin)
library(bigstatsr)
library(bigsnpr)
##################





##########################

##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################


set.seed(42)
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 6) {stop("not enough Arguments received")}

plinkGenotypesLoc=args[1]
phenoLabelsLoc=args[2]
outLoc=args[3]
simFunctLoc=args[4]
isPopKin=F
if (args[5] == "1") { isPopKin = T}

pops = c()
for (i in 6:length(args)) { # loop through the rest of the arguments, where each is supposed to be a pop label
  pops = c(pops,args[i] )
}

source(simFunctLoc)

# Constants
gctaoloc = "/home/mk907/software/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"

## DEBUG VARS
# plinkGenotypesLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/raw/1445_ALL"
# phenoLabelsLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/height_1445_pop"
# outLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED/BLV/height_res"
# pops = c("BELIEVE", "UK_SEA")

#  simFunctLoc="C:/Users/mk23/GoogleDrive_Cam/project/BangladeshExome/scripts/R/sims/pheDecomp_functs.R"

#  simFunctLoc="C:/Users/mk23/GoogleDrive_Cam/project/BangladeshExome/scripts/R/sims/pheDecomp_functs.R"


# load in and process all data
loadedData = loadInData (plinkGenotypesLoc,phenoLabelsLoc,outLoc,pops, isPopKin = isPopKin)
fam_pheno = loadedData$fam_pheno
X_all = loadedData$X_all
X_pops_indices = loadedData$X_pops_indices
K_standard = loadedData$K_standard
K_subpops = loadedData$K_subpops
famData = loadedData$famData

####################a

# create Env indicator, as last pop's indices
envIndicator  = rep(0,length(fam_pheno$pheno))
envIndicator[X_pops_indices[[length(X_pops_indices)]] ] = 1 # as EnvFixedBeta =1, this is the same as Envterm



  # perform all 3 models on all 3 truth scenarios
print(paste0("Starting performAllTests, K_standard: ", dim(K_standard), " / K_subpops[[1]]: ", dim(K_subpops[[1]])) )
  truthVg = performAllTests(K_standard, famData, outLoc, envIndicator, fam_pheno$pheno)
  print("Finished performAllTests")
  # store results
  
  
  sink(paste0(outLoc, "_results.txt") )
  
  print("VG model:")
  print(paste0("VG: ",truthVg$modelVG[1] ) )
  print(paste0("Ve: ",truthVg$modelVG[2])  )
  print("______________________")
  
  
  print("VG+Env model:")
  print(paste0("VG: ",truthVg$modelVGEnv[1])  )
  print(paste0("Ve: ",truthVg$modelVGEnv[2] ) )
  print(paste0("Env: ",truthVg$modelVGEnv[3])  )
  print("______________________")
  
  print("VG+Env+GxE model:")
  print(paste0("VG: ",truthVg$modelVGEnvGxE[1])  )
  print(paste0("Ve: ",truthVg$modelVGEnvGxE[2] ) )
  print(paste0("Env: ",truthVg$modelVGEnvGxE[3] ) )
  print(paste0("GxE: ",truthVg$modelVGEnvGxE[4])  )
  print("______________________")
  
  

  # also perform basic Vg/Ve test on each sub populations
  k=1
  for (k in 1:length(X_pops_indices)) {
    popIndicecs = X_pops_indices[[k]]
    
 
   # Xz = scale(X_all[popIndicecs,])
    # K_sub = Xz %*% t(Xz) / ncol(X_all)   
    K_sub = K_subpops[[k]]
    #K_PopKin = popkin::popkin(X_all, loci_on_cols = TRUE)
    
    # write out common data for GCTA
    iddata= paste0("indi_",1:length(popIndicecs)) 
    famData_sub = cbind.data.frame(iddata,iddata)
    colnames(famData_sub) = c("fam","id")
    outLoc_sub = paste0(outLoc,"_sub",k)
    genio::write_grm(outLoc_sub,K_sub, fam = famData_sub)
    
     
    truthVg_subpop = estimate_VgVe(K_sub, famData_sub, outLoc_sub, fam_pheno$pheno[popIndicecs])

    # store results
    print(paste0("pop ",pops[k]," VG: ",truthVg_subpop[1] ) )
    print("______________________")
  }
  sink()

print(paste0("Analysis finished and results written to ", outLoc))



  
