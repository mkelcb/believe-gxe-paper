# Performs simulation with a h2 of 0.25 on the BELIEVE datasets, 2 populations
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
initResultDF = function(){
  DF <- matrix(NA, numTests, 2+3+4) 
  colnames(DF) <- c("VG","Ve","VG","Ve","FEnv","VG","Ve","FEnv","VGxE") # add colnames
  return(DF)
}

initResultDF_subpop = function(pops){
  DF <- matrix(NA, numTests, length(pops)) 
  colnames(DF) <-pops 
  return(DF)
}

############

a = c(1,2,3)
b = c(2,4,6)
c = cbind(a,b)
x = c(2,2,2)
c*x

n=10
p=3
set.seed(42)
X = matrix(rnorm(n*p),n,p)
EnvTerm=rnorm(n)
GxETerms = matrix(NA,n,p)
print("GxETerms = matrix(NA,n,p)")
print(gc())
# X = scale(X)
for (k in 1:p) {
  GxETerms[,k] = X[,k] * EnvTerm 
}
GxETerms2 = X*EnvTerm


# X = X_all
# Vg = h2
# VGxE_true
simPheno_GxE = function(X, EnvTerm, Vg, VGxE_true){  # X: genotype matrix, EnvTerm: environtment matrix, h2: target heritability, h2_env: the % of the total phenotypic variance coming from the Env, and h2_GxE, the h2 due to GxE
  n = nrow(X)
  p = ncol(X)
  
  beta = rnorm(p, mean = 0, sd = 1) 
  GV = X%*%beta
  scale = as.numeric( sqrt( Vg / var(GV) ) )
  g = GV * scale # scale the pure genetic values 
  gc()
  print("X%*%beta")
  # calculate env component
  # 2. via linear model of GxE= X1 * Env + X2*Env, IE env interacts with every causal variant separately
  EnvTerm_z = scale(EnvTerm)[,1]
 # GxETerms = matrix(NA,n,p)
 # print("GxETerms = matrix(NA,n,p)")
 # 
  # X = scale(X)
  # for (k in 1:p) {
  #   GxETerms[,k] = X[,k] * EnvTerm 
  # }
  GxETerms = X*EnvTerm # this is the same as above
  print(gc())
  print("GxETerms = X*EnvTerm")
  # REML will underestimate VGxE if cor(GxEV,EnvTerm) is too large, so we keep regenerating GxEV until it is not
  GxECor=1
  while (abs(GxECor) > 0.1) {
    beta_gxe = rnorm(p, mean = 0, sd = 1)
    GxEV = GxETerms%*%beta_gxe
    GxECor = cor(GxEV,EnvTerm) 
    #cor(GxEV,EnvTerm) 
   # cor(gxe,EnvTerm)
    gc()
  }
  remove(GxETerms)
  print("remove(GxETerms)")
  print(gc())
  gc()

  
  scale_GxE = as.numeric( sqrt( VGxE_true / var(GxEV) ) )
  gxe = GxEV * scale_GxE # scale the pure genetic values
  
  noise = rnorm(n, mean = 0, sd = sqrt(1-VGxE_true- Vg ))  # sqrt(Ve) # 
  
  y = g  +EnvTerm+ noise + gxe #  
  
  #results <- list("y" = y, "beta" = beta, "beta_env" = beta_env, "beta_gxe"= beta_gxe)
  #results <- list("y" = y, "GxECor" = GxECor)
  results <- list("g" = g, "noise" = noise, "gxe" = gxe)
  
  remove(X)
  gc()
  return(results)
}

##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################



options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 8) {stop("not enough Arguments received")}

plinkGenotypesLoc=args[1]
phenoLabelsLoc=args[2]
outLoc=args[3]
simFunctLoc=args[4]
isPopKin=F
if (args[5] == "1") { isPopKin = T}
print(paste0("isPopKin: ", isPopKin))
numbStart=as.numeric(args[6])
numbEnd=as.numeric(args[7])
set.seed(numbStart)
precalculatedKinshipLoc=F
if (args[8] != "0") { precalculatedKinshipLoc = args[8]}

pops = c()
for (i in 9:length(args)) { # loop through the rest of the arguments, where each is supposed to be a pop label
  pops = c(pops,args[i] )
}

source(simFunctLoc)

# Constants
h2= 0.333
VGxE_true= 0.333 # how much variance GxE explains
EnvFixedBeta= 1# by using 1, the Envterm is simply the same as the EnvIndicator 0.5 # make sure this is on the same scale as h2/VGxE so that we can display them in the same plot
numTests=20
LDAKEnabled=F
VG_ =1
VG_Env = 2
VE_ENV_GXE = 3

gctaoloc = "/home/mk907/software/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"

## DEBUG VARS
# plinkGenotypesLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/raw/1445_ALL"
# phenoLabelsLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/data/height_1445_pop"
# outLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED/BLV/height_res"
# outLoc="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/believe_pca/results/UKBB_BELIEVE_MATCHED/BLV/sims_res"
# pops = c("BELIEVE", "UK_SEA")

#  simFunctLoc="C:/Users/mk23/GoogleDrive_Cam/project/BangladeshExome/scripts/R/sims/pheDecomp_functs.R"
#  simFunctLoc="C:/Users/mk23/GoogleDrive_Cam/project/BangladeshExome/scripts/R/sims/pheDecomp_functs.R"

# load in and process all data
loadedData = loadInData (plinkGenotypesLoc,phenoLabelsLoc,outLoc,pops, isPopKin = isPopKin, precalculatedKinshipLoc = precalculatedKinshipLoc)
fam_pheno = loadedData$fam_pheno
X_all = loadedData$X_all
X_pops_indices = loadedData$X_pops_indices
K_standard = loadedData$K_standard
K_subpops = loadedData$K_subpops
famData = loadedData$famData

###########################################
# 3 scenarios to test:
# 1. just same baseline h2
# 2. the last subpop has an Env Fixed effect
# 3. the last subpop has an Env and GxE Effect

# plot 3 facets for the 3 true scenarios:
# 1 scenario has Title: "truth: Vg: 0.25 / Env: 0, GxE:0", in different colors (green, blue, orange)
# has 3x3 groups of results per plot for the 3 models, each with 3 estimates (grey "N/A" for models that did not have that)
# at the level of the true value, also add a black horizontal bar

# 1 data frame per model: rows: 20 replicates, 3 columns each variance estimates
# so 3 dataframes per scenario (3 models)
# and we have a list of 3 for the 3 scenarios (ie 3x3 Dfs)

scenarios=list() # this will have 3 DFs
scenarios[[VG_]] =initResultDF() # truth: VG 
scenarios[[VG_Env]] = initResultDF() # truth: VG+Env 
scenarios[[VE_ENV_GXE]] = initResultDF() # truth: VG+Env+GxE 
# and each data frame will have 3 columns, one for each estimate



# per subpop plots (these tell us how h2 looks, if we only had one population)
# 3 true scenarios as before
# within, list of 3 dataframes, 1 for each truth scenario, each with as many columns as pops
scenarios_subpop=list()
scenarios_subpop[[VG_]] = initResultDF_subpop(pops)
scenarios_subpop[[VG_Env]] = initResultDF_subpop(pops)
scenarios_subpop[[VE_ENV_GXE]] = initResultDF_subpop(pops)


# single dataframe for the mean genetic values for each pop
GV_subpop <- matrix(NA, numTests, length(pops)) 
colnames(GV_subpop) <- pops



i=1
gc()

for (i in numbStart:numbEnd) { 
  print(paste0("NUM TEST: ", i))
  print(gc())
  # always perform same 3+1 tests, just Vg, VG+Env and VG+Env+GxE plus the individual subpops on  their own
  # also export out GV for each group (and calculate Fst, to check if pairwise Fst and GV has some relationship)

  # generate a fixed effect Env Term
  EnvTerm = rep(0,nrow(X_all))
  EnvTerm[X_pops_indices[[length(X_pops_indices)]] ] = EnvFixedBeta # add the fixed effect to the last pop's indices
  # EnvTerm  =rbinom(N*2, 1, 0.5) 

  envIndicator  = rep(0,nrow(X_all))
  envIndicator[X_pops_indices[[length(X_pops_indices)]] ] = 1 # as EnvFixedBeta =1, this is the same as Envterm
  
  # generate phenotypes
  simData =simPheno_GxE (X_all, EnvTerm, h2, VGxE_true)
  gc()
  print(gc())
  print("Finished simPheno_GxE")
  #I) simulate base pheno for all, of the same overall Vg
  VGnoise = rnorm(nrow(X_all), mean = 0, sd = sqrt(1-h2)) # add noise with a variance, which is always everything other than the other VCs 
  y_Vg = simData$g + VGnoise
  
  # II) Add fixed Effect Env to the second population
  y_VgEnv = simData$g + VGnoise + EnvTerm 
  
  # III) Vg + Fixed Env + GxE (same Vg, ie 'underlying' h2, but we have a fixed effect which then interacts with the Genes
  y_VgEnvGxE = simData$g + simData$noise + simData$gxe  + EnvTerm
  

  # perform all 3 models on all 3 truth scenarios
  truthVg = performAllTests(K_standard, famData, outLoc, envIndicator, y_Vg)
  truthVgEnv = performAllTests(K_standard, famData, outLoc, envIndicator, y_VgEnv)
  truthVgEnvGxE =  performAllTests(K_standard, famData, outLoc, envIndicator, y_VgEnvGxE)

  # store results
  # 1. truth scenario only VG, model only VG
  scenarios[[VG_]][i,1] = truthVg$modelVG[1] 
  scenarios[[VG_]][i,2] = truthVg$modelVG[2]
  
  # 1. truth scenario only VG, model VG+Env
  scenarios[[VG_]][i,3] = truthVg$modelVGEnv[1] 
  scenarios[[VG_]][i,4] = truthVg$modelVGEnv[2]
  scenarios[[VG_]][i,5] = truthVg$modelVGEnv[3]
  
  # 1. truth scenario only VG, model VG+Env+GxE
  scenarios[[VG_]][i,6] = truthVg$modelVGEnvGxE[1] 
  scenarios[[VG_]][i,7] = truthVg$modelVGEnvGxE[2]
  scenarios[[VG_]][i,8] = truthVg$modelVGEnvGxE[3]
  scenarios[[VG_]][i,9] = truthVg$modelVGEnvGxE[4]
  
  
  # 2. truth scenario VG +Env, model only VG
  scenarios[[VG_Env]][i,1] = truthVgEnv$modelVG[1] 
  scenarios[[VG_Env]][i,2] = truthVgEnv$modelVG[2]
  
  # 2. truth scenario VG+Env, model VG+Env
  scenarios[[VG_Env]][i,3] = truthVgEnv$modelVGEnv[1] 
  scenarios[[VG_Env]][i,4] = truthVgEnv$modelVGEnv[2]
  scenarios[[VG_Env]][i,5] = truthVgEnv$modelVGEnv[3]
  
  # 2. truth scenario VG+Env, model VG+Env+GxE
  scenarios[[VG_Env]][i,6] = truthVgEnv$modelVGEnvGxE[1] 
  scenarios[[VG_Env]][i,7] = truthVgEnv$modelVGEnvGxE[2]
  scenarios[[VG_Env]][i,8] = truthVgEnv$modelVGEnvGxE[3]
  scenarios[[VG_Env]][i,9] = truthVgEnv$modelVGEnvGxE[4]

  # 3. truth scenario VG +Env +Gxe, model only VG
  scenarios[[VE_ENV_GXE]][i,1] = truthVgEnvGxE$modelVG[1] 
  scenarios[[VE_ENV_GXE]][i,2] = truthVgEnvGxE$modelVG[2]
  
  # 3. truth scenario VG+Env +Gxe, model VG+Env
  scenarios[[VE_ENV_GXE]][i,3] = truthVgEnvGxE$modelVGEnv[1] 
  scenarios[[VE_ENV_GXE]][i,4] = truthVgEnvGxE$modelVGEnv[2]
  scenarios[[VE_ENV_GXE]][i,5] = truthVgEnvGxE$modelVGEnv[3]
  
  # 3. truth scenario VG+Env +Gxe, model VG+Env+GxE
  scenarios[[VE_ENV_GXE]][i,6] = truthVgEnvGxE$modelVGEnvGxE[1] 
  scenarios[[VE_ENV_GXE]][i,7] = truthVgEnvGxE$modelVGEnvGxE[2]
  scenarios[[VE_ENV_GXE]][i,8] = truthVgEnvGxE$modelVGEnvGxE[3]
  scenarios[[VE_ENV_GXE]][i,9] = truthVgEnvGxE$modelVGEnvGxE[4]
  


  # also perform basic Vg/Ve test on each sub populations
  k=1
  remove(K_standard)
  gc()
  for (k in 1:length(X_pops_indices)) {
    popIndicecs = X_pops_indices[[k]]
    
    GV_subpop[i,k] =  mean( simData$g[popIndicecs] ) # also save the mean breeding values for each pop

    
   # Xz = scale(X_all[popIndicecs,])
    # K_sub = Xz %*% t(Xz) / ncol(X_all)   
    K_sub = K_subpops[[k]]
    #K_PopKin = 2 * popkin::popkin(X_all, loci_on_cols = TRUE)
    
    # write out common data for GCTA
    iddata= paste0("indi_",1:length(popIndicecs)) 
    famData_sub = cbind.data.frame(iddata,iddata)
    colnames(famData_sub) = c("fam","id")
    outLoc_sub = paste0(outLoc,"_sub",k)
    genio::write_grm(outLoc_sub,K_sub, fam = famData_sub)
    
     
    truthVg_subpop = estimate_VgVe(K_sub, famData_sub, outLoc_sub, y_Vg[popIndicecs])
    truthVgEnv_subpop = estimate_VgVe(K_sub, famData_sub, outLoc_sub, y_VgEnv[popIndicecs])
    truthVgEnvGxE_subpop = estimate_VgVe(K_sub, famData_sub, outLoc_sub, y_VgEnvGxE[popIndicecs])
    
    # store results
    scenarios_subpop[[VG_]][i,k] = truthVg_subpop[1] # 1. truth scenario only VG
    scenarios_subpop[[VG_Env]][i,k] = truthVgEnv_subpop[1]  # 2. truth scenario VG+Env
    scenarios_subpop[[VE_ENV_GXE]][i,k] = truthVgEnvGxE_subpop[1]  # 3. truth scenario VG+Env+GxE
    remove(K_sub)
    remove(truthVg_subpop)
    remove(truthVgEnv_subpop)
    remove(truthVgEnvGxE_subpop)
    gc()
  }
  remove(  truthVg)
  remove(truthVgEnv)
  remove(truthVgEnvGxE)
  remove(simData)
  gc()
}


# write results to disk
write.table(scenarios[[VG_]], paste0(outLoc,"_VG.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)
write.table(scenarios[[VE_ENV_GXE]], paste0(outLoc,"_VG_Env.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)
write.table(scenarios[[VG_Env]], paste0(outLoc,"_VE_ENV_GXE.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)

write.table(scenarios_subpop[[VG_]], paste0(outLoc,"_VG_sub.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)
write.table(scenarios_subpop[[VG_Env]], paste0(outLoc,"_VG_Env_sub.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)
write.table(scenarios_subpop[[VE_ENV_GXE]], paste0(outLoc,"_VE_ENV_GXE_sub.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)

write.table(GV_subpop, paste0(outLoc,"_GV.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)


print(paste0("Analysis finished and results written to ", outLoc))


# 
# A = rnorm(1000)
# B = rnorm(1000)
# cor(A,B)
# cor(A,B+10)
# cor(A,B*2)
# 
# cor(A,B) == cor(A,B+10)


