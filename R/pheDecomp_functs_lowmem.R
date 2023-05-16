# the reusable functions for the pheno decomposition
library(stringr)

loadInData = function(plinkGenotypesLoc,phenoLabelsLoc,outLoc,pops, isPopKin = F, precalculatedKinshipLoc = F, parts = 4) { # the reusable functions to load in the data
  #library(data.table)
  print (paste0("V2 of loadInData, precalculatedKinshipLoc: ", precalculatedKinshipLoc))
  if(precalculatedKinshipLoc != F) { kinLoc = precalculatedKinshipLoc } else {kinLoc = outLoc}
  
  phenoLabels=read.table(phenoLabelsLoc, header=T)
  #phenoLabels$indi_index=1:nrow(phenoLabels)
  tempDirLoc=paste0(outLoc,"temp/")
  unlink(tempDirLoc, recursive = TRUE)
  
  (obj.bed <- bed(paste0(plinkGenotypesLoc,".bed") )) # this is only needed for the Fst calculation
  
  
  snpData = snp_readBed2(paste0(plinkGenotypesLoc,".bed"), # ind.row = ind.row,
                         backingfile = tempDirLoc, # sub_bed(bedfile)
                         ncores = 1) # if 16, then this causes segmentation errors
  
  snpAttached <- snp_attach(snpData)				   
 # X_all = snpAttached$genotypes
  
  # need to load/convert it into integers in parts, otherwise I get OOM
  numIndisTotal = nrow(snpAttached$genotypes) #  17895 786150

    
 
  
  chunker <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  indiIndices = 1:numIndisTotal
  partIndices = chunker(indiIndices,parts)
  #X_all = NULL
  X_all = matrix(0L, nrow(snpAttached$genotypes), ncol(snpAttached$genotypes))
  gc()
  for(i in 1:length(partIndices)) {
    print(paste0("processing part",i,"/",parts ))
    X_part = snpAttached$genotypes[partIndices[[i]],] # 37515078216 bytes -> to 18757539216 bytes
    gc()
    mode(X_part) <- "integer"
    gc()
    X_all[partIndices[[i]],] = X_part
    # 
    # if(i == 1) {
    #   X_all = X_part
    # } else { 
    #   #X_all = rbind(X_all,X_part) # thisi s slow
    #   X_all <- rbindlist(list(X_all,X_part))
    #   }

    remove(X_part)
    gc()
  }


  #X_all = as.integer(X_all[]) # flatten 112545234216 bytes 
  #mode(X_all) <- "integer"
  #object.size(X_all) # 56272617216 bytes

  
  fam2=snpAttached$fam
  fam2$indi_index=1:nrow(fam2) # get the original indices the are aligned to the genotypes
  
  # subset the fam list to those groups we want to keep for the analysis
  keep_indices = c()
  pop_indices = list()
  for (i in 1:length(pops)) {
    pop_indices[[i]] = which(phenoLabels$pop == pops[i])
    keep_indices = c(keep_indices, pop_indices[[i]]  )
  }

  
  # get the IIDs of those we want to keep
  phenoLabels = phenoLabels[keep_indices,] 
  head(phenoLabels)
  
  #merge out the fam indis with the above
  fam_pheno = merge (fam2, phenoLabels,by.x="sample.ID", by.y="IID" )
  fam2_all = fam2[fam_pheno$indi_index,] # create a new fam file that is aligned to the subset X_all
  fam2_all$new_index = 1:nrow(fam2_all)
  
  # now subset the genotype data too
  X_all_map = 1:nrow(X_all)
  X_all_map_selected = X_all_map[fam_pheno$indi_index] # this will have a map of {original_index, new_index }
  # so this will have the indices of the ORIGINAL
  # but what we want is the subpop indices mapped to the 'X_all', which is already SUBSET
  
  
  # this creates an OOM, as this temporarily has to hold both arrays in RAM (before X_all would be overwritten by the smaller array)
  X_all = X_all[fam_pheno$indi_index,] 
  gc()
  print("finished assembling genotype matrix")
  
  # also subset 
  X_pops_indices = list()
  for (i in 1:length(pops)) {
    fam_sub = fam_pheno[which(fam_pheno$pop == pops[i]),] # subset the pheno file to subpop, that has the indices aligned to X_ordered
    # now map fam_sub$indi_index, to X_all
    # this is done via fam2_all$new_index, as this matched to X_all, but has also the same as fam_sub$indi_index
    # so we subset fam2_all to be the same people as fam_sub, and then use its fam2_all$new_index to get the map
    
    #               select X_all indices,  where     the subset pops orig indices, match the indices of X_all
    X_all_sub_indices = fam2_all$new_index[ match(fam_sub$indi_index, fam2_all$indi_index )   ]
    
    
    X_pops_indices[[i]] = X_all_sub_indices  # need to find the indices of the subpop, mapped to the X_all, ie the subset 
  }
  # need to map the subpop's indices onto the NEW X_all, as that is against what we aligned the phenotypes
  
  

  
  
  
  i=1
  
  # generate Kinship matrices, standard and PopKin way
  #Xz = scale(X_all)
  # K_standard = Xz %*% t(Xz) / ncol(X_all)   
  #K_PopKin = popkin::popkin(X_all, loci_on_cols = TRUE)
  
  dim(X_all)
  
  K_standard = calcKinship(snpAttached, X_all, fam2_all, kinLoc, isPopKin = isPopKin)
  gc()
  print("finished Kinship Calc")
  # write out common data for GCTA
  iddata= paste0("indi_",1:nrow(X_all)) 
  famData = cbind.data.frame(iddata,iddata)
  colnames(famData) = c("fam","id")
  genio::write_grm(outLoc,K_standard, fam = famData)
  #remove(K_standard)
  gc()
  
  # Fst: calculate this AFTER we have the kinship, as we can reuse that for popkin
  if(isPopKin ==T) {
    calcFst_popKin(K_standard /2, fam_pheno, outLoc)  # for the Fst calculations we need to undo the 2 * correction...
  
  } else {
    calcFst_WC(obj.bed, fam_pheno, outLoc, pops)
  }
  # calcFst(obj.bed, fam_pheno, outLoc, pops, isPopKin = isPopKin )
  gc()

  
  # this is an expensive operation, let's  cache the subpop kinships too
  K_subpops = list()
  for (k in 1:length(X_pops_indices)) {
    popIndicecs = X_pops_indices[[k]]
    
    #Xz = scale(X_all[popIndicecs,])
    # K_sub = Xz %*% t(Xz) / ncol(X_all)  
    fam_current = fam2_all[popIndicecs,]
    K_sub = calcKinship(snpAttached, X_all[popIndicecs,], fam_current, paste0(kinLoc,"_",pops[k]), isPopKin = isPopKin ) 
    gc()
    K_subpops[[k]] = K_sub
  }
  print("finished subpop Kinship Calc")

  
  # get rid of all missing data, as phenotype simulation will choke on it
  print("will look for missing genotypes")
  # OLD code for getting rid of missing genotype data: causes OOM
  # sum(is.na(X_all))
  # where_genotype_is_missing = which(is.na(X_all))
  # gc()
  #X_all[where_genotype_is_missing] = 0L # mean impute data, take care to NOT to convert to double, by using integer for 0   #typeof(X_all)
  #########################
  
  # do it in parts too to lower RAM requirements
  chunker <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  indiIndices = 1:nrow(X_all)
  partIndices = chunker(indiIndices,parts)
  X_part = NULL
  for(i in 1:length(partIndices)) {
    print(paste0("processing missing genotypes part",i,"/",parts ))
    inde = partIndices[[i]]
    
    X_part = X_all[inde,]
    where_genotype_is_missing = which(is.na(X_part))
    X_part[where_genotype_is_missing] = 0L
    X_all[inde,] = X_part
  }
  remove(X_part)
  gc()
  
  
  print("Loaded all Data")
  return( list( "fam_pheno" = fam_pheno , "X_all" = X_all, "X_pops_indices" = X_pops_indices , "K_standard" = K_standard, "K_subpops" = K_subpops, "famData" = famData, "fam2_all" = fam2_all) )
}




##########################
# estimate h2 via GCTA
# K = K_standard
# y = fam_pheno$pheno
estimate_VgVe = function(K, famData, outLoc, y ){ # basic h2 estimation, no Env or GxE
  famData$pheno = y # scale(y)[,1], DO NOT z-score!
  write.table(famData,paste0(outLoc,".phen"),col.names = F,row.names = F, quote=F)
  resultloc = paste0(outLoc,".hsq")
  if(file.exists(resultloc) ) {file.remove(resultloc) } # delete any previous runs
  # genio::write_grm(outLoc,K, fam = famData)
  cmd=paste0(gctaoloc,"  --thread-num 8 --reml-est-fix --reml-alg 1  --reml --grm ",outLoc," --pheno ",outLoc,".phen --out ",outLoc)
  system(cmd) 
  
  # check if file exist, and if it does not then return dummy values
  if(file.exists(resultloc) == F ) {return(c(-1,-1)) }
 
  Vg = readLines(file(resultloc, "r"), 5)[2] # read up until in the 5th line, but then return the 2nd line, which is Vg
  Vg = as.numeric(strsplit(Vg, "\t")[[1]] [2] ) # clean it
  
  Ve = readLines(file(paste0(outLoc,".hsq"), "r"), 5)[3] # read up until in the 5th line, but then return the 3rd line, which is Ve
  Ve = as.numeric(strsplit(Ve, "\t")[[1]] [2] ) # clean it
  
  
  # Perform corrections: scale the variance components back to 100% (as the fixed effect estimation will cause var(y) < 1, but to recover the original scale VCs we need to scale them by the ratio of each VC against the sum of the total random effects)
  Vp= Vg+Ve
  ratio = Vp /1 # correct based on the standardized variance of 1, ie just 1
  Vg = Vg/ratio
  Ve = Ve/ratio
  remove(K)
  gc()
  return(c(Vg,Ve))
}

# y = y_2pop
# envIndicator
estimate_VgVe_Env = function(K, famData, outLoc, y, envIndicator ){ # h2 estimation with a fixed effect Env
  famData$pheno = y # scale(y)[,1], DO NOT z-score!
  write.table(famData,paste0(outLoc,".phen"),col.names = F,row.names = F, quote=F)
  resultloc = paste0(outLoc,"_env.hsq")
  if(file.exists(resultloc) ) {file.remove(resultloc) } # delete any previous runs
  
  famData$pheno = scale( envIndicator)[,1] # envIndicator # EnvTerm # EnvTerm_z
  write.table(famData,paste0(outLoc,".qcovar"),col.names = F,row.names = F, quote=F)

  # use "covar" instead of "qcovar", as we are using a binary indicator env effect not a quantitative one
  cmd=paste0(gctaoloc," --thread-num 8 --reml-est-fix --reml-alg 1  --reml --grm ",outLoc," --pheno ",outLoc,".phen --qcovar ", outLoc, ".qcovar --out ",paste0(outLoc,"_env"))
  system(cmd) 
  
  
  # check if file exist, and if it does not then return dummy values
  if(file.exists(resultloc)  == F) {return(c(-1,-1, -1)) }
  
  Vg = readLines(file(resultloc, "r"), 5)[2] # read up until in the 5th line, but then return  2nd
  Vg = as.numeric(strsplit(Vg, "\t")[[1]] [2] ) # clean it
  
  Ve = readLines(file(resultloc, "r"), 5)[3]  # read up until in the 5th line, but then return 3rd
  Ve = as.numeric(strsplit(Ve, "\t")[[1]] [2] ) # clean it
  
  FEnv = readLines(file(resultloc, "r"), 15)[15]  # read up until in the 15th line, but then return 15th
  FEnv = as.numeric(strsplit(FEnv, "\t")[[1]] [1] ) # clean it
  FEnv = FEnv / sd(envIndicator) # restore Env effect to its original scale
  
  
  # Perform corrections: scale the variance components back to 100% (as the fixed effect estimation will cause var(y) < 1, but to recover the original scale VCs we need to scale them by the ratio of each VC against the sum of the total random effects)
  Vp= Vg+Ve
  ratio = Vp /1 # correct based on the standardized variance of 1, ie just 1
  Vg = Vg/ratio
  Ve = Ve/ratio
  remove(K)
  gc()
  return(c(Vg,Ve, FEnv))
}



# K = K_standard
# y = y_gxe
# envIndicator = EnvTerm
estimate_VgVe_Env_GxE = function(K, famData, outLoc, y, envIndicator ){ # h2 estimation with both fixed effect Env + GxE
  mean(envIndicator)
  var(envIndicator)
  mean(y)
  var(y)
  Fmodel = lm(y ~ envIndicator)
  summary(Fmodel) # this recovers the correct Env Beta, despite the interaction
  
  famData$pheno = y#   Fmodel$residuals # y # scale(y)[,1], DO NOT z-score!
  write.table(famData,paste0(outLoc,".phen"),col.names = F,row.names = F, quote=F)
  resultloc = paste0(outLoc,".hsq")
  if(file.exists(resultloc)  ) {file.remove(resultloc) } # delete any previous runs
  
  famData$pheno = envIndicator# scale( envIndicator)[,1] # envIndicator # EnvTerm # EnvTerm_z
  write.table(famData,paste0(outLoc,".qcovar"),col.names = F,row.names = F, quote=F)
  
  # remove: --qcovar ", outLoc, ".qcovar, otherwise will get error "Error: analysis stopped because more than half of the variance components are constrained. The result would be unreliable"
  cmd=paste0(gctaoloc," --thread-num 8 --reml-est-fix --reml-alg 1  --reml --reml-lrt 2 --grm ",outLoc," --pheno ",outLoc,".phen --gxe ",outLoc,".qcovar  --out ", paste0(outLoc) )
  
  system(cmd) # https://www.r-bloggers.com/2021/09/how-to-use-system-commands-in-your-r-script-or-package/
  
  
  # check if file exist, and if it does not then return dummy values
  if(file.exists(resultloc)== F ) {return(c(-1,-1,-1,-1)) }
  Vg = readLines(file(resultloc, "r"), 5)[2]  # read up until in the 5th line, but then return  2nd
  Vg = as.numeric(strsplit(Vg, "\t")[[1]] [2] ) # clean it
  
  VGxE = readLines(file(resultloc, "r"), 5)[3]  # read up until in the 5th line, but then return 3rd
  VGxE = as.numeric(strsplit(VGxE, "\t")[[1]] [2] ) # clean it
  
  Ve = readLines(file(resultloc, "r"), 5)[4]  # read up until in the 5th line, but then return 4th
  Ve = as.numeric(strsplit(Ve, "\t")[[1]] [2] ) # clean it
  
  FEnv = readLines(file(resultloc, "r"), 19)[19]  # read up until in the 19th line, but then return  19th
  FEnv = as.numeric(strsplit(FEnv, "\t")[[1]] [1] ) # clean it
  
  # Perform corrections: scale the variance components back to 100% (as the fixed effect estimation will cause var(y) < 1, but to recover the original scale VCs we need to scale them by the ratio of each VC against the sum of the total random effects)
  Vp= Vg+Ve+VGxE
  ratio = Vp /1 # correct based on the standardized variance of 1, ie just 1
  Vg = Vg/ratio
  Ve = Ve/ratio
  VGxE = VGxE/ratio
  remove(K)
  gc()
  return(c(Vg,Ve, FEnv, VGxE))
}

# y_current = y_Vg
# y_current = y_VgEnv
# y_current = y_VgEnvGxE
performAllTests = function(K_standard, famData, outLoc, envIndicator, y_current) {
  gc()
  modelVG = estimate_VgVe(K_standard, famData, outLoc, y_current)
  gc()
  modelVGEnv = estimate_VgVe_Env(K_standard, famData, outLoc, y_current, envIndicator)
  gc()
  modelVGEnvGxE = estimate_VgVe_Env_GxE(K_standard, famData, outLoc, y_current, envIndicator)
  remove(K_standard)
  gc()
  
  return (list ("modelVG" = modelVG, "modelVGEnv" = modelVGEnv, "modelVGEnvGxE" = modelVGEnvGxE ))
}

###################

# estimate h2 via GCTA
# K = K_standard
# y = fam_pheno$pheno
estimate_VgVe_bin = function(K, famData, outLoc, fam_pheno, prev ){ # basic h2 estimation, no Env or GxE
  y = fam_pheno$pheno
  
  famData$pheno = y # scale(y)[,1], DO NOT z-score!
  write.table(famData,paste0(outLoc,".phen"),col.names = F,row.names = F, quote=F)
  
  # print(paste0("estimate_VgVe_bin: nrow(famData): ",nrow(famData) , " / nrow(fam_pheno): ",nrow(fam_pheno) ))
  # print(str(famData))
  # print("______________")
  # 
  # print(str(fam_pheno))
  # print("______________")
  
  # write covariates into files 
  write.table(cbind.data.frame(famData[,1], famData[,2],fam_pheno$age, fam_pheno$age2),paste0(outLoc,".qcovar"),col.names = F,row.names = F, quote=F)
  write.table(cbind.data.frame(famData[,1], famData[,2],fam_pheno$sex.y),paste0(outLoc,".covar"),col.names = F,row.names = F, quote=F)

  resultloc = paste0(outLoc,".hsq")
  if(file.exists(resultloc) ) {file.remove(resultloc) } # delete any previous runs
  # genio::write_grm(outLoc,K, fam = famData)
  cmd=paste0(gctaoloc,"  --thread-num 8 --reml-est-fix --reml-alg 1  --reml --prevalence ",prev," --grm ",outLoc," --pheno ",outLoc,".phen --covar ",outLoc,".covar --qcovar ",outLoc,".qcovar --out ",outLoc)
  system(cmd) 
  
  # check if file exist, and if it does not then return dummy values
  if(file.exists(resultloc) == F ) {return(c(-1,-1)) }
  
  # as binary outcomes are transferred back to the liability scale, we read those directly
  h2 = readLines(file(resultloc, "r"), 8)[8] # read in the 9th line which contains the h2 estimate
  h2 = as.numeric(strsplit(h2, "\t")[[1]] [2] ) # clean it
  # reverse engineer the liability scale Ve, which is not reported by GCTA:
  Ve=1-h2 # as we dont have GxE here

  remove(K)
  gc()
  return(c(h2,Ve))
}





# y = y_2pop
# envIndicator
estimate_VgVe_Env_bin = function(K, famData, outLoc, fam_pheno, envIndicator, prev ){ # h2 estimation with a fixed effect Env
  y = fam_pheno$pheno
  
  famData$pheno = y # scale(y)[,1], DO NOT z-score!
  write.table(famData,paste0(outLoc,".phen"),col.names = F,row.names = F, quote=F)
  

  resultloc = paste0(outLoc,"_env.hsq")
  if(file.exists(resultloc) ) {file.remove(resultloc) } # delete any previous runs
  
  famData$pheno = scale( envIndicator)[,1] # envIndicator # EnvTerm # EnvTerm_z
  #write.table(famData,paste0(outLoc,".qcovar"),col.names = F,row.names = F, quote=F)
  
  # write covariates into files 
  write.table(cbind.data.frame(famData[,1], famData[,2], famData$pheno,fam_pheno$age, fam_pheno$age2),paste0(outLoc,".qcovar"),col.names = F,row.names = F, quote=F)
  write.table(cbind.data.frame(famData[,1], famData[,2],fam_pheno$sex.y),paste0(outLoc,".covar"),col.names = F,row.names = F, quote=F)
  
  
  
  # use "covar" instead of "qcovar", as we are using a binary indicator env effect not a quantitative one
  cmd=paste0(gctaoloc," --thread-num 8 --reml-est-fix --reml-alg 1  --reml --prevalence ",prev," --grm ",outLoc," --pheno ",outLoc,".phen --qcovar ", outLoc, ".qcovar --covar ",outLoc,".covar --out ",paste0(outLoc,"_env") )
  system(cmd) 
  
  
  # check if file exist, and if it does not then return dummy values
  if(file.exists(resultloc)  == F) {return(c(-1,-1, -1)) }
  
  # as binary outcomes are transferred back to the liability scale, we read those directly
  h2 = readLines(file(resultloc, "r"), 8)[8] # read in the 9th line which contains the h2 estimate
  h2 = as.numeric(strsplit(h2, "\t")[[1]] [2] ) # clean it
  # reverse engineer the liability scale Ve, which is not reported by GCTA:
  Ve=1-h2 # as we dont have GxE here
  
  FEnv = readLines(file(resultloc, "r"), 18)[18] # read in the 5th line which contains the h2 estimate
  FEnv = as.numeric(strsplit(FEnv, "\t")[[1]] [1] ) # clean it
  FEnv = FEnv / sd(envIndicator) # restore Env effect to its original scale
  
  
  remove(K)
  gc()
  return(c(h2,Ve, FEnv))
}



# K = K_standard
# y = y_gxe
# envIndicator = EnvTerm
estimate_VgVe_Env_GxE_bin = function(K, famData, outLoc, fam_pheno, envIndicator, prev ){ # h2 estimation with both fixed effect Env + GxE
  y = fam_pheno$pheno
  
  Fmodel = lm(y ~ envIndicator)
  summary(Fmodel) # this recovers the correct Env Beta, despite the interaction
  
  famData$pheno = y#   Fmodel$residuals # y # scale(y)[,1], DO NOT z-score!
  write.table(famData,paste0(outLoc,".phen"),col.names = F,row.names = F, quote=F)
  resultloc = paste0(outLoc,".hsq")
  if(file.exists(resultloc)  ) {file.remove(resultloc) } # delete any previous runs
  
  famData$pheno = envIndicator# scale( envIndicator)[,1] # envIndicator # EnvTerm # EnvTerm_z
  write.table(famData,paste0(outLoc,".env"),col.names = F,row.names = F, quote=F) # make sure this goes int oa separate file, and the qcovar file does NOT have the gxe term too
  
  # write the other covariates out too
  write.table(cbind.data.frame(famData[,1], famData[,2],fam_pheno$age, fam_pheno$age2),paste0(outLoc,".qcovar"),col.names = F,row.names = F, quote=F)
  write.table(cbind.data.frame(famData[,1], famData[,2],fam_pheno$sex.y),paste0(outLoc,".covar"),col.names = F,row.names = F, quote=F)
  

  # remove: --qcovar ", outLoc, ".qcovar, otherwise will get error "Error: analysis stopped because more than half of the variance components are constrained. The result would be unreliable"
  cmd=paste0(gctaoloc," --thread-num 8 --reml-est-fix --reml-alg 1  --reml --prevalence ",prev," --reml-lrt 2 --grm ",outLoc," --pheno ",outLoc,".phen --gxe ",outLoc,".env --qcovar ", outLoc, ".qcovar --covar ",outLoc,".covar  --out ", paste0(outLoc) )
  
  system(cmd) # https://www.r-bloggers.com/2021/09/how-to-use-system-commands-in-your-r-script-or-package/
  
  
  # check if file exist, and if it does not then return dummy values
  if(file.exists(resultloc)== F ) {return(c(-1,-1,-1,-1)) }
 
  # as binary outcomes are transferred back to the liability scale, we read those directly
  h2 = readLines(file(resultloc, "r"), 12)[12] # read in the 12th line which contains the h2 estimate
  h2 = as.numeric(strsplit(h2, "\t")[[1]] [2] ) # clean it
  
  
  VGxE = readLines(file(resultloc, "r"), 13)[13] # read in the 13th line which contains the h2 estimate
  VGxE = as.numeric(strsplit(VGxE, "\t")[[1]] [2] ) # clean it
  
  # reverse engineer the liability scale Ve, which is not reported by GCTA:
  Ve=1-h2 - VGxE # 



  FEnv = readLines(file(resultloc, "r"), 25)[25] # read in the 5th line which contains the h2 estimate
  FEnv = as.numeric(strsplit(FEnv, "\t")[[1]] [1] ) # clean it
  

  remove(K)
  gc()
  return(c(h2,Ve, FEnv, VGxE))
}


performAllTests_bin = function(K_standard, famData, outLoc, envIndicator, fam_pheno, prev) {
  gc()
  modelVG = estimate_VgVe_bin(K_standard, famData, outLoc, fam_pheno, prev)
  gc()
  modelVGEnv = estimate_VgVe_Env_bin(K_standard, famData, outLoc, fam_pheno, envIndicator, prev)
  gc()
  modelVGEnvGxE = estimate_VgVe_Env_GxE_bin(K_standard, famData, outLoc, fam_pheno, envIndicator, prev)
  remove(K_standard)
  gc()
  
  return (list ("modelVG" = modelVG, "modelVGEnv" = modelVGEnv, "modelVGEnvGxE" = modelVGEnvGxE ))
}


calcKinship = function(snpAttached, X_current, fam_current, outputloc, isPopKin = F) {

  if(isPopKin == T) {
    return (calcKinshipViaPopKin(X_current,fam_current, outputloc))
  } else {
    return (calcKinshipViaGCTA(snpAttached, X_current, fam_current, outputloc))
  }

}

#X_current = X_all
#outputloc = paste0(outLoc, "_popkintest")
#fam_current <- cbind.data.frame(1:nrow(X_current), 1:nrow(X_current))
#K_PopKin_loaded = read_grm(KinLoc)$kinship
#K_PopKin_loaded2 = read_grm(KinLoc, size_bytes = 8 )$kinship

calcKinshipViaPopKin = function(X_current, fam_current, outputloc) {
  gc()
  KinLoc = paste0(outputloc,"_POPKIN")
  print("calcKinshipViaPopKin() called")
  
  if (file.exists(paste0(KinLoc,".grm.bin")) == F) {
    print("Kinship matrix did not exit, we create it now")
  K_PopKin =  2* popkin::popkin(X_current, loci_on_cols = TRUE, mem_lim = 40) # apply the coeficient to relatedness conversion
  fam = fam_current[,1:2]
  colnames(fam) = c("fam", "id")
  # write popkin to disk in GCTA format
  write_grm(KinLoc, K_PopKin, fam = fam) # , size_bytes = 8, popkin defaults to 8 bytes per entry, but its probably OK to store it at 4 bytes same as GCTA
 
  # free up RAM
  remove(K_PopKin)
  gc()
  } else {
    print("Kinship matrix already calculated, loading it from disk")
  }
  return( read_grm(KinLoc)$kinship )

}


calcKinshipViaGCTA = function(snpAttached, X_current, fam_current, outputloc) {
  print("calcKinshipViaGCTA() called")
  dummyPlinkFile = cbind.data.frame(fam_current$sample.ID,fam_current$sample.ID,0,0,0,-9)             
  colnames(dummyPlinkFile) = c("fam", "id", "pat", "mat", "sex", "pheno")   
  
  dummyBimFile= snpAttached$map
  colnames(dummyBimFile) = c("chr", "id", "posg", "pos", "ref", "alt")
  
  
  # to make sure pop labels are aligned to the randomly selected indis, we write this out too
  write.table(fam_current, paste0(outputloc,"labels_GCTA.txt"), sep="\t", row.names = F, col.names = T, quote = FALSE)
  
  # write the selected subset onto disk
  gc()
  plinkFile=paste0(outputloc,"_GCTAintermediate")
  GCTA_KinLoc = paste0(outputloc,"_GCTA")
  if (file.exists(paste0(GCTA_KinLoc,".grm.bin")) == F) {
    

  write_plink(
    plinkFile,
    t(X_current),  # genio assumes that SNPs are on rows...
    bim = dummyBimFile,
    fam = dummyPlinkFile
  )
  #  generate GRM via GCTA CLI
  
  # if there are say, more than 400K SNPs, then we will likely get an OOM, so we do it in parts then
  numSNP= ncol(X_current)
  if (numSNP > 400000) {
    gc()
    print("more  than 400K SNPs, doing it in 3 parts") # https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 1 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 2 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 3 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 4 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 5 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 6 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 7 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 8 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 9 --out ",GCTA_KinLoc); system(cmd) 
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-part 10 10 --out ",GCTA_KinLoc); system(cmd) 
    
    
    print("copying together the parts")
    cmd=paste0("cat ",GCTA_KinLoc,".part_10_*.grm.id > ",GCTA_KinLoc,".grm.id"); system(cmd) 
    cmd=paste0("cat ",GCTA_KinLoc,".part_10_*.grm.bin > ",GCTA_KinLoc,".grm.bin"); system(cmd) 
    cmd=paste0("cat ",GCTA_KinLoc,".part_10_*.grm.N.bin > ",GCTA_KinLoc,".grm.N.bin"); system(cmd) 
    print(paste0(file.exists(paste0(GCTA_KinLoc,".grm.bin")), ": GRM exists!") ) # delete any previous runs
  
    # cleanup intermediate files to save space
    for (i in 1:10) {
      unlink(paste0(GCTA_KinLoc,".part_10_",str_pad(i,2, pad = "0"),".grm.id"))
      unlink(paste0(GCTA_KinLoc,".part_10_",str_pad(i,2, pad = "0"),".grm.bin"))
      unlink(paste0(GCTA_KinLoc,".part_10_",str_pad(i,2, pad = "0"),".grm.N.bin"))
      unlink(paste0(GCTA_KinLoc,".part_10_",str_pad(i,2, pad = "0"),".log"))
    }
    unlink(paste0(plinkFile,".bed"))
    unlink(paste0(plinkFile,".fam"))
    unlink(paste0(plinkFile,".bim"))
    
    } else {
    cmd=paste0(gctaoloc," --thread-num 8 --bfile ",plinkFile," --make-grm-bin --out ",GCTA_KinLoc)
    system(cmd) 
  }

  gc()
  } else { print(paste0(file.exists(paste0(GCTA_KinLoc,".grm.bin")), ": ALREADY exists!") )}
  
  
  return( read_grm(GCTA_KinLoc)$kinship )
  
}


#fam2_ordered = fam_current
#fam2_ordered$pop = "pop1"
#fam2_ordered$pop[ (nrow(fam2_ordered)/2 + 1):nrow(fam2_ordered)] = "pop2"
#
#fam2_ordered$pop[ (nrow(fam2_ordered)/3 + 1):( 2* nrow(fam2_ordered)/3)] = "pop2"
#fam2_ordered$pop[ (2* nrow(fam2_ordered)/3 + 1):nrow(fam2_ordered)] = "pop3"
#X_ordered = X_all

# PopKin way of calculating Fst
calcFst_popKin = function(K_PopKin, fam2_ordered, outLoc){
  
  split_pop <- split(rows_along(fam2_ordered), fam2_ordered$pop) # this gets the indices of for the indis in each pop
  npop <- length(split_pop) # 79 how many (sub) populations in total
  dist_fst_popkin <- matrix(0, npop, npop) # create a blank distance matrix
  colnames(dist_fst_popkin) <- rownames(dist_fst_popkin) <- names(split_pop) # add colnames
  i=2
  j=1
  
  # double loop
  for (i in 2:npop) {
    print(paste0("processing pop ",i, " out of ", npop) )
    for (j in 1:(i - 1)) {
      ind <- unlist(split_pop[c(i, j)], use.names = FALSE) # grab the indis for a pair of population
      
      kinship_sub <- popkin::rescale_popkin(K_PopKin[ind, ind], fam2_ordered$pop[ ind ]) # must rescale the kinship matrix too by the same pops, otherwise we get inflation
      w <- popkin::weights_subpops( fam2_ordered$pop[ ind ] ) # get pairwise wieghts
      
      Fst_popkin <- popkin::fst(kinship_sub,w) # calculate pairise Fst
      Fst_popkin # 0.1038986, this is a lot closer
      dist_fst_popkin[i, j] <- dist_fst_popkin[j, i] <- print(Fst_popkin) # mirror the Fst into both places of the symmetric matrix
    }
  }
  #write.table(dist_fst_popkin, paste0(outLoc,"pairwiseFst_fixSNP.txt"), sep="\t", row.names = T, col.names = T, quote = FALSE)
  
  # use the 'wrong' filename for backward-compatibility, so that the same loops in bash will work
  write.table(dist_fst_popkin, paste0(outLoc,"_","popFSt_prive.txt"), sep="\t", row.names = T, col.names = T, quote = FALSE)
  
  
  # the below would produce individual-level pairwise Fst, whereas we want pop-level pairwise
  # kinship <- popkin::popkin(X_ordered, subpops = fam2_ordered$pop, loci_on_cols = TRUE)
  # 
  # # get Individual level Fst
  # pairwise_fst <- pwfst(kinship) # estimated matrix
  # mean(pairwise_fst) # 0.09634726, the mean pairwise Fst is the usual 10%
  # max(pairwise_fst) #  0.3797001, but goes up to 37%!
  # 
  # remove(kinship)
  # gc()
  # 
  # # write fst results
  # write.table(all_fst, paste0(outLoc,"_","popFSt_popKin.txt"), sep="\t", row.names = T, col.names = T, quote = FALSE)
}





calcFst_WC = function(obj.bed, fam_pheno, outLoc, pops){
  
  #####################
  # Prive's Fst 
  # https://github.com/privefl/freq-ancestry/blob/f5b03a2d639eb8244889f5f32a4c8fb2e1fcd2e9/code/1-prepare-1000G-freq.R
  
  all_freq = NULL
  for (i in 1:length(pops)) {
    fam_sub = fam_pheno[which(fam_pheno$pop == pops[i]),]
    popAFs = bed_MAF(obj.bed, ind.row = fam_sub$indi_index, ncores = nb_cores())$af
    all_freq = cbind(all_freq,popAFs)
  }
  
  # loop and calculate Fsts
  # https://github.com/privefl/freq-ancestry/blob/cba3ca3fa0b346a7e136b0374add563a46b288d2/code/4-filter-merge-freq.R
  
  all_N <- table(fam_pheno$pop)
  npop <- length(all_N)
  
  popFsts_c = c()
  #pops = c()
  all_fst <- matrix(0, npop, npop) # create a blank distance matrix
  colnames(all_fst) <- rownames(all_fst) <- pops # names(all_N) # add colnames names(all_N) would be wrong order!
  for (i in 2:npop) {
    print(paste0("processing pop ",i, " out of ", npop) )
    for (j in 1:(i - 1)) {
      
      fst <- snp_fst(list(
        data.frame(af = all_freq[, i], N = all_N[i]),
        data.frame(af = all_freq[, j], N = all_N[j])
      ), overall = TRUE)
      all_fst[i, j] <- all_fst[j, i] <- print(fst)
    }
  }
  

  write.table(all_fst, paste0(outLoc,"_","popFSt_prive.txt"), sep="\t", row.names = T, col.names = T, quote = FALSE)
}


