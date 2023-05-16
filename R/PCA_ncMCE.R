#############################
# Performs PCA and then ncMCE
#############################

# imports
# install.packages("hexbin")
# install.packages("igraph");
# install.packages("parallelDist")
# install.packages("rsvd")
# install.packages("bigsnpr")
# install.packages("bigstatsr")
# ##install.packages("ggplot2")

library(bigstatsr)
library(bigsnpr)
#library(ggplot2)
library(hexbin)
library(igraph)
#library(parallelDist)
library(rsvd)


options(error=traceback)
# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")}


# arguments
inputBedLoc=args[1] # location of the .bed file (should already have basic QC)
outLoc=args[2] # the folder and filename stem where the results are to be saved
numCores=as.numeric(args[3]) # how many cores to use
numDims=as.numeric(args[4]) # the number of dimensions to extract
runncMCE = F
if ( is.na(args[5]) == F && args[5] == "1") { runncMCE = T}
centring= F

#if ( is.na(args[5]) == F && args[5] == "1") { centring = T}


########################
# I PCA via bigSNPR
########################
set.seed(42)
outname=paste0(outLoc,"_PC")
PCAresName = paste0(outname,"_PCA")
usedIndisLoc = paste0(outLoc,"usedIndis.txt")
usedIndis_idxLoc = paste0(outLoc,"usedIndis_idx.txt")
usedSNPs_idxLoc =  paste0(outLoc,"usedSNPs_idx.txt")
initialSNPfilterLoc=paste0(outLoc,"initialSNPFilter_idx.txt")
initialIndifilterLoc=paste0(outLoc,"initialIndiFilter_idx.txt")
pcaSVDLoc=paste0(outLoc,"_svd")
bigparallelr::set_blas_ncores(numCores)

#tempDirLoc_PCA=paste0(outLoc,"temp_PCA/")
tempDirLoc=paste0(outLoc,"temp/")



# 1. PLINK2
plink2 <- download_plink2(outLoc)


(obj.bed <- bed(inputBedLoc)) # this needs to be loaded for both PCa and ncMCE

if( file.exists(PCAresName) == F) {
  print("Starting PCA")

# 2. initial QC
#tmp_PCA <- tempfile(tmpdir = tempDirLoc_PCA)
#bedfile <- snp_readBed(inputBedLoc, backingfile = tmp_PCA) # , backingfile = sub_bed(bedfile)

# check if initial QC was not already run
if( file.exists(initialSNPfilterLoc) == F) {
  print("Performing initial QC step")

start_time <- Sys.time()
# First, let us detect all pairs of related individuals.
rel <- snp_plinkKINGQC(
  plink2.path = plink2,
  bedfile.in = inputBedLoc,
  thr.king = 2^-4.5,
  make.bed = FALSE,
  ncores = numCores #nb_cores()
)


# get rid of too relateds
ind.rel <- match(c(rel$IID1, rel$IID2), obj.bed$fam$sample.ID)  # /!\ use $ID1 instead with old PLINK
ind.norel <- rows_along(obj.bed)[-ind.rel]


# 2. initial PCA (which performs QC)
obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.norel, k = numDims,ncores = numCores) # nb_cores()

# look for outliers
prob <- bigutilsr::prob_dist(obj.svd$u, ncores = numCores) # nb_cores()
S <- prob$dist.self / sqrt(prob$dist.nn)


# 3. rerun PCA without outliers
ind.row <- ind.norel[S < 0.5]
ind.col <- attr(obj.svd, "subset")

# write partial results
write.table(ind.row, initialIndifilterLoc, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ind.col, initialSNPfilterLoc, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
end_time <- Sys.time(); print(end_time - start_time)
} else {
  print("Loading previous run of QC SNP filter")
  ind.row = unlist(read.table(initialIndifilterLoc,  header = FALSE) )
  ind.col = unlist(read.table(initialSNPfilterLoc,  header = FALSE) )
}

print("Starting final PCA")
start_time <- Sys.time()

obj.svd2 <- bed_autoSVD(obj.bed, ind.row = ind.row,ind.col = ind.col, thr.r2 = NA,k = numDims, ncores = numCores)  # nb_cores()
ind.col2 <- attr(obj.svd2, "subset") # get the final SNP list that excludes all long range LD

# get PCs
usedIndis=obj.bed$.fam$sample.ID[ind.row]
usedSNPs =obj.bed$.map$marker.ID[ind.col2]
PCs <- predict(obj.svd2, ind.row = ind.row) # get the actual PCs via predict
PCmat_indis = cbind(usedIndis,PCs) # concat it with the indi IDs

# here I should probably also save the loadings, so that we could generate PCs for the full dataset via OADP
# write the entire obj.svd2 to disk
saveRDS(obj.svd2, file = pcaSVDLoc ) 

# write PCs and final keeplists for both indis and SNPs to disk
write.table(PCmat_indis,PCAresName , sep = "\t", row.names = F, col.names = F, quote = FALSE)
print(paste("written PCs for ",nrow(PCmat_indis)," indis to", PCAresName))

write.table(usedIndis, usedIndisLoc, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(usedSNPs, paste0(outLoc,"usedSNPs.txt"), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(ind.row, usedIndis_idxLoc, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ind.col2,usedSNPs_idxLoc, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
end_time <- Sys.time(); print(end_time - start_time)


# delete temp files
#unlink(tempDirLoc_PCA, recursive = TRUE)

# plot results as a sanity check
png(filename=paste0(outname,"PCA.png")  , width=1280, height=720, res=128)
plot(PCs[,1],PCs[,2], main="PCA", xlab="PC1", ylab="PC2")
dev.off()

if(numDims >= 4) {
  png(filename=paste0(outname,"PCA_3_4.png")  , width=1280, height=720, res=128)
  plot(PCs[,3],PCs[,4], main="PCA", xlab="PC3", ylab="PC4")
  dev.off()
}

} else {
  print("PCA already run, we load final keeplist results from disk")
  usedIndis=unlist(read.table(usedIndisLoc,  header = FALSE))
  ind.row = unlist(read.table(usedIndis_idxLoc,  header = FALSE) )
  ind.col2 = unlist(read.table(usedSNPs_idxLoc,  header = FALSE) )
}


########################
# II ncMCE
########################
print("Starting ncMCE")

if(runncMCE) {
  
  big_tcrossprodSelf2 <- function (X,tmp_K, fun.scaling = big_scale(center = FALSE, scale = FALSE),
                                   ind.row = rows_along(X), ind.col = cols_along(X), block.size = block_size(nrow(X)))
  {
    bigstatsr:::check_args()
    n <- length(ind.row)
    K <- FBM(n, n, init = 0, backingfile = tmp_K)
    m <- length(ind.col)
    means <- numeric(m)
    sds <- numeric(m)
    intervals <- bigstatsr:::CutBySize(m, block.size)
    X_part_temp <- matrix(0, n, max(intervals[, "size"]))
    for (j in rows_along(intervals)) {
      ind <- bigstatsr:::seq2(intervals[j, ])
      ind.col.ind <- ind.col[ind]
      ms <- fun.scaling(X, ind.row = ind.row, ind.col = ind.col.ind)
      means[ind] <- ms$center
      sds[ind] <- ms$scale
      bigstatsr:::increment_scaled_tcrossprod(K, X_part_temp, X, ind.row,
                                              ind.col.ind, ms$center, ms$scale)
    }
    structure(K, center = means, scale = sds)
  }
  
  
  
  distantMatLoc=paste0(outLoc,"distanceMat")
  if( file.exists(distantMatLoc) == F) {
    print("calculating Distance Matrix")
    
    #set.seed(42) # setting the seed does not affect tempdir creation
    # remove any previously failed runs
    unlink(tempDirLoc, recursive = TRUE)
    tmp <- tempfile(tmpdir = tempDirLoc)
    tmp_fbm <- tempfile(tmpdir = tempDirLoc)
    
    tmp_K <- tempfile(tmpdir = tempDirLoc)
    
    # generate the euclidean distance matrix, I used to do this via parDist, but that does not work for large matrices
    # bed_tcrossprodSelf also does not work, as it is either very slow and/or generates NA results (probably has to do with the fun.scaling)
    # so I decided to copy the SNP data into a regular FBM and do it via big_tcrossprodSelf, as that has no scaling by default, and is also 3x faster
    # the cross product based formula for the euclidean distance matrix is from : https://hal.archives-ouvertes.fr/hal-02047514/document
    # read the SNP data in
    snpData = snp_readBed2(inputBedLoc,
                           backingfile = tmp, # sub_bed(bedfile)
                           ind.row = ind.row,
                           ind.col = ind.col2,
                           ncores = numCores)
    
    # Loading the data from backing files
    snpAttached <- snp_attach(snpData)
    
    # need to make sure there aren't missing values
    imputedgenotypes = snp_fastImputeSimple(snpAttached$genotypes)  # https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html   # https://github.com/privefl/bigsnpr/issues/124
    
    # convert to fbm format
    fmbSNP = big_copy(imputedgenotypes, backingfile = tmp_fbm) # https://privefl.github.io/bigstatsr/reference/big_copy.html 
    
    
    start_time <- Sys.time()
    bigparallelr::set_blas_ncores(numCores)
    #K <- big_tcrossprodSelf(fmbSNP) # https://www.rdocumentation.org/packages/bigstatsr/versions/1.5.0/topics/big_tcrossprodSelf
    K <- big_tcrossprodSelf2(fmbSNP, tmp_K)
    
    
    x = K[1:nrow(K),1:ncol(K)] # flatten
    e = rep(1,nrow(x)) # just a column of 1s
    s = diag(x) # a column of the square sums: which is the diagonals of the Gram are the ||x||^2 of each column
    x = s %*%t(e) + e %*% t(s) -2*x
    x = sqrt(x) # get the square root, as the above gives the squared dists
    
    end_time <- Sys.time(); print(end_time - start_time)
    
    
    # save this to disk as an R object
    saveRDS(x, file = distantMatLoc ) 
    
    # delete temp files
    unlink(tempDirLoc, recursive = TRUE)
    
  } else {
    print("loading distance matrix from disk")
    # load distance matrix from disk into RAM:
    x <- readRDS(distantMatLoc)
  }
  
  
  # perform ncMCE Kernel generation,  each step if not already run
  #Make sure the matrix is symmetric
  x <- pmax(x, t(x));
  
  graphLoc=paste0(outLoc,"_g")
  if( file.exists(graphLoc) == F) {
    print("calculating graph")
    start_time <- Sys.time()
    
    #Create a graph object out of the adjacency matrix x
    g <- graph.adjacency(x, mode = "undirected", weighted = TRUE);
    
    #sparse_matrix <- as(x, "sparseMatrix")
    
    #g <- graph_from_adjacency_matrix(sparse_matrix, mode = "undirected", weighted = TRUE)
    
    saveRDS(g, file = graphLoc )  # save to disk
    end_time <- Sys.time(); print(end_time - start_time)
  } else {
    print("loading graph")
    g <- readRDS(graphLoc)
  }
  
  mstLoc=paste0(outLoc,"_mst")
  if( file.exists(mstLoc) == F) {
    print("calculating mst")
    start_time <- Sys.time()
    
    #MC-kernel computation
    mst <- minimum.spanning.tree(g);
    saveRDS(mst, file = mstLoc )  # save to disk
    end_time <- Sys.time(); print(end_time - start_time)
  } else {
    print("loading mst")
    mst <- readRDS(mstLoc)
  }
  
  kernelLoc=paste0(outLoc,"_kernel")
  if( file.exists(kernelLoc) == F) {
    print("calculating kernel")
    start_time <- Sys.time()
    
    kernel <- shortest.paths(mst);
    #Kernel centring
    if(centring == T){
      N <- nrow(kernel);
      J <- diag(N) - (1/N)*matrix(1, N, N); #Form the centring matrix J
      kernel <- (-0.5)*(J %*% kernel^2 %*% J);
    }
    saveRDS(kernel, file = kernelLoc )  # save to disk
    end_time <- Sys.time(); print(end_time - start_time)
  } else {
    print("loading kernel")
    kernel <- readRDS(kernelLoc)
  }
  
  
  # perform randomized SVD, if not already run
  resLoc=paste0(outLoc,"_res")
  if( file.exists(resLoc) == F) {
    print("calculating ncMCE SVD")
    start_time <- Sys.time()
    #SVD-based Embedding
    #res <- svd(kernel)
    res = rsvd(kernel, k=numDims, nv = numDims, nu = 0) # dont need left singular vector
    saveRDS(res, file = resLoc )  # save to disk
    end_time <- Sys.time(); print(end_time - start_time)
  } else {
    print("loading ncMCE SVD")
    res <- readRDS(resLoc)
  }
  
  L <- diag(res$d) # singular values
  V <- res$v # right singular vectors 
  # u would be left singular vectors
  
  sqrtL <- sqrt(L[1:numDims, 1:numDims]) # take sqrt of them to get eigen values
  V <- V[, 1:numDims] # truncate this to length
  
  s <- t(sqrtL %*% t(V)) # this is numIndis x numPCs i dont understand as he gets his PCs as eigenVal * right Singular vector, whereas usually you get this from Left_Singular_Vector * singularVal 
  
  
  # write results to disk
  MCE_indis = cbind(usedIndis,s) # concat it with the indi IDs
  
  MCEresName = paste0(outname,"_ncMCE")
  # write PCs and final keeplists for both indis and SNPs to disk
  write.table(MCE_indis,MCEresName , sep = "\t", row.names = F, col.names = F, quote = FALSE)
  print(paste("written ",nrow(MCE_indis)," ncMCEs to", MCEresName))
  
  # plot results as a sanity check
  png(filename=paste0(outname,"MCE.png")  , width=1280, height=720, res=128)
  plot(s[,1],s[,2], main="ncMCE", xlab="dim1", ylab="dim2")
  dev.off()
  
  if(numDims >= 4) {
    png(filename=paste0(outname,"MCE_3_4.png")  , width=1280, height=720, res=128)
    plot(s[,3],s[,4], main="ncMCE", xlab="dim3", ylab="dim4")
    dev.off()
  }
  
  # need to figure out what is the 'projection matrix' or the equivalent to 'SNP loadings' for ncMCE
  
} else {print("ncMCE skipped") }

