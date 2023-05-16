#############################
# Performs PCA / ncMCE visualisation
#############################



#options(error=traceback)
# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 5) {stop("not enough Arguments received")}


# arguments
inputPCsLoc=args[1]
PCsToVisualiseLoc=args[2] # list of predictors who match the colnames of the first arg
groupLabelsLoc=args[3] # each signature IID, group, color (eg HS0001, "case", "red")
outPutLoc=args[4]
dimName=args[5] # the name of the dimension ie "PC" or "Dim"
plotTitle=""
if(length(args) == 6){ # title is optional
  plotTitle=args[6]
  plotTitle = gsub("_", " ", plotTitle) #remove underscores
  
}


# load PCs (or MCEs)
PCs = read.table(inputPCsLoc, header=F) 

# load list to visualise
PCsToVisualise= read.table(PCsToVisualiseLoc, header=F) 
if(nrow(PCsToVisualise) < 2) {
  tryCatch(stop("not enough PCs to visualise"), finally = {cat("")})
}


# match the column names to the PCs
indicesToPick = match(PCsToVisualise$V2, colnames(PCs))

# subset the PCs to the 
PCs_needed = PCs[,c(1,indicesToPick)] # keep the first one too, as that has the IIDs

# rename the colnames of the PC dataframe to be '1', '2', ... instead of V2, so that we could easily construct the axis labels

renamedCols = colnames(PCs_needed)
renamedCols = as.numeric( gsub("V", "", renamedCols) ) -1
#renamedCols
colnames(PCs_needed) = c("IID", renamedCols[2:length(renamedCols)])
#colnames(PCs_needed) 
# load group labels
groupLabels = read.table(groupLabelsLoc, sep="\t", header=T, comment.char = "") # need this otherwise the colors with #123456 wouldn't be able to be read

# merge group labels with the PCs
labels_PCs = merge(groupLabels, PCs_needed, by = "IID")


# need to create legend names / colors from the final data table

legendNames = names( table(labels_PCs$group) ) # problem is that table does not return the colors/groups matched, so I cannot just do legendCols = names( table(labels_PCs$color) )
legendCols = c()
# so I will just use the group names as a lookup to get the colors
for (i in 1:length(legendNames)) {
  groupColors = labels_PCs$color[labels_PCs$group == legendNames[i]]
  legendCols = c(legendCols, groupColors[1]) # just grab the first one
}

legendNcol = round(length(legendCols) / 2)


# loop in 2s of PCs to visualise (if odd number then we visualise last-1 and last)
for(i in seq(from=4, to=ncol(labels_PCs), by=2)){
  if(i == ncol(labels_PCs)) { # if we are at the last predictor to visualise, then we take one step back so that we could have two
    #print("last PC!")
    i = i-1
  }
  cname1=colnames(labels_PCs)[i]
  cname2=colnames(labels_PCs)[(i+1)]
  print(paste0(i, ": ",cname1, " / ", cname2 ))
  xName = paste0(dimName,cname1)
  yName = paste0(dimName,cname2)
  
  filen = paste0(outPutLoc,xName,"_", yName,".png")
  png(filename=filen  , width=1280, height=720, res=128)
  par(mar=c(4.1, 4.1, 7, 4.1), xpd=TRUE)
  plot(labels_PCs[,i],labels_PCs[, (i+1) ], col=labels_PCs$color, xlab=xName, ylab=yName)

  title(plotTitle, line = 4.5) # , main=plotTitle #  line = 3.5
  legend(x="topright", legend =legendNames, inset=c(0,-0.22), col=legendCols,pch=16, ncol = legendNcol)
  dev.off()
  print(paste0("written plot to: ", filen) )
}


# 
# 
# # Random data to plot:
# A <- data.frame(x=rnorm(100, 20, 2), y=rnorm(100, 20, 2))
# B <- data.frame(x=rnorm(100, 21, 1), y=rnorm(100, 21, 1))
# 
# # Add extra space to right of plot area; change clipping to figure
# par(mar=c(4.1, 4.1, 6.1, 4.1), xpd=TRUE)
# 
# # Plot both groups
# plot(y ~ x, A, ylim=range(c(A$y, B$y)), xlim=range(c(A$x, B$x)), pch=1,
#      main="Scatter plot of two groups")
# points(y ~ x, B, pch=3)
# 
# # Add legend to top right, outside plot region
# legend("topright", inset=c(0,-0.3), legend=c("A","B"), pch=c(1,3), ncol=2)
# 
# 
