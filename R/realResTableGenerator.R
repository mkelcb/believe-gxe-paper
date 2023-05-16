# Aggregates results into a table (hardcoded 10 traits and assumes 23 lines per results file!)
##################
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")}

baseLoc=args[1]
outLoc=args[2]

# DEBUG VARS
#baseLoc="C:/0LocalHPC/believe/0QC/realResults/"
#outLoc="C:/0LocalHPC/believe/0QC/realResults/table"

myarray=c( "height", "bmi" ,"highCholesterol", "edu", "MI", "stroke", "t2d", "angina", "sbp", "dbp" )


textProcess = function(lines, varCompLine,pValLine = F) {
  h2 = lines[varCompLine]
  h2 =unlist( strsplit(h2, " "))
  h2 = h2[length(h2)]
  h2 = strsplit(h2, "\"")
  h2 = as.numeric(h2)
  h2 = round(h2,3)
  h2 = as.character(h2)
  
  # add asterisk to significant results
  extraText= "" # assume that it is NOT significant by default
  if(pValLine != F) { # if there was a line for the pvalue specified we load it and check if its significant
    pval= lines[pValLine] 
    pval =unlist( strsplit(pval, " "))[2]
    pval = as.numeric(pval)
    if(pval <=0.05) {extraText = "*"}      
  }
  
  h2 = paste0(h2,extraText)   
  return(h2)
}

GxEResults = NULL
SeparatePopResults = NULL
GxEResults = cbind.data.frame(c("pheno", "h2", "Country(in UK)","GxE"))
SeparatePopResults = cbind.data.frame(c("pheno", "UK SAS", "BELIEVE"))

i=1
for (i in 1:length(myarray)) {
  fileLoc = paste0(baseLoc,myarray[i],"_res_results.txt")
  lines = readLines(file(fileLoc, "r"), 23) # read up until in the 5th line, but then return the 3rd line, which is Ve
 
  # create data frame for GxE results
  h2 = textProcess(lines,11,21)
  Env = textProcess(lines,13,F) 
  GxE = textProcess(lines,14,20) 
  GxEResults = cbind.data.frame(GxEResults, c(myarray[i],h2,Env,GxE))
    
  # create data frame for separate populations h2 results
  BELIEVE_h2 = textProcess(lines,16,22) 
  UK_SAS_h2 = textProcess(lines,18,23) 
  SeparatePopResults = cbind.data.frame(SeparatePopResults, c(myarray[i],UK_SAS_h2,BELIEVE_h2))
  
  
}


write.table(GxEResults,paste0(outLoc,"_gxe.csv"), sep=",",col.names = F,quote = F,row.names = F)

write.table(SeparatePopResults,paste0(outLoc,"_pop.csv"), sep=",",col.names = F,quote = F,row.names = F)

print(paste0("Written results to ",outLoc))
