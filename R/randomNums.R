# gets a seeded random list of numbers

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 



numItems= as.numeric(args[1])
mySeed = as.numeric(args[2])
outputLoc = args[3]

#numItems=100
#mySeed=42
#outputLoc="C:/0Datasets/shaPRS/0shaPRS/0ukbb/results_190821/sims_eval_ukbb_rnd"

set.seed(mySeed)

myArray = rep(1:numItems)
myArray = sample(myArray)


write.table(myArray,outputLoc, col.names = F, row.names = F, quote = F)

print(paste0("written ", numItems, " randomised numbers to: ",outputLoc ))
