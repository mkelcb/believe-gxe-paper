
# Code repository for "Investigating the effect of gene-country interactions on health and anthropometric traits in South Asian populations"

DOI: XXX

This respository represents the latest snapshot of the bash and R scripts used to generate our results and is provided as-is. As the analysis involved a lot of input/output operations on very large files, these were generated asynchronously on a cluster. Thus these scripts are meant to be executed on the command line manually, block-by-block, waiting for the remote jobs to finish and verifying the integrity of the resulting files at each step. To reduce code duplication, certaint functions that were reused multiple times are defined only once across all files, however, they may be called from different scripts.

Code for the paper can be found under bash/: 

1. BELIEVE_functions.sh: variables and reusable functions
2. BELIEVE_pub3.sh: data pre-processing and analyses


Auxilliary R scripts used by the above scripts can be found under R/.

