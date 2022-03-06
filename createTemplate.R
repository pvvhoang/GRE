library("optparse")
library(MASS)

# Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
# Output: 
# outFile

option_list = list(
  make_option(c("-i", "--inFile"), type = "character"),
  make_option(c("-o", "--outFile"), type = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #=================================================
# opt$inFile="brain_coregreml.out2"
# opt$outFile="brain_coregreml.do"
# #=================================================
  
sink(opt$outFile)
cat(opt$inFile, "\n") # line 1: specify the file with parameter estimates
cat(4, "\n") # line 2: tot. # of variance & covariance components in the file
cat("R 2 1 3 4 4", "\n") # line 3: compute prop. of variance due to genetics ('R' is to get the ratio of variance #2, i.e. #2 / (#1 + #2 + #3 + #4 + #4))
cat("R 3 1 2 4 4", "\n") # line 4: compute prop. of variance due to transcriptome
cat("C 4 2 3", "\n") # line 5: compute correlation between g & b ('C' is to get the correlation of covariance #4, i.e. #4 / sqrt(#2 * #3))
sink()
  