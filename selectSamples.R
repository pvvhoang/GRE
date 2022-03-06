library("optparse")

# Rscript ../selectSamples.R --fam_file=${fam_file} --tra_file="expression_id.txt"  --number=50000 --outFile1="selectedFam.fam" --outFile2="selectedExpressionId.txt"
# Output: 
# selectedFam.fam
# selectedExpressionId.txt

option_list = list(
  make_option(c("-f", "--fam_file"), type = "character"),
  make_option(c("-t", "--tra_file"), type = "character"),
  make_option(c("-n", "--number"), type = "numeric"),
  make_option(c("-s", "--outFile1"), type = "character"),
  make_option(c("-e", "--outFile2"), type = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #==================================================
# opt$fam_file="qced_rdm_005_ukbb3.fam"
# opt$tra_file="expression_id.txt"
# opt$number=50000
# opt$outFile1="selectedFam.fam"
# opt$outFile2="selectedExpressionId.txt"
# #==================================================

#-------------------------------
# Select samples

fam <- read.table(opt$fam_file, stringsAsFactors=F, header=F)
tra <- read.table(opt$tra_file, stringsAsFactors=F, header=T)
int <- intersect(fam$V2, tra$IID)
total <- length(int)

set.seed(2)
id <- sample(1:total, opt$number)
idv = int[id]
idv <- sort(idv)

filename = opt$outFile1
writeLines(paste("Writing file: ", filename, sep = ""))
sink(filename)
write.table(fam[which(fam$V2 %in% idv),], row.names = F, quote = F, col.names = F)
sink()

filename = opt$outFile2
writeLines(paste("Writing file: ", filename, sep = ""))
sink(filename)
write.table(tra[which(tra$IID %in% idv),], row.names = F, quote = F, col.names = F)
sink()
#-------------------------------