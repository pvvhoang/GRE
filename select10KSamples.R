library("optparse")

# Rscript ../select10KSamples.R --fam_file=${fam_file} --tra_file="expression_id.txt"
# Output: 
# fam_10K.fam
# expression_id_10K.txt

option_list = list(
  make_option(c("-f", "--fam_file"), type = "character"),
  make_option(c("-t", "--tra_file"), type = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #==================================================
# opt$fam_file="qced_rdm_005_ukbb3.fam"
# opt$tra_file="expression_id.txt"
# #==================================================

#-------------------------------
# Select 10K samples

fam <- read.table(opt$fam_file, stringsAsFactors=F, header=F)
tra <- read.table(opt$tra_file, stringsAsFactors=F, header=T)
int <- intersect(fam$V2, tra$IID)
total <- length(int)

set.seed(2)
id <- sample(1:total, 10000)
idv = int[id]
idv <- sort(idv)

filename = "fam_10K.fam"
writeLines(paste("Writing file: ", filename, sep = ""))
sink(filename)
write.table(fam[which(fam$V2 %in% idv),], row.names = F, quote = F, col.names = F)
sink()

filename = "expression_id_10K.txt"
writeLines(paste("Writing file: ", filename, sep = ""))
sink(filename)
write.table(tra[which(tra$IID %in% idv),], row.names = F, quote = F, col.names = F)
sink()
#-------------------------------