library("optparse")

#Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --selected_id=${selected_id} --out_file=${out_file} --pre=${pre}

# Out
# out_file

option_list = list(
  make_option(c("--pheno_file"), type = "character"),
  make_option(c("--selected_id"), type = "character"),
  make_option(c("--out_file"), type = "character"),
  make_option(c("--pre"), type = "character", default="")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #==================================================
# opt$pheno_file="edu_trait.txt"
# opt$selected_id="selectedFam.fam"
# opt$out_file="edu_random_50K.dat"
# opt$pre="edu"
# #==================================================

#-------------------------------
# Read data
writeLines("Reading data ...")

fam=read.table(opt$selected_id)
head(fam)
nrow(fam)
# > head(fam)
# V1      V2 V3 V4 V5 V6
# 1 1000105 1000105  0  0  1 -9
# 2 1000137 1000137  0  0  2 -9
# 3 1000269 1000269  0  0  1 -9
# 4 1000331 1000331  0  0  1 -9
# 5 1000362 1000362  0  0  1 -9
# 6 1000443 1000443  0  0  1 -9
# > nrow(fam)
# [1] 50000

phe=read.table(opt$pheno_file)
colnames(phe) <- c("id", "trait")
head(phe)
nrow(phe)
# > head(phe)
# id trait
# 1      id   edu
# 2 2080762    10
# 3 3165676    20
# 4 1605581    10
# 5 5336523    20
# 6 2211334    15
# > nrow(phe)
# [1] 502648

#-------------------------------

#-------------------------------
# Check with fam file

# trait
m1 <- match(fam$V1, phe$id)
dat0 <- phe[m1,]
dat0 <- cbind(dat0[,1], dat0)
# > head(dat0)
# dat0[, 1]      id trait
# 460023   1000105 1000105    19
# 51797    1000137 1000137    20
# 129545   1000269 1000269    20
# 294888   1000331 1000331    19
# 31297    1000362 1000362    13
# 104249   1000443 1000443    15

# trait <- as.numeric(as.character(dat0[,3]))
# # > head(trait)
# # [1] 19 20 20  7 19 19
# # > length(trait)
# # [1] 50000

#-------------------------------

#-------------------------------
# Writing files

writeLines("Writing files ...")

pre <- opt$pre
if(pre != "") {
  pre <- paste(pre, "_", sep = "")
}

writeLines(paste("Writing file: ", opt$out_file, sep = ""))
sink(opt$out_file)
write.table(dat0,quote=F,col.name=F,row.name=F, sep = "\t")
sink()
#-------------------------------
