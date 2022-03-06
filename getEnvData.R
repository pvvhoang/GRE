library("optparse")

#Rscript ../getEnvData.R --other_pheno_file=${other_pheno_file} --covariate_file=${covariate_file} --out_file=${out_file}

# Out
# out_file

option_list = list(
  make_option(c("--other_pheno_file"), type = "character"),
  make_option(c("--covariate_file"), type = "character"),
  make_option(c("--out_file"), type = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #==================================================
# opt$other_pheno_file="ukb6247.csv"
# opt$covariate_file="my_sample_qcinfo.txt"
# opt$out_file="env_data.txt"
# #==================================================

#-------------------------------
# Read data
writeLines("Reading data ...")

cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,$4,$3,$2353,$170,$18,$2354}' ", opt$other_pheno_file, sep = "") # eid, year of birth, sex, age, Townsend deprivation index at recruitment (TDI), centre, batch
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
# > head(dat)
# V1     V2     V3        V4        V5     V6        V7
# 1     eid 34-0.0 31-0.0 21022-0.0   189-0.0 54-0.0 22000-0.0
# 2 2080762   1958      1        49  -3.62553  11009      2000
# 3 3165676   1952      1        57 -0.402275  11020        15
# 4 1605581   1946      0        63  -4.00346  11017      2000
# 5 5336523   1953      1        54  -5.38678  11008      2000
# 6 2211334   1957      1        52 -0.248725  11009         3

cmd = paste("awk '{print $1,$2,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41}' ", opt$covariate_file, sep = "") # id & principal components
dat2 = read.table(pipe(cmd), header=T)
# > head(dat2)
# fid     iid      pc1     pc2       pc3       pc4      pc5       pc6
# 1 5438743 5438743 -12.1725 5.39163 -1.281030  0.841765 -5.26521 -1.786570
# 2 2890570 2890570 -13.0245 6.41514 -0.183365  2.927610 -5.88964  0.940534
# 3 2689958 2689958 -11.4712 3.48383 -1.154580  3.083830  7.65160 -0.913399
# 4 5800351 5800351 -12.1327 4.02976 -0.988080  0.750294 -2.36431  0.431658
# 5 5644632 5644632 -12.2171 3.50821 -1.625990 -1.226800 -5.34580  3.816790
# 6 5121015 5121015 -11.5847 5.41131  0.785815 -0.940348 -6.08543  1.465410
# pc7        pc8      pc9      pc10
# 1  3.109920 -2.6308500  2.39288  0.307537
# 2  1.141060 -1.9821300 -2.70226  2.507750
# 3 -1.548790  1.4789300 -1.20895  1.064900
# 4 -0.534071 -0.6543670 -6.59593 -1.533560
# 5  0.579155 -1.1695800  1.01855 -0.795227
# 6  2.698600  0.0793596 -2.85344  3.514010

#-------------------------------

#-------------------------------

m1 <- match(dat2$iid, dat$V1)

year_of_birth <- as.numeric(as.character(dat[m1,2]))
sex <- as.numeric(as.character(dat[m1,3]))
age <- as.numeric(as.character(dat[m1,4]))
tdi <- as.numeric(as.character(dat[m1,5]))
centre <- as.numeric(as.character(dat[m1,6]))
batch <- as.numeric(as.character(dat[m1,7]))

# sex <- paste("sex", sex, sep = "")
# centre <- paste("centre", centre, sep = "")
# batch <- paste("batch", batch, sep = "")

pc1 <- dat2[,3]
pc2 <- dat2[,4]
pc3 <- dat2[,5]
pc4 <- dat2[,6]
pc5 <- dat2[,7]
pc6 <- dat2[,8]
pc7 <- dat2[,9]
pc8 <- dat2[,10]
pc9 <- dat2[,11]
pc10 <- dat2[,12]

#-------------------------------

#-------------------------------

covar = data.frame(dat2[, c(1,2)], year_of_birth, sex, age, tdi, centre, batch, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)
rownames(covar) = dat2$iid
# > head(covar)
# fid     iid year_of_birth  sex age      tdi      centre  batch
# 5438743 5438743 5438743          1967 sex1  41 -3.66848 centre11011 batch1
# 2890570 2890570 2890570          1962 sex0  46 -2.61299 centre11011 batch1
# 2689958 2689958 2689958          1957 sex0  52 -5.86031 centre11011 batch1
# 5800351 5800351 5800351          1943 sex0  65 -3.97211 centre11011 batch1
# 5644632 5644632 5644632          1952 sex0  56 -1.15201 centre11011 batch1
# 5121015 5121015 5121015          1959 sex1  50 -2.76298 centre11011 batch1
# pc1     pc2       pc3       pc4      pc5       pc6       pc7
# 5438743 -12.1725 5.39163 -1.281030  0.841765 -5.26521 -1.786570  3.109920
# 2890570 -13.0245 6.41514 -0.183365  2.927610 -5.88964  0.940534  1.141060
# 2689958 -11.4712 3.48383 -1.154580  3.083830  7.65160 -0.913399 -1.548790
# 5800351 -12.1327 4.02976 -0.988080  0.750294 -2.36431  0.431658 -0.534071
# 5644632 -12.2171 3.50821 -1.625990 -1.226800 -5.34580  3.816790  0.579155
# 5121015 -11.5847 5.41131  0.785815 -0.940348 -6.08543  1.465410  2.698600
# pc8      pc9      pc10
# 5438743 -2.6308500  2.39288  0.307537
# 2890570 -1.9821300 -2.70226  2.507750
# 2689958  1.4789300 -1.20895  1.064900
# 5800351 -0.6543670 -6.59593 -1.533560
# 5644632 -1.1695800  1.01855 -0.795227
# 5121015  0.0793596 -2.85344  3.514010
# > nrow(covar)
# [1] 488377

covar = subset(covar, complete.cases(covar))
# > nrow(covar)
# [1] 487765

colnames(covar) = c("FID", "IID", "year", "sex", "age", "tdi", "centre", "batch","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")

#-------------------------------

#-------------------------------
# Writing files

writeLines("Writing files ...")

writeLines(paste("Writing file: ", opt$out_file, sep = ""))
sink(opt$out_file)
write.table(covar,quote=F,col.name=T,row.name=F, sep = "\t")
sink()
#-------------------------------
