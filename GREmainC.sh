#===============================
# For 50K samples
# 1) Sampling distributions of model parameters for the genome-transcriptome partitioning model.
#===============================

# Input data
# .fam						fam file
# .txt						data of transcriptomic covariates (imputed data)
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# Run on statgen server
export DATA=/data/alh-admvhp1/GRE/Data
export OUTPUT=/data/alh-admvhp1/GRE/Data

cd $OUTPUT

# # tango
# export DATA=/home/users/phamvv/GRE/Data
# export OUTPUT=/home/users/phamvv/GRE/Data
# 
# cd $OUTPUT

# 1.1) fam file (should remove samples who have withdrawn)
# qced_rdm_005_ukbb3.fam

# 1.2) .txt: data of transcriptomic covariates (imputed data)
# ln -s /mnt/rdfs/WholeGenomeApproach/transcriptome/imputed_expression/Brain_Substantia_nigra_predicted_expression.txt ./Brain_Substantia_nigra_predicted_expression.txt

# 1.1&2) Select 50K samples
fam_file="./qced_rdm_005_ukbb3.fam"
tra_file="./Brain_Substantia_nigra_predicted_expression.txt"
awk '{print $2}' $tra_file > "expression_id.txt"
Rscript ../selectSamples.R --fam_file=${fam_file} --tra_file="expression_id.txt" --number=50000 --outFile1="selectedFam.fam" --outFile2="selectedExpressionId.txt"
grep -Fwf "selectedExpressionId.txt" "./Brain_Substantia_nigra_predicted_expression.txt" > "expression_selected.txt"
# Remove genes with 0 for all samples
### R ###
R

library(dplyr)

setwd("/data/alh-admvhp1/GRE/Data")

# expression file
exp=read.table("expression_selected.txt", header=F)
exp2=exp[,-c(1,2)] # exclude 1st & 2 nd col, ids
ncov=dim(exp2)[2] # no. of covariates (i.e. genes)
# > ncov
# [1] 2041

new_df <- colSums(exp2)
idx <- which(new_df == 0)
# > idx
#  V145  V686  V952 V1009 V1018 V1420 V1433 V1435 V1501 V1530 V1600 V1674 V1688
#   143   684   950  1007  1016  1418  1431  1433  1499  1528  1598  1672  1686
l <- length(idx)
# > l
# [1] 13
exp3 <- exp2[,-idx]
ncov=dim(exp3)[2]
# > ncov
# [1] 2028

# add in id
exp4=cbind(exp$V1, exp$V2, exp3)
# export dat
write.table(exp4, "expression_selected_processed.txt", col.names=F, row.names=F,quote=F)

q() # exit R
### R ###

# 1.3) Genomic ralationship matrix
../Tool/plink --bfile qced_rdm_005_ukbb3 --keep selectedFam.fam --make-grm-gz --out grm_qced_rdm_005_ukbb3_selected
gunzip grm_qced_rdm_005_ukbb3_selected.grm.gz
# remove the third column
awk '{print $1,$2,$4}' grm_qced_rdm_005_ukbb3_selected.grm > grm_qced_rdm_005_ukbb3_selected_processed.grm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.4) (**) Compute kernel matrix K

# i. compute B/sqrt(m)

### R ###
R

setwd("/data/alh-admvhp1/GRE/Data")

# covariate file
cov=read.table("expression_selected_processed.txt", header=F)
# > nrow(cov)
# [1] 50000 # number of patients
cov2=cov[,-c(1,2)] # exclude 1st & 2 nd col, ids
ncov=dim(cov2)[2] # no. of covariates (i.e. genes)
# > ncov
# [1] 2028

# standardization
cov3=data.frame(matrix(nrow=dim(cov2)[1], ncol=dim(cov2)[2]))
for (i in 1:ncov){
   sel=cov2[,i]
   ave=mean(sel, na.rm=T)
   std=sd(sel, na.rm=T)
   # apply standization
   cov3[,i]=(sel-ave)/(std*sqrt(ncov))
  }

# add in id
cov4=cbind(cov$V1, cov$V2, cov3)
# export dat
write.table(cov4, "brain.stdcov", col.names=F, row.names=F,quote=F)

q() # exit R
### R ###

# ii. compute  K=BB'/m
../Tool/mtg2 -p selectedFam.fam -pdmx brain.stdcov -thread 80 -out brain.bmat

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1.5 (***) 
# 2.2. Compute kernel matrix Q

# i. ensure positive definite kernel matrices
# .	for A
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -thread 80 -bend 2

# .	for K
../Tool/mtg2 -p selectedFam.fam -g brain.bmat -thread 80 -bend 2

# ii. Cholesky decomposition
# .	for A
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend -thread 80 -chol 1

# .	for K
../Tool/mtg2 -p selectedFam.fam -g brain.bmat.bend -thread 80 -chol 1

# iii. compute Q
../Tool/mtg2 -p selectedFam.fam -mg brain_selected_grm_bmat.chol -thread 80 -matmat 2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
# C1) Continuous data
#===============================

# Input data
# .fam						fam file
# .txt						data of transcriptomic covariates (imputed data)
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# # Run on statgen server
# export DATA=/data/alh-admvhp1/GRE/Data
# export OUTPUT=/data/alh-admvhp1/GRE/Data
# 
# cd $OUTPUT

# # tango
# export DATA=/home/users/phamvv/GRE/Data
# export OUTPUT=/home/users/phamvv/GRE/Data
# 
# cd $OUTPUT

# gadi
export DATA=/scratch/eu82/vp8928/GRE/Data
export OUTPUT=/scratch/eu82/vp8928/GRE/Data

cd $OUTPUT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C1.1) Simulate phenotype data & fit
# .dat: phenotypic data of the main trait
# cov = 0
cov=0
if [ -e res.greml.out ]
then
    rm res.greml.out
fi
if [ -e res.coregreml.out ]
then
    rm res.coregreml.out
fi
if [ -e res.coregreml.out.converted ]
then
    rm res.coregreml.out.converted
fi
for ((i=1;i<=50;i++))
    do
    
    echo "cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # brain_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> res.coregreml.out
    #++++++++++++++++++++++++
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> res.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> res.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> res.coregreml.out.converted
    #++++++++++++++++++++++++

  done

# cov = 0.2
cov=0.2
if [ -e resp.greml.out ]
then
    rm resp.greml.out
fi
if [ -e resp.coregreml.out ]
then
    rm resp.coregreml.out
fi
if [ -e resp.coregreml.out.converted ]
then
    rm resp.coregreml.out.converted
fi
for ((i=1;i<=50;i++))
    do
    
    echo "cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # brain_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> resp.coregreml.out
    #++++++++++++++++++++++++
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> resp.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> resp.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> resp.coregreml.out.converted
    #++++++++++++++++++++++++
  done

# cov = -0.2
cov=-0.2
if [ -e resn.greml.out ]
then
    rm resn.greml.out
fi
if [ -e resn.coregreml.out ]
then
    rm resn.coregreml.out
fi
if [ -e resn.coregreml.out.converted ]
then
    rm resn.coregreml.out.converted
fi
for ((i=1;i<=50;i++))
    do
    
    echo "cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # brain_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> resn.coregreml.out
    #++++++++++++++++++++++++
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> resn.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> resn.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> resn.coregreml.out.converted
    #++++++++++++++++++++++++
  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# C2) For binary data, correct the transformation, scale
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C2.1) Simulate phenotype data & fit
# .dat: phenotypic data of the main trait
# cov = 0
cov=0
if [ -e bi.res.greml.out ]
then
    rm bi.res.greml.out
fi
if [ -e bi.res.coregreml.out ]
then
    rm bi.res.coregreml.out
fi
if [ -e bi.res.coregreml.out.converted ]
then
    rm bi.res.coregreml.out.converted
fi
if [ -e bi.res.greml.out.transformed ]
then
    rm bi.res.greml.out.transformed
fi
if [ -e bi.res.coregreml.out.transformed ]
then
    rm bi.res.coregreml.out.transformed
fi
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # brain_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    
    #++++++++++++++++++++++++
    # estimate h2 (refer to iii. estimate functions of model parameters)
    # Estimate proporitions of phenotypic variance explained by genetics and transcriptome,
    # correlation between genetic and transcriptomic effects on phenotypes,
    # and standard errors of these estimates
    
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.res.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.res.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out.converted
    #++++++++++++++++++++++++

  done

#-----------------------
# cov = 0, transformation

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.res.greml.out.transformed
  done < bi.res.greml.out
 
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.coregreml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[10]} ${line[11]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.coregreml.out.transformed
    
    # for cov
    ../Tool/mtg2 -trf_h2 ${line[12]} ${line[13]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.res.coregreml.out.transformed
  done < bi.res.coregreml.out.converted
  
#-----------------------

# cov = 0.2
cov=0.2
if [ -e bi.resp.greml.out ]
then
    rm bi.resp.greml.out
fi
if [ -e bi.resp.coregreml.out ]
then
    rm bi.resp.coregreml.out
fi
if [ -e bi.resp.coregreml.out.converted ]
then
    rm bi.resp.coregreml.out.converted
fi
if [ -e bi.resp.greml.out.transformed ]
then
    rm bi.resp.greml.out.transformed
fi
if [ -e bi.resp.coregreml.out.transformed ]
then
    rm bi.resp.coregreml.out.transformed
fi
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # brain_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    
    #++++++++++++++++++++++++
    # estimate h2 (refer to iii. estimate functions of model parameters)
    # Estimate proporitions of phenotypic variance explained by genetics and transcriptome,
    # correlation between genetic and transcriptomic effects on phenotypes,
    # and standard errors of these estimates
    
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resp.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resp.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out.converted
    #++++++++++++++++++++++++

  done
  
#-----------------------
# cov = 0.2, transformation

# Transform file bi.resp.greml.out

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resp.greml.out.transformed
  done < bi.resp.greml.out

# Transform file bi.resp.coregreml.out.converted

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.coregreml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[10]} ${line[11]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.coregreml.out.transformed
    
    # for cov
    ../Tool/mtg2 -trf_h2 ${line[12]} ${line[13]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resp.coregreml.out.transformed
  done < bi.resp.coregreml.out.converted
#-----------------------

# cov = -0.2
cov=-0.2
if [ -e bi.resn.greml.out ]
then
    rm bi.resn.greml.out
fi
if [ -e bi.resn.coregreml.out ]
then
    rm bi.resn.coregreml.out
fi
if [ -e bi.resn.coregreml.out.converted ]
then
    rm bi.resn.coregreml.out.converted
fi
if [ -e bi.resn.greml.out.transformed ]
then
    rm bi.resn.greml.out.transformed
fi
if [ -e bi.resn.coregreml.out.transformed ]
then
    rm bi.resn.coregreml.out.transformed
fi
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # brain.bmat
    # brain_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    
    #++++++++++++++++++++++++
    # estimate h2 (refer to iii. estimate functions of model parameters)
    # Estimate proporitions of phenotypic variance explained by genetics and transcriptome,
    # correlation between genetic and transcriptomic effects on phenotypes,
    # and standard errors of these estimates
    
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resn.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resn.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out.converted
    #++++++++++++++++++++++++

  done
  
#-----------------------
# cov = -0.2, transformation

# Transform file bi.resn.greml.out

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resn.greml.out.transformed
  done < bi.resn.greml.out

# Transform file bi.resn.coregreml.out.converted

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.coregreml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[10]} ${line[11]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.coregreml.out.transformed
    
    # for cov
    ../Tool/mtg2 -trf_h2 ${line[12]} ${line[13]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resn.coregreml.out.transformed
  done < bi.resn.coregreml.out.converted
#-----------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
# F
#===============================

#===============================
# For 50K samples
# 1) Sampling distributions of model parameters for the genome-transcriptome partitioning model.
#===============================

# Input data
# .fam						fam file
# .txt						data of transcriptomic covariates (imputed data)
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# Run on statgen server
export DATA=/data/alh-admvhp1/GRE/gEnv
export OUTPUT=/data/alh-admvhp1/GRE/gEnv

cd $OUTPUT

# # tango
# export DATA=/home/users/phamvv/GRE/Data
# export OUTPUT=/home/users/phamvv/GRE/Data
# 
# cd $OUTPUT

# 1.1) fam file
# qced_rdm_005_ukbb3.fam

# 1.2) .txt: environment data
# extract environment data
other_pheno_file="ukb6247.csv" # year of birth, sex, age, etc.
covariate_file="my_sample_qcinfo.txt" # PCs
out_file="env_data.txt"
Rscript ../getEnvData.R --other_pheno_file=${other_pheno_file} --covariate_file=${covariate_file} --out_file=${out_file}

# 1.1&2) Select 50K samples
fam_file="./qced_rdm_005_ukbb3.fam"
tra_file="./env_data.txt"
awk '{print $2}' $tra_file > "env_id.txt"
Rscript ../selectSamples.R --fam_file=${fam_file} --tra_file="env_id.txt" --number=50000 --outFile1="selectedFam.fam" --outFile2="selectedEnvId.txt"
grep -Fwf "selectedEnvId.txt" "./env_data.txt" > "env_selected.txt"

# 1.3) Genomic ralationship matrix
../Tool/plink --bfile qced_rdm_005_ukbb3 --keep selectedFam.fam --make-grm-gz --out grm_qced_rdm_005_ukbb3_selected
gunzip grm_qced_rdm_005_ukbb3_selected.grm.gz
# remove the third column
awk '{print $1,$2,$4}' grm_qced_rdm_005_ukbb3_selected.grm > grm_qced_rdm_005_ukbb3_selected_processed.grm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.4) (**) Compute kernel matrix K

# i. compute B/sqrt(m)

### R ###
R

setwd("/data/alh-admvhp1/GRE/gEnv")

# covariate file
cov=read.table("env_selected.txt", header=F)
# > nrow(cov)
# [1] 50000 # number of patients
cov2=cov[,-c(1,2)] # exclude 1st & 2 nd col, ids
ncov=dim(cov2)[2] # no. of covariates (i.e. genes)
# > ncov
# [1] 16

# standardization
cov3=data.frame(matrix(nrow=dim(cov2)[1], ncol=dim(cov2)[2]))
for (i in 1:ncov){
   sel=cov2[,i]
   ave=mean(sel, na.rm=T)
   std=sd(sel, na.rm=T)
   # apply standization
   cov3[,i]=(sel-ave)/(std*sqrt(ncov))
  }

# add in id
cov4=cbind(cov$V1, cov$V2, cov3)
cov5 <- cov4[order(cov4[,1], decreasing=FALSE),]
# export dat
write.table(cov5, "env.stdcov", col.names=F, row.names=F,quote=F)

q() # exit R
### R ###

# ii. compute  K=BB'/m
../Tool/mtg2 -p selectedFam.fam -pdmx env.stdcov -thread 80 -out env.bmat > out.log

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1.5 (***) 
# 2.2. Compute kernel matrix Q

# i. ensure positive definite kernel matrices
# .	for A
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -thread 80 -bend 2 > out22i.log

# .	for K
../Tool/mtg2 -p selectedFam.fam -g env.bmat -thread 80 -bend 2 > out.log

# ii. Cholesky decomposition
# .	for A
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend -thread 80 -chol 1 > out22iiA.log

# .	for K
../Tool/mtg2 -p selectedFam.fam -g env.bmat.bend -thread 80 -chol 1 > out22iiK.log

# iii. compute Q
../Tool/mtg2 -p selectedFam.fam -mg env_selected_grm_bmat.chol -thread 80 -matmat 2 > out22iii.log

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
# F1) Continuous data
#===============================

# Input data
# .fam						fam file
# .txt						data of transcriptomic covariates (imputed data)
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# Run on statgen server
export DATA=/data/alh-admvhp1/GRE/gEnv
export OUTPUT=/data/alh-admvhp1/GRE/gEnv

cd $OUTPUT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# F1-1) Standing height
# .dat: phenotypic data of the main trait

#sample_pheno.dat: phenotype

# [alh-admvhp1@hscpl-statgen Data]$ wc -l selectedFam.fam
# 50000 selectedFam.fam

# R ------------------
# for education
R

dat=read.table("edu_corrected")
head(dat)
# > head(dat)
#        id edu
# 1 2080762  10
# 2 3165676  20
# 3 1605581  10
# 4 5336523  20
# 5 2211334  15
# 6 5391704  19

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("edu_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=T)
sink()

q()
# R ------------------

# extract phenotype
#pre="edu"
# pre="sh" # ***update this*** Standing height
#pre="sith" # sitting height
# pre="heel" # E6-1.3) Get phenotype - Heel bone mineral density 3148 3148-0.0 680
# pre="heelauto" # E6-1.3*) Get phenotype - Heel bone mineral density (BMD) T-score, automated 78 78-0.0 (22)
# pre="wei" # E6-1.4) Get phenotype - Weight 21002 21002-0.0	21002-1.0	21002-2.0 (2348 2349 2350)
# pre="fluid" # E6-1.5) Get phenotype - Fluid intelligence 20016 20016-0.0	20016-1.0	20016-2.0 (2306 2307 2308)
# pre="bmi" # E6-1.6) Get phenotype - BMI 21001 21001-0.0	21001-1.0 (2346 2347)
pre="hip" # E6-1.7) Get phenotype - Hip circumference 49 49-0.0	49-1.0	49-2.0 (8 9 10)
# pre="wai" # E6-1.8) Get phenotype - Waist circumference 48 48-0.0	48-1.0	48-2.0 (5 6 7)
#pre="blood" # E6-1.9) Get phenotype - Diastolic blood pressure 4079 4079-0.0	4079-0.1	4079-1.0	4079-1.1	4079-2.0	4079-2.1 (774 775 776 777 778 779)
#pre="10body" # Body fat percentage (23099) ? (ukb29455) 23099-0.0	23099-1.0 (1437 1438)
#pre="14blood" # E6-1.14) # 14 Systolic blood pressure, automated reading (4080) 4080-0.0 -> 4080-2.1 (780 -> 785)
#pre="16age" # E6-1.16) # 16 Age when periods started (menarche) (2714) 2714-0.0	2714-1.0	2714-2.0 (618 619 620) , useSex=F
#pre="17age" # E6-1.17) # 17 Age at first live birth (2754) (2754-0.0	2754-1.0	2754-2.0) (627 628 629)  , useSex=F
#pre="20score" # E6-1.20) # 20 Neuroticism score (20127) ??? 20127-0.0 (2339)
# pre="21rate" # E6-1.21) # 21 Pulse rate (4194) ? (ukb43545) 4194-0.0	4194-1.0	4194-2.0	4194-3.0 (14 15 16 17)
pheno_file="${pre}_trait.txt"
selected_id="selectedFam.fam"
out_file="${pre}_random_50K.dat"
Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --selected_id=${selected_id} --out_file=${out_file} --pre=${pre}
# edu_random_50K.dat

# model fitting
# o	inside env_greml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# env.bmat
# .	fit GREML:
# -d sample_pheno.dat: phenotype
../Tool/mtg2 -p selectedFam.fam -mg env_greml.matlist -d ${pre}_random_50K.dat -mod 1 -thread 20 -out ${pre}_greml.out > ${pre}_greml.log
 
# .	fit CORE GREML:
# o	inside env_coregreml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# env.bmat
# env_selected_grm_bmat.chol.matmat2
../Tool/mtg2 -p selectedFam.fam -mg env_coregreml.matlist -d ${pre}_random_50K.dat -mod 1 -thread 20 -out ${pre}_coregreml.out > ${pre}_coregreml.log

#++++++++++++++++++++++++
# o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
# These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
grep -vwE '(LKH|h2)' ${pre}_coregreml.out > ${pre}_coregreml.out2

# o	Prepare a do file for MTG2 to construct functions of model parameters.
Rscript ../createTemplate.R --inFile="${pre}_coregreml.out2" --outFile="${pre}_coregreml.do"

# .	Estimation:
../Tool/mtg2 -delta2 ${pre}_coregreml.do > ${pre}_coregreml.out3

#++++++++++++++++++++++++
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# F1-2) For binary - traits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.12) (Categorical) Get phenotype - # Vascular/heart problems diagnosed by doctor: High blood pressure (6150_4) 6150-0.0	6150-0.1	6150-0.2	6150-0.3	6150-1.0	6150-1.1	6150-1.2	6150-1.3	6150-2.0	6150-2.1	6150-2.2	6150-2.3 (1712 -> 1723)

# R ------------------
# get only related columns
R

cols <- c("6150-0.0", "6150-0.1", "6150-0.2", "6150-0.3", "6150-1.0", "6150-1.1", "6150-1.2", "6150-1.3", "6150-2.0", "6150-2.1", "6150-2.2", "6150-2.3")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
#  [1] 1712 1713 1714 1715 1716 1717 1718 1719 1720 1721 1722 1723
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4       V5       V6       V7       V8
# 1     eid 6150-0.0 6150-0.1 6150-0.2 6150-0.3 6150-1.0 6150-1.1 6150-1.2
# 2 2080762        4
# 3 3165676       -7
# 4 1605581        4
# 5 5336523       -7
# 6 2211334        3

# Coding	Meaning
# 1	Heart attack
# 2	Angina
# 3	Stroke
# 4	High blood pressure
# -7	None of the above
# -3	Prefer not to answer
dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-3,NA,dat[,2])
dat <- dat[-1,]
dat[,2] <- ifelse(dat[,2]==1,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==2,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==3,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==-7,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==4,1,dat[,2])

writeLines("Getting traits ...")
sink("12heart_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.18) (Categorical) # 18 Alcohol intake frequency (1558) ??? 1558-0.0	1558-1.0	1558-2.0 (441 442 443)

# Coding	Meaning
# 1	Daily or almost daily
# 2	Three or four times a week
# 3	Once or twice a week
# 4	One to three times a month
# 5	Special occasions only
# 6	Never
# -3	Prefer not to answer

# R ------------------
# get only related columns
R

cols <- c("1558-0.0", "1558-1.0", "1558-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 441 442 443

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4
# 1     eid 1558-0.0 1558-1.0 1558-2.0
# 2 2080762        2
# 3 3165676        2
# 4 1605581        2
# 5 5336523        2
# 6 2211334        1

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-3,NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==6,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==0,0,1)

writeLines("Getting traits ...")
sink("18alcohol_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.22) (Categorical) # 22 Happiness (4526) ??? (4526-0.0	4526-1.0	4526-2.0) (1418 1419 1420)

# Coding	Meaning
# 1	Extremely happy
# 2	Very happy
# 3	Moderately happy
# 4	Moderately unhappy
# 5	Very unhappy
# 6	Extremely unhappy
# -1	Do not know
# -3	Prefer not to answer

# R ------------------
# get only related columns
R

cols <- c("4526-0.0", "4526-1.0", "4526-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 1418 1419 1420

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4
# 1     eid 4526-0.0 4526-1.0 4526-2.0
# 2 2080762
# 3 3165676        2
# 4 1605581        2
# 5 5336523
# 6 2211334

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==4,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==5,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==6,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==-1,NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-3,NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==2,1,dat[,2])
dat[,2] <- ifelse(dat[,2]==3,1,dat[,2])

writeLines("Getting traits ...")
sink("21happiness_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# F1-2) For binary
# .dat: phenotypic data of the main trait

#sample_pheno.dat: phenotype

# [alh-admvhp1@hscpl-statgen Data]$ wc -l selectedFam.fam
# 50000 selectedFam.fam

# extract phenotype
# pre="12heart" # E6-1.12) (Categorical) Get phenotype - # Vascular/heart problems diagnosed by doctor: High blood pressure (6150_4) 6150-0.0	6150-0.1	6150-0.2	6150-0.3	6150-1.0	6150-1.1	6150-1.2	6150-1.3	6150-2.0	6150-2.1	6150-2.2	6150-2.3 (1712 -> 1723)
# pre="13non" # E6-1.13) (Categorical) 13 Non-cancer illness code, self-reported: hypertension (20002_1065) 87 from 20002-0.0 to 20002-2.28 (1901 -> 1987)
# pre="15ever" # E6-1.15) (Categorical) # 15 Ever smoked (20160) 20160-0.0 (2340)
# pre="18alcohol" # E6-1.18) (Categorical) # 18 Alcohol intake frequency (1558) ??? 1558-0.0	1558-1.0	1558-2.0 (441 442 443)
#pre="19mood" # E6-1.19) (Categorical) # 19 Mood swings (1920) (1920-0.0	1920-1.0	1920-2.0) 495 496 497
pre="21happiness" # E6-1.22) (Categorical) # 22 Happiness (4526) ??? (4526-0.0	4526-1.0	4526-2.0) (1418 1419 1420)
pheno_file="${pre}_trait.txt"
# ln -s /data/alh-admhl/sfile_ukbb_general/info/ukb6247.csv /data/alh-admvhp1/GRERes/Data/ukb6247.csv
# ln -s /mnt/rdfs/WholeGenomeApproach/ukbb3/my_sample_qcinfo.txt /data/alh-admvhp1/GRERes/Data/my_sample_qcinfo.txt
selected_id="selectedFam.fam"
out_file="${pre}_random_50K_bi.dat"
Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --selected_id=${selected_id} --out_file=${out_file}  --pre=${pre}
# sh_random_50K_adjusted.dat

# model fitting
# o	inside env_greml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# env.bmat
# .	fit GREML:
# -d sample_pheno.dat: phenotype
../Tool/mtg2 -p selectedFam.fam -mg env_greml.matlist -d ${pre}_random_50K_bi.dat -mod 1 -thread 20 -out ${pre}_greml.out > ${pre}_greml.log
awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' ${pre}_greml.out > ${pre}.bi.res.greml.out
 
# .	fit CORE GREML:
# o	inside brain_coregreml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# brain.bmat
# brain_selected_grm_bmat.chol.matmat2
../Tool/mtg2 -p selectedFam.fam -mg env_coregreml.matlist -d ${pre}_random_50K_bi.dat -mod 1 -thread 20 -out ${pre}_coregreml.out > ${pre}_coregreml.log
awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' ${pre}_coregreml.out > ${pre}.bi.res.coregreml.out

#++++++++++++++++++++++++
# o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
# These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
grep -vwE '(LKH|h2)' ${pre}_coregreml.out > ${pre}_coregreml.out2

# o	Prepare a do file for MTG2 to construct functions of model parameters.
Rscript ../createTemplate.R --inFile="${pre}_coregreml.out2" --outFile="${pre}_coregreml.do"

# .	Estimation:
../Tool/mtg2 -delta2 ${pre}_coregreml.do > ${pre}_coregreml.out3

# Write h2 and se to file
awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' ${pre}_coregreml.out > ${pre}.bi.res.coregreml.out.converted
awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' ${pre}_coregreml.out3 >> ${pre}.bi.res.coregreml.out.converted
awk '$1=="LKH" {printf $2"\n"}' ${pre}_coregreml.out >> ${pre}.bi.res.coregreml.out.converted
#++++++++++++++++++++++++

#-----------------------
# transformation

#------------
# R
R

dat <- read.table("12heart_random_50K_bi.dat")
nrow(dat[which(dat[,3]==1),])/50000

q()
#------------

#k : population prevalence
#k=0.24 # pre="12heart"
#k=0.22718 # pre="13non"
#k=0.5452 # pre="15ever"
#k=0.93516 # pre="18alcohol"
#k=0.43608 # pre="19mood"
k=0.31674 # pre="21happiness"
#p : proportion of cases in training sample
#p=0.24 # pre="12heart"
#p=0.22718 # pre="13non"
#p=0.5452 # pre="15ever"
#p=0.93516 # pre="18alcohol"
#p=0.43608 # pre="19mood"
p=0.31674 # pre="21happiness"

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} ${k} ${p} obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out > ${pre}.bi.res.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} ${k} ${p} obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> ${pre}.bi.res.greml.out.transformed
  done < ${pre}.bi.res.greml.out
 
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} ${k} ${p} obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out > ${pre}.bi.res.coregreml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[10]} ${line[11]} ${k} ${p} obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> ${pre}.bi.res.coregreml.out.transformed
    
    # for cov
    ../Tool/mtg2 -trf_h2 ${line[12]} ${line[13]} ${k} ${p} obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> ${pre}.bi.res.coregreml.out.transformed
  done < ${pre}.bi.res.coregreml.out.converted
  
#-----------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


