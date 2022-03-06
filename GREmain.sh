#===============================
# 1) Sampling distributions of model parameters for the genome-transcriptome partitioning model.
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

# tango
export DATA=/home/users/phamvv/GRE/Data
export OUTPUT=/home/users/phamvv/GRE/Data

cd $OUTPUT

# 1.1) fam file
# qced_rdm_005_ukbb3.fam

# 1.2) .txt: data of transcriptomic covariates (imputed data)
# ln -s /mnt/rdfs/WholeGenomeApproach/transcriptome/imputed_expression/Brain_Substantia_nigra_predicted_expression.txt ./Brain_Substantia_nigra_predicted_expression.txt

# 1.1&2) Select 10K samples
fam_file="./qced_rdm_005_ukbb3.fam"
tra_file="./Brain_Substantia_nigra_predicted_expression.txt"
awk '{print $2}' $tra_file > "expression_id.txt"
Rscript ../select10KSamples.R --fam_file=${fam_file} --tra_file="expression_id.txt"
grep -Fwf "expression_id_10K.txt" "./Brain_Substantia_nigra_predicted_expression.txt" > "expression_10K.txt"
# Remove genes with 0 for all samples
### R ###
R

library(dplyr)

setwd("/data/alh-admvhp1/GRE/Data")

# expression file
exp=read.table("expression_10K.txt", header=F)
# > exp[1:5,1:6]
#        V1      V2          V3        V4         V5          V6
# 1 1000584 1000584 -0.10675104 0.2017973 -0.5075195 -0.06146026
# 2 1001316 1001316 -0.07897366 0.2909065 -0.2537598 -0.45022481
# 3 1001492 1001492  0.93160198 0.5210385 -0.2537598  0.19850169
# 4 1002199 1002199  0.63067778 0.3965852 -0.5075195 -0.43226838
# 5 1002524 1002524  0.62830556 0.8557313 -0.5075195 -0.15294678
# > nrow(exp)
# [1] 10000 # number of patients
exp2=exp[,-c(1,2)] # exclude 1st & 2 nd col, ids
# > exp2[1:5, 1:6]
#            V3        V4         V5          V6        V7         V8
# 1 -0.10675104 0.2017973 -0.5075195 -0.06146026 0.4615964 -0.1125180
# 2 -0.07897366 0.2909065 -0.2537598 -0.45022481 0.6444207  0.1329098
# 3  0.93160198 0.5210385 -0.2537598  0.19850169 0.7515206  0.6400939
# 4  0.63067778 0.3965852 -0.5075195 -0.43226838 0.2434703  0.4377443
# 5  0.62830556 0.8557313 -0.5075195 -0.15294678 0.5063792  0.4107797
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
# > exp4[1:5, 1:6]
#    exp$V1  exp$V2          V3        V4         V5          V6
# 1 1000584 1000584 -0.10675104 0.2017973 -0.5075195 -0.06146026
# 2 1001316 1001316 -0.07897366 0.2909065 -0.2537598 -0.45022481
# 3 1001492 1001492  0.93160198 0.5210385 -0.2537598  0.19850169
# 4 1002199 1002199  0.63067778 0.3965852 -0.5075195 -0.43226838
# 5 1002524 1002524  0.62830556 0.8557313 -0.5075195 -0.15294678
# export dat
write.table(exp4, "expression_10K_processed.txt", col.names=F, row.names=F,quote=F)

q() # exit R
### R ###

# 1.3) Genomic ralationship matrix
../Tool/plink --bfile qced_rdm_005_ukbb3 --keep fam_10K.fam --make-grm-gz --out grm_qced_rdm_005_ukbb3_10K
gunzip grm_qced_rdm_005_ukbb3_10K.grm.gz
# remove the third column
awk '{print $1,$2,$4}' grm_qced_rdm_005_ukbb3_10K.grm > grm_qced_rdm_005_ukbb3_10K_processed.grm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.4) (**) Compute kernel matrix K

# i. compute B/sqrt(m)

### R ###
R

setwd("/data/alh-admvhp1/GRE/Data")

# covariate file
cov=read.table("expression_10K_processed.txt", header=F)
# > cov[1:5,1:6]
#        V1      V2          V3        V4         V5          V6
# 1 1000584 1000584 -0.10675104 0.2017973 -0.5075195 -0.06146026
# 2 1001316 1001316 -0.07897366 0.2909065 -0.2537598 -0.45022481
# 3 1001492 1001492  0.93160198 0.5210385 -0.2537598  0.19850169
# 4 1002199 1002199  0.63067778 0.3965852 -0.5075195 -0.43226838
# 5 1002524 1002524  0.62830556 0.8557313 -0.5075195 -0.15294678
# > nrow(cov)
# [1] 10000 # number of patients
cov2=cov[,-c(1,2)] # exclude 1st & 2 nd col, ids
# > cov2[1:5, 1:6]
#            V3        V4         V5          V6        V7         V8
# 1 -0.10675104 0.2017973 -0.5075195 -0.06146026 0.4615964 -0.1125180
# 2 -0.07897366 0.2909065 -0.2537598 -0.45022481 0.6444207  0.1329098
# 3  0.93160198 0.5210385 -0.2537598  0.19850169 0.7515206  0.6400939
# 4  0.63067778 0.3965852 -0.5075195 -0.43226838 0.2434703  0.4377443
# 5  0.62830556 0.8557313 -0.5075195 -0.15294678 0.5063792  0.4107797
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
# > cov4[1:5, 1:6]
#    cov$V1  cov$V2          X1           X2           X3           X4
# 1 1000584 1000584 -0.07639415 -0.009382701 -0.025584088 -0.008018851
# 2 1001316 1001316 -0.07366890 -0.004020889  0.006587869 -0.034241803
# 3 1001492 1001492  0.02547924  0.009826429  0.006587869  0.009516106
# 4 1002199 1002199 -0.00404460  0.002337926 -0.025584088 -0.033030607
# 5 1002524 1002524 -0.00427734  0.029965303 -0.025584088 -0.014189802
# export dat
write.table(cov4, "brain.stdcov", col.names=F, row.names=F,quote=F)

q() # exit R
### R ###

# ii. compute  K=BB'/m
../Tool/mtg2 -p fam_10K.fam -pdmx brain.stdcov -thread 80 -out brain.bmat

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1.5 (***) 
# 2.2. Compute kernel matrix Q

# i. ensure positive definite kernel matrices
# .	for A
../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm -thread 80 -bend 2

# .	for K
../Tool/mtg2 -p fam_10K.fam -g brain.bmat -thread 80 -bend 2

# ii. Cholesky decomposition
# .	for A
../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend -thread 80 -chol 1 

# .	for K
../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend -thread 80 -chol 1 

# iii. compute Q
../Tool/mtg2 -p fam_10K.fam -mg brain_grm_bmat.chol -thread 80 -matmat 2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.6) Simulate phenotype data & fit
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
for ((i=1;i<=100;i++))
    do
    
    echo "cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> res.coregreml.out

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
for ((i=1;i<=100;i++))
    do
    
    echo "cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> resp.coregreml.out

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
for ((i=1;i<=100;i++))
    do
    
    echo "cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> resn.coregreml.out

  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
# 2) For binary data
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1) Simulate phenotype data & fit
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
for ((i=1;i<=100;i++))
    do
    
    echo "binary - cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out

  done

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
for ((i=1;i<=100;i++))
    do
    
    echo "binary - cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out

  done

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
for ((i=1;i<=100;i++))
    do
    
    echo "binary - cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out

  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.2) Transform

# cov = 0

# Transform file bi.res.greml.out
# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.greml.out
# 0.0658 0.0030   0.0111 0.0030   0.0116 0.0009   0.1260 0.0340     0.1306 0.0097       7310.1967
# 0.0721 0.0031   0.0065 0.0030   0.0116 0.0009   0.0722 0.0334     0.1290 0.0096       7197.3073
# 0.0672 0.0032   0.0117 0.0032   0.0115 0.0009   0.1295 0.0354     0.1273 0.0097       7194.1964
# 0.0657 0.0031   0.0105 0.0031   0.0114 0.0009   0.1196 0.0348     0.1303 0.0097       7361.4565
# 0.0667 0.0031   0.0107 0.0031   0.0132 0.0010   0.1182 0.0344     0.1458 0.0100       7227.4160
# 0.0661 0.0030   0.0078 0.0030   0.0125 0.0009   0.0901 0.0344     0.1448 0.0099       7462.7374
# 0.0645 0.0032   0.0140 0.0032   0.0115 0.0009   0.1553 0.0351     0.1281 0.0096       7228.4475
# 0.0692 0.0032   0.0085 0.0031   0.0138 0.0010   0.0926 0.0339     0.1510 0.0100       7185.9310
# 0.0648 0.0031   0.0116 0.0031   0.0122 0.0010   0.1310 0.0344     0.1382 0.0099       7322.3072
# 0.0678 0.0032   0.0101 0.0031   0.0123 0.0010   0.1115 0.0348     0.1369 0.0098       7229.7157

if [ -e bi.res.greml.out.transformed ]
then
    rm bi.res.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.res.greml.out.transformed
  done < bi.res.greml.out

# Transform file bi.res.coregreml.out

# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.coregreml.out
# 0.0657 0.0030   0.0115 0.0030   0.0122 0.0010   -0.0017 0.0012   0.1313 0.0346     0.1390 0.0115     -0.0199 0.0140       7311.2162
# 0.0721 0.0031   0.0067 0.0030   0.0120 0.0010   -0.0009 0.0012   0.0747 0.0338     0.1331 0.0112     -0.0104 0.0134       7197.6217
# 0.0672 0.0032   0.0114 0.0032   0.0110 0.0010   0.0015 0.0012   0.1253 0.0353     0.1208 0.0109     0.0162 0.0131       7194.9544
# 0.0654 0.0031   0.0112 0.0031   0.0123 0.0010   -0.0025 0.0012   0.1301 0.0356     0.1424 0.0118     -0.0284 0.0144       7363.3587
# 0.0667 0.0031   0.0108 0.0031   0.0133 0.0011   -0.0001 0.0013   0.1187 0.0347     0.1465 0.0116     -0.0016 0.0142       7227.4225
# 0.0661 0.0030   0.0076 0.0030   0.0122 0.0010   0.0009 0.0012   0.0871 0.0344     0.1403 0.0113     0.0105 0.0134       7463.0482
# 0.0645 0.0032   0.0145 0.0032   0.0126 0.0011   -0.0032 0.0012   0.1643 0.0360     0.1422 0.0116     -0.0357 0.0142       7231.7636
# 0.0694 0.0032   0.0079 0.0031   0.0131 0.0011   0.0021 0.0012   0.0859 0.0337     0.1417 0.0111     0.0225 0.0131       7187.3881
# 0.0647 0.0031   0.0117 0.0031   0.0123 0.0010   -0.0002 0.0012   0.1323 0.0348     0.1393 0.0114     -0.0027 0.0134       7322.3283
# 0.0678 0.0032   0.0101 0.0032   0.0124 0.0011   -0.0002 0.0012   0.1121 0.0350     0.1378 0.0114     -0.0023 0.0138       7229.7291

if [ -e bi.res.coregreml.out.transformed ]
then
    rm bi.res.coregreml.out.transformed
fi
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
  done < bi.res.coregreml.out
  
# cov = 0.2

# Transform file bi.resp.greml.out
if [ -e bi.resp.greml.out.transformed ]
then
    rm bi.resp.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resp.greml.out.transformed
  done < bi.resp.greml.out

# Transform file bi.resp.coregreml.out
if [ -e bi.resp.coregreml.out.transformed ]
then
    rm bi.resp.coregreml.out.transformed
fi
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
  done < bi.resp.coregreml.out

# cov = -0.2

# Transform file bi.resn.greml.out
if [ -e bi.resn.greml.out.transformed ]
then
    rm bi.resn.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resn.greml.out.transformed
  done < bi.resn.greml.out

# Transform file bi.resn.coregreml.out
if [ -e bi.resn.coregreml.out.transformed ]
then
    rm bi.resn.coregreml.out.transformed
fi
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
  done < bi.resn.coregreml.out

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# 3) Check the results
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.0) check the results

#******************************
# likelihood ratio test (LRT)

# run local
### R ###
#R

# setwd("C:/Users/phamvv/MyDoc/31SNP/GRE/Data")
setwd("C:/Users/phamvv/MyDoc/GRE/Data")

# continuous data, cov = 0

res_greml=read.table("res.greml.out", header=F)
# > head(res_greml)
#       V1     V2     V3     V4     V5     V6     V7     V8     V9    V10       V11
# 1 0.2526 0.0249 0.3472 0.0278 0.4030 0.0179 0.3462 0.0273 0.4019 0.0118 -3647.753
# 2 0.2190 0.0253 0.3901 0.0286 0.4008 0.0179 0.3863 0.0277 0.3968 0.0119 -3682.751
# 3 0.2264 0.0251 0.3722 0.0283 0.3872 0.0173 0.3776 0.0281 0.3928 0.0119 -3591.361
# 4 0.2171 0.0249 0.3855 0.0282 0.3882 0.0174 0.3890 0.0279 0.3918 0.0119 -3613.823
# 5 0.2308 0.0253 0.3684 0.0284 0.4145 0.0182 0.3634 0.0276 0.4089 0.0118 -3647.917
# 6 0.2707 0.0248 0.3179 0.0274 0.3961 0.0174 0.3228 0.0274 0.4023 0.0117 -3570.107

res_coregreml=read.table("res.coregreml.out", header=F)
# > head(res_coregreml)
#       V1     V2     V3     V4     V5     V6      V7     V8     V9    V10    V11    V12     V13    V14       V15
# 1 0.2528 0.0249 0.3484 0.0279 0.4060 0.0189 -0.0078 0.0155 0.3486 0.0278 0.4063 0.0148 -0.0078 0.0155 -3647.628
# 2 0.2208 0.0253 0.3929 0.0287 0.4116 0.0192 -0.0278 0.0155 0.3939 0.0285 0.4126 0.0150 -0.0279 0.0157 -3681.108
# 3 0.2264 0.0251 0.3674 0.0283 0.3769 0.0179  0.0274 0.0152 0.3681 0.0281 0.3776 0.0142  0.0274 0.0150 -3589.799
# 4 0.2171 0.0249 0.3862 0.0283 0.3899 0.0183 -0.0046 0.0152 0.3906 0.0284 0.3944 0.0146 -0.0046 0.0154 -3613.778
# 5 0.2307 0.0253 0.3674 0.0285 0.4119 0.0191  0.0069 0.0157 0.3613 0.0279 0.4051 0.0146  0.0068 0.0153 -3647.823
# 6 0.2705 0.0248 0.3163 0.0275 0.3925 0.0182  0.0099 0.0149 0.3197 0.0276 0.3968 0.0143  0.0100 0.0150 -3569.885

n <- 100
type1_rate <- c()
for (i in 1:n) {
  chi = -2*(res_greml[i,11] - res_coregreml[i,15])
  p = pchisq(chi,1,lower.tail=F)
  type1_rate <- c(type1_rate, p)
}
idx <- which(type1_rate < 0.05)
# percentage of type 1 error rate
per <- length(idx)/length(type1_rate)*100
# > per
# [1] 5

# binary data, cov = 0

res_greml=read.table("bi.res.greml.out", header=F)
# > head(res_greml)
#       V1     V2     V3     V4     V5    V6     V7     V8     V9    V10      V11
# 1 0.0658 0.0030 0.0111 0.0030 0.0116 9e-04 0.1260 0.0340 0.1306 0.0097 7310.197
# 2 0.0721 0.0031 0.0065 0.0030 0.0116 9e-04 0.0722 0.0334 0.1290 0.0096 7197.307
# 3 0.0672 0.0032 0.0117 0.0032 0.0115 9e-04 0.1295 0.0354 0.1273 0.0097 7194.196
# 4 0.0657 0.0031 0.0105 0.0031 0.0114 9e-04 0.1196 0.0348 0.1303 0.0097 7361.457
# 5 0.0667 0.0031 0.0107 0.0031 0.0132 1e-03 0.1182 0.0344 0.1458 0.0100 7227.416
# 6 0.0661 0.0030 0.0078 0.0030 0.0125 9e-04 0.0901 0.0344 0.1448 0.0099 7462.737

res_coregreml=read.table("bi.res.coregreml.out", header=F)
# > head(res_coregreml)
#       V1     V2     V3     V4     V5     V6      V7     V8     V9    V10    V11    V12     V13    V14      V15
# 1 0.0657 0.0030 0.0115 0.0030 0.0122 0.0010 -0.0017 0.0012 0.1313 0.0346 0.1390 0.0115 -0.0199 0.0140 7311.216
# 2 0.0721 0.0031 0.0067 0.0030 0.0120 0.0010 -0.0009 0.0012 0.0747 0.0338 0.1331 0.0112 -0.0104 0.0134 7197.622
# 3 0.0672 0.0032 0.0114 0.0032 0.0110 0.0010  0.0015 0.0012 0.1253 0.0353 0.1208 0.0109  0.0162 0.0131 7194.954
# 4 0.0654 0.0031 0.0112 0.0031 0.0123 0.0010 -0.0025 0.0012 0.1301 0.0356 0.1424 0.0118 -0.0284 0.0144 7363.359
# 5 0.0667 0.0031 0.0108 0.0031 0.0133 0.0011 -0.0001 0.0013 0.1187 0.0347 0.1465 0.0116 -0.0016 0.0142 7227.422
# 6 0.0661 0.0030 0.0076 0.0030 0.0122 0.0010  0.0009 0.0012 0.0871 0.0344 0.1403 0.0113  0.0105 0.0134 7463.048

n <- 100
type1_rate <- c()
for (i in 1:n) {
  chi = -2*(res_greml[i,11] - res_coregreml[i,15])
  p = pchisq(chi,1,lower.tail=F)
  type1_rate <- c(type1_rate, p)
}
idx <- which(type1_rate < 0.05)
# percentage of type 1 error rate
per <- length(idx)/length(type1_rate)*100
# > per
# [1] 16

# q() # exit R
### R ###
#******************************

#******************************
# mean, sd, number of cases of pheno

# cov = -0.2
cov=-0.2
i=1
echo "binary - cov = -0.2 - i = " $i

#i=1
Rscript ../random.R --cov=${cov} --seed=${i}

../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
# Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
# Out: brain.bmat.bend.chol.matvec

R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
R CMD BATCH --no-save ../combine_ge_cor.r     #combine

awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
# [alh-admvhp1@hscpl-statgen Data]$ head sample_pheno.dat
# 1000584 1000584 0.399299052480251
# 1001316 1001316 0.592648203065912
# 1001492 1001492 0.916271042033746
# 1002199 1002199 1.60385095053335
# 1002524 1002524 -0.151439901729948
# 1002548 1002548 0.408030673567955
# 1004952 1004952 0.883277887280014
# 1005476 1005476 -0.647416350209848
# 1006220 1006220 -0.685912998195519
# 1006545 1006545 -1.17972574499003

### R ###
R

dat <- read.table("sample_pheno.dat", header=F)
# > head(dat)
#        V1      V2         V3
# 1 1000584 1000584  0.3992991
# 2 1001316 1001316  0.5926482
# 3 1001492 1001492  0.9162710
# 4 1002199 1002199  1.6038510
# 5 1002524 1002524 -0.1514399
# 6 1002548 1002548  0.4080307
mean(dat[,3])
sd(dat[,3])
# > mean(dat[,3])
# [1] 3.442258e-17
# > sd(dat[,3])
# [1] 0.9426564

q() # exit R
### R ###

# Convert the phenotype data to binary
Rscript ../toBinary.R --phe="sample_pheno.dat"

### R ###
R

dat <- read.table("sample_pheno.dat", header=F)
# > head(dat)
#        V1      V2 V3
# 1 1000584 1000584  0
# 2 1001316 1001316  0
# 3 1001492 1001492  0
# 4 1002199 1002199  1
# 5 1002524 1002524  0
# 6 1002548 1002548  0
# number of cases
ind <- which(dat[,3] == 1)
length(ind)
# > length(ind)
# [1] 892

q() # exit R
### R ###

#******************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# 4) For binary data, correct the transformation (compute heritability before the transformation)
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.1) Simulate phenotype data & fit
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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.res.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.res.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    #awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resp.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resp.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resn.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resn.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    #++++++++++++++++++++++++

  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.2) Transform

# cov = 0

# Transform file bi.res.greml.out
# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.greml.out
# 0.0658 0.0030   0.0111 0.0030   0.0116 0.0009   0.1260 0.0340     0.1306 0.0097       7310.1967
# 0.0721 0.0031   0.0065 0.0030   0.0116 0.0009   0.0722 0.0334     0.1290 0.0096       7197.3073
# 0.0672 0.0032   0.0117 0.0032   0.0115 0.0009   0.1295 0.0354     0.1273 0.0097       7194.1964
# 0.0657 0.0031   0.0105 0.0031   0.0114 0.0009   0.1196 0.0348     0.1303 0.0097       7361.4565
# 0.0667 0.0031   0.0107 0.0031   0.0132 0.0010   0.1182 0.0344     0.1458 0.0100       7227.4160
# 0.0661 0.0030   0.0078 0.0030   0.0125 0.0009   0.0901 0.0344     0.1448 0.0099       7462.7374
# 0.0645 0.0032   0.0140 0.0032   0.0115 0.0009   0.1553 0.0351     0.1281 0.0096       7228.4475
# 0.0692 0.0032   0.0085 0.0031   0.0138 0.0010   0.0926 0.0339     0.1510 0.0100       7185.9310
# 0.0648 0.0031   0.0116 0.0031   0.0122 0.0010   0.1310 0.0344     0.1382 0.0099       7322.3072
# 0.0678 0.0032   0.0101 0.0031   0.0123 0.0010   0.1115 0.0348     0.1369 0.0098       7229.7157

if [ -e bi.res.greml.out.transformed ]
then
    rm bi.res.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.res.greml.out.transformed
  done < bi.res.greml.out
  
# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.greml.out.transformed
# 0.3681858 9.9351741E-02 0.3816276 2.8344465E-02
# 0.2109763 9.7598463E-02 0.3769521 2.8052256E-02
# 0.3784132 0.1034427 0.3719846 2.8344465E-02
# 0.3494843 0.1016894 0.3807509 2.8344465E-02
# 0.3453934 0.1005206 0.4260436 2.9221097E-02
# 0.2632821 0.1005206 0.4231215 2.8928887E-02
# 0.4538037 0.1025661 0.3743222 2.8052256E-02
# 0.2705874 9.9059522E-02 0.4412386 2.9221097E-02
# 0.3827964 0.1005206 0.4038356 2.8928887E-02
# 0.3258153 0.1016894 0.4000368 2.8636677E-02

# Transform file bi.res.coregreml.out

# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.coregreml.out
# 0.0657 0.0030   0.0115 0.0030   0.0122 0.0010   -0.0017 0.0012   0.1337209 3.5622749E-02 0.1418605 1.2971596E-02 -0.1435225 9.9912226E-02 7311.2162
# 0.0721 0.0031   0.0067 0.0030   0.0120 0.0010   -0.0009 0.0012   7.5280897E-02 3.4235440E-02 0.1348315 1.2312169E-02 -0.1003724 0.1316023 7197.6217
# 0.0672 0.0032   0.0114 0.0032   0.0110 0.0010   0.0015 0.0012   0.1231102 3.4903757E-02 0.1187905 1.1513434E-02 0.1339499 0.1129285 7194.9544
# 0.0654 0.0031   0.0112 0.0031   0.0123 0.0010   -0.0025 0.0012   0.1334922 3.7069891E-02 0.1466031 1.3464548E-02 -0.2129994 0.1010828 7363.3587
# 0.0667 0.0031   0.0108 0.0031   0.0133 0.0011   -0.0001 0.0013   0.1192053 3.4951422E-02 0.1467991 1.2826558E-02 -8.3437692E-03 0.1069064 7227.4225
# 0.0661 0.0030   0.0076 0.0030   0.0122 0.0010   0.0009 0.0012   8.6659066E-02 3.4170900E-02 0.1391106 1.2166515E-02 9.3466461E-02 0.1256599 7463.0482
# 0.0645 0.0032   0.0145 0.0032   0.0126 0.0011   -0.0032 0.0012   0.1701878 3.7838388E-02 0.1478873 1.3316347E-02 -0.2367449 8.7957807E-02 7231.7636
# 0.0694 0.0032   0.0079 0.0031   0.0131 0.0011   0.0021 0.0012   8.3509512E-02 3.3109974E-02 0.1384778 1.1837375E-02 0.2064287 0.1335897 7187.3881
# 0.0647 0.0031   0.0117 0.0031   0.0123 0.0010   -0.0002 0.0012   0.1325028 3.5200141E-02 0.1392978 1.2438717E-02 -1.6671877E-02 9.8236732E-02 7322.3283
# 0.0678 0.0032   0.0101 0.0032   0.0124 0.0011   -0.0002 0.0012   0.1123471 3.5295583E-02 0.1379310 1.2481890E-02 -1.7871395E-02 0.1103516 7229.7291

if [ -e bi.res.coregreml.out.transformed ]
then
    rm bi.res.coregreml.out.transformed
fi
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
  done < bi.res.coregreml.out
  
#   [alh-admvhp1@hscpl-statgen Data]$ head bi.res.coregreml.out.transformed
# 0.3907472 0.1040936 0.4145320 3.7904426E-02 -0.4193885 0.2919545
# 0.2199790 0.1000397 0.3939925 3.5977513E-02 -0.2932992 0.3845564
# 0.3597415 0.1019926 0.3471189 3.3643521E-02 0.3914163 0.3299895
# 0.3900789 0.1083223 0.4283904 3.9344888E-02 -0.6224077 0.2953750
# 0.3483310 0.1021319 0.4289631 3.7480611E-02 -2.4381410E-02 0.3123922
# 0.2532273 9.9851124E-02 0.4064964 3.5551894E-02 0.2731193 0.3671920
# 0.4973074 0.1105679 0.4321429 3.8911831E-02 -0.6917946 0.2570224
# 0.2440240 9.6750982E-02 0.4046474 3.4590110E-02 0.6032073 0.3903638
# 0.3871877 0.1028587 0.4070435 3.6347300E-02 -4.8717055E-02 0.2870585
# 0.3282906 0.1031376 0.4030496 3.6473453E-02 -5.2222177E-02 0.3224595

# cov = 0.2

# Transform file bi.resp.greml.out
if [ -e bi.resp.greml.out.transformed ]
then
    rm bi.resp.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resp.greml.out.transformed
  done < bi.resp.greml.out

# Transform file bi.resp.coregreml.out
if [ -e bi.resp.coregreml.out.transformed ]
then
    rm bi.resp.coregreml.out.transformed
fi
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
  done < bi.resp.coregreml.out

# > res_coregreml=read.table("bi.resp.coregreml.out", header=F)
# > res_coregreml_transformed=read.table("bi.resp.coregreml.out.transformed", header=F)
# > head(res_coregreml)
#       V1     V2     V3     V4     V5     V6     V7     V8         V9        V10       V11        V12       V13        V14      V15
# 1 0.0731 0.0033 0.0095 0.0033 0.0123 0.0011 0.0080 0.0013 0.08566276 0.02977700 0.1109107 0.01023556 0.7400748 0.18715860 6811.479
# 2 0.0706 0.0034 0.0124 0.0034 0.0143 0.0012 0.0050 0.0014 0.11556380 0.03134278 0.1332712 0.01133990 0.3754838 0.12198210 6785.320
# 3 0.0687 0.0034 0.0131 0.0034 0.0123 0.0011 0.0081 0.0013 0.11876700 0.03083862 0.1115141 0.01008273 0.6381118 0.14207930 6860.251
# 4 0.0618 0.0034 0.0214 0.0035 0.0113 0.0011 0.0066 0.0013 0.19870010 0.03251858 0.1049211 0.01010849 0.4244219 0.09675611 6860.354
# 5 0.0711 0.0036 0.0145 0.0036 0.0137 0.0012 0.0085 0.0014 0.12467760 0.03100586 0.1177988 0.01028974 0.6030796 0.13333590 6617.617
# 6 0.0680 0.0034 0.0152 0.0034 0.0130 0.0011 0.0077 0.0013 0.13620070 0.03059928 0.1164875 0.01033382 0.5477688 0.12169220 6773.480
# > head(res_coregreml_transformed)
#          V1         V2        V3         V4       V5        V6
# 1 0.2503160 0.08701166 0.3240933 0.02990944 2.162580 0.5468980
# 2 0.3376901 0.09158705 0.3894331 0.03313643 1.097205 0.3564451
# 3 0.3470502 0.09011383 0.3258564 0.02946285 1.864633 0.4151713
# 4 0.5806235 0.09502285 0.3065910 0.02953810 1.240207 0.2827320
# 5 0.3643216 0.09060253 0.3442210 0.03006776 1.762265 0.3896222
# 6 0.3979934 0.08941446 0.3403893 0.03019654 1.600641 0.3555980

# cov = -0.2

# Transform file bi.resn.greml.out
if [ -e bi.resn.greml.out.transformed ]
then
    rm bi.resn.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resn.greml.out.transformed
  done < bi.resn.greml.out

# Transform file bi.resn.coregreml.out
if [ -e bi.resn.coregreml.out.transformed ]
then
    rm bi.resn.coregreml.out.transformed
fi
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
  done < bi.resn.coregreml.out

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# 5) For binary data, correct the transformation, scale
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.1) Simulate phenotype data & fit
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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.res.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.res.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    #awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resp.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resp.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resn.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resn.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    #++++++++++++++++++++++++

  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.2) Transform

# cov = 0

# Transform file bi.res.greml.out
# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.greml.out
# 0.0658 0.0030   0.0105 0.0030   0.0114 0.0009   0.1193 0.0340     0.1303 0.0096       7355.0683
# 0.0730 0.0031   0.0062 0.0030   0.0117 0.0009   0.0682 0.0333     0.1285 0.0096       7160.8152
# 0.0684 0.0033   0.0117 0.0033   0.0118 0.0010   0.1272 0.0353     0.1288 0.0098       7113.5014
# 0.0673 0.0031   0.0099 0.0031   0.0119 0.0009   0.1109 0.0347     0.1338 0.0098       7281.5921
# 0.0656 0.0031   0.0112 0.0031   0.0129 0.0010   0.1250 0.0346     0.1440 0.0099       7272.7813
# 0.0657 0.0030   0.0080 0.0030   0.0125 0.0009   0.0932 0.0345     0.1447 0.0099       7472.1840
# 0.0646 0.0031   0.0134 0.0032   0.0115 0.0009   0.1493 0.0350     0.1283 0.0096       7254.0113
# 0.0693 0.0032   0.0089 0.0031   0.0140 0.0010   0.0967 0.0340     0.1520 0.0100       7151.3153
# 0.0645 0.0031   0.0116 0.0030   0.0122 0.0009   0.1315 0.0344     0.1383 0.0099       7336.7993
# 0.0676 0.0032   0.0096 0.0031   0.0123 0.0010   0.1071 0.0348     0.1374 0.0098       7265.4862

if [ -e bi.res.greml.out.transformed ]
then
    rm bi.res.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.res.greml.out.transformed
  done < bi.res.greml.out
  
# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.greml.out.transformed
# 0.3486077 9.9351741E-02 0.3807509 2.8052256E-02
# 0.1992879 9.7306259E-02 0.3754911 2.8052256E-02
# 0.3716924 0.1031505 0.3763677 2.8636677E-02
# 0.3240620 0.1013972 0.3909783 2.8636677E-02
# 0.3652637 0.1011050 0.4207838 2.8928887E-02
# 0.2723406 0.1008128 0.4228293 2.8928887E-02
# 0.4362710 0.1022738 0.3749067 2.8052256E-02
# 0.2825680 9.9351741E-02 0.4441607 2.9221097E-02
# 0.3842575 0.1005206 0.4041278 2.8928887E-02
# 0.3129580 0.1016894 0.4014979 2.8636677E-02

# Transform file bi.res.coregreml.out

# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.coregreml.out
# 0.0657 0.0030   0.0108 0.0030   0.0120 0.0010   -0.0016 0.0012   0.1266120 3.5455160E-02 0.1406800 1.2896303E-02 -0.1405457 0.1028758 7355.8954
# 0.0730 0.0031   0.0064 0.0031   0.0121 0.0010   -0.0011 0.0012   7.1668535E-02 3.4285117E-02 0.1354983 1.2369169E-02 -0.1250000 0.1350913 7161.2496
# 0.0685 0.0033   0.0113 0.0033   0.0112 0.0010   0.0019 0.0012   0.1191983 3.4617614E-02 0.1181435 1.1397997E-02 0.1688906 0.1161139 7114.7358
# 0.0669 0.0031   0.0107 0.0031   0.0128 0.0011   -0.0025 0.0013   0.1252927 3.6933146E-02 0.1498829 1.3531250E-02 -0.2136206 0.1032312 7283.4724
# 0.0656 0.0031   0.0112 0.0031   0.0129 0.0011   -0.0000 0.0013   0.1248606 3.5186116E-02 0.1438127 1.2742001E-02 0.0000000E+00 NaN 7272.7818
# 0.0657 0.0030   0.0078 0.0030   0.0121 0.0010   0.0009 0.0012   8.9244850E-02 3.4286886E-02 0.1384439 1.2177112E-02 9.2640847E-02 0.1241829 7472.5088
# 0.0647 0.0031   0.0139 0.0032   0.0125 0.0010   -0.0031 0.0012   0.1637220 3.7640661E-02 0.1472320 1.3299043E-02 -0.2351794 8.9887783E-02 7257.2303
# 0.0694 0.0032   0.0084 0.0032   0.0133 0.0011   0.0020 0.0012   8.8328078E-02 3.3246342E-02 0.1398528 1.1919702E-02 0.1892189 0.1282998 7152.6277
# 0.0645 0.0031   0.0117 0.0031   0.0123 0.0010   -0.0003 0.0012   0.1331058 3.5283826E-02 0.1399317 1.2466669E-02 -2.5007816E-02 9.7674049E-02 7336.8300
# 0.0676 0.0032   0.0096 0.0031   0.0124 0.0010   -0.0002 0.0012   0.1076233 3.5312459E-02 0.1390135 1.2509165E-02 -1.8330889E-02 0.1123090 7265.4976

if [ -e bi.res.coregreml.out.transformed ]
then
    rm bi.res.coregreml.out.transformed
fi
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
  done < bi.res.coregreml.out
  
# [alh-admvhp1@hscpl-statgen Data]$ head bi.res.coregreml.out.transformed
# 0.3699741 0.1036039 0.4110824 3.7684415E-02 -0.4106899 0.3006144
# 0.2094233 0.1001849 0.3959409 3.6144070E-02 -0.3652637 0.3947516
# 0.3483105 0.1011565 0.3452283 3.3306200E-02 0.4935169 0.3392976
# 0.3661191 0.1079227 0.4379743 3.9539799E-02 -0.6242229 0.3016529
# 0.3648564 0.1028177 0.4202365 3.7233524E-02 0.0000000E+00 NaN
# 0.2607833 0.1001901 0.4045483 3.5582859E-02 0.2707067 0.3628761
# 0.4784136 0.1099901 0.4302281 3.8861264E-02 -0.6872200 0.2626620
# 0.2581044 9.7149462E-02 0.4086653 3.4830678E-02 0.5529184 0.3749061
# 0.3889498 0.1031032 0.4088958 3.6428977E-02 -7.3075585E-02 0.2854143
# 0.3144871 0.1031869 0.4062127 3.6553156E-02 -5.3564869E-02 0.3281792

# cov = 0.2

# Transform file bi.resp.greml.out
if [ -e bi.resp.greml.out.transformed ]
then
    rm bi.resp.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resp.greml.out.transformed
  done < bi.resp.greml.out

# Transform file bi.resp.coregreml.out
if [ -e bi.resp.coregreml.out.transformed ]
then
    rm bi.resp.coregreml.out.transformed
fi
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
  done < bi.resp.coregreml.out

# cov = -0.2

# Transform file bi.resn.greml.out
if [ -e bi.resn.greml.out.transformed ]
then
    rm bi.resn.greml.out.transformed
fi
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.greml.out.transformed
    
    # for transcriptome
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resn.greml.out.transformed
  done < bi.resn.greml.out

# Transform file bi.resn.coregreml.out
if [ -e bi.resn.coregreml.out.transformed ]
then
    rm bi.resn.coregreml.out.transformed
fi
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
  done < bi.resn.coregreml.out

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# 6) For binary data, correct the transformation, scale, k = 0.01
#===============================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6.1) Simulate phenotype data & fit
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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat" --kval=0.01
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.res.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.res.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat" --kval=0.01
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    #awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resp.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resp.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat" --kval=0.01
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resn.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resn.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    #++++++++++++++++++++++++

  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# 7) For binary data, correct the transformation, scale, k = 0.5
#===============================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7.1) Simulate phenotype data & fit
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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat" --kval=0.5
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.res.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.res.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat" --kval=0.5
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    #awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resp.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resp.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    #++++++++++++++++++++++++

  done

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
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i}
    
    ../Tool/mtg2 -p fam_10K.fam -g grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p fam_10K.fam -g brain.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: brain.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat" --kval=0.5
    
    # model fitting
    # o	inside brain_greml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p fam_10K.fam -mg brain_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    
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
    # o	Inside brain_coregreml.out3
     # ******************************************************************
     # MTG2 version 2.21 (Apr2021)
     # ******************************************************************
     # delta2    : brain_coregreml.do
     # 
     # var-cov-info matrix file: brain_coregreml.out2
     # Ratio:  0.1337209     SE:  3.5622749E-02  p-value:  1.7416898E-04
     # Ratio:  0.1418605     SE:  1.2971596E-02  p-value:  7.7769562E-28
     # Cor  : -0.1435225     SE:  9.9912226E-02  p-value:  0.1508640 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resn.coregreml.out
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resn.coregreml.out
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    #++++++++++++++++++++++++

  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================
