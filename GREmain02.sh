#===============================
# B1) Continuous data
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B1.1) Simulate phenotype data & fit
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
# B2) For binary data, correct the transformation, scale
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B2.1) Simulate phenotype data & fit
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
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside brain_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_10K_processed.grm
    # brain.bmat
    # brain_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p fam_10K.fam -mg brain_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
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
#===============================

#===============================
#===============================

# Correlation between grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol and brain.bmat.bend.chol

# o	ouput: grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol
# [alh-admvhp1@hscpl-statgen Data]$ wc -l grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol
# 50005000 grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol
# [alh-admvhp1@hscpl-statgen Data]$ head grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol
#                      1                     1   1.006697
#                      2                     1  2.7675326E-03
#                      2                     2   1.003686
#                      3                     1  6.4125904E-03
#                      3                     2  6.6622100E-03
#                      3                     3   1.006144
#                      4                     1 -1.8581041E-03
#                      4                     2  1.0506219E-03
#                      4                     3 -1.9141618E-03
#                      4                     4   1.001036

# o	output: brain.bmat.bend.chol
# [alh-admvhp1@hscpl-statgen Data]$ wc -l brain.bmat.bend.chol
# 50005000 brain.bmat.bend.chol
# [alh-admvhp1@hscpl-statgen Data]$ head brain.bmat.bend.chol
#                      1                     1   1.012592
#                      2                     1 -1.7754473E-02
#                      2                     2   1.011695
#                      3                     1 -3.1189149E-02
#                      3                     2 -6.9634723E-03
#                      3                     3   1.038747
#                      4                     1 -2.0474184E-03
#                      4                     2 -4.0480113E-03
#                      4                     3  2.4522025E-02
#                      4                     4  0.9794471

awk '($1 != $2) {printf $3"\n"}' grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol > f1
# [alh-admvhp1@hscpl-statgen Data]$ wc -l f1
# 49995000 f1
# [alh-admvhp1@hscpl-statgen Data]$ head f1
# 2.7675326E-03
# 6.4125904E-03
# 6.6622100E-03
# -1.8581041E-03
# 1.0506219E-03
# -1.9141618E-03
# -2.5497659E-03
# -3.3871667E-04
# 5.0336082E-04
# -1.0990744E-02

awk '($1 != $2) {printf $3"\n"}' brain.bmat.bend.chol > f2
# [alh-admvhp1@hscpl-statgen Data]$ wc -l f2
# 49995000 f2
# [alh-admvhp1@hscpl-statgen Data]$ head f2
# -1.7754473E-02
# -3.1189149E-02
# -6.9634723E-03
# -2.0474184E-03
# -4.0480113E-03
# 2.4522025E-02
# -3.0175031E-03
# 5.6893909E-03
# -2.2660891E-02
# -4.3839455E-02

### R ###
#R

# setwd("C:/Users/phamvv/MyDoc/31SNP/GRE/Data")
setwd("/data/alh-admvhp1/GRE/Data")

f1=read.table("f1", header=F)
f2=read.table("f2", header=F)
# > head(f1)
#             V1
# 1  0.002767533
# 2  0.006412590
# 3  0.006662210
# 4 -0.001858104
# 5  0.001050622
# 6 -0.001914162
# > head(f2)
#             V1
# 1 -0.017754473
# 2 -0.031189149
# 3 -0.006963472
# 4 -0.002047418
# 5 -0.004048011
# 6  0.024522025

cor(f1, f2)
# > cor(f1, f2)
#           V1
# V1 0.1002223

# > cov(f1,f2)
#              V1
# V1 5.179013e-06

# q() # exit R
### R ###

#===============================
#===============================