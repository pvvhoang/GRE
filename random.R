library("optparse")
library(MASS)

# Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
# Output: 
# rnd.v
# rnd.v2
# rnd.v3

option_list = list(
  make_option(c("-c", "--cov"), type = "numeric"),
  make_option(c("-s", "--seed"), type = "numeric"),
  make_option(c("-n", "--number"), type = "numeric", default = 10000)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #=================================================
# opt$cov=0
# opt$seed=1
# opt$number=10000
# #=================================================
  
  set.seed(opt$seed)

  tn=opt$number
  #tn
  gec <- opt$cov
  ge_cor=matrix(0,3,3)
  ge_cor[1,1]=1
  ge_cor[2,2]=1
  ge_cor[3,3]=1
  ge_cor[1,2]=gec/0.4
  ge_cor[2,1]=gec/0.4
  
  v1=mvrnorm(tn,c(0,0,0),ge_cor, empirical=TRUE)
  
  # > tn=10000
  # > gec <- 0.2
  # > ge_cor=matrix(0,3,3)
  # > ge_cor[1,1]=1
  # > ge_cor[2,2]=1
  # > ge_cor[3,3]=1
  # > ge_cor[1,2]=gec/0.4
  # > ge_cor[2,1]=gec/0.4
  # > v1=mvrnorm(tn,c(0,0,0),ge_cor, empirical=TRUE)
  # > cov(v1[,1],v1[,2])
  # [1] 0.5
  # > cov(v1[,1],v1[,3])
  # [1] -2.175926e-15
  # > cov(v1[,2],v1[,3])
  # [1] -4.165496e-16
  # > cov(v1[,1],v1[,1])
  # [1] 1
  # > cov(v1[,2],v1[,2])
  # [1] 1
  # > cov(v1[,3],v1[,3])
  # [1] 1
  # > v1new <- v1[,1]*0.4^.5
  # > v2new <- v1[,2]*0.4^.5
  # > v3new <- v1[,3]*0.2^.5
  # > cov(v1new,v2new)
  # [1] 0.2
  # > cov(v1new,v3new)
  # [1] -6.154706e-16
  # > cov(v2new,v3new)
  # [1] -1.177306e-16
  # > cov(v1new,v1new)
  # [1] 0.4
  # > cov(v2new,v2new)
  # [1] 0.4
  # > cov(v3new,v3new)
  # [1] 0.2
  
  sink("rnd.v")
  write.table(v1[,1],quote=F,col.name=F,row.name=F)
  sink()
  sink("rnd.v2")
  write.table(v1[,2],quote=F,col.name=F,row.name=F)
  sink()
  sink("rnd.v3")
  write.table(v1[,3],quote=F,col.name=F,row.name=F)
  sink()
  