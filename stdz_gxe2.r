

  cmat=as.matrix(read.table("legp.vcv"))

  chol=read.table("grm_qced_rdm_005_ukbb3_10K_processed.grm.bend.chol.matvec")
  yv1=scale(chol$V1)
  yv1=yv1*cmat[1,1]^.5
  
  chol=read.table("brain.bmat.bend.chol.matvec")
  yv2=scale(chol$V1)
  yv2=yv2*cmat[2,2]^.5

  yv=read.table("rnd.v3")
  yv3=scale(yv$V1)
  yv3=yv3*cmat[3,3]^.5

  sink("mvnorm.out")
  write.table (cbind(yv1,yv2,yv3),row.names=F,col.names=F,quote=F)
  sink()


