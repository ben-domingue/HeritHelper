make_grm<-function(ids, #this is a nx2 (family, individual) matrix of ids that you want to use for creation of grm
                   plink.file.path="/hrsshare/analytic_sample/hrs_geno_final", #stem of the plink file name you want to use
                   out.name="hrs_tmp",
                   np=10
                   ) {
  write.table(ids,file="/tmp/ids_tmp.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  system("cd /tmp")
  paste("gcta64 --bfile",plink.file.path," --keep /tmp/ids_tmp.txt --autosome --maf 0.01 --make-grm --out ",out.name," --thread-num",np)->cmd
  system(cmd)
}
  
