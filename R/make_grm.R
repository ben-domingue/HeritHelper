make_grm<-function(ids, #this is a nx2 (family, individual) matrix of ids that you want to use for creation of grm
                   #you probably don't need to change this next option
                   plink.file.path="/hrsshare/analytic_sample/hrs_geno_final", #path/stem of the plink file name you want to use
                   out.name="/tmp/hrs_tmp", #path/stem to where you want output
                   extra.txt=NULL,np=10
                   )
{
  write.table(ids,file="/tmp/ids_tmp.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  paste("gcta64 --bfile",plink.file.path," --keep /tmp/ids_tmp.txt --autosome --maf 0.01 --make-grm --out ",out.name," --thread-num",np)->cmd
  if (!is.null(extra.txt)) paste(cmd,extra.txt)->cmd
  system(cmd)
}
  
