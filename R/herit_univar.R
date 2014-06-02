herit_univar<-function(df, #data frame with all the variables.
                       pheno, #name of outcome
                       grm="/tmp/hrs_tmp", #default matches make_grm default.
                       id.cols=c("family","subjectID"), #pass a character vector with their names
                       qcovar=NULL, #names of quantitative covariates
                       covar=NULL, #names of dichotomous covariates
                       cutoff=0.025,
                       out.name="/tmp/hrs_tmp",
                       np=10,
                       extra.txt=NULL,
                       return.call=FALSE
                       )
{
  df[,c(id.cols,pheno)]->pheno
  paste0(out.name,".pheno")->fn
  write.table(pheno,file=fn,quote=FALSE,row.names=FALSE,col.names=FALSE)
  paste("--pheno ",fn)->pheno.txt
  if (!is.null(qcovar)) {
    df[,c(id.cols,qcovar)]->qcovar
    fn<-paste0(out.name,".qcovar")
    write.table(qcovar,file=fn,quote=FALSE,row.names=FALSE,col.names=FALSE)
    paste("--qcovar ",fn)->qcovar.txt
  } else ""->qcovar.txt
  if (!is.null(covar)) {
    df[,c(id.cols,covar)]->covar
    paste0(out.name,".covar")->fn
    write.table(covar,file=fn,quote=FALSE,row.names=FALSE,col.names=FALSE)
    paste("--covar ",fn)->covar.txt
  } else ""->covar.txt
  ifelse(""==covar.txt & ""==qcovar.txt,"","--reml-est-fix")->fe.txt
  paste("gcta64 --reml --grm-gz ",grm," --grm-cutoff",cutoff,pheno.txt,qcovar.txt,covar.txt,fe.txt,"--out",out.name,"--thread-num",np)->cmd
  if (!is.null(extra.txt)) paste(cmd,extra.txt)->cmd
  if (!return.call) {
    system(cmd,intern=TRUE)->txt
    #
    tr<-list()
    grep("V(G)/Vp",txt,fixed=TRUE)->index
    tr$h2<-strsplit(txt[index],"\t")[[1]][2:3]
    #
    grep("Estimatesof fixed effects:",txt)->i1
    grep("Summary result of REML analysis has been saved in the file",txt)->i2
    if (length(i1)>0 & length(i2)>0) txt[i1:(i2-1)]->tr$fe
    #
    tr$all<-txt
    tr
  } else {
    cmd
  }
}
