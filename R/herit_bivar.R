herit_bivar<-function(df, #data frame with all the variables.
                      pheno, #names of outcome
                      grm="/tmp/hrs_tmp", #default matches make_grm default.
                      id.cols=c("family","subjectID"), #pass a character vector with their names
                      qcovar=NULL, #names of quantitative covariates
                      covar=NULL, #names of dichotomous covariates
                      cutoff=0.025,
                      out.name="/tmp/hrs_tmp",
                      np=10,
                      extra.txt=NULL
                       )
{
  if (length(pheno)!=2) stop("need two phenotypes")
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
  paste("gcta64 --reml-bivar 1 2 --grm-gz ",grm," --grm-cutoff",cutoff,pheno.txt,qcovar.txt,covar.txt,fe.txt,"--out",out.name,"--thread-num",np)->cmd
  if (!is.null(extra.txt)) 
    cmd <- paste(cmd, extra.txt)
  system(cmd,intern=TRUE)->txt
  txt
}
