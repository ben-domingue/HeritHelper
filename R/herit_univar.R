herit_univar<-function(df, #data frame with all the variables.
                       pheno,
                       id.cols=c("family","subjectID"), #pass a character vector with their names
                       qcovar=NULL,
                       covar=NULL,
                       cutoff=0.025
                       )
{
  #gcta64 --reml --reml-est-fix --grm hrs --grm-cutoff 0.025  --pheno hrs_graduate.txt --qcovar hrs_covars.txt  --out hrs_graduate_PCs --thread-num 20                         
}


