this is meant to make some of the work for the SSM paper easier.

@@@begin code

library(HeritHelper)

load("/hrsshare/hrs_linked.Rdata")

x[x$raracem=="1.white/caucasian",]->x

x[x$rahispan=="0. not hispanic",]->x

x[x$rabyear>=1930 & x$rabyear<1951,]->x

read.table("/hrsshare/cleaned/hrs_mds.mds",header=TRUE)->mds

NULL->mds$SOL

names(mds)[1:2]<-c("family","subjectID")

merge(x,mds)->x


make_grm(ids=x[,c("family","subjectID")],np=20)

herit_univar(x,pheno="r8bmi",np=20)->h1

herit_univar(x,pheno="r8bmi",qcovar=c("C1","C2"),np=20)->h2

herit_univar(x,pheno="r8bmi",qcovar=c("raedyrs","C1","C2"),np=20)->h3

@@@end code

