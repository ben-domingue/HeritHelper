this is meant to make some of the work for the SSM paper easier.

@@@begin code
library(HeritHelper)

load("/hrsshare/hrs_linked.Rdata")
x[x$raracem=="1.white/caucasian",]->x
x[x$rahispan=="0. not hispanic",]->x
x[x$rabyear>=1930 & x$rabyear<1951,]->x

x[,c("family","subjectID")]->ids
make_grm(ids,np=20)

@@@end code

