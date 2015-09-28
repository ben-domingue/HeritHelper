## plink.file="/hrsshare/cleaned/v2_hrs_geno_final_translated"
## gwas.file="/tmp/GWAS.result"
## wd="/tmp/grs/"
## out.name<-"smoke"

## system("awk '{print $2, $4, $5, $11, $9}' ~/gwas_results/tag.evrsmk.tbl > /tmp/GWAS.result") 
## #make_pgs(out.name="eversmoke") 
## setwd("/tmp/smoke")

make_pgs<-function(plink.file="hrs_geno_final_translated",gwas.file="/tmp/GWAS.result",wd="/tmp/grs/",out.name,clump=TRUE) {
    tr<-list()
    #################################
    getwd()->orig.dir
    system(paste("mkdir ",wd))
    setwd(wd)
    #################################
    #i'm just creating symbolic links here so that i can work in one directory. not really imperative that you do something like this, but if you don't youll need to change directory stuff downstream.
    system(paste("ln -s ",plink.file,".bed ./gen.bed",sep=""))
    system(paste("ln -s ",plink.file,".bim ./gen.bim",sep=""))
    system(paste("ln -s ",plink.file,".fam ./gen.fam",sep=""))
    #################################
    #get agct snps from gwas file
    #system(paste("awk '{print $2}' gen.bim > snps.txt"))
    #system("plink --bfile gen --extract snps.txt --make-bed --out gen2")
    #remove ambiguous snps (strand issues). remember that we're only going to use SNPs which have a quickly identifiable strand.
    system(paste("awk '{ print $1, $2 $3, $4, $5}' ",gwas.file," > temp"))
    read.table("temp")->tmp
    toupper(tmp[,2])->tmp[,2]
    write.table(tmp,file="temp",quote=FALSE,row.names=FALSE,col.names=FALSE)
    #
    "awk '{ if ($2 == ffACff || $2 == ffAGff || $2 == ffCAff || $2 == ffCTff || $2 == ffGAff || $2== ffGTff || $2 == ffTCff || $2 == ffTGff ) print $0}' temp > GWAS.noambig"->txt
    gsub('ff','"',txt)->txt
    system(txt)
    #
    system('echo "SNP Allele1Allele2 P W" > head.txt')
    system("cat head.txt GWAS.noambig > GWAS2.noambig")
    #################################
    #get snps from plink files
    #system("awk '{print $2}' gen.bim > available.snps")
    system("awk '{print $1}' GWAS.noambig > available.snps")
    #system("plink --bfile gen --extract available.snps --make-bed --out gen --silent")
    #system("wc -l gen.bim",intern=TRUE)->tr$common
    #now get those from gwas
    read.table("available.snps")->gwas
    read.table("gen.bim",header=TRUE)->data
    intersect(gwas[,1],data[,2])->common
    length(common) -> tr$common
    #write.table(common,file="common.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
    #################################
    if (clump) {
        #Clump data in 2 rounds using plink2
        #1st clumping & extract tops snps for 2nd round
        fun<-function(i) {
            paste("plink --bfile gen --chr ",i,"  --clump GWAS2.noambig  --clump-p1 1 --clump-p2 1 --clump-r2 .5 --clump-kb 250 --out traitX",i,".round1 --silent",sep="")->cmd
            system(cmd)
            system(paste("awk '{print $3, $5}' traitX",i,".round1.clumped > traitX",i,".round2.input",sep=""))
            system(paste("awk '{print $3}' traitX",i,".round1.clumped > traitX",i,".extract2",sep=""))
            cmd
        }
        library(parallel)
        makeCluster(22)->cl #this was based on a big machine with 22 cores. you could not do this on a smaller machine and just run the below
        clusterApply(cl,1:22,fun)->garbage
        garbage[[1]]->tr$clump1
        #2nd clumping & extract tops snps for profile
        fun<-function(i) {
            paste("plink  --bfile gen --chr ",i," --extract traitX",i,".extract2 --clump traitX",i,".round2.input --clump-p1 1 --clump-p2 1 --clump-r2 .2 --clump-kb 5000 --out traitX",i,".round2 --silent",sep="")->cmd
            system(cmd)
            system(paste("awk '{print $3}' traitX",i,".round2.clumped > traitX",i,".selected",sep=""))
            cmd
        }
        clusterApply(cl,1:22,fun)->garbage
        garbage[[1]]->tr$clump2
        stopCluster(cl)
        system("cat traitX1.selected traitX2.selected traitX3.selected traitX4.selected traitX5.selected traitX6.selected traitX7.selected traitX8.selected traitX9.selected traitX10.selected traitX11.selected traitX12.selected traitX13.selected traitX14.selected traitX15.selected traitX16.selected traitX17.selected traitX18.selected traitX19.selected traitX20.selected traitX21.selected traitX22.selected > traitX.selected")
    } else {
        system("echo 'SNP' > /tmp/head.txt")
        system("cat /tmp/head.txt available.snps > traitX.selected")
    }
    #################################
    # The traitX"$i".selected files will contain the lists of top snps
    # Merge the alleles, effect & P values onto these files
    #R
    read.table("traitX.selected",header=TRUE)->selected
    read.table(gwas.file,header=TRUE)->effects
    #
    effects[!duplicated(effects),]->effects
    #
    effects[,1] %in% selected$SNP -> index
    effects[index,]->effects
    effects[,2]<-toupper(effects[,2])
    effects[,3]<-toupper(effects[,3])
    #################################
    #make sure strands are aligned
    bim <- read.table(file="gen.bim",header=FALSE)
    names(bim)<-c("chr","snp","a","b","a1.bim","a2.bim")
    names(effects)<-c("snp","a1.eff","a2.eff","pv","beta")
    NULL->bim$a->bim$b
    merge(bim,effects,by="snp")->test
    #table(test$a1.bim,test$a1.eff)
    #get rid of ambig strands from bim (already yanked from gwas data)
    test$a1.bim=="T" & test$a2.bim=="A" -> i1
    test$a1.bim=="A" & test$a2.bim=="T" -> i2
    test[!(i1 | i2),]->test
    test$a1.bim=="C" & test$a2.bim=="G" -> i1
    test$a1.bim=="G" & test$a2.bim=="C" -> i2
    test[!(i1 | i2),]->test
    #now make everything a/c
    for (nm in c("a1.bim","a2.bim","a1.eff","a2.eff")) {
        ifelse(test[[nm]]=="T","A",test[[nm]])->test[[nm]]
        ifelse(test[[nm]]=="G","C",test[[nm]])->test[[nm]]
    }
    ifelse(test$a1.bim!=test$a1.eff,-1*test$beta,test$beta)->test$beta
    #old
    ## ifelse(test$a1.eff==test$a1.bim | test$a1.eff==test$a2.bim,1,0)->flip
    ## ifelse(flip==0 & test$a1.eff=="A","T",test$a1.eff)->test$a1.eff
    ## ifelse(flip==0 & test$a1.eff=="T","A",test$a1.eff)->test$a1.eff
    ## ifelse(flip==0 & test$a1.eff=="C","G",test$a1.eff)->test$a1.eff
    ## ifelse(flip==0 & test$a1.eff=="G","C",test$a1.eff)->test$a1.eff
    test[,c("snp","a1.eff","beta")]->z
    write.table(z,file="score_file.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
    nrow(z)->tr$final.n
    #################################
    #create score!
    setwd(orig.dir)
    system(paste("plink --bfile ",wd,"gen --score ",wd,"score_file.txt --silent --out ",out.name,sep=""))
    dump("tr",file=paste(out.name,".metadata",sep=""))
    system(paste("rm -r ",wd))
    tr
}




## read.table("/tmp/smoke/smoke.profile",header=TRUE)->grs.new
## read.table("~/hrs/smoking/grs/hrs_eversmoke_score.profile",header=TRUE)->grs.old
## plot(grs.old[,6],grs.new[,6]); abline(0,1)




