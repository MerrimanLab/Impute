#########################################################################################
# This gets the rsids from a plink .bim file & compares to the rsIDS from imputed 1000G
# - see how many are in common
#######################################################################################

rev.comp<-function(x,rev=TRUE) {
  x<-toupper(x)
  y<-x
        y[x=="A"]<-"T"
        y[x=="C"]<-"G"
        y[x=="G"]<-"C"
        y[x=="T"]<-"A"
  y
}

args <- commandArgs(T)

i <- args[1]
bim.file <- args[2]
dir<-args[3]
path<-args[4]

setwd(paste(dir, "/original_data",sep=""))
print(getwd())
out.dir<-paste(dir,"/fixstrand/", sep="")

bim <- read.delim(bim.file, header=F,sep="\t",skip=0,fill=FALSE,stringsAsFactors=FALSE)
bim.num <- 0
bim.tot <- dim(bim)[1]

### 1000Genomes merged file -  have to do it per xsome
genomes.file <- paste(path,"/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3_chr",i,".legend", sep="")

out.delete <- paste(out.dir,"_chr",i,".delete", sep="")
out.flip <- paste(out.dir,"_chr",i,".flip", sep="")
out.pos <- paste(out.dir,"_chr",i,".pos", sep="")
out.chr <- paste(out.dir,"_chr",i,".chr", sep="")

genomes <-read.table(genomes.file, header=T,stringsAsFactors=F)
genomes_rs <- genomes[,1]
f <- function(s) strsplit(s, ":")[[1]][1]
g_rs <-sapply(genomes_rs, f)
posns <- match(bim[,2],g_rs)
missing<-is.na(posns)
sum(!missing)
rs.match <- bim[!missing,]
match<-genomes[posns[!missing],]

dim(rs.match)
dim(match)

print(rs.match[1:10,])
print(match[1:10,])

flip.me <- {}
flip.me2 <- {}
rsRemove <- {}

## Alleles need to be flipped by plink as on wrong strand
to.flip <- (rs.match[,5] != match[,3])  & (rs.match[,5] != match[,4]) & (rs.match[,6] != match[,3]) & (rs.match[,6] != match[,4]) #{ ## FLIP as allles match the wrong way around

## THIS IS WHAT WE WILL BE using plink to enforce the REF ALLELE
flip.me <- rs.match[to.flip,2]
flip.me[1:5]

## CHECK THE A-T & C-G SNPs
to.freq <- ( (rs.match[,5] == "A")  & (rs.match[,6] == "T")) |  ( (rs.match[,5] == "T")  & (rs.match[,6] == "A")) |  ( (rs.match[,5] == "C")  & (rs.match[,6] == "G")) |  ( (rs.match[,5] == "G")  & (rs.match[,6] == "C"))

match[to.freq,][1:10,]
rs.match[to.freq,][1:10,]

atcg.match <- match[to.freq,]
atcg.rs.match <- rs.match[to.freq,]

atcg.match[1:10,]
atcg.rs.match[1:10,]

minthr <- .45
maxthr <- .55
posns <- (atcg.match[,"eur.aaf"] >=  minthr & atcg.match[,"eur.aaf"]  <= maxthr)

missing<-is.na(posns)
atcg.match.new <- atcg.match[!posns,]
atcg.rs.match.new <- atcg.rs.match[!posns,]

## These A-T & C-G the MAF is too ambigous
rsRemove <-  atcg.match[posns,1]
dim(atcg.match.new)
dim(atcg.rs.match.new)

atcg.match <- atcg.match.new
atcg.rs.match <- atcg.rs.match.new

to.flip.atcg <- ((atcg.rs.match[,5] == atcg.match[,"a1"]) &  (as.numeric(atcg.match[,8]) > 0.5) ) |  ((atcg.rs.match[,5] == atcg.match[,"a0"]) &  (as.numeric(atcg.match[,8]) < 0.5) )  #|  ((atcg.rs.match[,6] == atcg.match[,"a1"]) &  (as.numeric(atcg.match[,8]) < 0.5) )

atcg.match[to.flip.atcg,][1:10,]
atcg.rs.match[to.flip.atcg,][1:10,]

dim(atcg.match[to.flip.atcg,])
flip.me2 <- atcg.rs.match[to.flip.atcg,2]

## Rev comp allelese to flip & A-T & C-G ones as well
write.table(flip.me,out.flip,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
write.table(flip.me2,out.flip,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=TRUE)

### A-T & C-G allelese to delete - as the MAF is too close to tell. So these will be imputed
write.table(rsRemove,out.delete,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

## WRITE the positions that have changed
to.write <- rs.match[,4] != match[,2]
print.pos <- cbind(rs.match[to.write,2], match[to.write,2])
write.table(print.pos,out.pos,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

###Write SNPs to be kept
print.chr <- cbind(rs.match[,2],i)
write.table(print.chr,out.chr,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
