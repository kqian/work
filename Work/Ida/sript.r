#This is Ida Episperm 1 sRNA-seq data
#working path

#input bam
bamFls = list.files(path="./data/2013-10-11/bowtie/",pattern=".bam$",full.names=TRUE)
#input bedFls
#check the chromosome name are the same as bam file head
bedFls = list.files(path="./annotation/",pattern="hsa",full.names=TRUE)
#counting  
source("scripts/dataAnalysis/snowCount.r")
data = snowCount(bamFls,bedFls,cpus=4)

#get annotation
#ensembl
gene.anno <- import("annotation/hsa_ensembl.gtf")
#miRBase v20
miRNA.anno<- import("annotation/hsa.gff3")
#piRNA database http://www.ibab.ac.in/pirna/Human.tar.gz
piRNA.anno<-  import("annotation/hsa_piRNA.gtf")
#tiRNA
#from ucsc Table 
tRNA.anno<-  import("annotation/hsa19_tRNA.bed")
#lincRNA
#from ucsc Table 
lincRNA.anno<-  import("annotation/hsa19_lincRNA.bed")
#Tandem Repeats Finder
#from ucsc Table 
trf<-import("annotation/hsa19_TRF.bed")
#Repeatmasker
#from ucsc Table 
repeatmask<- import("annotation/hsa19_repeatmasker.bed")

elementMetadata(miRNA.anno)[,5]->rownames(data$counts$hsa.gff3)

#piRNA http://www.ibab.ac.in/pirna/Human.tar.gz

save.image("countData.RData")


library(Rsamtools)
#keep all alignment into a list of GRanges
aln<-list()
for (i in 1:length(bamFls)) aln<-c(aln,granges(readGAlignmentsFromBam(bamFls[i])))
save(aln,file="aln.RData")

#
setwd("Ida/")
load("countData.RData")
#load phenodata
pd<-read.table(file="expSpec",sep="\t",header=T,stringsAsFactors=F)
#order 'pd' as sample order
pd<-pd[c(10:19,1,20:23,2:9),]
library(edgeR)

count.miRNA<-data$counts$hsa.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA"),]
count.pri.miRNA<-data$counts$hsa.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA_primary_transcript"),]

hist(rowSums(count.miRNA),xlim=c(0,49),1000000)
hist(rowSums(count.pri.miRNA),xlim=c(0,49),2000000)

grp <- factor(pd[,4])


de.edgeR(count.miRNA,grp,"miRNAtest.csv")

t1<-de.edgeR(count.miRNA,grp,"miRNA_obese_vs_lean_Tag_Filter.csv","lean","obese",T,T)
t2<-de.edgeR(count.miRNA,grp,"miRNA_obese_vs_lean_Trend_Filter.csv","lean","obese",F,T)
t3<-de.edgeR(count.miRNA,grp,"miRNA_obese_vs_lean_Tag_All.csv","lean","obese",T,F)
t4<-de.edgeR(count.miRNA,grp,"miRNA_obese_vs_lean_Trend_All.csv","lean","obese",F,F)

#keep t2 as final list
t2.2<-de.edgeR(count.pri.miRNA,grp,"Pri_miRNA_obese_vs_lean_Trend_Filter.csv","lean","obese",F,T)
  