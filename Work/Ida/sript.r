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
#piRNA database http://www.ibab.ac.in/pirna/Human.tar.gz are NCBI 36.
#use piRNA cluster anno instead
#http://www.uni-mainz.de/FB/Biologie/Anthropologie/493_ENG_HTML.php
piRNA.anno<-  import("annotation/hsa_piRNA.gtf")
#piRNA anno 
count.piRNA<-matrix(0,length(piRNA.anno),23)
#aln from aln.RData
#only count reads with at least 10nt overlap
for (i in 1:23) count.piRNA[,i]<- countOverlaps(piRNA.anno, aln[[i]],,minoverlap=15L)
data$counts$hsa_piRNA.gtf<-count.piRNA
colnames(count.piRNA)<-colnames(data$counts$hsa19_tRNA.bed)

rownames(count.piRNA) <-paste(elementMetadata(piRNA.anno)[,5],seqnames(piRNA.anno),start(piRNA.anno),end(piRNA.anno),strand(piRNA.anno),sep=";")

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

#intron from ucsc gene
intron<-import("annotation/intron_ucsc.bed")

save.image("countData.RData")

unanno<-gaps(c(reduce(gene.anno),reduce(miRNA.anno),reduce(piRNA.anno),reduce(trf),reduce(repeatmask),reduce(tRNA.anno),reduce(lincRNA.anno),reduce(intron)))
count.unRNA<-matrix(0,length(unanno),23)
for (i in 1:23) count.unRNA[,i]<- countOverlaps(unanno, aln[[i]],type= "within")

#count all sRNA loci without annotation
sapply(aln,reduce)->aln.re
for (i in 1:length(aln.re)) aln.re[[i]]<-aln.re[[i]][which(countOverlaps(aln.re[[i]],aln[[i]])>1)]
aln.re2<-aln.re[[1]]
for (i in 2:length(aln)) aln.re2<-c(aln.re2,aln.re[[i]])
aln.re3<-reduce(aln.re2)
rm(aln.re,aln.re2)

count.all<-matrix(0,length(aln.re3),length(aln))
for (i in 1:length(aln)) count.all[,i]<- countOverlaps(aln.re3, aln[[i]])
dge <- DGEList(count.all, group=grp)
keep<- rowSums(cpm(dge) > 1) >= (length(grp)/2)
dge <- dge[keep,]
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTrendedDisp(dge)
et <- exactTest(dge, pair=c("lean","obese"))
et.adj.p.val <- p.adjust(et$table$PValue, method = "BH")
aln.re4<-aln.re3[keep]
aln.re4[which(et$table$PValue<0.05)]->aln.re5
et2<-et[which(et$table$PValue<0.05),]$table
et2<-cbind(et2,et.adj.p.val[which(et$table$PValue<0.05)])
et2<-cbind(et2,countOverlaps(aln.re5,gene.anno))
et2<-cbind(et2,countOverlaps(aln.re5,miRNA.anno))
et2<-cbind(et2,countOverlaps(aln.re5,piRNA.anno))
et2<-cbind(et2,countOverlaps(aln.re5,tRNA.anno))
et2<-cbind(et2,countOverlaps(aln.re5,lincRNA.anno))
et2<-cbind(et2,countOverlaps(aln.re5,trf))
et2<-cbind(et2,countOverlaps(aln.re5,repeatmask))
et2<-cbind(et2,countOverlaps(aln.re5,intron))
et2<-cbind(et2,rowSums(et2[,5:12]))
et2<-cbind(et2,width(aln.re5))
colnames(et2)<-c("logFC","logCPM","PValue","FDR","Ensembl","miRNA","piRNA","tRNA","lincRNA","trf","repeatmark","intron","anno_sum","width")

count.all2<-count.all[keep,]
paste(grp,sapply(strsplit(bamFls,"\\//"),tail,1),sep=".")->colnames(count.all2)
count.all2<-count.all2[,order(colnames(count.all2))]

et2<-cbind(et2,count.all2[which(et$table$PValue<0.05),])
write.csv(et2[order(et2$PValue),],"sRNA_obese_vs_lean2.csv")
DEsiRNA.seq<-fetch.region.si(aln,aln.re5)

rm(aln)
save.image("countData.RData")


library(Rsamtools)
#keep all alignment into a list of GRanges
aln<-list()
for (i in 1:length(bamFls)) aln<-c(aln,granges(readGAlignmentsFromBam(bamFls[i])))
save(aln,file="aln.RData")

library(BSgenome.Hsapiens.UCSC.hg19)
library(ggbio)

temp<-sort(subsetByOverlaps(aln[[6]],temp.anno,minoverlap=1L))
p1<-
  autoplot(temp)
p2<-autoplot(Hsapiens, which =reduce(temp))
p3<-autoplot(Hsapiens, which =temp.anno)
tracks(p1,p2,p3,heights=c(5,1,1))

setwd("Ida/")
load("countData.RData")
#load phenodata
pd<-read.table(file="expSpec",sep="\t",header=T,stringsAsFactors=F)
#order 'pd' as sample order as bamFls
pd<-pd[c(10:19,1,20:23,2:9),]
grp <- factor(pd[,4])

library(edgeR)
elementMetadata(miRNA.anno)[,5]->rownames(data$counts$hsa.gff3)
count.miRNA<-data$counts$hsa.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA"),]
#count.pri.miRNA<-data$counts$hsa.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA_primary_transcript"),]
hist(rowSums(count.miRNA),xlim=c(0,49),1000000)
hist(rowSums(count.pri.miRNA),xlim=c(0,49),2000000)

#miRNA DE
DEmiRNA<-de.edgeR(count.miRNA,grp,"miRNA_obese_vs_lean_Trend_Filter.csv","lean","obese",F,T)
DEmiRNA.seq<-fetch.region(aln,DEmiRNA,miRNA.anno[which(elementMetadata(miRNA.anno)[,5] %in% rownames(DEmiRNA))])
write.csv(cbind(DEmiRNA[,1:8],DEmiRNA.seq,DEmiRNA[,9:32]),"miRNA_obese_vs_lean_seq.csv")

#piRNA DE
#count.piRNA<-data$counts$hsa_piRNA.gtf
#give location info in the rownames
#rownames(count.piRNA) <- elementMetadata(piRNA.anno)[,5]
DEpiRNA <- de.edgeR.pi(count.piRNA,grp,"piRNA_obese_vs_lean.csv","lean","obese",F,T)

#get piRNA seq
#build an anno same as rownames
piRNA.anno2<-piRNA.anno
elementMetadata(piRNA.anno2)[,5]<-paste(elementMetadata(piRNA.anno)[,5],seqnames(piRNA.anno),start(piRNA.anno),end(piRNA.anno),strand(piRNA.anno),sep=";")
DEpiRNA.seq<-fetch.region.pi(aln,DEpiRNA,piRNA.anno2[which(elementMetadata(piRNA.anno2)[,5] %in% rownames(DEpiRNA))])
DEpiRNA.c <- countOverlaps(piRNA.anno2[which(elementMetadata(piRNA.anno2)[,5] %in% rownames(DEpiRNA))],piRNA.anno2)
names(DEpiRNA.c)<- elementMetadata(piRNA.anno2[which(elementMetadata(piRNA.anno2)[,5] %in% rownames(DEpiRNA))])[,5]
DEpiRNA.c<-DEpiRNA.c[rownames(DEpiRNA)]
write.csv(cbind(DEpiRNA[,1:4],DEpiRNA.c,DEpiRNA.seq,DEpiRNA[,5:27]),"piRNA_obese_vs_lean_seq.csv")

#DE tRFs/tiRNA
#count.tRNA<-data$counts$hsa19_tRNA.bed
#build new tiRNA anno by seperate know tRNA to 3 region: 5'(25nt)/ M /3'(25nt) 
tRNA.5<-resize(tRNA.anno[1:623],width=25)
tRNA.3<-resize(tRNA.anno[1:623],width=25,fix="end")
tRNA.M<-narrow(tRNA.anno[1:623],start=26,end=width(tRNA.anno[1:623])-25)
elementMetadata(tRNA.5)[,1]<-paste(elementMetadata(tRNA.5)[,1],"5",sep=".")
elementMetadata(tRNA.3)[,1]<-paste(elementMetadata(tRNA.3)[,1],"3",sep=".")
elementMetadata(tRNA.M)[,1]<-paste(elementMetadata(tRNA.M)[,1],"M",sep=".")
tiRNA.anno<-c(tRNA.5,tRNA.3,tRNA.M)
count.tiRNA<-matrix(0,1869,23)
#aln from aln.RData
#only count reads with at least 10nt overlap
for (i in 1:23) count.tiRNA[,i]<- countOverlaps(tiRNA.anno, aln[[i]], maxgap=0L, minoverlap=10L)
colnames(count.tiRNA)<-colnames(data$counts$hsa19_tRNA.bed)
rownames(count.tiRNA) <- elementMetadata(tiRNA.anno)[,1]
#de.edgeR.pi == de.edgeR.ti
DEtiRNA <- de.edgeR.pi(count.tiRNA,grp,"DEtiRNA_obese_vs_lean.csv","lean","obese",F,T)
DEtiRNA.seq<-fetch.region.ti(aln,DEtiRNA,tiRNA.anno[which(elementMetadata(tiRNA.anno)[,1] %in% rownames(DEtiRNA))])
write.csv(cbind(DEtiRNA[,1:4],DEtiRNA.seq,DEtiRNA[,5:27]),"DEtiRNA_obese_vs_lean_seq.csv")
####Auto find tRNAf region by filter out 10% background reads
auto.region.ti<-function(alignlist,anno) {
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  new.range<-GRanges()
  
  for (i in 1:length(anno)) {
    range.read<-GRanges()
    print(i)
    temp.anno<-anno[i]
    temp<-sapply(temp.aln, subsetByOverlaps, temp.anno, minoverlap=10L)   
    temp<-temp[lapply(temp,length)>0]
    if (length(temp) == 0) next
    #get the tRNAf ranges
      for (j in 1:length(temp)) {
      #count each region
      newgr<-unique(temp[[j]])
      newgr.c<-countOverlaps(newgr,temp[[j]],type="equal")
      newgr<-newgr[which(newgr.c > (length(temp[[j]])*0.1))]
      range.read<-reduce(c(newgr,range.read))
        }  
    
    new.range<-c(new.range,range.read)
    }       
  return(new.range)
}

tRNAf.anno<-auto.region.ti(aln,tRNA.anno)

#0.05
de.edgeR<-function (counts,group,output,conA,conB,Tagwise=T,filter=T) {
  library(edgeR)
  library(GenomicRanges)
  counts=counts[,which(group %in% c(conA,conB))]
  group=factor(group[which(group %in% c(conA,conB))])
  
  dge <- DGEList(counts, group=group)
  
  #filter, as least expressed in half the libraries
  if (filter==T) {
    keep<- rowSums(cpm(dge) > 1) >= (length(group)/2)
    dge <- dge[keep,]
  }
  
  dge <- calcNormFactors(dge)
  
  dge <- estimateCommonDisp(dge)
  if (Tagwise==T) dge <- estimateTagwiseDisp(dge) else   dge <- estimateTrendedDisp(dge)
  
  et <- exactTest(dge, pair=c(conA,conB))
  npvalue<-length(which(et$table$PValue<0.05))
  
  tt<-topTags(et,n=npvalue)
  
  #get miRNA anno
  tab2<-as.data.frame(elementMetadata(miRNA.anno)[which( elementMetadata(miRNA.anno) [,5] %in% rownames(data.frame(tt))), c(2,5:8)])
  
  #get read count
  tab3<-counts[c(rownames(data.frame(tt))),]
  colnames(tab3)<-paste(as.character(group),colnames(tab3),sep=".")
  
  tab1<-tt[order(row.names(tt)),]
  tab2<-tab2[order(tab2$ID),]
  tab3<-tab3[order(row.names(tab3)),order(colnames(tab3))]
  result.tab<-cbind(tab1,tab2,tab3)
  result.tab<-result.tab[order(result.tab$PValue),]
  write.csv(result.tab,file=output)
  return(result.tab)
}

fetch.region<-function(alignlist, result,anno) {
  seqa<-c()
  seqb<-c()
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  for (i in 1:length(result[,1])) {
    print(i)
    temp.anno<-result$ID[i]
    temp.anno<-anno[elementMetadata(anno)$ID==temp.anno]
    seq1<-toString(getSeq(Hsapiens,temp.anno))
    seq2<-"NA"
    temp<-sapply(temp.aln,subsetByOverlaps,temp.anno)
    temp.s<-sapply(temp,start)
    temp.s<-sapply(temp.s, table)
    temp.s<-temp.s[lapply(temp.s,length)>0]
    temp.s<-table(names(sapply(temp.s,which.max)))
    start(temp.anno)<-as.numeric(names(which.max(temp.s)))
    
    temp.e<-sapply(temp,end)
    temp.e<-sapply(temp.e, table)
    temp.e<-temp.e[lapply(temp.e,length)>0]
    temp.e<-table(names(sapply(temp.e,which.max)))
    if (start(temp.anno)< as.numeric(names(which.max(temp.e)))) {
      end(temp.anno)<-as.numeric(names(which.max(temp.e)))
      seq2<-toString(getSeq(Hsapiens,temp.anno)) 
    }
    seqa<-c(seqa,seq1)
    seqb<-c(seqb,seq2)
  }        
  return(cbind(seqa,seqb))
}

de.edgeR.pi<-function (counts,group,output,conA,conB,Tagwise=T,filter=T) {
  library(edgeR)
  library(GenomicRanges)
  counts=counts[,which(group %in% c(conA,conB))]
  group=factor(group[which(group %in% c(conA,conB))])
  
  dge <- DGEList(counts, group=group)
  #dge<-dge[which(rowSums(dge$counts)>0),]
  #filter, as least expressed in half the libraries
  if (filter==T) {
    keep<- rowSums(cpm(dge) > 1) >= (length(group)/2)
    dge <- dge[keep,]
  }
  
  dge <- calcNormFactors(dge)
  
  dge <- estimateCommonDisp(dge)
  if (Tagwise==T) dge <- estimateTagwiseDisp(dge) else   dge <- estimateTrendedDisp(dge)
  
  et <- exactTest(dge, pair=c(conA,conB))
  npvalue<-length(which(et$table$PValue<0.05))
  
  tt<-topTags(et,n=npvalue)
  
  tt<-data.frame(tt)  
  #get read count
  tab3<-counts[row.names(tt),]
  colnames(tab3)<-paste(as.character(group),colnames(tab3),sep=".")
  
  tab1<-tt[order(row.names(tt)),]
  tab3<-tab3[order(row.names(tab3)),order(colnames(tab3))]
  result.tab<-cbind(tab1,tab3)
  result.tab<-result.tab[order(result.tab$PValue),]
  write.csv(result.tab,file=output)
  return(result.tab) 
}

fetch.region.pi<-function(alignlist, result,anno) {
  seqa<-c()
  seqb<-c()
  seqc<-c()
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  for (i in 1:length(result[,1])) {
    print(i)
    temp.anno<-anno[which(elementMetadata(anno)[,5]==rownames(result[i,]))]
    seq1<-toString(getSeq(Hsapiens,temp.anno))
    seq2<-"NA"
    seq3<-"NA"
       
    temp<-sapply(temp.aln,subsetByOverlaps,temp.anno, minoverlap=15L) 
    temp<-temp[lapply(temp,length)>0]
    temp.s<-sapply(temp,start)
    temp.s<-sapply(temp.s, table)
    temp.s<-temp.s[lapply(temp.s,length)>0]
    temp.s<-table(names(sapply(temp.s,which.max)))
    start(temp.anno)<-as.numeric(names(which.max(temp.s)))
    
    temp.e<-sapply(temp,end)
    temp.e<-sapply(temp.e, table)
    temp.e<-temp.e[lapply(temp.e,length)>0]
    temp.e<-table(names(sapply(temp.e,which.max)))
    if (start(temp.anno)< as.numeric(names(which.max(temp.e)))) {
      end(temp.anno)<-as.numeric(names(which.max(temp.e)))
      seq2<-toString(getSeq(Hsapiens,temp.anno)) 
    }
    
    #get the most abundent reads ranges
    range.read<-GRanges()
    for (j in 1:length(temp)) {
      anno.j<-unique(temp[[j]])[which.max(countOverlaps(unique(temp[[j]]),temp[[j]],type="equal"))]
      range.read<-c(range.read,anno.j)
    }
    if (length(unique(range.read)) == 1) {
      seq3<-toString(getSeq(Hsapiens,unique(range.read)))
    } else {
      range.read<-unique(range.read)[which.max(countOverlaps(unique(range.read),range.read,type="equal"))]
      seq3<-toString(getSeq(Hsapiens,unique(range.read)))
    }
    
    seqa<-c(seqa,seq1)
    seqb<-c(seqb,seq2)
    seqc<-c(seqc,seq3)
  }        
  return(cbind(seqa,seqb,seqc))
}

fetch.region.ti<-function(alignlist, result,anno) {
  seqa<-c()
  seqb<-c()
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  for (i in 1:length(result[,1])) {
    print(i)
    temp.anno<-anno[which(elementMetadata(anno)[,1]==rownames(result[i,]))]
    seq1<-"NA"
    seq2<-"NA"
    temp<-sapply(temp.aln, subsetByOverlaps, temp.anno, minoverlap=10L)   
    temp<-temp[lapply(temp,length)>0]
    #get the most abundent reads ranges
    range.read<-GRanges()
    for (j in 1:length(temp)) {
      anno.j<-unique(temp[[j]])[which.max(countOverlaps(unique(temp[[j]]),temp[[j]],type="equal"))]
      range.read<-c(range.read,anno.j)
    }
    if (length(unique(range.read)) == 1) {
      seq1<-toString(getSeq(Hsapiens,unique(range.read)))
    } else {
      range.read<-unique(range.read)[which.max(countOverlaps(unique(range.read),range.read,type="equal"))]
      seq1<-toString(getSeq(Hsapiens,unique(range.read)))
    }
    temp.s<-sapply(temp,start)
    temp.s<-sapply(temp.s, table)
    # temp.s<-temp.s[lapply(temp.s,length)>0]
    temp.s<-table(names(sapply(temp.s,which.max)))
    start(temp.anno)<-as.numeric(names(which.max(temp.s)))
    
    temp.e<-sapply(temp,end)
    temp.e<-sapply(temp.e, table)
    # temp.e<-temp.e[lapply(temp.e,length)>0]
    temp.e<-table(names(sapply(temp.e,which.max)))
    if (start(temp.anno)< as.numeric(names(which.max(temp.e)))) {
      end(temp.anno)<-as.numeric(names(which.max(temp.e)))
      seq2<-toString(getSeq(Hsapiens,temp.anno)) 
    }
    seqa<-c(seqa,seq1)
    seqb<-c(seqb,seq2)
  }        
  return(cbind(seqa,seqb))
}

fetch.region.si<-function(alignlist,anno) {
  seqa<-c()
  seqb<-c()
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  for (i in 1:length(anno)) {
    print(i)
    temp.anno<-anno[i]
    seq1<-"NA"
    seq2<-"NA"
    temp<-sapply(temp.aln, subsetByOverlaps, temp.anno, minoverlap=15L)   
    temp<-temp[lapply(temp,length)>0]
    if (length(temp) == 0) next
    #get the most abundent reads ranges
    range.read<-GRanges()
    for (j in 1:length(temp)) {
      anno.j<-unique(temp[[j]])[which.max(countOverlaps(unique(temp[[j]]),temp[[j]],type="equal"))]
      range.read<-c(range.read,anno.j)
    }
    if (length(unique(range.read)) == 1) {
      seq1<-toString(getSeq(Hsapiens,unique(range.read)))
    } else {
      range.read<-unique(range.read)[which.max(countOverlaps(unique(range.read),range.read,type="equal"))]
      seq1<-toString(getSeq(Hsapiens,unique(range.read)))
    }
    temp.s<-sapply(temp,start)
    temp.s<-sapply(temp.s, table)
    # temp.s<-temp.s[lapply(temp.s,length)>0]
    temp.s<-table(names(sapply(temp.s,which.max)))
    start(temp.anno)<-as.numeric(names(which.max(temp.s)))
    
    temp.e<-sapply(temp,end)
    temp.e<-sapply(temp.e, table)
    # temp.e<-temp.e[lapply(temp.e,length)>0]
    temp.e<-table(names(sapply(temp.e,which.max)))
    if (start(temp.anno)< as.numeric(names(which.max(temp.e)))) {
      end(temp.anno)<-as.numeric(names(which.max(temp.e)))
      seq2<-toString(getSeq(Hsapiens,temp.anno)) 
    }
    seqa<-c(seqa,seq1)
    seqb<-c(seqb,seq2)
  }        
  return(cbind(seqa,seqb))
}
