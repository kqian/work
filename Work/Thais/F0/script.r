bamFls = list.files(path="../../../Lars/Thais/results/2013-10-06/",pattern=".bam$",full.names=TRUE)

bedFls = list.files(path="../F1_sperm/annotation/",pattern="rn",full.names=TRUE)
  
source("../../scripts/dataAnalysis/snowCount.r")

data = snowCount(bamFls,bedFls,cpus=3)
  

gene.anno <- import("../F1_sperm/annotation/rn4.gtf")

miRNA.anno<- import("../F1_sperm/annotation/rno2.gff3")



#piRNA database http://www.ibab.ac.in/pirna/Rat.tar.gz
piRNA.anno<-  import("../F1_sperm/annotation/rn4_piRNA.gtf")


#tiRNA
library(Rsamtools)
library(rtracklayer)
library(AnnotationHub)
hub <- AnnotationHub()
filters(hub) <- list(Species="Rattus norvegicus")  
tRNA<-hub$goldenpath.rn4.database.tRNAs_0.0.1.RData
# tRNA.count<- summarizeOverlaps(tRNA,BamFileList(bamFls, index=character()))

#Tandem Repeats Finder
trf<-hub$goldenpath.rn4.database.simpleRepeat_0.0.1.RData
trf.count<- summarizeOverlaps(trf,BamFileList(bamFls, index=character()))
#Repeatmasker
#from ucsc
repeatmask<- import("../F1_sperm/annotation/rn4_repeatmasker.bed")
library(rtracklayer)

elementMetadata(miRNA.anno)[,5]->rownames(data$counts$rno2.gff3)

pd<-read.table(file="./expSpec.txt",sep="\t",header=T,stringsAsFactors=F)
pd<-pd[c(1,8,2:7,9,18,10:17),]
#miRNA DE
count.miRNA<-data$counts$rno2.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA"),]
f0.miRNA<-de.edgeR(count.miRNA,grp,"F0_miRNA_hfat_vs_ctrl_Trend_Filter.csv","ctrl","hfat",F,T)
f0.miRNA.seq<-fetch.region(aln,f0.miRNA,miRNA.anno)

#piRNA DE
count.piRNA<-data$counts$rn4_piRNA.gtf
rownames(count.piRNA) <- paste(elementMetadata(piRNA.anno)[,5],seqnames(piRNA.anno),start(piRNA.anno),end(piRNA.anno),strand(piRNA.anno),sep=";")
F0.piRNA <- de.edgeR.pi(count.piRNA,grp,"F0_piRNA_hfat_vs_ctrl_Trend_Filter.csv","ctrl","hfat",F,T)

#get seq
library(BSgenome.Rnorvegicus.UCSC.rn4)
piRNA.anno2<-piRNA.anno
elementMetadata(piRNA.anno2)[,5]<-paste(elementMetadata(piRNA.anno)[,5],seqnames(piRNA.anno),start(piRNA.anno),end(piRNA.anno),strand(piRNA.anno),sep=";")
F0_hfat_vs_ctrl.piRNA.seq<-fetch.region.pi(aln,F0.piRNA,piRNA.anno2)
write.csv(cbind(F0.piRNA[,1:4],F0_hfat_vs_ctrl.piRNA.seq,F0.piRNA[,5:22]),"F0_piRNA_hfat_vs_ctrl_seq.csv")

#DE tiRNA
count.tRNA<-data$counts$rn4_tRNA.bed
#check the order first, make sure the name are same
rownames(count.tRNA) <- elementMetadata(tRNA)[,1]
tRNA.5<-resize(tRNA,width=25)
tRNA.3<-resize(tRNA,width=25,fix="end")
tRNA.M<-narrow(tRNA,start=26,end=width(tRNA)-25)
elementMetadata(tRNA.5)[,1]<-paste(elementMetadata(tRNA.5)[,1],"5",sep=".")
elementMetadata(tRNA.3)[,1]<-paste(elementMetadata(tRNA.3)[,1],"3",sep=".")
elementMetadata(tRNA.M)[,1]<-paste(elementMetadata(tRNA.M)[,1],"M",sep=".")
tiRNA.anno<-c(tRNA.5,tRNA.3,tRNA.M)
count.tiRNA<-matrix(0,1332,18)
for (i in 1:18) count.tiRNA[,i]<- countOverlaps(tiRNA.anno, aln[[i]], maxgap=0L, minoverlap=10L)
colnames(count.tRNA)->colnames(count.tiRNA)
rownames(count.tiRNA) <- elementMetadata(tiRNA.anno)[,1]
#de.edgeR.pi == de.edgeR.ti
F0.tiRNA <- de.edgeR.pi(count.tiRNA,grp,"F0_tiRNA_hfat_vs_ctrl.csv","ctrl","hfat",F,T)
F0.tiRNA.seq<-fetch.region.ti(aln,F0.tiRNA,tiRNA.anno)
write.csv(cbind(F0.tiRNA[,1:4],F0.tiRNA.seq,F0.tiRNA[,5:22]),"F0_tiRNA_hfat_vs_ctrl_seq.csv")

#ggbio vis
library(ggbio)

temp.anno<-piRNA.anno[which(elementMetadata(piRNA.anno)[,5]==strsplit(row.names(F0.piRNA)[1],";")[[1]][1])]
temp.anno<-piRNA.anno[grep("DQ738740",elementMetadata(piRNA.anno)[,5])]
temp.anno<-tiRNA.anno[grep(strsplit(row.names(F0.tiRNA)[19],"\\.")[[1]][2],elementMetadata(tiRNA.anno)[,1])]
temp.anno<-tiRNA.anno[which(elementMetadata(tiRNA.anno)[,1]==row.names(F0.tiRNA)[19])]

temp.anno<-miRNA.anno[which(elementMetadata(miRNA.anno)[,5]==row.names(f0.miRNA)[12])]

for (i in 1:length(temp.anno)) {
  temp<-lapply(aln,subsetByOverlaps,temp.anno[i])
  temp<-temp[lapply(temp,length)>0]
  print(temp)
}

temp<-sort(subsetByOverlaps(aln[[1]],temp.anno,minoverlap=1L))
p1<-
  autoplot(temp)
p2<-autoplot(Rnorvegicus, which =reduce(temp))
p3<-autoplot(Rnorvegicus, which =temp.anno)
tracks(p1,p2,p3)

autoplot(Rnorvegicus, which =temp.anno[1])


sessionInfo()->sess.info
save.image("countData.RData")

write.csv(cbind(f0.miRNA[,1:8],f0.miRNA.seq,f0.miRNA[,9:27]),"f0.miRNA_seq.csv")
#check reads anno distribution
fetch.region<-function(alignlist, result,anno) {
  seqa<-c()
  seqb<-c()
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  for (i in 1:length(result[,1])) {
    print(i)
    temp.anno<-result$ID[i]
    temp.anno<-anno[elementMetadata(anno)$ID==temp.anno]
    seq1<-toString(getSeq(Rnorvegicus,temp.anno))
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
      seq2<-toString(getSeq(Rnorvegicus,temp.anno)) 
    }
    seqa<-c(seqa,seq1)
    seqb<-c(seqb,seq2)
  }        
  return(cbind(seqa,seqb))
}

fetch.region.pi<-function(alignlist, result,anno) {
  seqa<-c()
  seqb<-c()
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  for (i in 1:length(result[,1])) {
    print(i)
    temp.anno<-anno[which(elementMetadata(anno)[,5]==rownames(result[i,]))]
    seq1<-toString(getSeq(Rnorvegicus,temp.anno))
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
      seq2<-toString(getSeq(Rnorvegicus,temp.anno)) 
    }
    seqa<-c(seqa,seq1)
    seqb<-c(seqb,seq2)
  }        
  return(cbind(seqa,seqb))
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
      seq1<-toString(getSeq(Rnorvegicus,unique(range.read)))
    } else {
      range.read<-unique(range.read)[which.max(countOverlaps(unique(range.read),range.read,type="equal"))]
      seq1<-toString(getSeq(Rnorvegicus,unique(range.read)))
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
      seq2<-toString(getSeq(Rnorvegicus,temp.anno)) 
    }
    seqa<-c(seqa,seq1)
    seqb<-c(seqb,seq2)
  }        
  return(cbind(seqa,seqb))
}
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
#0.05
de.edgeR.pi<-function (counts,group,output,conA,conB,Tagwise=T,filter=T) {
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
