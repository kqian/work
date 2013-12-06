#This is Thais Rat sperm F1 sRNA-seq data
#working path
#larsroed@padawan.cbs.dtu.dk:/home/people/larsroed/projects/Kui/Thais/F1_sperm
#input files, DTU server
bamFls = list.files(path="./data/2013-10-18/",pattern=".bam$",full.names=TRUE)
#annotation files, DTU server
#check the chromosome name are the same as bam file head
bedFls = list.files(path="./annotation/",pattern="rn",full.names=TRUE)
source("scripts/dataAnalysis/snowCount.r")
data = snowCount(bamFls,bedFls,cpus=3)

#get annotation
#ensembl
#Rattus_norvegicus.RGSC3.4.69.gtf.gz
gene.anno <- import("annotation/rn4.gtf")
#miRbase  v20
miRNA.anno<- import("annotation/rno2.gff3")
#piRNA database http://www.ibab.ac.in/pirna/Rat.tar.gz
piRNA.anno<-  import("annotation/rn4_piRNA.gtf")

#Tandem Repeats Finder
trf<-hub$goldenpath.rn4.database.simpleRepeat_0.0.1.RData
trf.count<- summarizeOverlaps(trf,BamFileList(bamFls, index=character()))
#Repeatmasker
#from ucsc
repeatmask<- import("annotation/rn4_repeatmasker.bed")

library(Rsamtools)
#keep all alignment into a list of GRanges
aln<-list()
for (i in 1:length(bamFls)) aln<-c(aln,granges(readGAlignmentsFromBam(bamFls[i])))
save(aln,file="align.RData")

library(rtracklayer)
library(BSgenome.Rnorvegicus.UCSC.rn4)

#following code is executed in iMac xpm526@sund-it-mac-4:/Users/xpm526/Work/Thais/F0
#get tRNA annotation
#becareful about the order of this vs the one from UCSC, make sure the tRNA ID be consistent

#tiRNA
library(AnnotationHub)
hub <- AnnotationHub()
filters(hub) <- list(Species="Rattus norvegicus")  
tRNA<-hub$goldenpath.rn4.database.tRNAs_0.0.1.RData
# tRNA.count<- summarizeOverlaps(tRNA,BamFileList(bamFls, index=character()))


#vis % of data
boxplot( sapply(data$counts, colSums) *100/data$reads)

#ensem<-split(data$counts$rn4_repeatmasker.bed, as.factor(elementMetadata(repeatmask)[,1]))
ensem<-split(data$counts$rn4.gtf,elementMetadata(gene.anno)[,1])
ensem<-sapply(ensem, matrix,  ncol = 38)
ensem<-sapply(ensem,colSums)
#ensem*100/data$reads
boxplot( ensem*100/data$reads) 

boxplot( ensem[,apply(ensem * 100/data$reads, 2, median)>0.3]*100/data$reads)

grp <- factor(pd[,4])
design <- model.matrix(~0+grp)
colnames(design)<-levels(grp)
filter<-T


#miRNA DE
pd<-read.table(file="result/2013-10-18/expSpec",header=T,stringsAsFactors=F)
elementMetadata(miRNA.anno)[,5]->rownames(data$counts$rno2.gff3)
count.miRNA<-data$counts$rno2.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA"),]
dge <- DGEList(count.miRNA, group=grp)
#piRNA DE
count.piRNA<-data$counts$rn4_piRNA.gtf
rownames(count.piRNA) <- paste(elementMetadata(piRNA.anno)[,5],seqnames(piRNA.anno),start(piRNA.anno),end(piRNA.anno),strand(piRNA.anno),sep=";")
dge <- DGEList(count.piRNA, group=grp)
#tRNA DE
count.tRNA<-data$counts$rn4_tRNA.bed
rownames(count.tRNA) <- elementMetadata(tRNA)[,1]
tRNA.5<-resize(tRNA,width=25)
tRNA.3<-resize(tRNA,width=25,fix="end")
tRNA.M<-narrow(tRNA,start=26,end=width(tRNA)-25)
elementMetadata(tRNA.5)[,1]<-paste(elementMetadata(tRNA.5)[,1],"5",sep=".")
elementMetadata(tRNA.3)[,1]<-paste(elementMetadata(tRNA.3)[,1],"3",sep=".")
elementMetadata(tRNA.M)[,1]<-paste(elementMetadata(tRNA.M)[,1],"M",sep=".")
tiRNA.anno<-c(tRNA.5,tRNA.3,tRNA.M)
count.tiRNA<-matrix(0,1332,38)
for (i in 1:38) count.tiRNA[,i]<- countOverlaps(tiRNA.anno, aln[[i]], maxgap=0L, minoverlap=10L)
colnames(data$counts$rn4.gtf)->colnames(count.tiRNA)
rownames(count.tiRNA) <- elementMetadata(tiRNA.anno)[,1]
dge <- DGEList(count.tiRNA, group=grp)

if (filter==T) {
  keep<- rowSums(cpm(dge) > 1) >= (length(grp)/2)
  dge <- dge[keep,]
}
#mds plot
plotMDS(dge,col=as.numeric(grp),label=grp)

dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
#lrt <- glmLRT(fit)
#
lrt <- glmLRT(fit,contrast=c(-1,0,1,0))


#miRNA
F1_GLM_PHC_vs_PCC.miRNA<-lrt.output.miRNA(glmLRT(fit,contrast=c(-1,0,1,0)),count.miRNA,grp,"F1_miRNA_PHC_vs_PCC_GLM.csv")
F1_GLM_PCH_vs_PCC.miRNA<-lrt.output.miRNA(glmLRT(fit,contrast=c(-1,1,0,0)),count.miRNA,grp,"F1_miRNA_PCH_vs_PCC_GLM.csv")
F1_GLM_PHH_vs_PHC.miRNA<-lrt.output.miRNA(glmLRT(fit,contrast=c(0,0,-1,1)),count.miRNA,grp,"F1_miRNA_PHH_vs_PHC_GLM.csv")
F1_GLM_PHH_vs_PCH.miRNA<-lrt.output.miRNA(glmLRT(fit,contrast=c(0,-1,0,1)),count.miRNA,grp,"F1_miRNA_PHH_vs_PCH_GLM.csv")
F1_GLM_PHC_vs_PCC.miRNA.seq<-fetch.region(aln,F1_GLM_PHC_vs_PCC.miRNA,miRNA.anno)
F1_GLM_PCH_vs_PCC.miRNA.seq<-fetch.region(aln,F1_GLM_PCH_vs_PCC.miRNA,miRNA.anno)
F1_GLM_PHH_vs_PHC.miRNA.seq<-fetch.region(aln,F1_GLM_PHH_vs_PHC.miRNA,miRNA.anno)
F1_GLM_PHH_vs_PCH.miRNA.seq<-fetch.region(aln,F1_GLM_PHH_vs_PCH.miRNA,miRNA.anno)
write.csv(cbind(F1_GLM_PHC_vs_PCC.miRNA[,1:8],F1_GLM_PHC_vs_PCC.miRNA.seq,F1_GLM_PHC_vs_PCC.miRNA[,9:48]),"F1_miRNA_PHC_vs_PCC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PCH_vs_PCC.miRNA[,1:8],F1_GLM_PCH_vs_PCC.miRNA.seq,F1_GLM_PCH_vs_PCC.miRNA[,9:48]),"F1_miRNA_PCH_vs_PCC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PHH_vs_PHC.miRNA[,1:8],F1_GLM_PHH_vs_PHC.miRNA.seq,F1_GLM_PHH_vs_PHC.miRNA[,9:48]),"F1_miRNA_PHH_vs_PHC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PHH_vs_PCH.miRNA[,1:8],F1_GLM_PHH_vs_PCH.miRNA.seq,F1_GLM_PHH_vs_PCH.miRNA[,9:48]),"F1_miRNA_PHH_vs_PCH_GLM_seq.csv")

#piRNA
F1_GLM_PHC_vs_PCC.piRNA<-lrt.output.piRNA(glmLRT(fit,contrast=c(-1,0,1,0)),count.piRNA,grp,"F1_piRNA_PHC_vs_PCC_GLM.csv")
F1_GLM_PCH_vs_PCC.piRNA<-lrt.output.piRNA(glmLRT(fit,contrast=c(-1,1,0,0)),count.piRNA,grp,"F1_piRNA_PCH_vs_PCC_GLM.csv")
F1_GLM_PHH_vs_PHC.piRNA<-lrt.output.piRNA(glmLRT(fit,contrast=c(0,0,-1,1)),count.piRNA,grp,"F1_piRNA_PHH_vs_PHC_GLM.csv")
F1_GLM_PHH_vs_PCH.piRNA<-lrt.output.piRNA(glmLRT(fit,contrast=c(0,-1,0,1)),count.piRNA,grp,"F1_piRNA_PHH_vs_PCH_GLM.csv")
piRNA.anno2<-piRNA.anno
elementMetadata(piRNA.anno2)[,5]<-paste(elementMetadata(piRNA.anno)[,5],seqnames(piRNA.anno),start(piRNA.anno),end(piRNA.anno),strand(piRNA.anno),sep=";")

F1_GLM_PHC_vs_PCC.piRNA.seq<-fetch.region.pi(aln,F1_GLM_PHC_vs_PCC.piRNA,piRNA.anno2)
F1_GLM_PCH_vs_PCC.piRNA.seq<-fetch.region.pi(aln,F1_GLM_PCH_vs_PCC.piRNA,piRNA.anno2)
F1_GLM_PHH_vs_PHC.piRNA.seq<-fetch.region.pi(aln,F1_GLM_PHH_vs_PHC.piRNA,piRNA.anno2)
F1_GLM_PHH_vs_PCH.piRNA.seq<-fetch.region.pi(aln,F1_GLM_PHH_vs_PCH.piRNA,piRNA.anno2)
write.csv(cbind(F1_GLM_PHC_vs_PCC.piRNA[,1:5],F1_GLM_PHC_vs_PCC.piRNA.seq,F1_GLM_PHC_vs_PCC.piRNA[,6:43]),"F1_piRNA_PHC_vs_PCC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PCH_vs_PCC.piRNA[,1:5],F1_GLM_PCH_vs_PCC.piRNA.seq,F1_GLM_PCH_vs_PCC.piRNA[,6:43]),"F1_piRNA_PCH_vs_PCC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PHH_vs_PHC.piRNA[,1:5],F1_GLM_PHH_vs_PHC.piRNA.seq,F1_GLM_PHH_vs_PHC.piRNA[,6:43]),"F1_piRNA_PHH_vs_PHC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PHH_vs_PCH.piRNA[,1:5],F1_GLM_PHH_vs_PCH.piRNA.seq,F1_GLM_PHH_vs_PCH.piRNA[,6:43]),"F1_piRNA_PHH_vs_PCH_GLM_seq.csv")
rm(aln)
save.image("countData.RData")
 
#tiRNA
F1_GLM_PHC_vs_PCC.tiRNA<-lrt.output.tiRNA(glmLRT(fit,contrast=c(-1,0,1,0)),count.tiRNA,grp,"F1_tiRNA_PHC_vs_PCC_GLM.csv")
F1_GLM_PCH_vs_PCC.tiRNA<-lrt.output.tiRNA(glmLRT(fit,contrast=c(-1,1,0,0)),count.tiRNA,grp,"F1_tiRNA_PCH_vs_PCC_GLM.csv")
F1_GLM_PHH_vs_PHC.tiRNA<-lrt.output.tiRNA(glmLRT(fit,contrast=c(0,0,-1,1)),count.tiRNA,grp,"F1_tiRNA_PHH_vs_PHC_GLM.csv")
F1_GLM_PHH_vs_PCH.tiRNA<-lrt.output.tiRNA(glmLRT(fit,contrast=c(0,-1,0,1)),count.tiRNA,grp,"F1_tiRNA_PHH_vs_PCH_GLM.csv")
F1_GLM_PHC_vs_PCC.tiRNA.seq<-fetch.region.ti(aln,F1_GLM_PHC_vs_PCC.tiRNA,tiRNA.anno)
F1_GLM_PCH_vs_PCC.tiRNA.seq<-fetch.region.ti(aln,F1_GLM_PCH_vs_PCC.tiRNA,tiRNA.anno)
F1_GLM_PHH_vs_PHC.tiRNA.seq<-fetch.region.ti(aln,F1_GLM_PHH_vs_PHC.tiRNA,tiRNA.anno)
F1_GLM_PHH_vs_PCH.tiRNA.seq<-fetch.region.ti(aln,F1_GLM_PHH_vs_PCH.tiRNA,tiRNA.anno)
write.csv(cbind(F1_GLM_PHC_vs_PCC.tiRNA[,1:5],F1_GLM_PHC_vs_PCC.tiRNA.seq,F1_GLM_PHC_vs_PCC.tiRNA[,6:43]),"F1_tiRNA_PHC_vs_PCC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PCH_vs_PCC.tiRNA[,1:5],F1_GLM_PCH_vs_PCC.tiRNA.seq,F1_GLM_PCH_vs_PCC.tiRNA[,6:43]),"F1_tiRNA_PCH_vs_PCC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PHH_vs_PHC.tiRNA[,1:5],F1_GLM_PHH_vs_PHC.tiRNA.seq,F1_GLM_PHH_vs_PHC.tiRNA[,6:43]),"F1_tiRNA_PHH_vs_PHC_GLM_seq.csv")
write.csv(cbind(F1_GLM_PHH_vs_PCH.tiRNA[,1:5],F1_GLM_PHH_vs_PCH.tiRNA.seq,F1_GLM_PHH_vs_PCH.tiRNA[,6:43]),"F1_tiRNA_PHH_vs_PCH_GLM_seq.csv")
rm(aln)
save.image("countData.RData")

library(ggbio)
#vcf call in DTU server
samtools mpileup -uf bowtieIndexes/rn4.fa data/2013-10-18/*.bam | ../../Soft/samtools-0.1.19/bcftools/bcftools view -bvcg - > data/2013-10-18/var.raw.bcf 
../../Soft/samtools-0.1.19/bcftools/bcftools view data/2013-10-18/var.raw.bcf | ../../Soft/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D100 > data/2013-10-18/var.flt.vcf
			
library(VariantAnnotation)
vcf<-readVcf("var.flt.vcf","rn4")
autoplot(vcf,which=tRNA[1])

library(BSgenome.Rnorvegicus.UCSC.rn4)
p1<-autoplot(sort(subsetByOverlaps(aln[[1]],miRNA.anno[2])))
p2<-autoplot(Rnorvegicus, which =miRNA.anno[2])
tracks(p1,p2)

#region info

save.image("countData.RData")
#pvalue 0.05
F1_GLM_PHC_vs_PCC.miRNA.seq<-fetch.region(aln,F1_GLM_PHC_vs_PCC.miRNA,miRNA.anno)

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
  seqb<-c()
  temp.aln<-sapply(alignlist,subsetByOverlaps,anno)
  for (i in 1:length(result[,1])) {
    print(i)
    temp.anno<-anno[which(elementMetadata(anno)[,1]==rownames(result[i,]))]
    seq2<-"NA"
    temp<-sapply(temp.aln, subsetByOverlaps, temp.anno, minoverlap=10L)   
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
    seqb<-c(seqb,seq2)
  }        
  return(seqb)
}

lrt.output.miRNA <- function (lrt,counts,group,output) {
  npvalue<-length(which(lrt$table$PValue<0.05))
  tt<-topTags(lrt,n=npvalue)
  
  #get miRNA anno
  tab2<-as.data.frame(elementMetadata(miRNA.anno)[which( elementMetadata(miRNA.anno) [,5] %in% rownames(data.frame(tt))), c(2,5:8)])
  #get read count
  tab3<-counts[rownames(data.frame(tt)),]
  colnames(tab3)<-paste(as.character(group),colnames(tab3),sep=".")
  
  tab1<-tt[order(row.names(tt)),]
  tab2<-tab2[order(tab2$ID),]
  tab3<-tab3[order(row.names(tab3)),order(colnames(tab3))]
  result.tab<-cbind(tab1,tab2,tab3)
  result.tab<-result.tab[order(result.tab$PValue),]
  write.csv(result.tab,file=output)
  return(result.tab) 
}
#pvalue 0.01
lrt.output.piRNA <- function (lrt,counts,group,output) {
  npvalue<-length(which(lrt$table$PValue<0.01))
  tt<-topTags(lrt,n=npvalue)
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
#
lrt.output.tiRNA <- function (lrt,counts,group,output) {
  npvalue<-length(which(lrt$table$PValue<0.01))
  tt<-topTags(lrt,n=npvalue)
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





