
#
#set dir
setwd("../Thais/")
load("F0/countData.RData")
data.F0<-data

load("F1/countData.RData")
data.F1<-data

#merge two data
data$reads<-c(data.F0$reads,data.F1$reads)
data$counts$rn4.gtf<-cbind(data.F0$counts$rn4.gtf,data.F1$counts$rn4.gtf)
data$counts$rno2.gff3<-cbind(data.F0$counts$rno2.gff3,data.F1$counts$rno2.gff3)

#load phenodata, check sample order
pd0<-read.table(file="F0/expSpec.txt",sep="\t",header=T,stringsAsFactors=F)
pd1<-read.table(file="F1/expSpec",sep="\t",header=T,stringsAsFactors=F)
#check pd order and change
pd0<-pd0[c(1,8,2:7,9,18,10:17),]
pd1[,1]#pd<-read.table(file="result/2013-10-18/expSpec",header=T,stringsAsFactors=F)
pd<-rbind(pd0,pd1)
#check 'pd' as sample order
library(edgeR)
count.miRNA<-data$counts$rno2.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA"),]
count.pri.miRNA<-data$counts$rno2.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA_primary_transcript"),]

count.miRNA.F1<-data.F1$counts$rno2.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA"),]
count.pri.miRNA.F1<-data.F1$counts$rno2.gff3[which(elementMetadata(miRNA.anno)[,2]=="miRNA_primary_transcript"),]

#hist(rowSums(count.miRNA),xlim=c(0,49),1000000)
#hist(rowSums(count.pri.miRNA),xlim=c(0,49),2000000)

grp <- factor(pd[,4])


# Two pair compair
# F0 
f0.miRNA<-de.edgeR(count.miRNA,grp,"F0_miRNA_hfat_vs_ctrl_Trend_Filter.csv","ctrl","hfat",F,T)

f1_PHC_vs_PCC.miRNA<-de.edgeR(count.miRNA,grp,"F1_miRNA_PHC_vs_PCC_Trend_Filter.csv","PCC","PHC",F,T)

# mix model compair
#build target
#target <- pd[,c(1,2,4)]
#colnames(target)[2:3]<-c("filial","Diet")
#target[c(1:18),2] <- "F0"
#target[c(19:56),2] <- "F1"
#target[which(target$Diet %in% c("PCC","PCH")),2] <- "F1C"
#target[which(target$Diet %in% c("PHC","PHH")),2] <- "F1H"
#target[c(1:8),3] <- "Chow"
#target[c(9:18),3] <- "HFD"
#target[which(target$Diet %in% c("PCC","PHC")),3] <- "Chow"
#target[which(target$Diet %in% c("PCH","PHH")),3] <- "HFD"
#grp2 <- factor(paste(target$filial,target$Diet,sep="."))
grp1 <- factor(pd1[,4])
dge1 <- DGEList(count.miRNA.F1, group=grp1)
design <- model.matrix(~0+grp1)
colnames(design)<-levels(grp1)
filter<-T

if (filter==T) {
  keep<- rowSums(cpm(dge1) > 1) >= (length(grp1)/2)
  dge1 <- dge1[keep,]
}
#mds plot
plotMDS(dge1,col=as.numeric(grp1),label=grp1)

dge1 <- calcNormFactors(dge1)
dge1 <- estimateGLMCommonDisp(dge1, design)
dge1 <- estimateGLMTrendedDisp(dge1, design)
dge1 <- estimateGLMTagwiseDisp(dge1, design)
fit <- glmFit(dge1, design)
#lrt <- glmLRT(fit)
lrt <- glmLRT(fit,contrast=c(-1,0,1,0))

lrt.output <- function (lrt,counts,group,output) {
  npvalue<-length(which(lrt$table$PValue<0.05))
  tt<-topTags(lrt,n=npvalue)
  
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

F1_GLM_PHC_vs_PCC.miRNA<-lrt.output(glmLRT(fit,contrast=c(-1,0,1,0)),count.miRNA.F1,grp1,"F1_miRNA_PHC_vs_PCC_GLM.csv")
F1_GLM_PCH_vs_PCC.miRNA<-lrt.output(glmLRT(fit,contrast=c(-1,1,0,0)),count.miRNA.F1,grp1,"F1_miRNA_PCH_vs_PCC_GLM.csv")
F1_GLM_PHH_vs_PHC.miRNA<-lrt.output(glmLRT(fit,contrast=c(0,0,-1,1)),count.miRNA.F1,grp1,"F1_miRNA_PHH_vs_PHC_GLM.csv")
F1_GLM_PHH_vs_PCH.miRNA<-lrt.output(glmLRT(fit,contrast=c(0,-1,0,1)),count.miRNA.F1,grp1,"F1_miRNA_PHH_vs_PCH_GLM.csv")

###just for  pictures!!!
f0.lars<-read.csv("F0/miRNA_finalized/significantmirnas/sigCountsMirBase.csv")
rownames(f0.lars)[11] <- "rno-miR-128-3p"
#compare based on IDs
all.miRNA.ID <- unique(c(f0.miRNA[,6], F1_GLM_PHC_vs_PCC.miRNA[,7]))
a1<- ( all.miRNA.ID %in% f0.miRNA[,6] )
# a2<-  ( all.miRNA.ID %in% rownames(f0.lars) )
a3<- ( all.miRNA.ID %in% F1_GLM_PHC_vs_PCC.miRNA[,7] )
a3<- ( all.miRNA.ID %in% F1_GLM_PCH_vs_PCC.miRNA[,7] )
a3<- ( all.miRNA.ID %in% F1_GLM_PHH_vs_PHC.miRNA[,7] )
a3<- ( all.miRNA.ID %in% F1_GLM_PHH_vs_PCH.miRNA[,7] )
aa<- cbind(a1,a3)
aa <- vennCounts(aa)  
vennDiagram(aa,names=c("F0","F1"))
library(venneuler)
plot(venneuler(c(F0=15, F1=83, "F0&F1"=3)))

#build f0.dge, which is the one produced in de.edgeR
f0.dge$pseudo.counts[rownames(f0.miRNA),]->ht.f0
sapply(strsplit(colnames(ht.f0),"[.]"),"[[",1)->colnames(ht.f0)
ht.f0<-as.matrix(ht.f0)

heatmap(ht.f0,Colv=NA,labRow=NA,ylab="miRNA",xlab="sample")
heatmap(ht.f0[intersect(rownames(f0.miRNA), rownames(F1_GLM_PHC_vs_PCC.miRNA)),],Colv=NA,labRow=NA,ylab="miRNA",xlab="sample")
heatmap(ht.f0[setdiff(rownames(f0.miRNA), rownames(F1_GLM_PHC_vs_PCC.miRNA)),],Colv=NA,labRow=NA,ylab="miRNA",xlab="sample")
heatmap(as.matrix(pirna[,8:25]),Colv=NA,labRow=NA,ylab="piRNA",xlab="sample")

save.image("test.RData")