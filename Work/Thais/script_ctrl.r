library(BSgenome.Rnorvegicus.UCSC.rn4)
readDNAStringSet("rat_srna.fasta","fasta")->rat.srna
runAnalysis1(rat.srna, outfile="rat_srna_hit.txt")

GRanges(seqnames=hits1[,1],ranges=IRanges(start=hits1[,2],end=hits1[,3]),strand=hits1[,4],ID=hits1[,5])->ctrl.sRNA
count.gene<-data$counts$rn4.gtf
rownames(count.gene) <- paste(elementMetadata(gene.anno)[,8],elementMetadata(gene.anno)[,9],elementMetadata(gene.anno)[,6],sep=";")
count.snRNA<-count.gene[grep("snRNA",row.names(count.gene)),]

library(edgeR)
  library(GenomicRanges)
  
  counts=counts[,which(group %in% c(conA,conB))]
  counts<-count.snRNA
  counts<-count.miRNA
  group=factor(group[which(group %in% c(conA,conB))])
  
  dge <- DGEList(counts, group=group)
  
  #filter, as least expressed in half the libraries
    keep<- rowSums(cpm(dge) > 10) >= (length(group)/2)
    dge <- dge[keep,]
   
  dge <- calcNormFactors(dge)
  
  dge <- estimateCommonDisp(dge)
     dge <- estimateTagwiseDisp(dge) 
#  dge <- estimateTrendedDisp(dge)
  
   et <- exactTest(dge, pair=c("ctrl","hfat"))
  
  apply(dge$pseudo.counts,1,sd)->sd.dge
  cbind(et$table,sd.dge)->var.dge
  summary(var.dge)
f0.sn.var[which(f0.sn.var[,2]>median(f0.sn.var[,2]) & f0.sn.var[,3]>median(f0.sn.var[,3]) & f0.sn.var[,4]<median(f0.sn.var[,4]) ),]
  var.dge->f0.sn.var
  

  var.dge->f0.mi.var
f0.mi.var[which(f0.mi.var[,2]>5 & f0.mi.var[,3]>0.7 & f0.mi.var[,4]<5 ),]
  
 dge <- DGEList(count.miRNA, group=grp)
dge <- DGEList(count.snRNA, group=grp)
keep<- rowSums(cpm(dge) > 10) >= (length(grp)/2)
dge <- dge[keep,]
dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit,contrast=c(-1,0,1,0))
apply(lrt$fitted.values,1,sd)->sd.dge
cbind(lrt$table,sd.dge)->var.dge
var.dge->f1.mi.var
var.dge->f1.sn.var

f1.mi.var[which(f1.mi.var[,2]>median(f1.mi.var[,2]) & f1.mi.var[,4]>median(f1.mi.var[,4]) & f1.mi.var[,5]<median(f1.mi.var[,5]) ),]

f1.mi.var[which(f1.mi.var[,2]> & f1.mi.var[,4]>0.6 & f1.mi.var[,5]<30 ),]

intersect(  row.names(    f1.mi.var[which(f1.mi.var[,2]>median(f1.mi.var[,2])
                                          & f1.mi.var[,5]<mean(f1.mi.var[,5])     
                                          ),] ),
  row.names(f0.mi.var[which(f0.mi.var[,2]>median(f0.mi.var[,2])
                            & f0.mi.var[,4]<median(f0.mi.var[,4])
                            
                            ),])) -> idlist

intersect(  row.names(    f1.sn.var[which(f1.sn.var[,2]>median(f1.sn.var[,2])
                                          & f1.sn.var[,5]<mean(f1.sn.var[,5])     
),] ),
            row.names(f0.sn.var[which(f0.sn.var[,2]>median(f0.sn.var[,2])
                                      & f0.sn.var[,4]<median(f0.sn.var[,4])
                                      
            ),])) -> idlist

writeHits <- function(seqname, matches, strand, file="", append=FALSE)
{
  if (file.exists(file) && !append)
    warning("existing file ", file, " will be overwritten with 'append=FALSE'")
  if (!file.exists(file) && append)
    warning("new file ", file, " will have no header with 'append=TRUE'")
  hits <- data.frame(seqname=rep.int(seqname, length(matches)),
                     start=start(matches),
                     end=end(matches),
                     strand=rep.int(strand, length(matches)),
                     patternID=names(matches),
                     check.names=FALSE)
  write.table(hits, file=file, append=append, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=!append)
}

runAnalysis1 <- function(dict0, outfile="")
{
  library(BSgenome.Rnorvegicus.UCSC.rn4)
  seqnames <- seqnames(Rnorvegicus)
  seqnames_in1string <- paste(seqnames, collapse=", ")
  cat("Target:", providerVersion(Rnorvegicus),
      "chromosomes", seqnames_in1string, "\n")
  append <- FALSE
  for (seqname in seqnames) {
    subject <- Rnorvegicus[[seqname]]
    cat(">>> Finding all hits in chromosome", seqname, "...\n")
    for (i in seq_len(length(dict0))) {
      patternID <- names(dict0)[i]
      pattern <- dict0[[i]]
      plus_matches <- matchPattern(pattern, subject)
      names(plus_matches) <- rep.int(patternID, length(plus_matches))
      writeHits(seqname, plus_matches, "+", file=outfile, append=append)
      append <- TRUE
      rcpattern <- reverseComplement(pattern)
      minus_matches <- matchPattern(rcpattern, subject)
      names(minus_matches) <- rep.int(patternID, length(minus_matches))
      writeHits(seqname, minus_matches, "-", file=outfile, append=append)
    }
    cat(">>> DONE\n")
  }
}
