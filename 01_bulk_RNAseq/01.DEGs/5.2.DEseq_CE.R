library(DESeq2)
mycounts<-read.csv("CE02_counts_for_DEseq2.txt",sep='\t')
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
condition<-factor(c(rep("tumbler",6),rep("others",12)), levels = c("tumbler","others"))
breed<-factor(c(rep("tumbler",6),rep("racing",6),rep("utility",6)), levels = c("tumbler","racing","utility"))
colData<-data.frame(row.names=colnames(mycounts), condition, breed)
mycounts <- mycounts[apply(mycounts,1, sum) >1, ]
dds <- DESeqDataSetFromMatrix(mycounts,colData, design= ~ condition)
dds$condition <- relevel(dds$condition, ref = "others")
dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "tumbler", "others"))
res = res[order(res$pvalue),]
summary(res)
write.csv(res,file="All_results_CE02.csv")
diff_gene_deseq2<-subset(res, padj <0.05 & abs(log2FoldChange)>1)
write.csv(diff_gene_deseq2, file="DEG_CE02_tumbler_vs_others.csv")

library(ggplot2)
vsd=vst(dds)
#rld <- rlog(dds)
pcaData <- plotPCA(vsd,intgroup=c("breed"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
width <- 8
height <- 8
p=ggplot(pcaData, aes(x=PC1, y=PC2, color=breed, label=name)) +
  geom_text(size=3,hjust=0.5, vjust=2) +
  geom_point(size=5) +
  xlab(paste("PC1:",percentVar[1],"% variance")) +
  ylab(paste("PC2:",percentVar[2],"% variance")) +
  stat_ellipse(level = 0.7)
#  coord_fixed()
ggsave(p, filename="PCA_CE.pdf", width=width, height=height)

sampleDists<-dist(t(assay(vsd)))
library(Cairo)
library(RColorBrewer)
library(pheatmap)
sampleDistsMatrix <- as.matrix(sampleDists)
rownames(sampleDistsMatrix) <- paste(vsd$breed)
colnames(sampleDistsMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
CairoPDF(file="heatmap_sampleToSample_CE.pdf")
pheatmap(sampleDistsMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off()
