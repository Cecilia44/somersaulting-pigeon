library(GenomicFeatures)
txdb<-makeTxDbFromGFF("/project-whj/sis/2.tumbler/2.count/Pigeon_scaffolds.FINAL.final_annotation.gff", format = "gff")
exon_gene <- exonsBy(txdb, by="gene")
exons_gene_lens <- lapply(exon_gene, function(x){sum(width(reduce(x)))})
class(exons_gene_lens)
length(exons_gene_lens)
exons_gene_lens <- as.data.frame(exons_gene_lens)
dim(exons_gene_lens)
#exons_gene_lens <- t(exons_gene_lens)
exons_gene_lens <- t(exons_gene_lens)
dim(exons_gene_lens)
write.table(exons_gene_lens, file = "exons_gene_lens", sep = " ")

exons_gene_lens <- read.table("exons_gene_lens", row.names = 1, sep=" ")
colnames(exons_gene_lens) = "length"
counts_CB02<-read.table("/project-whj/sis/2.tumbler/3.DEseq/CB02_counts_for_DEseq2.txt",header = T)
counts_CE02<-read.table("/project-whj/sis/2.tumbler/3.DEseq/CE02_counts_for_DEseq2.txt",header = T)
counts_OLL01<-read.table("/project-whj/sis/2.tumbler/3.DEseq/OLL01_counts_for_DEseq2.txt",header = T)
count_with_length_CB02 <- cbind(counts_CB02, exons_gene_lens)
count_with_length_CB02 <- count_with_length_CB02[,-1]
count_with_length_CE02 <- cbind(counts_CE02, exons_gene_lens)
count_with_length_CE02 <- count_with_length_CE02[,-1]
count_with_length_OLL01 <- cbind(counts_OLL01, exons_gene_lens)
count_with_length_OLL01 <- count_with_length_OLL01[,-1]
kb <- count_with_length_CB02$length/1000
countdata_CB02<-count_with_length_CB02[,1:18]
countdata_CE02<-count_with_length_CE02[,1:18]
countdata_OLL01<-count_with_length_OLL01[,1:18]
rpk_CB02<-countdata_CB02/kb
rpk_CE02<-countdata_CE02/kb
rpk_OLL01<-countdata_OLL01/kb
tpm_CB02 <- t(t(rpk_CB02)/colSums(rpk_CB02)*10^6)
tpm_CE02 <- t(t(rpk_CE02)/colSums(rpk_CE02)*10^6)
tpm_OLL01 <- t(t(rpk_OLL01)/colSums(rpk_OLL01)*10^6)
write.csv(tpm_CB02, file = "tpm_CB02.csv", quote = F)
write.csv(tpm_CE02, file = "tpm_CE02.csv", quote = F)
write.csv(tpm_OLL01, file = "tpm_OLL01.csv", quote = F)
