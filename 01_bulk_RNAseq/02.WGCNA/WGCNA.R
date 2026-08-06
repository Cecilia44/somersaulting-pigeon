setwd("/project-whj/sis/2.tumbler/5.WGCNA_tpm/3.2_3brainRegion")
library(WGCNA)
options(stringsAsFactors = FALSE)

femData = read.csv("../tpm_3brainRegion.csv")
dim(femData)
datExpr = as.data.frame(t(femData))
names(datExpr) = femData$X
datExpr <- datExpr[-1, ]
datExpr[] <- lapply(datExpr, function(x) as.numeric(as.character(x)))
datExpr0 <- datExpr[, !apply(datExpr, 2, function(x){sum(x==0)>1}),]
dim(datExpr0)
datExpr1 <- datExpr0[, !apply(datExpr0, 2, function(x) { sum(x == 0) > 1 })]
count_less_than_1 <- apply(datExpr1, 2, function(x) sum(x < 1))
cols_to_keep <- count_less_than_1 < 1
datExpr1 <- datExpr1[, cols_to_keep]



gsg = goodSamplesGenes(datExpr1, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr1 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr1), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf("sampleTree.pdf", width = 8, height = 6, pointsize = 12)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

pdf("sampleTree_remove_outliers.pdf", width = 8, height = 6, pointsize = 12)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 16000, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
dev.off()
table(clust)   
keepSamples = (clust==1)
datExpr1 = datExpr1[keepSamples, ] 
nGenes = ncol(datExpr1)
nSamples = nrow(datExpr1)

traitData <- read.csv("trait.csv")

Samples = rownames(datExpr1)
traitRows = match(Samples, traitData$ID)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()

pdf("sampleTree2.pdf", width = 8, height = 6, pointsize = 12)
sampleTree2 = hclust(dist(datExpr1), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(datExpr1, datTraits, file = "tumbler.RData")

getwd()
workingDir = "."
setwd(workingDir)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads() 
lnames = load(file = "tumbler.RData")
lnames

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr1, powerVector = powers, verbose = 5)
pdf("sft.pdf", width = 12, height = 8, pointsize = 12)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red") 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft$powerEstimate

net = blockwiseModules(datExpr1, power = 12,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "tumblerTOM",
                       verbose = 3)
table(net$colors)

pdf("cluster_dendrograms.pdf", width = 12, height = 9, pointsize = 12)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "tumbler-01-networkConstruction-auto.RData")

nGenes = ncol(datExpr1);
nSamples = nrow(datExpr1);
MEs0 = moduleEigengenes(datExpr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
pdf("labeledHeatmap.pdf", width = 10, height = 8)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

tumbling = as.data.frame(datTraits$phenotype)
names(tumbling) = "tumbling"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr1,tumbling, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(tumbling), sep="");
names(GSPvalue) = paste("p.GS.", names(tumbling), sep="");

module = "darkred"
column = match(module, modNames);
moduleGenes = moduleColors==module;
pdf("MM-GS_darkred.pdf", width = 7, height = 7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for tumbling",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

names(datExpr1)
names(datExpr1)[moduleColors=="darkred"]
annot = read.csv(file = "Pigeon_scaffolds.FINAL.final_annotation.gene.combine_and_delet.sort.list.csv")
dim(annot)
names(annot)
probes = names(datExpr1)
probes2annot = match(probes, annot$geneID)
sum(is.na(probes2annot))
geneInfo0 = data.frame(geneID = probes,
                       geneSymbol = annot$symbol[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
modOrder = order(-abs(cor(MEs, tumbling, use = "p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.tumbling))
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")

annot = read.csv(file = "Pigeon_scaffolds.FINAL.final_annotation.gene.combine_and_delet.sort.list.csv")
probes = names(datExpr1)
probes2annot = match(probes, annot$geneID)
allLLIDs = annot$symbol[probes2annot]
intModules = c("darkred")
for (module in intModules)
{
  modGenes = (moduleColors==module);
  modSymbols = allLLIDs[modGenes];
  fileName = paste("GeneSymbol-", module, ".txt", sep="");
  write.table(as.data.frame(modSymbols), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

