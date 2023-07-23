install.packages("BiocManager") 
BiocManager::install("WGCNA") 
library(WGCNA)
library(doParallel)

data = read.table('skin_adjusted_filtered_vst.txt',sep="\t",header=T,row.names=1)
sample_info = read.table("skin_sample_info_filtered.txt",sep="\t", header=T, row.names=1) #/PBMC
data = data[,colnames(data) %in% sample_info$skin_id] 

rv = rowVars(as.matrix(data))
select_gene = rownames(data[order(rv, decreasing=T)[seq_len(10000)],])
data2 = data[rownames(data) %in% select_gene, ] 
data2 = as.data.frame(t(data2))

enableWGCNAThreads(nThreads = NULL)
cpucores=7
options("mc.cores"=cpucores)
cores = makeCluster(detectCores(),type='PSOCK')
registerDoParallel(cores)

#soft threshold; 
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(v2, powerVector = powers, verbose = 5)

png("skin_WGCNA_soft_threshold.png", width = 300, height = 300)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9)
abline(h=0.80,col="red")
dev.off()

png("skin_WGCNA_mean_connectivity.png", width = 300, height = 300)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1)
dev.off()

#Adjacency Matrix , Topological Overlap Matrix（TOM）
adjacency = adjacency(data2, power = 6) 
TOM = TOMsimilarity(adjacency,TOMType = "unsigned")
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, 
                             distM = dissTOM,
                             deepSplit = 2, 
                             pamStage = FALSE,
                             pamRespectsDendro = TRUE,
                             minClusterSize = 50) 
dynamicColors <- labels2colors(dynamicMods)
colorOrder <- c("grey", standardColors(max(dynamicMods)))
MEList <- moduleEigengenes(data2, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.1

merge <- mergeCloseModules(v2, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

png("skin_WGCNA_Dendro_merge.png", width = 500, height = 500)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors <- dynamicColors
names(dynamicMods) = colnames(v2)
mModu <- match(merge$colors, colorOrder) 
module = data.frame(gene=colnames(data2),dModu=dynamicMods,mModu=mModu)
write.table(module, file="skin_WGCNA_gene_cluster.txt",sep='\t',quote=F)


#module score
module_scoring = function(x){
  gene1 = subset(module,dModu==x)
  data2 = data[rownames(data) %in% gene1$gene,] 
  pca = prcomp(t(data2), scale=F) 
  return(pca$x[,1,drop=F])
}
score = NULL
for (i in 2:max(module$dModu)){
  score = cbind(score,module_scoring(i))  
}
score2 = scale(score) 
write.table(score, file="skin_WGCNA_module_score.txt",sep='\t',quote=F)


#GO
module$gene_sym = gsub("\\.","-", module$gene)
module$entrezid = mapIds(org.Hs.eg.db, keys=as.character(module$gene_sym), keytype = "SYMBOL", column="ENTREZID") #ENSEMBL
module_2 = subset(module, !is.na(module$entrezid))

res1 = NULL
for (i in 1:max(module_2$dModu)){ 
  gene1 = subset(module_2, dModu==i)
  GO1 = as.data.frame(enrichPathway(gene=gene1$entrezid, readable=T))
  if (is.na(GO1[1,1])==T){
    GO1 = as.data.frame(matrix(rep(NA),nrow=1,ncol=10))
    colnames(GO1) = c("module","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
    GO1[1,1]=i
    res1 = rbind(res1, GO1)
  } else
    res1 = rbind(res1, data.frame(module=i, GO1))
}
write.table(res1, file="skin_WGCNA_module_GO.txt",sep='\t',quote=F,col.names=NA)
