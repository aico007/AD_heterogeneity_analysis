library(DESeq2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db) 
library(ReactomePA) 

data = read.table('skin_count_data_adjusted_filtered.txt',sep="\t",header=T,row.names=1) #/PBMC
sample_info = read.table("skin_cross-sectional_sample_info.txt",sep="\t", header=T, row.names=1) #/PBMC
data = as.matrix(data[, colnames(data) %in% sample_info$sample_name])

dds <- DESeqDataSetFromMatrix(countData = data, colData=sample_info, design= ~ disease)
dds <- DESeq(dds)
res <- results(dds,contrast=c("disease","AD","healthy"))
res1 <- data.frame(gene=rownames(data),logFC=res$log2FoldChange, FDR=res$padj)
write.table(res1, file="skin_logFC_AD_vs_healthy.txt",sep="\t",quote=F)

#skin
DEG_skin = subset(res1, FDR<0.01&abs(logFC)>2)
#PBMC
DEG_PBMC = subset(res1, FDR<0.05&abs(logFC)>1)

#Volcano plot
res1 = subset(res1,!is.na(res1$FDR))
res1$FDR_log10 = -log10(res1$FDR)
x=1
y=0.05
res1$thcolor = ifelse(res1$FDR >= y , 0,
                      ifelse(res1$logFC >=x, 1,
                             ifelse(res1$logFC <=-x, 2, 0)))
res2 = subset(res1, FDR<y&abs(logFC)>x|FDR<y&logFC<=(-x))

g = ggplot() +
  geom_point(aes(x=res1$logFC, y=-log10(res1$FDR), color=as.factor(res1$thcolor)), size = 0.3) +
  labs(x="",y="")+
  scale_color_manual(values = c("black", "red2", "blue2")) +
  geom_hline(yintercept=-log10(y), color = "darkgreen",size=0.2) + 
  geom_vline(xintercept=c(-x,x), color = "darkgreen",size=0.2) +
  scale_x_continuous(breaks = seq(-10, 10, by=2), limits=c(-11.1,11.1))+ #for skin
  scale_y_continuous(breaks = seq(0, 50, by=10), limits=c(0,47))+ #for skin
  theme_bw(base_size = 8) +
  theme(panel.background = element_blank(),panel.grid=element_blank())+ #(size=0.5,colour="black",fill="white"))+
  theme(legend.position ="none") +
  theme(axis.text.x = element_text(size = 8, colour = "black"),axis.text.y = element_text(size = 8, colour = "black"))
show(g)
ggsave("skin_AD_vs_healthy_DEG_volcano_plot.pdf",dpi=300, height=4,width=4) 


#GSEA
gene1 = subset(res1, FDR<0.05)
gene1$entrezid = mapIds(org.Hs.eg.db, keys=as.character(gene1$gene), keytype = "SYMBOL", column="ENTREZID") 
gene1 = subset(gene1, !is.na(gene1$entrezid)&!is.na(gene1$logFC))
gene2 = gene1$logFC
names(gene2) = gene1$entrezid
gene2 = gene2[rev(order(gene2))] 
gsea1 = gsePathway(gene2, organism="human", nPerm=10000, pvalueCutoff=0.1, pAdjustMethod="BH", verbose=FALSE)
gsea1 = as.data.frame(gsea1)
gsea2 = subset(gsea1, p.adjust<0.05&enrichmentScore>0.5) 
gsea2 = gsea2[order(gsea2$enrichmentScore,decreasing=T),]
write.table(gsea2, file="Skin_DEG_AD_vs_Healthy_GSEA.txt",sep="\t",quote=F, row.names=F) 

g = ggplot(gsea2, aes(x = Description, y = as.numeric(enrichmentScore))) 
g = g + geom_bar(stat = "identity", position = "dodge",fill="darkgreen")
g = g + coord_flip()
g = g+theme_classic() #bw(base_size = 8) #+theme(panel.background = element_rect(size=1,colour="black",fill="white"))
g = g + theme(axis.text.y = element_text(size = 12, angle = 0, hjust = 0, colour = "black"),axis.text.x = element_text(size=12, angle = 0, hjust = 0, colour = "black"),
               axis.title.x = element_blank(),axis.title.y = element_blank(),plot.margin = unit(c(0.1,1.2,0.1,0.1), "cm"))
g = g+ylim(0,0.85)
plot(g)
ggsave(g, file = "Skin_DEG_GSEA_bargraph.pdf",dpi=300,height=4,width=4)


