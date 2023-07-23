library("circlize")
library(ComplexHeatmap)
library(dplyr)

KEGG_data = read.table("KEGG_ligand_receptor_pair_inflammatory_list.txt",sep="\t",header=T)
sample_info = read.table("skin_pbmc_cross-sectional_sample_info.txt",sep="\t", header=T, row.names=1) 
skin_data = read.table("skin_cross-sectional_adjusted_filtered_vst.txt",sep="\t", header=T, row.names=1) 
pbmc_data = read.table("PBMC_cross-sectional_adjusted_filtered_vst.txt",sep="\t", header=T, row.names=1) 

count_link = function(x){
  patient1 = sample_info[x,] 
  skin_data2 = skin_data[, colnames(skin_data) %in% patient1$skin_id, drop=F]
  d2 = left_join(KEGG_data,transform(skin_data2, Cytokine=rownames(skin_data2)),by="Cytokine")
  d2 = left_join(d2,transform(skin_data2,Receptor=rownames(skin_data2)),by="Receptor")
  colnames(d2)[c(3,4)]=c("cytokine","receptor")#skin.cytokine - skin.receptor
  d2.1 = d2[!duplicated(d2$Receptor),]

  pbmc_data2 = pbmc_data[, colnames(pbmc_data) %in% patient1$pbmc_id, drop=F]
  d3 = left_join(KEGG_data,transform(pbmc_data2,Cytokine=rownames(pbmc_data2)),by="Cytokine")
  d3 = left_join(d3,transform(pbmc_data2,Receptor=rownames(pbmc_data2)),by="Receptor")
  colnames(d3)[c(3,4)]=c("cytokine","receptor")#pbmc.cytokine - pbmc.receptor
  d3$Cytokine = paste0("P_",d3$Cytokine)
  d3.1 = d3[!duplicated(d3$Receptor),]

  d2.3 = left_join(d2[,1:3],d3.1[,c("Receptor","receptor")],by="Receptor") #skin.cytokine - pbmc.receptor
  d2.3$Receptor = paste0("P_",d2.3$Receptor)
  d3.3 = left_join(d3[,1:3],d2.1[,c("Receptor","receptor")],by="Receptor") #pbmc.cytokine - skin.receptor
  d3$Receptor = paste0("P_",d3$Receptor)#pbmc.cytokine - pbmc.receptor
  
  link_all = rbind(d2,d3,d2.3,d3.3)
  link1 = subset(link_all, receptor>0&cytokine>0.5)
  link1$cytokine_type = ifelse(grepl("P_",link1$Cytokine)==T,"P","S")
  link1$receptor_type = ifelse(grepl("P_",link1$Receptor)==T,"P","S")
  link1$link_type = paste(link1$cytokine_type, link1$receptor_type, sep="_")
  return(cbind(patient1, 
               num_link=ifelse(is.na(link1[1,1])==T,0,nrow(link1)),
               S_S=ifelse(is.na(subset(link1,link_type=="S_S")[1,1])==T,0,nrow(subset(link1,link_type=="S_S"))), 
               S_P=ifelse(is.na(subset(link1,link_type=="S_P")[1,1])==T,0,nrow(subset(link1,link_type=="S_P"))),
               P_S=ifelse(is.na(subset(link1,link_type=="P_S")[1,1])==T,0,nrow(subset(link1,link_type=="P_S"))),
               P_P=ifelse(is.na(subset(link1,link_type=="P_P")[1,1])==T,0,nrow(subset(link1,link_type=="P_P")))))
}

res = NULL
for (i in 1:nrow(sample_info)){
  res = rbind(res, count_link(i))
}
write.table(res, file="cytokine_receptor_link_number.txt",sep='\t',quote=F,col.names=NA)


#statistics
library(lawstat)
res$disease = as.factor(res$disease)
colnames(res)[11:15]
#"num_link","S_S","S_P","P_S","P_P"
res.AD = subset(res,disease=="AD")
res.healthy = subset(res,disease=="Healthy")

link_num_stat = function(x){
  return(cbind(class=colnames(L2.1)[x], 
               BrunnerMunzel = brunner.munzel.test(res.AD[,x],res.healthy[,x])$p.value))
}
L3 = NULL
for (i in 11:15){
  L3 = rbind(L3,link_num_stat(i))
}
#write.table(L3, file="cytokine_receptor_link_number_stats.txt",sep='\t',quote=F,col.names=NA)


#Circos plot for individual patient 
heatmap_data1 = read.table('circos_cell_type_heatmap_reference_data.txt',sep="\t",header=T,row.names=1)
link0 = read.table('circos_cytokine_receptor_alignment_tissue_gene_order.txt',sep="\t",header=T,row.names=1)

split1 = factor(heatmap_data1$a, levels=c("skin","pbmc"))
heatmap_data2 = heatmap_data1[,-1]
link0 = data.frame(gene=rownames(heatmap_data1), tissue=heatmap_data1$a, od=c(1:185,1:133))

#link with order in circos plot
colnames(link0)[1]="Cytokine"
link2 = left_join(link1[,1:2],link0,by="Cytokine")
colnames(link2)[4] = "cytokine_od"
colnames(link0)[1]="Receptor"
link2 = left_join(link2,link0,by="Receptor")
colnames(link2)[6] = "receptor_od"
link2 = na.omit(link2)
link2 = link2[order(link2$tissue.x,link2$cytokine_od),]

KEGG_data2 = data.frame(gene=c(as.character(KEGG_data$Cytokine),as.character(KEGG_data$Receptor)),
                     type=c(rep("Cytokine",nrow(KEGG_data)),rep("Receptor",nrow(KEGG_data))))
KEGG_data3 = rbind(KEGG_data2, data.frame(gene=paste0("P_",KEGG_data2$gene), type=KEGG_data2$type))

cytokine_type = left_join(data.frame(od=1:nrow(heatmap_data2),gene=rownames(heatmap_data2)),KEGG_data3, by="gene")
cytokine_type$type2 = ifelse(cytokine_type$type=="Cytotkine",1,-1) 
cytokine_type = left_join(cytokine_type, 
                                 data.frame(gene=unique(c(as.character(link2$Cytokine),as.character(link2$Receptor))),exist=rep("yes")),by="gene")      
cytokine_type$type3 = ifelse(is.na(cytokine_type$exist)==TRUE,0,cytokine_type$type2)
rownames(cytokine_type) = cytokine_type$gene
cytokine_type2 = cytokine_type[,6,drop=F] 
link3 = subset(link2,tissue.x=="skin"&tissue.y=="pbmc")
link4 = subset(link2,tissue.x=="pbmc"&tissue.y=="skin")


circos.clear()
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1),start.degree = 90)
circos.initialize(split1, xlim=cbind(rep(10,2),c(length(split1[split1=="skin"])-5,length(split1[split1=="pbmc"]))))
circos.track(split1,ylim=c(0,1),track.height = 0.05)
draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "skin"),
            get.cell.meta.data("cell.end.degree", sector.index = "skin"),
            rou1 = get.cell.meta.data("cell.top.radius", track.index = 1),
            rou2 = get.cell.meta.data("cell.bottom.radius", track.index = 1), 
            col = "#CC99CC")
draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "pbmc"),
            get.cell.meta.data("cell.end.degree", sector.index = "pbmc"),
            rou1 = get.cell.meta.data("cell.top.radius", track.index = 1),
            rou2 = get.cell.meta.data("cell.bottom.radius", track.index = 1), 
            col = "red3")
circos.clear()
par(new = TRUE) 
circos.par("canvas.xlim" = c(-1.1, 1.1), "canvas.ylim" = c(-1.1, 1.1), gap.after = c(2.5, 10),start.degree = 85)
col1 = colorRamp2(0:20, c("gray90","darkgreen","purple4","red3","goldenrod1","salmon","navy","turquoise1","steelblue2","coral4","black",
                          "deeppink3","purple","plum4","deeppink","skyblue3","green1","green4","darkolivegreen4","palegreen1","olivedrab1"))
circos.heatmap(heatmap_data2, col = col1, split = split1, track.height = 0.3, cluster = FALSE,
               cell_width = rep(1, nrow(heatmap_data2)),
               bg.border="black", bg.lty = 1, bg.lwd = 1) #,cell.border = 1
col2 = colorRamp2(c(-1, 0, 1), c("green", "white", "orange"))
circos.heatmap(cytokine_type2, col = col2, track.height = 0.05)

for (x in seq_len(nrow(link2))){
  circos.link(link2[x,3],c(link2[x,4],link2[x,4]),link2[x,5],c(link2[x,6],link2[x,6]))
}

for (x in seq_len(nrow(link2))){
  circos.link(link3[x,3],c(link3[x,4],link3[x,4]),link3[x,5],c(link3[x,6],link3[x,6]), col="red")
}

for (x in seq_len(nrow(link2))){
  circos.link(link4[x,3],c(link4[x,4],link4[x,4]),link4[x,5],c(link4[x,6],link4[x,6]), col="blue")
}

dev.copy2pdf(file=paste0("cytokine_receptor_circos_plot_",patient1$sample_name,".pdf"), height=5, width=5) 
dev.off()


