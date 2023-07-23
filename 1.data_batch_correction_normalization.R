if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("sva")
#or
#devtools::install_github("zhangyuqing/sva-devel",force = TRUE)

library(sva)
library(DESeq2)
library(matrixStats)
library(dplyr)

#For skin data and PBMC data respectively
sample_info = read.table("skin_sample_info.txt",sep="\t", header=T, row.names=1) #/PBMC
count_data = read.table("skin_count_data.txt",sep="\t", header=T, row.names=1) #/PBMC

#batch correction
sample_info$batch = as.factor(sample_info$batch)
adjusted = ComBat_seq(count_data, batch=sample_info$batch, group=NULL, full_mod=FALSE)

#filtering
adjusted = subset(adjusted, !(grepl("RPS",rownames(adjusted))|grepl("RPL",rownames(adjusted))))
adjusted = subset(adjusted, !(rownames(adjusted)=="HBA1"|rownames(adjusted)=="HBA2"|rownames(adjusted)=="HBB"))

cutoff1 = function(num){
  if(num >= 1){
    return(1)
  }else return(0)
}
adjusted_logi <- apply(adjusted, c(1,2), cutoff1)
filter1 = rownames(adjusted_logi[rowSums(adjusted_logi)>ncol(adjusted_logi)*0.05, ]) #genes that are expressed in more than 5% of patients 
filter2 = rownames(adjusted[rowMaxs(adjusted)>8, ])
filter3 = intersect(filter1, filter2)
data_filtered = adjusted[rownames(adjusted) %in% filter3,]
write.table(data_filtered, file="skin_count_data_adjusted_filtered.txt",sep='\t',quote=F) #/PBMC

# vst
dds = DESeqDataSetFromMatrix(countData = data_filtered, colData=sample_info, design= ~1)
dds = DESeq(dds)
vsd = vst(dds, blind=FALSE)
mat1 = assay(vsd)
mat2 = mat1 - rowMedians(mat1)
write.table(v1,file="skin_adjusted_filtered_vst.txt",sep="\t",quote=F) #/PBMC



