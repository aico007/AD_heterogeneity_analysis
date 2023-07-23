#feature data generated with tsfresh (python module)
data = read.table('longitudinal_set_timeseries_extract_feature.csv',sep=",",header=T)
data$feature=sapply(strsplit(rownames(data),"\\__"), "[", 2)
data = subset(data, feature %in% c("maximum","mean","minimum","root_mean_square",
                                  "approximate_entropy__m_2__r_0.7","mean_abs_change","cid_ce__normalize_False"))
dadta = t(scale(t(data[,-31])))

colnames(data)[1:2]
#"term" "feature"
hist(as.matrix(data[,3:32]))
cutoff2 <- function(num){
  if(is.na(num)){
    return(NA)
  }else if(num > 2){
    return(2)
  }else if(num < -2){
    return(-2)
  }else return(num)
}
data[,3:32] = apply(as.matrix(data[,3:32]) ,c(1,2), cutoff2) 
data2 = subset(data,term!="EASI")
data2 = as.matrix(data2[,3:32])


#k-means clustering
library(factoextra)
#identification of optimal k by Silhouette method
pdf("longitudinal_feature_kmeans_Silhouette_plot.pdf", width = 4, height = 4)
fviz_nbclust(t(data2), kmeans, method = "silhouette", k.max = 10, nboot = 100,iter.max=30)+
  labs(subtitle = "Silhouette method")
dev.off()

set.seed(40)
data.km = kmeans(t(data2),3)
data.km2 = data.frame(patient_id=names(data.km$cluster),cluster=data.km$cluster)

#pca
pca = prcomp(t(data2), scale=F) 
plot = as.data.frame(pca$x)[,1:4]
plot = left_join(transform(plot,patient_id=rownames(plot)), data.km2, by="patient_id")
plot = left_join(plot, transform(t(data2),patient_id=rownames(t(data2))), by="patient_id")

g = ggplot(data=plot2,aes(x=PC1,y=PC2))
g = g+geom_point(aes(x=PC1,y=PC2), alpha=0.8, size=2.8) + 
  scale_fill_manual(values=c("#F8766D", "#00BA38","#619CFF")) + 
  scale_shape_manual(values = c(24, 23, 21)) 
g = g+theme_bw(base_size = 8) +
  theme(panel.background = element_rect(size=0.5,colour="black",fill="white"),axis.title = element_text(size = 8)) 
ggsave("longitudinal_feature_PCA_plot.pdf",width=2.4,height=2.4) 


#multiple comparison test 
library(lawstat)
cluster1 = subset(plot, cluster==1)
cluster2 = subset(plot, cluster==2)
cluster3 = subset(plot, cluster==3)

L1=NULL
for (i in 7:ncol(plot)){ 
  L1=rbind(L1, data.frame(term=colnames(plot)[i], 
                    p.value.kruskal=kruskal.test(x=list(cluster1[,i],cluster2[,i],cluster3[,i]))$p.value,
                    p.brunner.C1vsC2=brunner.munzel.test(cluster1[,i],cluster2[,i])$p.value,
                    p.brunner.C1vsC3=brunner.munzel.test(cluster1[,i],cluster3[,i])$p.value,
                    p.brunner.C2vsC3=brunner.munzel.test(cluster2[,i],cluster3[,i])$p.value))
}

L2 = as.matrix(L1[,3:5])
rownames(L2) = L1$term
for (i in 1:nrow(L2)){
  L2[i,]=p.adjust(L2[i,], method='holm')
}
write.table(L2, file="longitudinal_timeseries_feature_cluster_statistics_p.value.adj.holm.txt",sep='\t',quote=F,col.names=NA)
