# ---------------------------------------------------------------------
## PBMC8K+panT4k
load("p8k_t4k.RData") # 8381+4538 = 12919, genes 33694
bs1 = 8381;bs2 = 4538 #PBMC8K+panT4k
## PBMC8K+PBMC4K - not necessary because from the same donor
load("p8k_t4k.RData") # 8381+4538 = 12919, genes 33694
m1 = m[1:8381,] #pbmc8k
load("p4k_t3k_self.RData")
m2 = m_pbmc2[1:4340,] #pbmc4k
m = rbind(m1,m2)
bs1 = 8381;bs2 = 4340 #PBMC8K+PBMC4K

min.cell = 3
id_used = which(colSums(m)>=min.cell) #18084
m_seurat = .normalize(m,n=min.cell,nc=10000,log.space=T)$m
df_HV<-.get_variable_gene_poisson(m_seurat,n=1700) # must use log data
# For compare, use HV genes
set.seed(0)
tsne<-Rtsne(as.matrix(m_seurat[,df_HV$used]),pca=F,is.distance=F)
set.seed(0)
cl <- kmeans(as.matrix(m_seurat[,df_HV$used]),6,iter.max=100,algorithm="MacQueen")
save(tsne, cl, file="P8k+P4k_HV.RData")

# auto correcting
load("P8k+T4k_HV.RData")
m_pca = prcomp(m_seurat[,df_HV$used],center = T,scale. = T) # need use cor matrix
png("PC.png",width = 600, height = 480)
plot(m_pca, type = "l",npcs=30)
dev.off()
m_pca_data = predict(m_pca, newdata=m_seurat[,df_HV$used])[,1:20]

cl = .major_cluster(m_pca_data,bs=c(bs1,bs2),Nmajor=4)
rm_id = .batch_genes_detect(m_seurat,bs=c(bs1,bs2),cl=cl)
length(rm_id)
id_final = setdiff(which(df_HV$used),rm_id)
length(id_final)

# tsne and clustering
ttp <- as.matrix(m_seurat[,id_final])
set.seed(0)
tsne<-Rtsne(ttp,pca=F,is.distance=F)
#set.seed(0)
#cl6 <- kmeans(ttp,kpp_init(ttp, 6),iter.max=100,algorithm="MacQueen")
# better to clustering with PCA data
m_pca = prcomp(ttp,center = T,scale. = T) # need use cor matrix
m_pca_data = predict(m_pca, newdata=ttp)[,1:30]
set.seed(0)
#cl <- kmeans(m_pca_data,kpp_init(m_pca_data, 6),iter.max=100,algorithm="MacQueen")
cl10 <- kmeans(m_pca_data,10,iter.max=100,algorithm="MacQueen")
save(tsne,cl,cl10,file="P8k+T4k_auto.RData")

# ---------------------------------------------------------------------
# PBMC8K+PBMC4K
# ---------------------------------------------------------------------

## visualization
png(filename="batchP8k+T4k_auto.png")
.batch_plot(tsne,c(bs1,bs2))
dev.off()

png(filename="clusterP8k+T4k_auto_cl10.png",width = 800, height = 800)
.cluster_plot(tsne,cl10)
dev.off()

marker = c("CD3D", "IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7")
png("biomarker_P8k+T4k_auto.png",height=1200,width=1200)
.biomarker_plot(m[,marker],tsne$Y)
dev.off()
