# *********************************************************************
# -------------------------------Part 1-for compare--------------------
# *********************************************************************
# ---------------------------------------------------------------------
# load pbmc_6k
# ---------------------------------------------------------------------
m6k = readMM("pbmc_6k/hg19/matrix.mtx") %>% t()  #5419 barcode and 32738 genes
barcode_6k = read.table("pbmc_6k/hg19/barcodes.tsv",sep='\t',stringsAsFactors=F)
gene_6k = read.table("pbmc_6k/hg19/genes.tsv",sep='\t',stringsAsFactors=F)
rownames(m6k) = barcode_6k[,1]
colnames(m6k) = gene_6k[,1]

# ---------------------------------------------------------------------
# load pbmc_4k
# ---------------------------------------------------------------------
m4k = readMM("pbmc_4k/GRCh38/matrix.mtx") %>% t()  #4340 barcode and 33694 genes
barcode_4k = read.table("pbmc_4k/GRCh38/barcodes.tsv",sep='\t',stringsAsFactors=F)
gene_4k = read.table("pbmc_4k/GRCh38/genes.tsv",sep='\t',stringsAsFactors=F)
rownames(m4k) = barcode_4k[,1]
colnames(m4k) = gene_4k[,1]

# ---------------------------------------------------------------------
# combine to pbmc_6k4k and save the data
# ---------------------------------------------------------------------
length(intersect(barcode_4k[,1], barcode_6k[,1])) # 0 that's good
common_genes = intersect(gene_4k[,1], gene_6k[,1]) # 31232 common genes
m6k4k = rbind(m6k[,common_genes],m4k[,common_genes])
rownames(gene_4k) = gene_4k[,1]
colnames(m6k4k) = gene_4k[common_genes,2]
save(m6k4k,file="pbmc6k4k.RData") # 9759 31232
load("pbmc6k4k.RData")

gene_4k[gene_4k[,1]=="ENSG00000244734",]
gene_4k[gene_4k[,1]=="ENSG00000206172",]


# ---------------------------------------------------------------------
# all preprocess - QC -> Normalization -> HV genes selection
# ---------------------------------------------------------------------
m1 = m6k4k[1:5419,]
m2 = m6k4k[5420:9759,]
rs = rowSums(m6k4k)
min.cell = 3
id_used = which(colSums(m6k4k)>=min.cell) #18084
id_used1 = which(colSums(m1)>=min.cell) #16707
id_used2 = which(colSums(m2)>=min.cell) #15656
m_seurat = .normalize(m6k4k,n=min.cell,nc=10000,log.space=T)$m
m_seurat1 = .normalize(m1,n=min.cell,nc=10000,log.space=T)$m
m_seurat2 = .normalize(m2,n=min.cell,nc=10000,log.space=T)$m
df_HV<-.get_variable_gene_poisson(m_seurat,n=1700) # must use log data
df_HVb1<-.get_variable_gene_poisson(m_seurat1,n=1500) # must use log data
df_HVb2<-.get_variable_gene_poisson(m_seurat2,n=1500) # must use log data
png("HV_poisson.png",height=600,width=600)
.plot_HV(df_HV)
dev.off()
# ---------------------------------------------------------------------
# tsne and clustering/PCA clustering
# ---------------------------------------------------------------------
set.seed(0)
tsne<-Rtsne(as.matrix(m_seurat[,df_HV$used]),pca=F,is.distance=F)
set.seed(0)
cl <- kmeans(as.matrix(m_seurat[,df_HV$used]),6,iter.max=100,algorithm="MacQueen")
#or pca clustring
m_pca = prcomp(m_seurat[,df_HV$used],center = T,scale. = T) # need use cor matrix
png("PC.png",width = 600, height = 480)
plot(m_pca, type = "l",npcs=50)
dev.off()
m_pca_data = predict(m_pca, newdata=m_seurat[,df_HV$used])[,1:30]
set.seed(0)
cl_pca <- kmeans(m_pca_data,6,iter.max=100)
save(tsne,cl,cl_pca,file="tsne_cl_clpca.RData")
load("tsne_cl_clpca.RData")
# ---------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------
#plot batch
png(filename="batch_p6k4k.png")
.batch_plot(tsne,c(5419,4340))
dev.off()
#plot clusting
png(filename="cluster_p6k4k_pca.png",width = 800, height = 800)
.cluster_plot(tsne,cl_pca)
dev.off()
#plot markers
marker = c("CD3D", "IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7","HBA1","HBA2","HBB")
#marker=c("IL4R","IL1A",        "IL1B",        "IL1R1",       "IL1R2",       "IL1RAP","IL1RL1",      "IL1RL2",      "IL1RN",       "IL2","IL10","IL10RA","IL10RB","IL10RB-AS1","IL11RA")
png("biomarker_p6k4k.png",height=1200,width=1200)
.biomarker_plot(m6k4k[,marker],tsne$Y)
dev.off()

# *********************************************************************
# ----------------------------Part 2-Auto correcting-------------------
# *********************************************************************
# ---------------------------------------------------------------------
# Auto remove batch effect
# ---------------------------------------------------------------------
#KS significant genes
p.ks1 = .batch_gene(as.matrix(m_seurat),as.matrix(m_seurat[1:5419,]),test="ks")
p.ks2 = .batch_gene(as.matrix(m_seurat),as.matrix(m_seurat[5420:9759,]),test="ks")
p.ks12 = .batch_gene(as.matrix(m_seurat[1:5419,]),as.matrix(m_seurat[5420:9759,]),test="ks")
#save(p.ks1,p.ks2,p.ks12,file="p6k4k_ks.RData")
load("p6k4k_ks.RData")

p.chisq = .batch_gene(as.matrix(m6k4k[1:5419,]),as.matrix(m6k4k[5420:9759,]),test="ks")
#save(p.chisq,file="Pchisq12.RData")
load("Pchisq12.RData")

id_p = which(-log(p.chisq)>100)
length(intersect(id_p,id_HV)) #158

pks12 = .batch_gene(m_seurat[1:5419,match(id_HV,id_used)],m_seurat[5420:9759,match(id_HV,id_used)],test="ks")
id2 = id_HV[which(pks12<1e-30)]
length(id2)


df_ks = data.frame(ks1=-log(p.ks1), ks2=-log(p.ks2))
#which(df_ks[,1]>20)
#which(df_ks[,2]>20)
id_ks = id_used[intersect(which(df_ks[,1]>20),which(df_ks[,2]>20))] #30,30: idks=860; 20,20: idks=1278; 15,15:1631 #10,10:2125
id_ks12 = id_used[which(-log(p.ks12)>30)] #30,30: idks=860; 20,20: idks=1278; 15,15:1631 #10,10:2125
length(id_ks)
length(id_ks12)
length(intersect(id_ks,id_HV)) #100
length(intersect(id_ks12,id_HV))#226



#retrieve genes
rank = data.frame(rank12=rep(30000,31232), rank1=rep(30000,31232), rank2=rep(30000,31232))
rank$rank12[id_used[order(df_HV$dispersion_norm,decreasing=T)]] = 1:length(id_used)
rank$rank1[id_used1[order(df_HVb1$dispersion_norm,decreasing=T)]] = 1:length(id_used1)
rank$rank2[id_used2[order(df_HVb2$dispersion_norm,decreasing=T)]] = 1:length(id_used2)
norm_rank = .norm_rank_diff(rank$rank1,rank$rank2)

id_retrieve = which(norm_rank>quantile(norm_rank, 0.99))
#id_retrieve = which(rank[,3]>10000&rank[,2]<2000| rank[,2]>10000&rank[,3]<2000)
intersect(id_retrieve, id_p)
#intersect(id_retrieve, id_HV)


#final id
id_HV = id_used[df_HV$used]
id_final = setdiff(id_HV, setdiff(id_p,id_retrieve))## 30,30:1634;  20,20: 1605
id_final = setdiff(id_HV,id_ks)## 20,20: 1600; 15,15: 1565 #10,10:1511
length(id_final)

#id1 = id_HV
#id2 = id_ks = id_used[union(which(df_ks[,1]<30),which(df_ks[,2]<30))]
#id3 = id_retrieve = which(norm_rank>quantile(norm_rank, 0.99))
#length(intersect(id1, union(id2,id3)))


#compare with BESG in p4k+t3k
load("id_rm_p4kt3k.RData")
id_rm_4k3k = id_rm #127
id_rm_6k4k = intersect(id_ks,id_HV) #100
common6k4k_4k3k_BESG_id = intersect(gene_4k[id_rm_4k3k,1],common_genes[id_rm_6k4k]) #40 common genes
save(common6k4k_4k3k_BESG_id,file="common6k4k_4k3k_BESG_id.RData")

#tsne
ttp = as.matrix(m_seurat[,match(id_final,id_used)])
set.seed(0)
tsne<-Rtsne(ttp,pca=F,is.distance=F)
set.seed(0)
cl <- kmeans(ttp,6,iter.max=100,algorithm="MacQueen")
#or pca clustring
m_pca = prcomp(ttp,center = T,scale. = T) # need use cor matrix
m_pca_data = predict(m_pca, newdata=ttp)[,1:30]
set.seed(0)
cl_pca <- kmeans(m_pca_data,6,iter.max=100)
save(tsne,cl,cl_pca,file="tsne_cl_clpca_6k4k_chisq.RData")

load("tsne_cl_clpca_final_15_onlyRM.RData")
# Visualization
#plot batch
png(filename="batch.png")
.batch_plot(tsne,c(5419,4340))
dev.off()
#plot clusting
png(filename="cluster_pca_2.png",width = 800, height = 800)
.cluster_plot(tsne,cl_pca)
dev.off()
png(filename="cluster.png",width = 800, height = 800)
.cluster_plot(tsne,cl)
dev.off()
#plot markers
marker = c("CD3D", "IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7","HBA1","HBA2","HBB")
#marker=c("IL4R","IL1A",        "IL1B",        "IL1R1",       "IL1R2",       "IL1RAP","IL1RL1",      "IL1RL2",      "IL1RN",       "IL2","IL10","IL10RA","IL10RB","IL10RB-AS1","IL11RA")
png("biomarker.png",height=1200,width=1200)
.biomarker_plot(m6k4k[,marker],tsne$Y)
dev.off()

# *********************************************************************
# -----------------------Part 3-semi-Auto-----------------------------
# *********************************************************************
# ---------------------------------------------------------------------
# manually remove batch effect
# ---------------------------------------------------------------------
# cl_pca: 1,3,4-T; 5-M; 6-B, 2-NK
summary(as.factor(cl_pca$cluster)) #1:1805, 2:478, 3:2430, 4:1054, 5:2637, 6:1355
length(id_used)
id_t = which(cl_pca$cluster %in% c(1,3,4))
id_t1 = id_t[id_t<=5419];id_t2 = id_t[id_t>5419] #2958+2331=5289
id_m = which(cl_pca$cluster ==5)
id_m1 = id_m[id_m<=5419];id_m2 = id_m[id_m>5419] #2637+1459=1178
id_b = which(cl_pca$cluster == 6)
id_b1 = id_b[id_b<=5419];id_b2 = id_b[id_b>5419] #648+707=1355
id_nk = which(cl_pca$cluster ==2)
id_nk1 = id_nk[id_nk<=5419];id_nk2 = id_nk[id_nk>5419] #183+295=478
length(id_nk)

p_t = .batch_gene(m_seurat[id_t1,],m_seurat[id_t2,],test="ks")
p_m = .batch_gene(m_seurat[id_m1,],m_seurat[id_m2,],test="ks")
p_b = .batch_gene(m_seurat[id_b1,],m_seurat[id_b2,],test="ks")
p_nk = .batch_gene(m_seurat[id_nk1,],m_seurat[id_nk2,],test="ks")
df_ks = data.frame(pt=-log(p_t),pm=-log(p_m),pb=-log(p_b),pnk=-log(p_nk))
save(df_ks, file="pks6k4k_sub.RData")
load("pks6k4k_sub.RData")

sum(df_ks$pnk>30)
id_rm3 = id_used[with(df_ks,pt>50|pm>50|pb>50|pnk>50)] #30:1987(t:1553,m:1245,b:466,nk:80), 50:1676
length(id_rm3)

id_final3 = setdiff(id_used[df_HV$used],id_rm3) #30:1495, 50:1533
length(id_final3)


ttp = as.matrix(m_seurat[,match(id_final3,id_used)])
set.seed(0)
tsne<-Rtsne(ttp,pca=F,is.distance=F)
set.seed(0)
cl <- kmeans(ttp,6,iter.max=100,algorithm="MacQueen")
#or pca clustring
m_pca = prcomp(ttp,center = T,scale. = T) # need use cor matrix
m_pca_data = predict(m_pca, newdata=ttp)[,1:30]
set.seed(0)
cl_pca <- kmeans(m_pca_data,6,iter.max=100)
save(tsne,cl,cl_pca,file="tsne_cl_clpca_part3_30.RData")

### Visualization
#plot batch
png(filename="batch_part3_30.png")
.batch_plot(tsne,c(5419,4340))
dev.off()
#plot clusting
png(filename="cluster_part3_30.png",width = 800, height = 800)
.cluster_plot(tsne,cl_pca)
dev.off()
#plot markers
marker = c("CD3D", "IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7","HBA1","HBA2","HBB")
#marker=c("IL4R","IL1A",        "IL1B",        "IL1R1",       "IL1R2",       "IL1RAP","IL1RL1",      "IL1RL2",      "IL1RN",       "IL2","IL10","IL10RA","IL10RB","IL10RB-AS1","IL11RA")
png("biomarker_part3_30.png",height=1200,width=1200)
.biomarker_plot(m6k4k[,marker],tsne$Y)
dev.off()

# *********************************************************************
# -----------------------Part 4-MNN ------------------------------
# *********************************************************************
# (on app.01)
require("scran")
load("pbmc6k4k.RData")
min.cell = 3
id_used = which(Matrix::colSums(m6k4k)>=min.cell) #18084
m_seurat = .normalize(m6k4k,n=min.cell,nc=10000,log.space=T)$m
df_HV<-.get_variable_gene_poisson(m_seurat,n=1700) # must use log data
result_mnn <- mnnCorrect(as.matrix(t(m_seurat[1:5419,])),as.matrix(t(m_seurat[5420:9759,])))
save(result_mnn,file="p6k3k_mnn_all_genes.RData")
result_mnn_hv <- mnnCorrect(as.matrix(t(m_seurat[1:5419,])),as.matrix(t(m_seurat[5420:9759,])),hvg.genes=which(df_HV$used))
save(result_mnn_hv,file="p6k3k_mnn_hv_genes.RData")

# back to mario
load("p6k3k_mnn_all_genes.RData")
load("p6k3k_mnn_hv_genes.RData")
mnn1 = result_mnn_hv$correct[[1]]
mnn2 = result_mnn_hv$correct[[2]]
mnn = t(cbind(mnn1,mnn2))
#gene_used = which(Matrix::colSums(mnn)>=min.cell)
df_pbmc<-.get_variable_gene_poisson(mnn)
set.seed(0)
tsne_pbmc<-Rtsne(mnn[,df_pbmc$used],pca=F,is.distance=F)
set.seed(0)
cl <- kmeans(mnn[,df_pbmc$used],6,iter.max=100,algorithm="MacQueen")
m_pca = prcomp(mnn[,df_pbmc$used],center = T,scale. = T) # need use cor matrix
png("PC.png",width = 600, height = 480)
plot(m_pca, type = "l",npcs=50)
dev.off()
m_pca_data = predict(m_pca, newdata=mnn[,df_pbmc$used])[,1:30]
set.seed(0)
cl_pca <- kmeans(m_pca_data,6,iter.max=100)
save(tsne,cl,cl_pca,file="tsne_cl_clpca_mnn_hv.RData")

### Visualization
#plot batch
png(filename="batch_part4_hv_genes.png")
.batch_plot(tsne,c(5419,4340))
dev.off()
#plot clusting
png(filename="cluster_part4_hv_genes.png",width = 800, height = 800)
.cluster_plot(tsne,cl)
dev.off()
#plot markers
marker = c("CD3D", "IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7","HBA1","HBA2","HBB")
#marker=c("IL4R","IL1A",        "IL1B",        "IL1R1",       "IL1R2",       "IL1RAP","IL1RL1",      "IL1RL2",      "IL1RN",       "IL2","IL10","IL10RA","IL10RB","IL10RB-AS1","IL11RA")
png("biomarker_part4_hv_genes.png",height=1200,width=1200)
.biomarker_plot(m6k4k[,marker],tsne$Y)
dev.off()
