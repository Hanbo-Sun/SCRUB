###/////////////////////////////////////////////////////////////////
# This R file "BECC_p4kt3k.R" is to combine various pvalue/distance,
# by which to to identify BESGs as accurate as possible.
###/////////////////////////////////////////////////////////////////
# ---------------------------------------------------------------------
# Prepare: 1.Obtain rm genes; 2.load p4k+t3k; 3.normalize & HV
# ---------------------------------------------------------------------
#load rm genes index
load("id_rm_p4kt3k.RData")
#load p4k+t3k
load("p4k_t3k_self.RData")
m = m_pbmc2
id_b1 = 1:4340;id_b2 = 4341:7895
min.cell = 3
id_used = which(colSums(m)>=min.cell) #18084
#normalize and HV
m_seurat = .normalize(m,n=min.cell,nc=10000,log.space=T)$m
df_HV<-.get_variable_gene_poisson(m_seurat[,id_used],n=1700) # must use log data
id_HV = id_used[df_HV$used]


# ---------------------------------------------------------------------
# Statistics: 1.ks 2.Bhattacharyya 3.Mann-Whitney
# ---------------------------------------------------------------------
#pv - ks exclude 0 (exclude 0 is better)
ks = .batch_gene(m_seurat[id_b1,df_HV$used],m_seurat[id_b2,df_HV$used],test="ks",zero="ignore")
id_ks = id_HV[which(ks<1e-6)]
length(id_ks)
intersect(id_ks, id_rm) #118|250
id_final = setdiff(id_HV, id_ks)
length(id_final)

#bc - Bhattacharyya dist, include 0 is better
bc = .batch_gene(m[id_b1,id_HV],m[id_b2,id_HV],test="bc",zero="use")
summary(bc)
id_bc = id_HV[which(bc>0.14)]
length(id_bc)
intersect(id_bc, id_rm) #116|369
id_final = setdiff(id_HV, id_bc)
length(id_final)

#Mann-Whitney’ test, exclude 0 is better
mw = .batch_gene(m[id_b1,id_HV],m[id_b2,id_HV],test="mw",zero="ignore")
summary(mw)
id_mw = id_HV[which(mw<1e-5)]
length(id_mw)
intersect(id_mw, id_rm) #115|280
id_final = setdiff(id_HV, id_mw)
length(id_final)

# combine - vote
# total:118/250,116/369,115/280,unique=464/899; intersect:106/183
# vote majority
tb = table(c(id_ks,id_bc,id_mw))
idd = names(tb[tb>=2]) %>% as.numeric()
intersect(idd, id_rm) #117/247
id_final = setdiff(id_HV,idd)


# ---------------------------------------------------------------------
# rank
# ---------------------------------------------------------------------
m1 = m[id_b1,];m2 = m[id_b2,]
id_used1 = which(colSums(m1)>=min.cell) #16707
id_used2 = which(colSums(m2)>=min.cell) #15656
m_seurat1 = .normalize(m1,n=min.cell,nc=10000,log.space=T)$m
m_seurat2 = .normalize(m2,n=min.cell,nc=10000,log.space=T)$m
df_HVb1<-.get_variable_gene_poisson(m_seurat1,n=1700) # must use log data
df_HVb2<-.get_variable_gene_poisson(m_seurat2,n=1700) # must use log data
#rank
rank = data.frame(rank12=rep(30000,33694), rank1=rep(30000,33694), rank2=rep(30000,33694))
rank$rank12[id_used[order(df_HV$dispersion_norm,decreasing=T)]] = 1:length(id_used)
rank$rank1[id_used1[order(df_HVb1$dispersion_norm,decreasing=T)]] = 1:length(id_used1)
rank$rank2[id_used2[order(df_HVb2$dispersion_norm,decreasing=T)]] = 1:length(id_used2)
norm_rank = .norm_rank_diff(rank$rank1,rank$rank2)
# retrieve
id_retrieve = which(norm_rank>quantile(norm_rank, 0.99)) # not sensitive to the quantile

# ---------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------
ttp = as.matrix(m_seurat[,match(id_final,id_used)])
#tsne
set.seed(0)
tsne<-Rtsne(ttp,pca=F,is.distance=F)
#clustering
m_pca = prcomp(ttp,center = T,scale. = T) # need use cor matrix
m_pca_data = predict(m_pca, newdata=ttp)[,1:30]
set.seed(4)
cl <- kmeans(m_pca_data,7,iter.max=100,algorithm="MacQueen")
save(tsne,cl,file="4k3k_ksNo0_106.RData")

#plot batch effect
png(filename="batch.png")
.batch_plot(tsne,c(4340,3555))
dev.off()
# plot clusting
png(filename="cluster.png",width = 800, height = 800)
.cluster_plot(tsne,cl)
dev.off()
#plot marker genes
marker = c("CD3D","IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7","HBA1","HBA2","HBB")
png("biomarker.png",height=1200,width=1200)
.biomarker_plot(m_pbmc2[,marker],tsne$Y)
dev.off()
