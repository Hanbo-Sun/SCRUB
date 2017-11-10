###/////////////////////////////////////////////////////////////////
# This R file "BEC_p4kt3k.R" is to correct batch effect by criterion
# by which to to identify BESGs as accurate as possible.
###/////////////////////////////////////////////////////////////////

# ---------------------------------------------------------------------
# get all BESG genes' absolute index (305 total)
# ---------------------------------------------------------------------
load("p4kTsub_t3k.RData") #2380+ 3555 =5935 # 55% T cell in PBMC oringinal data
m_sub = m_pbmc #5935x33694
id_used_sub = which(colSums(m_sub)>=3) #16925
m_sub = .normalize(m_sub,n=min.cell,nc=10000,log.space=T)$m
p.ks = .batch_gene(m_sub[1:2380,],m_sub[2381:5935,],test="ks")
#save(p.ks,file="pks_4K3K.RData")
id_besg = which(p.ks==0)
id305 = id_used_sub[id_besg]

# ---------------------------------------------------------------------
# get HV genes of batch 1+2, batch 1, batch 2
# ---------------------------------------------------------------------
load("p4k_t3k_self.RData")
m = m_pbmc2 #7895 33694 range:0-8029
rs = rowSums(m) #rs range:856 to 48443 1qt:3089 3qt:4754, median:3873
id_b1 = 1:4340;id_b2 = 4341:7895
#pbmc cells #PanT cells
m1 = m[id_b1,];m2 = m[id_b2,]
min.cell = 3
id_used = which(colSums(m)>=min.cell) #18084
id_used1 = which(colSums(m1)>=min.cell) #16707
id_used2 = which(colSums(m2)>=min.cell) #15656

# ---------------------------------------------------------------------
# normalization and HV
# ---------------------------------------------------------------------
m_seurat = .normalize(m,n=min.cell,nc=10000,log.space=T)$m
m_seurat1 = .normalize(m1,n=min.cell,nc=10000,log.space=T)$m
m_seurat2 = .normalize(m2,n=min.cell,nc=10000,log.space=T)$m
# HV
df_HV<-.get_variable_gene_poisson(m_seurat,n=1700) # must use log data
df_HVb1<-.get_variable_gene_poisson(m_seurat1,n=5000) # must use log data
df_HVb2<-.get_variable_gene_poisson(m_seurat2,n=5000) # must use log data
idHV = id_used[df_HV$used]


# ---------------------------------------------------------------------
# Get the 127 BESG absolute index
# ---------------------------------------------------------------------
id_rm2 = intersect(idHV,id305) #127
id_rm = intersect(idHV,id2) #127
length(id_rm)
length(intersect(id_rm2, id_rm))

save(id_rm,file="id_rm_p4kt3k.RData")
load("id_rm_p4kt3k.RData")
id_final = setdiff(idHV,id305) #1573
id_final = setdiff(idHV,id2) #1573

#
#id_final2 = setdiff(idHV,id305[150:305]) #1578 #1580 #1615 #150:305?
length(id_final)

#tsne
ttp = as.matrix(m_seurat[,match(id_final,id_used)])
ttp = as.matrix(m_seurat[,match(id_final2,id_used)])
ttp = as.matrix(m_seurat[,df_HV$used][,setdiff(1:1700,id)])

set.seed(0)
tsne<-Rtsne(ttp,pca=F,is.distance=F)
tsne<-Rtsne(m_pca_data,pca=F,is.distance=F)
#clustering
m_pca = prcomp(ttp,center = T,scale. = T) # need use cor matrix
m_pca_data = predict(m_pca, newdata=ttp)[,1:30]
set.seed(4)
cl <- kmeans(m_pca_data,kpp_init(m_pca_data, 6),iter.max=100,algorithm="MacQueen")
cl <- kmeans(m_pca_data,6,iter.max=100,algorithm="MacQueen")
set.seed(0)
cl <- kmeans(as.matrix(m_seurat[,match(id_final,id_used)]),6,iter.max=100,algorithm="MacQueen")
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

# ---------------------------------------------------------------------
# Give rank for each genes for batch 1+2, 1, 2
# ---------------------------------------------------------------------
rank = data.frame(rank12=rep(30000,33694), rank1=rep(30000,33694), rank2=rep(30000,33694))
rank$rank12[id_used[order(df_HV$dispersion_norm,decreasing=T)]] = 1:length(id_used)
rank$rank1[id_used1[order(df_HVb1$dispersion_norm,decreasing=T)]] = 1:length(id_used1)
rank$rank2[id_used2[order(df_HVb2$dispersion_norm,decreasing=T)]] = 1:length(id_used2)


plot_rank = rank[idHV,]
plot_rank$group = rep("Good gene",1700)
plot_rank$group[which(rownames(plot_rank) %in% as.character(id_rm))] = "Bad gene"
plot_rank$s1 = with(plot_rank,rank12-min(rank1,rank2))
plot_rank$s2 = apply(plot_rank[,c("rank1","rank2")], 1, max)
plot_rank$s3 = apply(plot_rank[,c("rank1","rank2")], 1, min)
plot_rank$s4 = with(plot_rank, ave(rank1,rank2))

summary(plot_rank$s1[which(plot_rank$group == "Bad gene")])
summary(plot_rank$s1[which(plot_rank$group == "Good gene")])

png("rank_plot2.png",height=800,width=800)
ggplot(plot_rank,aes(s1,s3,col=group))+geom_point()
dev.off()

png("rank_density2.png",height=800,width=800)
ggplot(plot_rank,aes(x=rank2, fill=group)) + geom_density(alpha=0.25)
dev.off()


# ---------------------------------------------------------------------
# try to find out decision to classify batch 1 and batch 2 - inspired by the view that 123/127 BESG still ks=0 for batch 1 v.s. batch 2
# ---------------------------------------------------------------------
p.ks = .batch_gene(as.matrix(m_seurat[1:4340,]),as.matrix(m_seurat[4341:7895,]),test="ks")
test = data.frame(ks=rep(0,18084), exp1=rep(0,18084), exp2=rep(0,18084))
test$ks = -log(p.ks)
test$ks[which(test$ks>30)] = 30
test$exp1 = colSums(m_pbmc2[1:4340,id_used])
test$exp2 = colSums(m_pbmc2[4341:7895,id_used])
test$group = rep("Good gene",18084)
test$group[match(id_rm,id_used)] = "Bad gene"

plot_test = test[df_HV$used,]
id30 = which(plot_test$ks==30)
idHV[id30]
plot30 = cbind(plot_test[id30,1:3],plot_rank[id30,])
plot30$group =as.factor(plot30$group)
colnames(plot30)
control <- trainControl(method="cv", number=5)
set.seed(0)
(lda_mod <- train(group~exp1+exp2+s1+s2+s4, data=plot30,method="C5.0", trControl=control, allowParallel=TRUE))
summary(plot30$group)

require("party")
iris_ctree <- ctree(group~exp1+exp2+s1+s2+s4, data=plot30)
print(iris_ctree)
plot(iris_ctree, cex=2)
dev.off()

id22 <- which((plot30$exp1<4000)&(plot30$exp2<600))
id_final = idHV
idHV[id30][id22]
id_final = setdiff()
id_final2 = union(setdiff(idHV,idHV[id30]),idHV[id30[id22]])
#1333
length(id_final2) #1478
setdiff(idHV,idHV[id30])
idHV[id30[id22]]

#plot_test = cbind(plot_test[,4],as.data.frame(scale(plot_test[,1:3])))
#colnames(plot_test)[1] = "group"
summary(as.factor(plot_test$group))

png("test_plot3.png",height=800,width=800)
ggplot(plot_test,aes(exp1,exp2,col=group))+geom_point()
dev.off()

sum(plot_test$ks==30)
id = which(plot_test$ks==30)

control <- trainControl(method="cv", number=5)
set.seed(0)
(rf_mod <- train(group~s1+s2, data=plot_rank,method="rf", trControl=control, allowParallel=TRUE))
set.seed(0)
(lda_mod <- train(group~ks+exp1+exp2, data=plot_test,method="lda", trControl=control, allowParallel=TRUE))
set.seed(0)
(lda_mod <- train(group~ks+exp1+exp2, data=plot_test,method="C5.0", trControl=control, allowParallel=TRUE))
set.seed(0)
(rf_mod <- train(group~ks+exp1+exp2, data=plot_test,method="rf", trControl=control, allowParallel=TRUE))


# ---------------------------------------------------------------------
# -----------------------continue batch effect genes 0918--------------
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# part 1 - ks test in batch 1+2 vs batch 1; and 1+2 v.s 2
# ---------------------------------------------------------------------
#p.ks1 = .batch_gene(as.matrix(m_seurat),as.matrix(m_seurat[1:4340,]),test="ks")
#p.ks2 = .batch_gene(as.matrix(m_seurat),as.matrix(m_seurat[4341:7895,]),test="ks")

p.chisq = .batch_gene(m_pbmc2[1:4340,idHV], m_pbmc2[4341:7895,idHV])
save(p.chisq, file = "p4kt3k_Pchisq.RData")
save(p.chisq, file = "p4kt3k_Pchisqnon0.RData")
summary(-log(p.chisq))

id_p = idHV[which(-log(p.chisq)>5)]
length(id_p)
length(intersect(id_p, id_rm))
intersect(setdiff(id_p,id_retrieve), id_rm)



rank = data.frame(rank12=rep(30000,33694), rank1=rep(30000,33694), rank2=rep(30000,33694))
rank$rank12[id_used[order(df_HV$dispersion_norm,decreasing=T)]] = 1:length(id_used)
rank$rank1[id_used1[order(df_HVb1$dispersion_norm,decreasing=T)]] = 1:length(id_used1)
rank$rank2[id_used2[order(df_HVb2$dispersion_norm,decreasing=T)]] = 1:length(id_used2)
norm_rank = .norm_rank_diff(rank$rank1,rank$rank2)


id_retrieve = which(norm_rank>quantile(norm_rank, 0.95)) # not sensitive to the quantile
#id_retrieve = which(rank[,3]>10000&rank[,2]<2000| rank[,2]>10000&rank[,3]<2000)
intersect(id_retrieve, id_p)
intersect(id_retrieve, id_rm2)


id_final = setdiff(idHV, setdiff(id3,id_retrieve))# Assume we can identify those 104 genes out of 368
length(id_final)

intersect(id_rm, idHV)

length(setdiff(id_p,id_retrieve))

intersect(id_final, id_rm)


df_ks = data.frame(ks1=-log(p.ks1), ks2=-log(p.ks2))
# consider truncated
#which(df_ks[,1]>30)
#which(df_ks[,2]>50)
id_ks = intersect(which(df_ks[,1]>30),which(df_ks[,2]>30)) #remove 217 from HV
id_ks = intersect(which(df_ks[,1]>20),which(df_ks[,2]>30)) # remove 247 from HV

length(id_ks)
#id_ks = union(which(df_ks[,1]>30),which(df_ks[,2]>30))
id_ks = id_used[id_ks]

rank[rank[,3]>10000,3]

tp = rank[setdiff(id_ks,id_rm),]
tp = tp[which(tp[,1]<=1700),]

tp = rank[id_ks,]
tp = tp[tp$rank2>10000,]
idPlus = rownames(tp)


sum(tp$rank2>(tp$rank12+tp$rank1))
tp = tp[which(tp$rank2>(tp$rank12+tp$rank1)),]
tp = tp[(tp$rank2>10000),]
dim(tp)
tp
tp$rank2>(tp$rank12+tp$rank1)


id_final = setdiff(idHV, id_ks)# rm all 217(from368,in HV 217) genes # remove 2
id_final = setdiff(idHV, setdiff(id_rm,id_retrieve))# Assume we can identify those 104 genes out of 368

length(id_retrieve)
length(id_final)
#tsne
ttp = as.matrix(m_seurat[,match(id_final,id_used)])
set.seed(0)
tsne<-Rtsne(ttp,pca=F,is.distance=F)
set.seed(0)
cl <- kmeans(ttp,6,iter.max=100,algorithm="MacQueen")
#or pca clustring
m_pca = prcomp(ttp,center = T,scale. = T) # need use cor matrix
png("PC.png",width = 600, height = 480)
plot(m_pca, type = "l",npcs=50)
dev.off()
m_pca_data = predict(m_pca, newdata=ttp)[,1:20]
set.seed(0)
cl_pca <- kmeans(m_pca_data,6,iter.max=100)
save(tsne,cl,cl_pca,file="tsne_cl_clpca_4k3k_kmmd_0.10non0.RData")

load("tsne_cl_clpca_final_15_onlyRM.RData")
# Visualization
#plot batch
png(filename="batch.png")
.batch_plot(tsne,c(4340,3555))
dev.off()
#plot clusting
png(filename="cluster_pca4.png",width = 800, height = 800)
.cluster_plot(tsne,cl_pca)
dev.off()
png(filename="cluster.png",width = 800, height = 800)
.cluster_plot(tsne,cl)
dev.off()
#plot markers
marker = c("CD3D", "IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7","HBA1","HBA2","HBB")
#marker=c("IL4R","IL1A",        "IL1B",        "IL1R1",       "IL1R2",       "IL1RAP","IL1RL1",      "IL1RL2",      "IL1RN",       "IL2","IL10","IL10RA","IL10RB","IL10RB-AS1","IL11RA")
png("biomarker.png",height=1200,width=1200)
.biomarker_plot(m_pbmc2[,marker],tsne$Y)
dev.off()
