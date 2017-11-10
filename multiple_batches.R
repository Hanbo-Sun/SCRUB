


# ---------------------------------------------------------------------
# Test on P4k + t3k
# ---------------------------------------------------------------------
#load rm genes index
load("id_rm_p4kt3k.RData")
#load p4k+t3k
load("p4k_t3k_self.RData")
min.cell = 3
m = m_pbmc2
id_used = which(Matrix::colSums(m)>=min.cell)
id_b1 = 1:4340;id_b2 = 4341:7895
vec_sig = .gene_test(list(m[id_b1,],m[id_b2,]))
vec_preserve = .gene_preserve(c(m[id_b1,],m[id_b2,]))
id_final = which(vec_sig - vec_preserve<=0)
m_seurat = .normalize(m,n=min.cell,nc=10000,log.space=T)$m #move our 116 out of 127 BESG




# ---------------------------------------------------------------------
# Test on p6k + P4k + t3k
# ---------------------------------------------------------------------
# load pbmc_6k
m6k = readMM("pbmc_6k/hg19/matrix.mtx") %>% t()  #5419 barcode and 32738 genes
barcode_6k = read.table("pbmc_6k/hg19/barcodes.tsv",sep='\t',stringsAsFactors=F)
gene_6k = read.table("pbmc_6k/hg19/genes.tsv",sep='\t',stringsAsFactors=F)
rownames(m6k) = barcode_6k[,1]
colnames(m6k) = gene_6k[,1]
# load pbmc_4k
m4k = readMM("pbmc_4k/GRCh38/matrix.mtx") %>% t()  #4340 barcode and 33694 genes
barcode_4k = read.table("pbmc_4k/GRCh38/barcodes.tsv",sep='\t',stringsAsFactors=F)
gene_4k = read.table("pbmc_4k/GRCh38/genes.tsv",sep='\t',stringsAsFactors=F)
rownames(m4k) = barcode_4k[,1]
colnames(m4k) = gene_4k[,1]
# load t_3k
t3k = readMM("t_3k/GRCh38/matrix.mtx") %>% t()  #4340 barcode and 33694 genes
barcode_3k = read.table("t_3k/GRCh38/barcodes.tsv",sep='\t',stringsAsFactors=F)
gene_3k = read.table("t_3k/GRCh38/genes.tsv",sep='\t',stringsAsFactors=F)
rownames(t3k) = barcode_3k[,1]
colnames(t3k) = gene_3k[,1]
# merge p6k p4k and t3k
length(intersect(barcode_4k[,1], barcode_6k[,1])) # 0 that's good
length(intersect(barcode_3k[,1], barcode_6k[,1])) # 0 that's good
length(intersect(barcode_3k[,1], barcode_4k[,1])) # 21 intersect
common_genes = intersect(gene_4k[,1], gene_6k[,1]) # 31232 common genes, known t3k and p4k same genes
mul = rbind(m6k[,common_genes],m4k[,common_genes],t3k[,common_genes])
rownames(gene_4k) = gene_4k[,1]
colnames(mul) = gene_4k[common_genes,2]
#save(mul,fil0e="p6k_p4k_t3k.RData") # 9759 31232
load("p6k_p4k_t3k.RData")
min.cell = 100
min.gene = 1000
m = mul

id_b1 = 1:5419; id_b2 = 5420:9759; id_b3 = 9760:13314

m1 = m[id_b1,]
m2 = m[id_b2,]
m3 = m[id_b3,]
rs1 = rowSums(m1)
rs2 = rowSums(m2)
rs3 = rowSums(m3)
m1 = m1[rs1>min.gene,]
m2 = m2[rs2>min.gene,]
m3 = m3[rs3>min.gene,]

m_qc = rbind(m1,m2,m3)
id_b1 = 1:nrow(m1)
id_b2 = (nrow(m1)+1):(nrow(m1)+nrow(m2))
id_b3 = (nrow(m1)+nrow(m2)+1):(nrow(m1)+nrow(m2)+nrow(m3))
id_used = which(Matrix::colSums(m_qc)>=min.cell)
length(id_used)

m1_qc = m_qc[id_b1,id_used]
m2_qc = m_qc[id_b2,id_used]
m3_qc = m_qc[id_b3,id_used]
dim(m1_qc)
# ---------------------------------------------------------------------
# combination method for  p6k + P4k + t3k
# ---------------------------------------------------------------------
vec_sig = .com(list(m1_qc,m2_qc,m3_qc),tks=0.08,tmw=0.08,tbc=0.04,min.cell=100)
m_seurat = vec_sig$m_seurat
id_HV = vec_sig$id_HV
vec_sig = vec_sig$vote
dim(m_seurat)
length(id_HV)

vec_preserve = .preserve(c(m1_qc,m2_qc,m3_qc),tpre=0.05)
vec_preserve = vec_preserve[id_HV]
length(vec_preserve)

id_final_com = which(vec_sig - vec_preserve<=0) #1536 used, 1771 hV;

id_final_com = setdiff(which(vec_sig - vec_preserve<=1),which(vec_sig==1&&vec_preserve==0)) #1536 used, 1771 hV;

#com_u
id_final_com = union(which(vec_sig - vec_preserve<=1), which(vec_sig==3))#1536 used, 1771 hV;
#com_u2
id_final_com = union(which(vec_sig - vec_preserve<=0),which(vec_sig==1&&vec_preserve==0)) #1536 used, 1771 hV;


id_final_com = which(vec_sig<=0) #1536 used, 1771 hV;
length(id_final_com)

table(vec_sig)
# ---------------------------------------------------------------------
# sequential method for  p6k + P4k + t3k
# ---------------------------------------------------------------------
vec_sig_seq = .seq(list(m[id_b1,],m[id_b2,],m[id_b3,]),tks=0.10,tmw=0.10,tbc=0.08)
id_final_seq = which(!vec_sig_seq)
length(id_final)


# ---------------------------------------------------------------------
# The final id p6k + P4k + t3k
# ---------------------------------------------------------------------
id_final = union(id_final_seq,id_final_com)
length(id_final)



# ---------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------
ttp = as.matrix(m_seurat[,id_HV[id_final_com]])

ttp = as.matrix(m_seurat)

dim(ttp)
#tsne
set.seed(0)
tsne<-Rtsne(ttp,pca=F,is.distance=F)
#clustering
m_pca = prcomp(ttp,center = T,scale. = T) # need use cor matrix
m_pca_data = predict(m_pca, newdata=ttp)[,1:30]
set.seed(0)
cl <- kmeans(m_pca_data,6,iter.max=100,algorithm="MacQueen")
save(tsne,cl,file="tsne_cl_com_more.RData")
#load("tsne_cl_com_seq_1644.RData")

#plot batch effect
png(filename="batch_com_more.png",width = 800, height = 800)
.batch_plot(tsne,c(length(id_b1),length(id_b2),length(id_b3)))
dev.off()

# plot clusting
png(filename="cluster_com_more.png",width = 800, height = 800)
.cluster_plot(tsne,cl)
dev.off()
#plot marker genes
marker = c("CD3D","IL7R", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP", "S100A4", "CCR7","HBA1","HBA2","HBB")
png("biomarker_com_more.png",height=1200,width=1200)
.biomarker_plot(mul[which(rowSums(mul)>min.gene),marker],tsne$Y)
dev.off()
