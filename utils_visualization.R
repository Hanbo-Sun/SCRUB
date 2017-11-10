# --------------------------------------------------------------------
# plot dispersion vs. expression and HV genes
# --------------------------------------------------------------------
.plot_HV = function(df){
  ggplot(df,aes(mean,dispersion,color=as.factor(used)))+geom_point(size=0.2)+
  scale_color_manual(values=c("black","red"))
}


#--------------------------------------------------------------------------
# Visualization biomarker plot
#--------------------------------------------------------------------------
.biomarker_plot = function(m,Y){
  nm = colnames(m)
  n = length(nm)
  p <- list()
  for (i in 1:n){
    p[[i]] = ggplot(data.frame(tSNE1=Y[,1],tSNE2=Y[,2],Expression=log(m[,i]+1)),aes(tSNE1,tSNE2,color=Expression))+
    geom_point(size=0.3)+scale_color_continuous(name='Expression in Log scale',low='grey', high='blue')+theme(legend.position="none")+coord_fixed(1)+
    ggtitle(nm[i])
  }
  grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]])
}

# --------------------------------------------------------------------
# Batch effect -
# --------------------------------------------------------------------
.batch_plot = function(tsne_pbmc,batchVec){
  nb = length(batchVec)
  myColors <- c("blue"," green1", "hotpink","mediumorchid2", "red", "yellow", "purple1", "darkgoldenrod2","deepskyblue","darkgreen")
  dff = data.frame(tsne1 = tsne_pbmc$Y[,1],tsne2 = tsne_pbmc$Y[,2])
  dff$batch = c(rep(1,batchVec[1]),rep(2,batchVec[2]),rep(3,batchVec[3])) %>% as.factor()
  ggplot(dff,aes(tsne1,tsne2,col=batch))+geom_point(size=0.4,alpha=0.6)+theme_classic()+
     scale_color_manual(values=myColors[1:nb])+coord_fixed(1) # "blue","yellow"; "yellow","blue"
}

# --------------------------------------------------------------------
# clustering - visualization
# --------------------------------------------------------------------
.cluster_plot = function(tsne_pbmc,cl_pbmc){
  myColors <- c("gray48"," green1", "hotpink","mediumorchid2", "red", "yellow", "purple1", "darkgoldenrod2","deepskyblue","darkgreen")
  dff = data.frame(tsne1 = tsne_pbmc$Y[,1],tsne2 = tsne_pbmc$Y[,2],cluster=cl_pbmc$cluster)
  k = length(unique(cl_pbmc$cluster))
  dff$cluster = as.factor(dff$cluster)
  ggplot(dff,aes(tsne1,tsne2,col=cluster))+geom_point(size=0.4,alpha=0.6)+
    scale_color_manual(values=myColors[1:k])+
    theme_classic()+coord_fixed(1)
}

# --------------------------------------------------------------------
# label - visualization
# --------------------------------------------------------------------
.label_plot <- function(Y, tittle = NULL,subset=NULL, ...) {
    myColors <- c("gray48"," green1", "hotpink","mediumorchid2", "red", "yellow", "purple1", "darkgoldenrod2","deepskyblue","darkgreen")
    k = length(unique(rownames(m)))
    dff = data.frame(tsne1 = Y[,1],tsne2 = Y[,2],type=rownames(m))
    ggplot(dff,aes(tsne1,tsne2,col=as.factor(type)))+geom_point(size=0.4,alpha=0.6)+
      scale_color_manual(values=myColors[1:k])+
      theme_classic()+coord_fixed(1)

}
