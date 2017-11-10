
.seq = function(m,tks=0.15,tmw=0.15,tbc=0.15,min.cell=3,np=0.1){
  nb = length(m) # No. Batches
  #check if g1==g2==...
  G = rep(0,nb)
  B = NULL
  for(i in 1:nb){
    G[i] = ncol(m[[i]])
    B = c(B,nrow(m[[i]]))
  }
  if(length(unique(G))!=1){print("Warning: different No. genes")}
  CSB = cumsum(B) # cumulative batch size
  # merge together -> normalize
  MM = m[[1]]
  if (nb!=1){
    for(i in 2:nb){
      MM = rbind(MM,m[[i]])
      print(dim(MM))
    }
  }
  m_seurat = .normalize(MM,n=min.cell,nc=10000,log.space=T)
  id_used = m_seurat$used_genes
  m_seurat = m_seurat$m
  g = np*length(id_used)
  df_HV<-.get_variable_gene_poisson(m_seurat,n=g) # must use log data
  id_HV = id_used[df_HV$used]
  # separate
  M = list()
  M[[1]] = as.matrix(m_seurat[1:CSB[1],df_HV$used])
  print(CSB)
  for (i in 2:nb){M[[i]] = as.matrix(m_seurat[(CSB[i-1]+1):CSB[i],df_HV$used])}
  # all possible combinations
  filter_vote = rep(0,g)
  for (i in 1:(nb-1)){
    if(i==1){m1 = M[[1]];m2 = rbind(M[[2]],M[[3]])}
    if(i==2){m1 = M[[2]];m2 = rbind(M[[1]],M[[3]])}
    if(i==3){m1 = M[[3]];m2 = rbind(M[[1]],M[[2]])}
    #if(i==1){m1 = M[[1]]}
    #else{m1=rbind(m1,M[[i]])}
    #m2 = M[[i+1]]
    # KS test
    pv = .batch_ks(m1,m2,g,mark="use") # !!!!
    vote_ks = pv<=quantile(pv,tks)
    print(quantile(pv,tks))
    print(sum(vote_ks))
    # MW test
    mw = .batch_mw(m1,m2,g,mark="ignore")
    mw[is.na(mw)]=1 # ALL 0 for both compared batch, regard as same
    vote_mw = mw<=quantile(mw,tmw,na.rm=T)
    print(quantile(mw,tmw))
    print(sum(vote_mw))
    # bc distance
    bc = .batch_mw(m1,m2,g,mark="use")
    print(sum(is.na(bc)))
    bc[is.na(bc)]=0
    vote_bc = bc>=quantile(bc,1-tbc)
    print(quantile(bc,1-tbc))
    print(sum(vote_bc))
    # filter votes:
    temp = (vote_ks + vote_mw + vote_bc)>=2
    filter_vote = filter_vote + temp
  }
  print(table(filter_vote))
  hv_sig = rep(TRUE,G[1]) #non highly variable genes remove anyway
  hv_sig[id_HV] = ((filter_vote)>0)
  hv_sig
}

.com = function(m,tks=0.15,tmw=0.15,tbc=0.15,min.cell=100,np=0.1,use=c(TRUE,FALSE,TRUE)){
  nb = length(m) # No. Batches
  #check if g1==g2==...
  G = rep(0,nb)
  B = NULL
  for(i in 1:nb){
    G[i] = ncol(m[[i]])
    B = c(B,nrow(m[[i]]))
  }
  if(length(unique(G))!=1){print("Warning: different No. genes")}
  CSB = cumsum(B) # cumulative batch size
  # merge together -> normalize
  MM = m[[1]]
  if (nb!=1){
    for(i in 2:nb){
      MM = rbind(MM,m[[i]])
      print(dim(MM))
    }
  }
  m_seurat = .normalize(MM,n=min.cell,nc=10000,log.space=T)
  id_used = m_seurat$used_genes
  m_seurat = m_seurat$m
  g = round(np*length(id_used))
  df_HV<-.get_variable_gene_poisson(m_seurat,g) # must use log data
  id_HV = id_used[df_HV$used]
  # separate
  M = list()
  M[[1]] = as.matrix(m_seurat[1:CSB[1],df_HV$used])
  print(CSB)
  for (i in 2:nb){M[[i]] = as.matrix(m_seurat[(CSB[i-1]+1):CSB[i],df_HV$used])}
  # all possible combinations
  cb = combn(nb,2)
  vote = matrix(0,nrow=ncol(cb),ncol=g)
  for (i in 1:ncol(cb)){
    m1 = M[[cb[1,i]]]
    m2 = M[[cb[2,i]]]
    # KS test
    pv = .batch_ks(m1,m2,g,mark=ifelse(use[1],"use","ignore")) # !!!!
    vote_ks = pv<=quantile(pv,tks)
    print(quantile(pv,tks))
    print(sum(vote_ks))
    # MW test
    mw = .batch_mw(m1,m2,g,mark=ifelse(use[2],"use","ignore"))
    mw[is.na(mw)]=1 # ALL 0 for both compared batch, regard as same
    vote_mw = mw<=quantile(mw,tmw,na.rm=T)
    print(quantile(mw,tmw))
    print(sum(vote_mw))
    # bc distance
    bc = .batch_mw(m1,m2,g,mark=ifelse(use[3],"use","ignore"))
    print(sum(is.na(bc)))
    bc[is.na(bc)]=0
    vote_bc = bc>=quantile(bc,1-tbc)
    print(quantile(bc,1-tbc))
    print(sum(vote_bc))
    # filter votes:
    temp = (vote_ks + vote_mw + vote_bc)>=2
    vote[i,] = temp
  }
  vote = colSums(vote)
  print(table(vote))
  list(id_HV = id_HV, m_seurat=m_seurat, vote = vote)
}


.preserve = function(m,tpre=0.05,np=0.1){
  m1 = m[[1]]
  m2 = m[[2]]
  m3 = m[[3]]
  m_seurat1 = .normalize(m1,n=min.cell,nc=10000,log.space=T)
  id1 = m_seurat1$used_genes
  m_seurat1 = m_seurat1$m
  m_seurat2 = .normalize(m2,n=min.cell,nc=10000,log.space=T)
  id2 = m_seurat2$used_genes
  m_seurat2 = m_seurat2$m
  m_seurat3 = .normalize(m3,n=min.cell,nc=10000,log.space=T)
  id3 = m_seurat3$used_genes
  m_seurat3 = m_seurat3$m
  d1 = .get_variable_gene_poisson(m_seurat1,n=round(ncol(m1)*np))$dispersion_norm
  d2 = .get_variable_gene_poisson(m_seurat2,n=round(ncol(m2)*np))$dispersion_norm
  d3 = .get_variable_gene_poisson(m_seurat3,n=round(ncol(m3)*np))$dispersion_norm
  rank1 = rep(ncol(m_seurat1)+1,ncol(m1))
  rank1[id1] = rank(-d1,ties.method = "first")
  rank2 = rep(ncol(m_seurat2)+1,ncol(m2))
  rank2[id2] = rank(-d2,ties.method = "first")
  rank3 = rep(ncol(m_seurat3)+1,ncol(m3))
  rank3[id3] = rank(-d3,ties.method = "first")
  rank12 = .norm_rank_diff(rank1,rank2)
  rank13 = .norm_rank_diff(rank1,rank3)
  rank23 = .norm_rank_diff(rank2,rank3)
  pre1 = rank12>quantile(rank12,1-tpre)
  pre2 = rank13>quantile(rank13,1-tpre)
  pre3 = rank23>quantile(rank23,1-tpre)
  res = pre1 + pre2 + pre3
  res
}
