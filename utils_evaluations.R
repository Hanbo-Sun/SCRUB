# ----------------------------------------------------------------------------
# calculate_between_within
# ----------------------------------------------------------------------------
#column 1 and 2 tSNE$Y column 3 clustering label
.calculate_between_within = function(x){
  n = nrow(x)
  tb <- table(x[,3])
  k=length(tb)
  ac_tb = rep(0,k)
  for (i in 1:k){
    ac_tb[i] = sum(tb[1:i])
  }

  x_sort = x[order(x[,3]),]
  #print(x_sort[1:10,1:2])
  dist = as.matrix(dist(x_sort[,1:2]))
  #print(dist[1:10,1:10])
  within = rep(0,k)
  for (i in 1:k){
    starti = 1
    endi = ac_tb[i]
    if(i!=1){starti = ac_tb[i-1]+1}
    #print(starti)
    #print(endi)
    within[i] = sum(sum(dist[starti:endi,starti:endi]))
  }
  within[within==0]=0

  #print(class(within))
  count_within = (sum(tb^2)-n)
  count_bet = n^2 - sum(tb^2)
  within_sum = sum(within)
  within_av = within_sum/count_within
  within = within/(tb*(tb-1))
  within[is.na(within)]=0
  #print(is.na(within))
  dist_sum = sum(sum(dist))
  #print(dist_sum)
  bet_av = (dist_sum - within_sum)/count_bet
  ratio = bet_av/within_av
  list(within=within,within_av=within_av,bet = bet_av,ratio=ratio)
}

# --------------------------------------------------------------------------
# adjust cluter
# --------------------------------------------------------------------------
.adjust_cluter = function(cl){
  tb <- sort(table(cl))
  for (i in 1:length(tb)){
    cl[cl==names(tb)[i]]=100*i
  }
  cl = cl/100
  cl
}

# --------------------------------------------------------------------------
# calculate trivial(Null) consistency
# --------------------------------------------------------------------------
.trivial_consist = function(cl1,cl2){
  k = length(table(cl1))
  x = matrix(0,nrow=k,ncol=k)
  freq1 = table(cl1)/length(cl1)
  freq2 = table(cl2)/length(cl2)
  for (i in 1:k){
    for (j in 1:k){
      x[i,j] = freq1[i]*freq2[j]
    }
  }
  list(x=x,trivial=sum(diag(x)))
}




# -----------------------------------------------------------------------------
# L2FC
# -----------------------------------------------------------------------------
.l2fc = function(data,cluster){
  p = nrow(data)
  n = ncol(data)
  tb = table(cluster)
  k = length(tb)
  #initialize
  counts = matrix(0,nrow=k,ncol=p)
  adj_counts = matrix(0,nrow=k,ncol=p)
  best_cl = rep(0,p)
  adj_best_cl = rep(0,p)
  l2fc = matrix(0,nrow=p,ncol=3)
  adj_l2fc = matrix(0,nrow=p,ncol=3)
  df = t(data)
  #df %>% group_by(cl) %>% summarise()
  for (j in 1:p){
    for (i in 1:n){
      counts[cluster[i],j] = counts[cluster[i],j] + data[j,i]
    }
  }

  adj_counts = counts/matrix(tb,nrow=k,ncol=p)
  best_cl = max.col(t(counts), 'first')
  adj_best_cl = max.col(t(adj_counts), 'first')
  for (i in 1:p){
    #print(i)
    #print(best_cl[i])
    l2fc[i,1] = counts[best_cl[i],i] #col1: best sum counts
    adj_l2fc[i,1] = counts[adj_best_cl[i],i]/tb[adj_best_cl[i]] # select per cell highest
    l2fc[i,2] = sum(counts[-best_cl[i],i])
    adj_l2fc[i,2] = sum(counts[-adj_best_cl[i],i])/sum(tb[-adj_best_cl[i]])
  }
  l2fc[,3] = log2(l2fc[,1]/l2fc[,2])
  adj_l2fc[,3] = log2(adj_l2fc[,1]/adj_l2fc[,2])
  list(best_cluster=best_cl,adj_best_cluster=adj_best_cl,l2fc = l2fc,adj_l2fc=adj_l2fc)
}


#------------------------------------------------------------------------------
# tsne - r square
#------------------------------------------------------------------------------
# return p value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

.tsne_exp_reg = function(data,tsne){
  # row genes(markers), col cells
  p = nrow(data)
  r2 = rep(0,p)
  pvalue = rep(0,p)
  ynames = rownames(data)
  df = cbind(t(data)%>% as.matrix(),tsne) %>% as.data.frame()
  colnames(df) = c(ynames,"tSNE1","tSNE2")
  for (i in 1:p){
    #print(i)
    lm_tp = lm(paste(ynames[i],'~tSNE1+tSNE2'),data = df)
    r2[i] = summary(lm_tp)$r.squared
    pvalue[i] = lmp(lm_tp)
  }
  list(r2,pvalue)
}
