
require("Matrix")
require("dplyr")
require("plyr")
require("Rtsne")
require("LaplacesDemon")
require("entropy")
require("ggplot2")
require("data.table")
require("matrixStats")
require("gridExtra")

# ---------------------------------------------------------------------
# normalize barcode matrix - by cell(cosine normalization), Macosko(cosine -> log)
# ---------------------------------------------------------------------
#x: matrix or dgTMatrix with row cell and column gene
#nc: scale factor, if nc==0, use the median of rowSums as scale factor(10x)
.normalize <- function(x,nc=10000,n=3,log.space=T){
  cs <- Matrix::colSums(x)
  used_genes <- which(cs >= n)
  x <- x[,used_genes]
  rs <- Matrix::rowSums(x)
  if(nc==0){nc = median(rs);log.space=F} #10x
  x = x/(rs/nc)
  if(log.space){x <- log(x+1)}
  list(m=x,used_genes=used_genes)
}


# --------------------------------------------------------------
# get variable genes from normalized UMI counts
# --------------------------------------------------------------
.get_variable_gene_poisson<-function(m,n=1700) {
  if(sum(Matrix::colSums(m)==0)>0){print("Warning: all-0 columns exist")}
  else{
    m = exp(m)-1
    df<-data.frame(mean=Matrix::colMeans(m),cv=apply(m,2,sd)/Matrix::colMeans(m),var=apply(m,2,var))
    df$dispersion<-with(df,var/mean)
    df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
    var_by_bin<-ddply(df,"mean_bin",function(x) {
      data.frame(bin_median=median(x$dispersion),
                 bin_mad=mad(x$dispersion))
    })
    df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
    df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
    df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
    # decide HV genes
    cut_off<-sort(df$dispersion_norm,decreasing=T)[n]
    df$used<-df$dispersion_norm >= cut_off
    df
  }
}



# --------------------------------------------------------------------------
# "normalize" the rank of HV genes, consider of both two batch rank
# --------------------------------------------------------------------------
# .norm_rank_diff = function(rank1,rank2){
#   q1 = rank1/max(rank1)
#   q2 = rank2/max(rank2)
#   abs(q1-q2)/min(q1,q2)
#
#   (max_rank-min_rank)/min_rank
# }

.norm_rank_diff = function(rank1,rank2){
  q1 = rank1/max(rank1)
  q2 = rank2/max(rank2)
  abs(q1-q2)
}
