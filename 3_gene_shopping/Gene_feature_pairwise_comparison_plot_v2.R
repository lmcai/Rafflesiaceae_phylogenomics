library(Hmisc)
library(psych)
library(corrplot)

#x=read.csv('/Users/limingcai/Dropbox/Rafflesiaceae_phylogeny_manuscript/TableS1_G2141_gene_features_full.csv')
#drops <- c("Gene_ID","stem_raff_branch_len","ancestor_of_stem_raff_bran_len","raff_placement_SHaLRT","raff_placement_UFBP","ID","dist2root")
#x=x[ , !(names(x) %in% drops)]
#x=x[x$stem_raff_branch_len.1<1,]
#x=x[x$ancestor_of_stem_raff_bran_len.1<0.1,]
#x=x[x$root2tip_variance<0.5,]

setwd('Dropbox/Rafflesiaceae_phylogeny/Figures/Fig1_Gene_feature_pairwise_comparison/')
x=read.csv('Gene_feature4plot_v2.csv')
x=x[ , !(names(x) %in% c("Gene_ID","Median_r2t_distance","align_len"))]
x=x[x$Median_r2t_distance<2,]
x=x[x$BrLen_ancestor_of_stem_Apodanthaceae<0.1,]
x=x[x$BrLen_ancestor_of_stem_Rafflesiaceae<0.1,]
x=x[x$Median_r2t_distance<0.5,]
x=x[x$r2t_variance<0.075,]
#pot scatterplot
pairs.panels(x[,1:10], method = "pearson", pch = 20,cex=0.1, ellipses=FALSE,smooth=FALSE,density=FALSE,)

res2 <- rcorr(as.matrix(x))

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

flattenCorrMatrix(res2$r, res2$P)
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

#heatmap
heatmap(x = res2$r, col = colorRampPalette(brewer.pal(8, "Oranges"))(25), symm = TRUE)
legend(x="bottomright", legend=c("min", "ave", "max"), fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))