########################################################################################
##Necessary libraries
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(glue)
library(dplyr)
library(psych)
library(PMA)
#######################################################################################
## Refer - 
##1.Daniela M Witten and Robert J Tibshirani. Extensions of sparse canonical correlation analysis with applications to genomic data. 
##Statistical applications in genetics and molecular biology, 8(1), 200
##2. Code was adpated from the work by - Nguyen, Quang P., et al. "Associations between the gut microbiome and metabolome in early life." 
##BMC microbiology 21.1 (2021): 1-19.
########################################################################################
data_cell <- read.csv('Filename.csv', header=TRUE, row.names = 2, check.names = FALSE) 
data_cell <- na.omit(data_cell)
data_metab <- read.csv('Filename.csv', header=TRUE, row.names = 2, check.names = FALSE)
data_metab <- na.omit(data_metab)
########################################################################################
## Use the CCA.permute function from the PMA package to perform a permutational test to determine 
##the best penalty values

perm.out <- CCA.permute(data_cell,data_metab,typex="ordered",typez="standard",nperms=50)

##The obtained penalty values are used to fit the main model. 
##The number of latent variable, K is set to 1

out <- CCA(data_cell,data_metab,typex="ordered",typez="standard",K=1,
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz)
########################################################################################
## Select the non-zero variables in each canonical variate
cell_idx <- which(out$u != 0)
met_idx <- which(out$v != 0)

## For all the variables find the correlation using the spearman correlation with a 
##Bejamini-Hochberg method for multiple testing, FDR <0.05
corr <- corr.test(data_cell,data_metab,use="pairwise",method="spearman",adjust="BH", alpha = 0.05,ci=TRUE)
sig <- t(apply(corr$p, 1, function(x){
  ifelse(x < 0.05, 1,0)
}))
cor_to_use <- corr$r
cor_to_use[which(sig == 0, arr.ind = 1)] <- NA
## Plot the observed correlation as a heatmap
reduced_heatmap <-  pheatmap(
  mat = t(cor_to_use),
  #na_col = "beige",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 315
)
###################################################################################################
## Plot the correlations for SCCA selected metabolites only, along with their loadings.
cca_mat <- cor_to_use[cell_idx, met_idx]
sig_mat <- sig[cell_idx, met_idx]
cca_mat[which(sig_mat == 0, arr.ind = 1)] <- NA
row <- data.frame("sCCA Loading" = as.factor(ifelse(out$u[cell_idx] > 0, "+","-")), check.names = F)
rownames(row) <- rownames(cca_mat)
col <- data.frame("sCCA Loading" = as.factor(ifelse(out$v[met_idx] > 0, "+","-")), check.names = F)
rownames(col) <- colnames(cca_mat)
ann_colors <- list("sCCA Loading" = c(cividis(10)[1], cividis(10)[10]))
names(ann_colors$"sCCA Loading") <- as.factor(c("+", "-")) 
reduced_heatmap <-  pheatmap(
  mat               = t(cca_mat),
  annotation_row    = col,
  annotation_colors = ann_colors,
  annotation_col    = row,
  annotation_names_row = T,
  annotation_names_col = T,
  cluster_rows = F,
  cluster_cols = F,
  show_colnames     = T,
  show_rownames     = T,
  drop_levels       = TRUE,
  fontsize          = 13,
  display_numbers = F,
  cex = 1, 
  legend = T,
  na_col = "beige",
  angle_col=315
  
)
###################################################################################################
## An additional permutational test is consucted on 5000 permutations to find the 
##significance of the observed the correlations
null_corr <- c()
n_perms = 5000

for (i in 1:n_perms){
  cell =cell
  met = met
  
  perm <- CCA.permute(x = cell, z = met, typex = "ordered",typez = "ordered", 
                      nperms = 50, niter = 25, standardize = T)
  cca_mod <- CCA(x = cell, z = met, typex = "ordered", typez = "ordered", 
                 niter = 25, penaltyx = perm$bestpenaltyx, 
                 penaltyz = perm$bestpenaltyz, standardize = T)
  null_corr[i] <- cca_mod$cors
} 

## p-value is determined 
raw_pval <- length(which(perm.out$cors >= out$cors))/length(perm.out$cor)
if(raw_pval == 0){
  raw_pval <- "0.001"
}
title <- glue("SCCA Correlation: {correlation}; 
               Permutational p-value: {pval}", correlation = round(out$cors, 3),
              pval = raw_pval)
###################################################################################################
