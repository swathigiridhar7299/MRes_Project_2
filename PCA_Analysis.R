#load libraries
library("factoextra")
library("FactoMineR")
library(ggpubr)
library(RColorBrewer)

data <- read.csv('file.csv') 
data <- na.omit(data)
pca.data <- PCA(data, scale.unit = TRUE, graph = FALSE)
fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))
b<-fviz_pca_var(pca.data, col.var = "cos2",
                gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
                repel = TRUE) 
a <- fviz_pca_ind(pca.data, col.ind = data$Treatment,
                  palette = "jco", addEllipses = FALSE, repel=FALSE, geom = "point")
PCA_Niraparib <- ggpar(a,
                       title = "Principal Component Analysis",
                       xlab = "PC1", ylab = "PC2",
                       legend.title = "Treatment", legend.position = "top",
                       ggtheme = theme_minimal())
PCA_Niraparib
