#######################################################################################################
#PCA analyses of the Microarray data
library(ggplot2)
#install.packages("ggfortify")
library(ggfortify)
#install.packages("devtools")
library('devtools')
#install_github("vqv/ggbiplot", force = TRUE)
library('ggbiplot')
library(readxl)
setwd("~/ownCloud/03_Data_ownCloud/01_Microarray/Expt2_2020/Expt2_Aug2021_Analysis/NEG_AND_STOP/")
#######################################################################################################
names <- read.table("unique.txt", header = T, sep = "\t")
#Expr<- t(qfilt$E)
aver_expr2 <- t(averaged$E)
str(aver_expr2)
# Expr <- read.table("~/Desktop/averaged.txt", header = T, sep = "\t")

#autoplot(prcomp(t(qfilt$E)))
data.pca <- prcomp(aver_expr2, scale. = T, center = T)
autoplot(data.pca, data = names, colour = 'Categories')
a <- 
  ggbiplot(data.pca, choices = c(1,2),
           obs.scale = 1,
           var.scale = 1,
           var.axes = F, 
           #groups = names$Treatment,
           #labels = names$Sample, labels.size = 3,
           ellipse = F, 
           circle = F) + 
  geom_point(aes(colour = names$Categories),
             size=3, alpha=1)+
  scale_colour_manual(values = c("#85305b", "#fb7363", "#000000"))+
  theme_bw(base_size = 14)+
  theme(axis.title = element_text(face = "bold", colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.text = element_text(colour ="black"),
        legend.title = element_blank())
a
# ggsave(path = ".", "PCA.pdf", device = "pdf", width = 7, height = 5)
#######################################################################################################
b <- 
  ggbiplot(data.pca, choices = c(1,2),
           obs.scale = 1,
           var.scale = 1,
           var.axes = F, 
           groups = names$Categories,
           labels = names$Sample, labels.size = 2.5,
           ellipse = F, 
           circle = F) +
 #ylim(c(-80,75))+
  scale_colour_manual(values = c("#85305b", "#fb7363", "#000000"))+
  theme_bw(base_size = 14)+
  theme(axis.title = element_text(face = "bold", colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(face = "bold"),
        #legend.text = element_text(colour ="black"),
        legend.title = element_blank())
b 
# ggsave(path = ".", "PCA_named.pdf", device = "pdf", width = 9, height = 6)
#######################################################################################################
screeplot(data.pca)
summary_pca <- summary(data.pca)
head(summary_pca)
summary_pca$sdev^2
screeplot(summary_pca, type = "lines")
#######################################################################################################
# Hex colors
# Neutrals - Induced = (#877458"), Uninduced = ("#cdc3a5")
# Conrol-pFLAG - Induced = ("#000000"), Uninduced = ("#8c8c8c")
# POS_Peps - Induced = ("#3b367e"), Uninduced = ("#9d9abe")
# Negs - induced = "#63123b", uninduced = "#b67f9b"
# NEG_STOPs - induced = "#fb7363"
