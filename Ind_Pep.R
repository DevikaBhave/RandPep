# Clear everything in environment
rm(list = ls())
######################################   PACKAGES   ###############################################
# For updating 
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install()
# Call packages from library for every new session
library("BiocManager")
# BiocManager::install("limma")
library(limma)
# install.packages("statmod")
library(statmod)
# install.packages("dplyr")
library(dplyr)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("gplots")
library(gplots)
# install.packages("openxlsx")
library(openxlsx)
library(cowplot)
# Set working directory 
setwd("~/ownCloud/03_Data_ownCloud/01_Microarray/Expt2_2020/Expt2_Aug2021_Analysis/NEG_AND_STOP/")

###################################################################################################
# Load text file with the corresponding names of files, samples, treatments, etc.
targets <- readTargets(file = "targets_ind_NEG.txt", sep = "\t", row.names = "No.")
# TARGET <- targets1[, "Target"]
# treatment <- targets[,"Treatment"]
# timepoints <- targets[,"Timepoint"]
# read the corresponding expression data containing (.txt) files 
z <- read.maimages(targets, 
                   path = "~/ownCloud/03_Data_ownCloud/01_Microarray/Expt2_2020/Target_Files/",
                   source = "agilent", green.only = TRUE, other.columns = "gIsWellAboveBG")
dim(z)
# ~~~~~ STEP 1 sanity check ~~~~~ extract file here ~~~~~
# write.xlsx(z, file = "Sanity_check_files/raw_export_sanitycheck1.xlsx")
# Layout <- getLayout(z$genes)
summary(z$E)
# test by plotting a scatter of the entire data
plotMD(z)
# plots the sample names to see array biases
plotMDS(z)
# Generates one scatter plot for each sample in a 3x2 matriz as a png file
# plotMA3by2(z)
plotDensities(z, legend = F)
# corrects the background in each array by normalizing the array to the dark well 
# normexp or edwards is recommended (see documentation)
q <- backgroundCorrect(z, method = "normexp")
# ~~~~~ STEP 2 sanity check ~~~~~ extract file here ~~~~~
# write.xlsx(q, file = "Sanity_check_files/bg_normexp_sanitycheck2.xlsx")
# Normalize expression intensities so the log-ratios have similar distributions across set of arrays
# Using the default for single channel "quantile"
q <- normalizeBetweenArrays(q, method = "quantile")
dim(q)
# ~~~~~ STEP 3 sanity check ~~~~~ extract file here ~~~~~
# write.xlsx(q, file = "Sanity_check_files/normbetarrays_sanitycheck3.xlsx")
# average computed from replicate probes in microarray (use 'GeneName' instead of 'ProbeName')
# q.avg <- avereps(q, ID = q$genes$ProbeName)
avg <- avereps(q, ID = q$genes$GeneName)
dim(avg)
# remove all the gene rows that have no value or are the controls
Control <- avg$genes$ControlType==1
Control1 <- avg$genes$ControlType==-1
# remove all the gene rows that have no gene name assigned
NoSymbol <- is.na(avg$genes$GeneName)
# keep all the gene rows that have an expression value >= 2 
IsExpr <- rowSums(avg$other$gIsWellAboveBG > 0) >=2
# new normalized elist
filt <- avg[!Control & !Control1 & !NoSymbol & IsExpr, ]
dim(filt)
head(filt)
names(filt$genes)
# ~~~~~ STEP 4 sanity check ~~~~~ extract file here ~~~~~
# write.xlsx(filt, file = "Sanity_check_files/sanitycheck4_Filter_data.xlsx")
# keep the columns you need and discard the rest
filt$genes <- filt$genes[, c("ProbeName", "GeneName", "ProbeUID", "Description")]
names(filt$genes)
head(filt$genes)
# ~~~~~ STEP 5 sanity check ~~~~~ extract file here ~~~~~
# write.xlsx(filt, file = "Sanity_check_files/Filter_cols_sanitycheck5.xlsx")
# average across technical replicates of samples -> only used for PCA plotting!
averaged <- avearrays(filt, ID = filt$targets$Levels)
avg_expr <- as.data.frame(averaged)
avg_expr <- select(avg_expr, -c(2,3))
# write.csv(avg_expr, "avg_expr_for_PCA.csv")
# scatter plot for normalized y and filtered yfilt
limma::plotMA(filt)
# plotMA3by2(qfilt)
plotDensities(filt, legend = F)
level<- c("pFLAG_i","NEG_Pep1_i","NEG_Pep1_STOP_i","NEG_Pep2_i","NEG_Pep2_STOP_i",
             "NEG_Pep3_i","NEG_Pep3_STOP_i","NEG_Pep4_i","NEG_Pep4_STOP_i","NEG_Pep5_i",
             "NEG_Pep5_STOP_i","NEG_Pep6_i","NEG_Pep6_STOP_i")
f <- factor(targets$Levels, levels = level)
design <- model.matrix(~f)
colnames(design) <- level
fit <- lmFit(filt, design)
Summary_table_fit <- summary(decideTests(fit))
# The first column is always the total no. of genes (input in the 'Up' row). 
# Subsequent columns contain the genes that are 'Down', 'NotSig' and 'Up' compared 
# to the given control
# write.xlsx(Summary_table_fit, file = "Summary.xlsx")
fit_ebayes <- eBayes(fit)
# ~~~~~ STEP 6 sanity check ~~~~~ extract file here ~~~~~
# write.xlsx(fit_ebayes, file = "fit_ebayes_sanitycheck5.xlsx")
Grand_table <- topTable(fit_ebayes, coef = NULL, 
                        sort.by = "F",
                        number = 50,
                        # number = nrow(fit_ebayes$genes), 
                        adjust.method = "BH")
# write.xlsx(Grand_table, file = "Grand_table_NEG_top50.xlsx")

###########################    Correction of batch effects    ###################################
# create batch_data file from Grand table file, such that samples as columns and logFC as rows
# delete rest of the columns including the (pFLAG) first sample column with avg expression
data_fc <- read.table("grand_batch.txt", header = T, sep = "\t")
batches <- read.table("batches.txt", header = T, sep = "\t")
batch_corr <- removeBatchEffect(data_fc, batch = batches$Batch)
# write.xlsx(batch_corr, file = "batch_corrected_NEG_i.xlsx")
# par(mfrow=c(1,1))
boxplot(as.data.frame(data_fc),main = "original")
boxplot(as.data.frame(batch_corr), main = "Batch corrected")

################################    Individual comparisons    ###################################
contrast_individual <- makeContrasts(
  "pFLAG_i",
  "NEG_Pep1_i",
  levels = design
)
fit_individual <- contrasts.fit(fit, contrast_individual)
fit_individual <- eBayes((fit_individual))
Toptable_individual <- topTable(fit_individual, coef = "NEG_Pep1_i", 
                                number = nrow(fit_individual$genes), 
                                adjust.method = "BH")
# write.xlsx(Toptable_individual, file = "Contrasts/topTable_POS_Pep3_STOP_i.xlsx")
# results_individual <- decideTests(fit_individual, method = "global", 
#                                   adjust.method = "BH", p.value = 0.05)
# summary(decideTests(fit_individual))
#vennDiagram(results_individual, circle.col = 2, cex = 0.9)
#limma::plotMA(fit_individual)
#################################################################################################

