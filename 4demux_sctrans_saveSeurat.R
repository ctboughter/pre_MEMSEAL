library(Seurat)
library(dplyr)
library(hdf5r)
library(glmGamPoi)
library(openxlsx)
library(HGNChelper)
library(openxlsx)
library(ggplot2)
library(SeuratDisk)
library(stringr)
library(future)
options(future.globals.maxSize = 10*1e9)

#rm(list=ls())

labl_hsh <- read.xlsx('4inputs/labels_hashtags.xlsx')
res <- 2 #define the resolution
PC <- 50

for (j in 61:76) {
  print(j)
  # Obviously, set your own directory here
  setwd('C:/Users/budha/Documents/R/YF_day.0.1.3.5.7')
  D.data <- Read10X_h5(paste('./', labl_hsh[j,1], '.h5', sep = ""))
  k1 <- str_split(labl_hsh[j,4], " ") # pick up the hashes
  k2 <- as.numeric(k1[[1]])
  
  rna.counts <- D.data$`Gene Expression` 
  hto.counts <- D.data$`Antibody Capture`[k2 + 192, ]
  adt.counts <- D.data$`Antibody Capture`[c(2, 3, 5, 9, 12, 13, 15, 16, 17, 18, 19,
                                            20, 21, 24, 25, 27, 29, 32, 37, 39, 40, 42, 
                                            43, 44, 46, 49, 50, 51, 52, 53, 54, 56, 57, 59,
                                            60, 62, 64, 70, 73, 74, 79, 80, 85, 86, 87, 91, 100, 
                                            104, 107, 117, 118, 126, 127, 135, 136, 137, 138, 140, 
                                            142, 144, 151, 153, 154, 156, 157, 158, 159, 160, 162, 
                                            175, 177, 178, 179, 182), ] # ADT
  
  D.data.joint <- intersect(colnames(rna.counts), colnames(hto.counts))
  #D.data.joint <- intersect(colnames(rna.counts), colnames(adt.counts))
  rna.counts <- rna.counts[, D.data.joint]
  hto.counts <- hto.counts[, D.data.joint]
  adt.counts <- adt.counts[, D.data.joint]
  
  D <- CreateSeuratObject(counts = rna.counts)
  D[["percent.mt"]] <- PercentageFeatureSet(D, pattern = "^mt-")
  D[["HTO"]] <- CreateAssayObject(counts = hto.counts)
  D[["ADT"]] <- CreateAssayObject(counts = adt.counts)
  DefaultAssay(D) <- "RNA"
  
  D <- subset(D, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 2.5 & nCount_RNA < 20000) 
  D <- NormalizeData(D, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
  D <- NormalizeData(D, assay = "ADT", normalization.method = "CLR", verbose = FALSE)
  
  D <- HTODemux(D, assay = "HTO", positive.quantile = 0.999, verbose = TRUE)
  
  Idents(D) <- "HTO_classification.global"
  D.singlet <- subset(D, idents = "Singlet")
  
  D.sct <- SCTransform(D.singlet, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = PC, seed.use = 42, verbose = FALSE)
  
  setwd('C:/Users/budha/Documents/R/YF_day.0.1.3.5.7_Sobjects_batchcorr')
  SaveH5Seurat(D.sct, filename = paste0(labl_hsh[j,1], ".h5seurat"), overwrite = TRUE, verbose = FALSE)
  rm(D, D.sct, D.data, HTO.list, rna.counts, adt.counts, D.singlet, hto.counts)
}
