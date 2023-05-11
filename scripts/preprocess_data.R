library(data.table)
library(dplyr)
library(Seurat)
library(SeuratWrappers) #remotes::install_github("satijalab/seurat-wrappers")
library(SeuratDisk)
save_loc <- "/Users/almas/Desktop/schwartz_lab/almas/artathon_analysis/"
ftp <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813_scc_scRNA_counts.txt.gz"
counts <- data.table::fread(ftp,data.table = F)

rownames(counts) <- counts$V1
counts$V1 <- NULL

obj <- Seurat::CreateSeuratObject(counts,row.names=rownames(counts))
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
print(unique(obj@meta.data$orig.ident)) # 11 different groups for batch effects
obj <- NormalizeData(obj)
# Let's try to account for batch effects 
obj <- ScaleData(obj, do.center = FALSE)

obj <- RunOptimizeALS(obj, k = 20, lambda = 5, split.by = "orig.ident")

obj <- SeuratWrappers::RunQuantileNorm(obj, split.by = "orig.ident") 
# Seurat Objects have rows as gene names and cells as columns
# When converting to an anndata it should flip the rows and columns using the Convert function
# Lets convert it to anndata and also save it as an h5Seurat if needed for the future
SaveH5Seurat(obj, filename = paste0(save_loc,"tcell_data.h5Seurat"))
Convert(paste0(save_loc,"tcell_data.h5Seurat"), dest = paste0(save_loc,"h5ad"))

