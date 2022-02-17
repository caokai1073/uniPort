library(Seurat)
library(SeuratDisk)
library(harmony)
library(rliger)
library(tidyverse)
library(SeuratWrappers)
library(EnsDb.Hsapiens.v75)
set.seed(1234)

file_path <- '/data/work/uniPort/spleen_ATAC/'
# load saved data from multiMAP tutorial
# RNA
spleen_count <- read.csv(paste0(file_path,'rna_count.txt'),header = FALSE)
name <- read.csv(paste0(file_path,'rna_count_name.txt'))
gene <- read.csv(paste0(file_path,'rna_count_gene.txt'))
spleen_count <- t(spleen_count)
rownames(spleen_count) <- gene[,1]
colnames(spleen_count) <- name[,1]
spleen_rna <- CreateSeuratObject(counts = spleen_count)

# ATAC
spleen_count <- read.csv(paste0(file_path,'atac_count.txt',header = FALSE))
name <- read.csv(paste0(file_path,'atac_count_name.txt'))
gene <- read.csv(paste0(file_path,'atac_count_gene.txt'))
spleen_count <- t(spleen_count)
rownames(spleen_count) <- gene[,1]
colnames(spleen_count) <- name[,1]
spleen_atac <- CreateSeuratObject(counts = spleen_count, assay = "activity")


# -----------------------------------------------
#       seurat integration
# -----------------------------------------------
spleen_rna <- NormalizeData(spleen_rna)
spleen_rna <- FindVariableFeatures(spleen_rna)

spleen_atac <- NormalizeData(spleen_atac)
spleen_atac <- ScaleData(spleen_atac, features = rownames(spleen_atac))

spleen_atac <- FindTopFeatures(spleen_atac, min.cutoff = 10)
spleen_atac <- RunTFIDF(spleen_atac)
spleen_atac <- RunSVD(spleen_atac)
spleen_atac <- RunUMAP(spleen_atac, reduction = "lsi", dims = 2:30)

transfer.anchors <- FindTransferAnchors(
  reference = spleen_rna, 
  query = spleen_atac, 
  features = VariableFeatures(spleen_rna),
  reference.assay = "RNA", 
  query.assay = "activity", 
  reduction = "cca")

genes.use <- VariableFeatures(spleen_rna)
refdata <- GetAssayData(spleen_rna,assay="RNA",slot="data")[genes.use, ]

imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = spleen_atac[["lsi"]],
  dims = 2:30)
spleen_atac[["RNA"]] <- imputation

meta <- read.csv(paste0(file_path,'rna_filter_meta.txt'),sep = '\t',header = FALSE)
spleen_rna$cluster <- meta$V1
meta <- read.csv(paste0(file_path,'atac_meta_combine.txt'),sep = '\t',header = FALSE)
spleen_atac$cluster <- meta$V1

spleen_atac$batch <- 'atac'
spleen_rna$batch <- 'rna'

coembed <- merge(x = spleen_rna, y = spleen_atac)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
DimPlot(coembed, group.by = c("batch", "cluster"))

saveRDS(coembed,file = paste0(file_path,'coembed_seurat.rds'))
write.table(coembed@reductions$umap@cell.embeddings,sep = '\t',quote = F,
            file = paste0(file_path,'seurat_integration_umap.txt'))


# -----------------------------------------------
#       harmony integration
# -----------------------------------------------
coembed_harmony <- coembed %>% 
  RunPCA(., features = genes.use, verbose = FALSE) %>% 
  RunHarmony(.,reduction = "pca",
             group.by.vars = "batch",
             plot_convergence = TRUE) %>% 
  RunUMAP(., reduction = "harmony", dims = 1:30)
DimPlot(coembed_harmony, group.by = c("batch", "cluster"),reduction = 'harmony')

write.table(coembed_harmony@reductions$harmony@cell.embeddings,sep = '\t',quote = F,
            file = paste0(file_path,'harmony_integration_umap.txt'))
            


# -----------------------------------------------
#       liger integration
# -----------------------------------------------
nFactors = 20   
coembed_liger <- RunOptimizeALS(coembed, k = nFactors, split.by = "batch")
coembed_liger <- RunQuantileNorm(coembed_liger, split.by = "batch")
coembed_liger <- FindNeighbors(coembed_liger, reduction = "iNMF", dims = 1:nFactors) %>% FindClusters()
coembed_liger <- RunUMAP(coembed_liger, dims = 1:nFactors, reduction = "iNMF")
DimPlot(coembed_liger, group.by = c("batch", "cluster"))

write.table(coembed_liger@reductions$iNMF@cell.embeddings, sep = '\t',quote = F,
            file = paste0(file_path,'~/kai/data/h5/liger_integration_umap.txt'))
            

