library(Signac)
library(Seurat)
library(SeuratWrappers)
library(EnsDb.Hsapiens.v75)
library(SeuratDisk)
set.seed(1234)

# download data
# # https://satijalab.org/signac/articles/pbmc_multiomic.html

file_path <- '/data/work/uniPort/signac_ATAC/'

# load the RNA and ATAC data
counts <- Read10X_h5(paste0(file_path,"pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"))
fragpath <- paste0(file_path,"pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg19"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# -------------------------------------------------------------
#                 Quality control
# -------------------------------------------------------------
DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc

gene.activities <- GeneActivity(pbmc)
write.table(as.data.frame(gene.activities),file = paste0(file_path,"gene_activity_signac.txt"), sep='\t',quote = F)

# -------------------------------------------------------------
#                    Peaks calling
# -------------------------------------------------------------
# call peaks using MACS2
DefaultAssay(pbmc) <- "ATAC"
peaks <- CallPeaks(pbmc, macs2.path = "~/anaconda3/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)
# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# -------------------------------------------------------------
#             Gene expression data processing
# -------------------------------------------------------------
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc,dims = 1:20)

# -------------------------------------------------------------
#             Gene expression data processing
# -------------------------------------------------------------
DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

# -------------------------------------------------------------
#                    Annotating cell types
# -------------------------------------------------------------
# load PBMC reference
reference <- LoadH5Seurat(paste0(file_path,"pbmc_multimodal.h5seurat"))
DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"
# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                  "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                  "NK Proliferating", "gdT",
                  "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                  "CD14 Mono", "CD16 Mono",
                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

# save file for uniPort input
meta <- pbmc@meta.data[,c("nCount_RNA","nFeature_RNA","nFeature_RNA","nFeature_ATAC",
                          "nucleosome_signal","nucleosome_percentile","predicted.id")]
write.table(meta,file = paste0(file_path,"meta_signac.txt"), sep='\t',quote = F)
write.table(pbmc@assays$RNA@counts,file = paste0(file_path,"RNA_Count_signac.txt"), sep='\t',quote = F)
saveRDS(pbmc, file = paste0(file_path,'pbmc_combined.rds'))


# -------------------------------------------------------------------------
#             ATAC analysis 
# -------------------------------------------------------------------------
pbmc <- readRDS(paste0(file_path,'pbmc_combined.rds'))
activity <- read.table(paste0(file_path,"gene_activity_signac.txt", sep='\t'))
colnames(activity) <- gsub('\\.','-',colnames(activity))
pbmc[["activity"]] <- CreateAssayObject(counts = activity)

DefaultAssay(pbmc) <- 'peaks'
pbmc_cd8 <- pbmc[,pbmc$cluster %in% c("CD8 Tmem","CD8 Naive","CD4 Naive","CD4 Tmem")]

pdf(paste0(file_path,'signac_atac_KLRC4-KLRK1_CoveragePlot.pdf'),width = 6,height = 5)
CoveragePlot(
  object = pbmc_cd8,
  region = c("chr12-10524952-10562745"),
  group.by = 'cluster',
  peaks = F,
  heights = 3.5)
dev.off()




# ---------------------------------------------------------------------------
#       Seurat, harmony, liger integration
# ---------------------------------------------------------------------------
pbmc <- readRDS(paste0(file_path,'pbmc_combined.rds'))

# seurat integration
pbmc_rna <- pbmc
DefaultAssay(pbmc_rna) <- "RNA"
pbmc_rna@assays$ATAC <- NULL
pbmc_rna <- RenameCells(pbmc_rna,new.names = unlist(lapply(colnames(pbmc_rna),function(x){paste0(x,'_RNA')})))

pbmc_atac <- pbmc
DefaultAssay(pbmc_atac) <- "ATAC"
pbmc_atac@assays$RNA <- NULL
pbmc_atac <- RenameCells(pbmc_atac,new.names = unlist(lapply(colnames(pbmc_atac),function(x){paste0(x,'_ATAC')})))

pbmc_rna <- NormalizeData(pbmc_rna)
pbmc_rna <- FindVariableFeatures(pbmc_rna)

gene.activities <- GeneActivity(pbmc_atac, features = VariableFeatures(pbmc_rna))
pbmc_atac[["activity"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(pbmc_atac) <- "activity"
pbmc_atac <- NormalizeData(pbmc_atac)
pbmc_atac <- ScaleData(pbmc_atac, features = rownames(pbmc_atac))
pbmc_atac <- FindTopFeatures(pbmc_atac, min.cutoff = 10)
pbmc_atac <- RunTFIDF(pbmc_atac)
pbmc_atac <- RunSVD(pbmc_atac)
pbmc_atac <- RunUMAP(pbmc_atac, reduction = "lsi", dims = 2:30)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna, 
  query = pbmc_atac, 
  features = VariableFeatures(pbmc_rna),
  reference.assay = "RNA", 
  query.assay = "activity", 
  reduction = "cca")
genes.use <- VariableFeatures(pbmc_rna)
refdata <- GetAssayData(pbmc_rna,assay="RNA",slot="data")[genes.use, ]
imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = pbmc_atac[["lsi"]],
  dims = 2:30)
pbmc_atac[["RNA"]] <- imputation

pbmc_atac$batch <- 'atac'
pbmc_rna$batch <- 'rna'

pbmc_coembed <- merge(x = pbmc_rna, y = pbmc_atac)
pbmc_coembed <- ScaleData(pbmc_coembed, features = genes.use, do.scale = FALSE)
pbmc_coembed <- RunPCA(pbmc_coembed, features = genes.use, verbose = FALSE)
pbmc_coembed <- RunUMAP(pbmc_coembed, dims = 1:30)

DimPlot(pbmc_coembed, group.by = c("batch", "cluster"))

saveRDS(pbmc_coembed,file = paste0(file_path,'pbmc_coembed_seurat.rds'))
write.table(pbmc_coembed@reductions$umap@cell.embeddings, sep = '\t',quote = F,
            file = paste0(file_path,'seurat_integration_umap.txt'))
            

# harmony integration
pbmc_coembed_harmony <- pbmc_coembed %>% 
  RunPCA(., features = genes.use, verbose = FALSE) %>% 
  RunHarmony(.,reduction = "pca",
             group.by.vars = "batch",
             plot_convergence = TRUE) %>% 
  RunUMAP(., reduction = "harmony", dims = 1:30)
DimPlot(pbmc_coembed_harmony, group.by = c("batch", "cluster"))

write.table(pbmc_coembed_harmony@reductions$harmony@cell.embeddings,sep = '\t',quote = F,
            file = paste0(file_path,'harmony_integration_umap.txt'))
            
# liger integration
nFactors = 20   
pbmc_coembed_liger <- RunOptimizeALS(pbmc_coembed, k = nFactors, split.by = "batch")
pbmc_coembed_liger <- RunQuantileNorm(pbmc_coembed_liger, split.by = "batch")
pbmc_coembed_liger <- FindNeighbors(pbmc_coembed_liger, reduction = "iNMF", dims = 1:nFactors) %>% FindClusters()
pbmc_coembed_liger <- RunUMAP(pbmc_coembed_liger, dims = 1:nFactors, reduction = "iNMF")
DimPlot(pbmc_coembed_liger, group.by = c("batch", "cluster"))

write.table(pbmc_coembed_liger@reductions$iNMF@cell.embeddings,sep = '\t',quote = F,
            file = paste0(file_path,'liger_integration_umap.txt'))



# ---------------------------------------------------------------------------
#       memory usage and run time
# ---------------------------------------------------------------------------
library(resample)
library(SeuratWrappers)
library(peakRAM)
set.seed(1234)

# 采样 5k, 10k, 20k, 40k, 80k, 160k, 320k
# load RNA and ACTIVITY data
rna <- read.table(paste0(file_path,'RNA_Count_signac.txt'), sep='\t')
pbmc.rna <- CreateSeuratObject(counts = rna,assay = "RNA")
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna,nfeatures = 10000)

atac <- read.table(paste0(file_path,"gene_activity_signac.txt"), sep='\t')
pbmc.atac <- CreateSeuratObject(counts = atac,assay = "ATAC")
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- FindVariableFeatures(pbmc.atac,nfeatures = 10000)

var.gene <- intersect(VariableFeatures(pbmc.rna),VariableFeatures(pbmc.atac))[1:2000]
VariableFeatures(pbmc.rna) <- var.gene
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna$batch <- 'RNA'

VariableFeatures(pbmc.atac) <- var.gene
pbmc.atac <- ScaleData(pbmc.atac)
pbmc.atac$batch <- 'ATAC'

# subset
k = 5000
rna.sub = pbmc.rna[,sample(colnames(pbmc.rna),size = k)]
print(object.size(rna.sub),units = "auto", standard = "IEC")
atac.sub = pbmc.atac[,sample(colnames(pbmc.atac),size = k)]
print(object.size(atac.sub),units = "auto", standard = "IEC")

# expand
k = 320000
rna.sub <- rna[,sample(col(rna),size = k,replace = TRUE)]
rna.sub <- CreateSeuratObject(counts = rna.sub,assay = "RNA")
rna.sub <- NormalizeData(rna.sub)
VariableFeatures(rna.sub) <- var.gene

atac.sub <- atac[,sample(col(atac),size = k,replace = TRUE)]
atac.sub <- CreateSeuratObject(counts = atac.sub,assay = "ATAC")
atac.sub <- NormalizeData(atac.sub)
VariableFeatures(atac.sub) <- var.gene

print(object.size(rna.sub),units = "auto", standard = "IEC")
print(object.size(atac.sub),units = "auto", standard = "IEC")

# seurat
peakRAM(transfer.anchors <- FindTransferAnchors(
  reference = rna.sub, query = atac.sub,
  reference.assay = "RNA", query.assay = "ATAC", 
  features = var.gene, reduction = "cca"))
print(object.size(transfer.anchors),units = "auto", standard = "IEC")


# harmony
atac.sub <- RenameAssays(atac.sub,ATAC='RNA')
rna.sub$batch <- 'RNA'
atac.sub$batch <- 'ATAC'
harmonyObj <- merge(rna.sub,atac.sub)
VariableFeatures(harmonyObj) <- var.gene
harmonyObj <- ScaleData(harmonyObj,split.by = "batch")
harmonyObj <- RunPCA(harmonyObj)
print(object.size(harmonyObj),units = "auto", standard = "IEC")
peakRAM(harmonyObj <- RunHarmony(harmonyObj,reduction = "pca",dims.use = 1:30,group.by.vars = "batch"))
print(object.size(harmonyObj),units = "auto", standard = "IEC")


# Liger
combined_liger <- harmonyObj
combined_liger <- ScaleData(combined_liger, split.by = "batch", do.center = FALSE)
print(object.size(combined_liger),units = "auto", standard = "IEC")
peakRAM(combined_liger <- RunOptimizeALS(combined_liger, k=20, split.by="batch"))
print(object.size(combined_liger),units = "auto", standard = "IEC")
pbmc.atac <- RenameAssays(pbmc.atac,ATAC='RNA')

nFactors=20   
combined_liger <- merge(pbmc.rna,pbmc.atac)
DefaultAssay(combined_liger) <- 'RNA'
combined_liger <- NormalizeData(combined_liger)
combined_liger <- ScaleData(combined_liger, split.by = "batch", do.center = FALSE)
combined_liger <- RunOptimizeALS(combined_liger, k=nFactors, split.by="batch",lambda = 10)
combined_liger <- RunQuantileNorm(combined_liger, split.by="batch")
combined_liger <- FindNeighbors(combined_liger, reduction = "iNMF", dims = 1:20)
combined_liger <- RunUMAP(combined_liger, dims=1:nFactors, reduction="iNMF")
DimPlot(combined_liger,label = T,split.by = 'batch')

