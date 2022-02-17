library(Seurat)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(uwot)
library(rliger)
library(harmony)
library(SeuratDisk)
library(SeuratWrappers)
options(stringsAsFactors = FALSE)
set.seed(1234)

# download data
# scRNA https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576
# MERFISH https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248

file_path <- '/data/work/uniPort/MERFISH/'

# load file 
# --------------------------------------------------------
#         scRNA data
# --------------------------------------------------------
scRNA <- Read10X(file_path) %>%
  CreateSeuratObject(.,min.cells = 2, min.features = 200) %>%
  magrittr::inset2("percent.mt",value=PercentageFeatureSet(.,pattern = "^MT-")) %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(.,selection.method = "vst", nfeatures = 2000)  %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) 

sc_cluster <- read.csv(paste0(file_path,'MERFISH_scRNA_Cluster.csv'))
colnames(sc_cluster) <- c("Cell.name","Sex","Replicate.number","Cell.class","Sub","Neuron")
sc_cluster$cluster_sub <- sc_cluster$Cell.class %>% {
  .[.== "Fibroblast"] <- sc_cluster$Sub[.== "Fibroblast"]
  .[.== "Mature oligodendrocyte"] <- sc_cluster$Sub[.== "Mature oligodendrocyte"]
  .[.== "Endothelial"] <- sc_cluster$Sub[.== "Endothelial"]
  .[.== "Microglia"] <- sc_cluster$Sub[.== "Microglia"]
  .[.== "Astrocytes"] <- sc_cluster$Sub[.== "Astrocytes"]
  .[.== "Ependymal"] <- sc_cluster$Sub[.== "Ependymal"]
  .[.== "Immature oligodendrocyte"] <- sc_cluster$Sub[.== "Immature oligodendrocyte"]
  .[.== "Newly formed oligodendrocyte"] <- sc_cluster$Sub[.== "Newly formed oligodendrocyte"]
  .[.== "Mural"] <- sc_cluster$Sub[.== "Mural"]
  ;.}

sc_cluster$cluster_main <- sc_cluster$Cell.class
sc_cluster$cluster_main <- gsub('Mature oligodendrocyte','OD Mature',sc_cluster$cluster_main)
sc_cluster$cluster_main <- gsub('Immature oligodendrocyte','OD Immature',sc_cluster$cluster_main)
sc_cluster$cluster_main <- gsub('Newly formed oligodendrocyte','Newly formed OD',sc_cluster$cluster_main)
sc_cluster$cluster_main <- gsub('Astrocytes','Astrocyte',sc_cluster$cluster_main)
sc_cluster$cluster_sub <- gsub('Mature oligodendrocyte','OD Mature',sc_cluster$cluster_sub)
sc_cluster$cluster_sub <- gsub('Immature oligodendrocyte','OD Immature',sc_cluster$cluster_sub)
sc_cluster$cluster_sub <- gsub('Newly formed oligodendrocyte','Newly formed OD',sc_cluster$cluster_sub)
sc_cluster$cluster_sub <- gsub('Astrocytes','Astrocyte',sc_cluster$cluster_sub)

rownames(sc_cluster) <- sc_cluster$Cell.name
sc_cluster <- sc_cluster[,c("cluster_main", "cluster_sub")]
meta <- scRNA@meta.data %>% cbind(.,sc_cluster)
scRNA@meta.data <- meta

# delete "Ambiguous" and "Unstable" cluster
scf <- scRNA[,!scRNA$cluster_main %in% c("Ambiguous","Unstable")]
scf_cluster <- sc_cluster[!sc_cluster$cluster_main %in% c("Ambiguous","Unstable"),]


# --------------------------------------------------------
#         MERFISH data
# --------------------------------------------------------
st_raw <- read.csv(paste0(file_path,'Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv'))
st_raw$cluster_main <- st_raw$Cell_class
st_raw$cluster_main[grep('OD Immature',st_raw$cluster_main)] <- 'OD Immature'
st_raw$cluster_main[grep('OD Mature',st_raw$cluster_main)] <- 'OD Mature'
st_raw$cluster_main[grep('Endothelial',st_raw$cluster_main)] <- 'Endothelial'
st_raw$cluster_main[grep('Endothelial',st_raw$cluster_main)] <- 'Endothelial'
st_raw$cluster_sub <- st_raw$Cell_class

# choose mouse1 expr as the training set
st <- st_raw[
  st_raw$Animal_ID == 1,
  c('cluster_main','cluster_sub',"Centroid_X","Centroid_Y",colnames(st_raw)[10:170])] %>% 
  .[,!colnames(.) %in% c('Blank_1','Blank_2','Blank_3','Blank_4','Blank_5')]
st <- st[,!colnames(st)=='Fos'] # expression is NaN 
rownames(st) <- paste0('cell',1:nrow(st))
st_cluster <- st[,1:2]

# delete "Ambiguous" cluster
stf <- st[!st$cluster_main == "Ambiguous",] 
stf_cluster <- stf[,1:2]
cord <- stf[,3:4] # coord of MERFISH cells
st <- st[,5:159]
stf <- stf[,5:159]
st <- as.data.frame(t(st))
stf <- as.data.frame(t(stf))

st_umap_all <- uwot::umap(t(stf))
plot(st_umap_all,pch=16,asp = 1,cex = 0.1,
     xlab = "UMAP_1",ylab = "UMAP_2")

# --------------------------------------------------------
#         save data for uniPort input
# --------------------------------------------------------
# save scRNA file
# filtered files used
write.table(as.data.frame(scRNA@assays$RNA@counts),sep = '\t',quote = F, file = paste0(file_path,'MERFISH_scRNA.txt'))
write.table(sc_cluster,sep = '\t',quote = F, file = paste0(file_path,'MERFISH_scRNA_cluster.txt'))
write.table(as.data.frame(scf@assays$RNA@counts),sep = '\t',quote = F,file = paste0(file_path,'MERFISH_scRNA_filter.txt'))
write.table(scf_cluster,quote = F,row.names = T, sep = '\t', file = paste0(file_path,'MERFISH_scRNA_filter_cluster.txt'))

# save ST file
# filtered files used
write.table(st, quote = F,row.names = T,sep = '\t', file = paste0(file_path,'HPRI/MERFISH_st.txt'))
write.table(st_cluster,quote = F,row.names = T,sep = '\t', file = paste0(file_path,'MERFISH_st_cluster.txt'))
write.table(stf, quote = F,row.names = T,sep = '\t', file = paste0(file_path,'MERFISH_st_filter.txt'))
write.table(stf_cluster,quote = F,row.names = T,sep = '\t', file = paste0(file_path,'MERFISH_st_filter_cluster.txt'))


# -----------------------------------------------------------------------
#    USE mouse2 MERFISH expr for prediction training 
# -----------------------------------------------------------------------
stm2 <- st_raw[
  st_raw$Animal_ID == 2,
  c('cluster_main','cluster_sub',"Centroid_X","Centroid_Y",colnames(st_raw)[10:170])] %>% 
  .[,!colnames(.) %in% c('Blank_1','Blank_2','Blank_3','Blank_4','Blank_5')]
stm2 <- stm2[,!colnames(stm2)=='Fos'] # NaN
rownames(stm2) <- paste0('cell',1:nrow(stm2))
stm2_cluster <- stm2[,1:2]
stm2f <- stm2[!stm2$cluster_main == "Ambiguous",] 
stm2f_cluster <- stm2f[,1:2]
stm2f <- stm2f[,5:159]
stm2f <- as.data.frame(t(stm2f))
# save stm2 file
write.table(stm2f, quote = F,row.names = T,sep = '\t', file = paste0(file_path,'MERFISH_mouse2_filter.txt'))
write.table(stm2f_cluster,quote = F,row.names = T,sep = '\t', file = paste0(file_path,'MERFISH_mouse2_filter_cluster.txt'))


# -----------------------------------------------------------------------
#   predict unmeasured gene expression mouse2 cells (outs)
# -----------------------------------------------------------------------
# load outs of predicited mouse2 expr
pred <- read.table(paste0(file_path,'Merfish_predict_RNA_mouse2.txt'))
umap <- read.table(paste0(file_path,'Merfish_umap_mouse2.txt')) %>% {colnames(.) <- c('X','Y');.}
umap$cell <- rownames(umap)
pred <- as.data.frame(t(pred))
st_cluster <- read.table(paste0(file_path,'MERFISH_mouse2_filter_cluster.txt',sep = '\t')) 
st_cluster$cell <- rownames(st_cluster)

# load mouse2 real expr
sc <- as.data.frame(read.table(paste0(file_path,'mouse2_scRNA_norm_count.txt',sep = ',')))
sc_gene <- read.table(paste0(file_path,'mouse2_genename.txt',header = T))
sc_cluster <- read.table(paste0(file_path,'mouse2_scRNA_cluster.txt',row.names = 1,header = T,sep = '\t'))
sc_cluster$cell <- rownames(sc_cluster)
colnames(sc) <- sc_gene[,1]
rownames(sc) <- rownames(sc_cluster)

cluster_list <- intersect(sc_cluster$cell_type,st_cluster$cluster_main)
sc <- sc[rownames(sc_cluster[sc_cluster$cell_type %in% cluster_list,]),]
st <- pred[rownames(st_cluster[st_cluster$cluster_main %in% cluster_list,]),]

df_st <- lapply(cluster_list,function(x){
  df_st <- as.data.frame(colMeans(st[rownames(st_cluster[st_cluster$cluster_main %in% x,]),])) %>% {colnames(.) <- x;.}
  return(df_st)
}) %>% do.call('cbind',.)
df_sc <- lapply(cluster_list,function(x){
  df_sc <- as.data.frame(colMeans(sc[rownames(sc_cluster[sc_cluster$cell_type %in% x,]),])) %>% {colnames(.) <- x;.}
  return(df_sc)
}) %>% do.call('cbind',.)
df_st$`All Gene Means` <- colMeans(st)
df_sc$`All Gene Means` <- colMeans(sc)

cor_plot_list <- lapply(colnames(df_st),function(x){
  dt <- data.frame('Real'= df_sc[,x],'Predict'= df_st[,x])
  fmt <- " ~ R^2 == %.3f * ',' ~ {P == %.3f}"
  sum <- summary(lm(dt[,'Real'] ~ dt[,'Predict'],dt))
  lab <- sprintf(fmt, sum$r.squared, coef(sum)[2, 4])
  p <- ggscatterhist(
    dt,  y = 'Predict', x = 'Real',  
    shape=19, color = '#9970ab', 
    size =2, alpha = 0.8,
    title = x,
    margin.plot.size = 1.2,
    margin.plot =  "histogram",
    margin.params = list(
      color = 'black',  
      fill = '#1f78b4'),
    xlim = c(min(dt[,'Real']),max(dt[,'Real'])),
    ylim = c(min(dt[,'Predict']),max(dt[,'Predict'])),
    margin.ggtheme = theme_void(),
    ggtheme = theme(
      axis.title = element_text(size = 15),
      panel.background = element_rect(fill = "transparent",size = 50), 
      plot.background = element_rect(fill = "transparent"),
      plot.title = element_text(hjust = .5, face = "bold", size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=0.8)
    )) %>% 
    {.$sp <- .$sp + 
      geom_hline(yintercept = median(dt[,2]), linetype = "dashed", color = "grey") +
      geom_vline(xintercept = median(dt[,1]), linetype = "dashed", color = "grey") +
      annotate("text", parse=TRUE, size = 4, label=lab,
               x= max(dt['Predict'])*0.2, 
               y= max(dt['Real'])*0.8) +
      stat_smooth(method=lm);.}
  return(p)
})

pdf(paste0(file_path,'mouse2_real_predict_RNA.pdf',width = 5,height = 5))
print(list(cor_plot_list[[1]],cor_plot_list[[2]],cor_plot_list[[3]],
           cor_plot_list[[4]],cor_plot_list[[5]],cor_plot_list[[6]],
           cor_plot_list[[7]],cor_plot_list[[8]],cor_plot_list[[9]]))
dev.off()



# -----------------------------------------------------------------------
#   seurat, harmony, liger integration (use mouse1)
# -----------------------------------------------------------------------
# filtered clusters used

# seurat CCA integration (all genes)
com_gene <- intersect(rownames(scf),rownames(stf)) 
sc <- as.data.frame(scf@assays$RNA@counts)
sc <- sc[com_gene,]
hp <- stf[com_gene,] 
co <- as.matrix(cbind(sc,hp))
co_cluster <- rbind(scf_cluster,stf_cluster)
co_cluster$batch <- ifelse(rownames(co_cluster) %in% colnames(sc),"SC",'ST')

hps <- CreateSeuratObject(counts = stf) 
hps <- NormalizeData(hps)
hps <- ScaleData(hps)
hps@meta.data <- stf_cluster
VariableFeatures(hps) <- rownames(hps)

cca <- RunCCA(scf,hps)
cca@meta.data <- co_cluster
cca <- RunUMAP(cca, reduction = "cca",dims = 1:30)
DimPlot(cca,group.by = 'cluster_main',split.by = 'batch',label = T)
DimPlot(cca, split.by = 'cluster_main',group.by = 'batch',label = T,ncol = 5)

write.table(
  as.data.frame(cca@reductions$cca@cell.embeddings),
  quote = F,row.names = T,sep = '\t',col.names = T,
  file = paste0(file_path,'MERFISH_CCA_embeddings_dim.txt'))

# seurat integration
ifnb.list <- list('SC'=scf,'ST'=hps)
features <- SelectIntegrationFeatures(object.list = ifnb.list)
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features) 
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

combined@meta.data <- co_cluster
DimPlot(combined,group.by = 'cluster_main',split.by = 'batch',label = T)

write.table(
  as.data.frame(combined@reductions$umap@cell.embeddings),
  quote = F,row.names = T,sep = '\t',col.names = T,
  file = paste0(file_path,'MERFISH_seurat_integration_embeddings.txt'))

# harmony all genes
harmonyObj_all <- FindVariableFeatures(cca)
harmonyObj_all <- RunPCA(harmonyObj_all)
harmonyObj_all <- RunHarmony(
  harmonyObj_all,
  reduction = "pca",
  dims.use = 1:30,
  group.by.vars = "batch",
  plot_convergence = TRUE)
harmonyObj_all <- RunUMAP(harmonyObj_all, reduction = "harmony", dims = 1:30)
DimPlot(harmonyObj_all,group.by = 'cluster_main',split.by = 'batch',label = T)

write.table(
  as.data.frame(harmonyObj@reductions$harmony@cell.embeddings),
  quote = F,row.names = T,sep = '\t',col.names = T,
  file = paste0(file_path,'MERFISH_harmony_embeddings_all_gene.txt'))


# liger integration
nFactors=20   
hps$batch <- 'ST'
scf$batch <- 'SC'
combined_liger <- merge(scf,hps)
DefaultAssay(combined_liger) <- 'RNA'
combined_liger <- NormalizeData(combined_liger)
combined_liger <- FindVariableFeatures(combined_liger)
combined_liger <- ScaleData(combined_liger, split.by = "batch", do.center = FALSE)
combined_liger <- RunOptimizeALS(combined_liger, k=nFactors, split.by="batch",lambda = 20)
combined_liger <- RunQuantileNorm(combined_liger, split.by="batch")
combined_liger <- FindNeighbors(combined_liger, reduction = "iNMF", dims = 1:20)
combined_liger <- RunUMAP(combined_liger, dims=1:nFactors, reduction="iNMF")
DimPlot(combined_liger, group.by = c("batch", "cluster_main"))

write.table(combined_liger@reductions$iNMF@cell.embeddings,sep = '\t',quote = F,
            file = paste0(file_path,'MERFISH_liger_iNMF.txt'))

