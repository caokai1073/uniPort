library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scatterpie)
library(RColorBrewer)
library(grDevices)
set.seed <- 1234

file_path <- '/data/work/uniPort/Brain_SPOT/'

#-----------------------------------------------------------
#    input data for uniPort
#-----------------------------------------------------------

# load st data from seurat
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# reference scRNA file from SPOTlight
path_to_data <- system.file(package = "SPOTlight")
cortex_sc <- readRDS(glue::glue("{path_to_data}/allen_cortex_dwn.rds"))

# save SPOT and RNA ref file for uniPort input
write.table(as.data.frame(brain@assays$Spatial@counts), 
            quote = F,row.names = T,sep = '\t',
            file = paste0(file_path,'SPOT_All_Count.txt'))

write.table(as.data.frame(cortex_sc@assays$RNA@counts), 
            quote = F,row.names = T,sep = '\t',
            file = paste0(file_path,'SPOT_Ref_RNA.txt'))
write.table(cortex_sc@meta.data['subclass'],quote = F,row.names = T,sep = '\t',
            file = paste0(file_path,'SPOT_Ref_RNA_cluster.txt'))


#-----------------------------------------------------------
#    analyze outputs of uniPort
#-----------------------------------------------------------
source(paste0(file_path,'spatial_function.R'))

ot <- read.table(paste0(file_path,'OT_Brain_SPOTRef_KL05_OT05.txt'),sep = '\t',header = T, row.names = 1)
ot_map <- mapCluster(ot,ref = cortex_sc, cluster = 'subclass')

# spatial scatter pie of cluster proportion
p <- stClusterPie(ot_map = ot_map, st = brain)

# single cluster proportion
stClusterExp(ot_map, brain, cluster = 'VLMC',cut = 0.25)

SpatialFeaturePlot(brain, features = c('Olig1','Atp1b2','Dcn','Olig2','Slc1a2','Osr1'),alpha = c(0.8,0.9)) 
















ot_map <- read.table('/data/yupines/kai/test/OT_Brain_LF_sub_KL0.5_OT0.5_T.txt',sep = '\t')
p <- stClusterPie(ot_map = ot_map, st = brain)

pdf('/data/yupines/kai/test/OT_Brain_LF_sub_KL0.5_OT0.5_T.pdf',width = 10,height = 7.5)
print(p)
dev.off()



