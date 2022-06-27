library(scatterpie)
library(RColorBrewer)
library(grDevices)
library(Seurat)
library(tidyverse)
library(reshape2)
set.seed <- 1234

file_path <- '/data/work/uniPort/BRCA_SPOT/'
setwd(file_path)

file_path <- '/data/yupines/kai/data/brca/'

#-----------------------------------------------------------
#    input data for uniPort
#-----------------------------------------------------------
# load st data
brca <- Load10X_Spatial(paste0(file_path,'st/'))
brca <- NormalizeData(brca)
brca <- ScaleData(brca)
# load scRNA data
brca_sc <- CreateSeuratObject(counts = Read10X(paste0(file_path,'sc/outs/')))
brca_sc <- NormalizeData(brca_sc)
brca_sc <- ScaleData(brca_sc)
# load scRNA cluster file
brca_cluster <- read.csv(paste0(file_path,'sc/Whole_miniatlas_meta.csv'), header = T,row.names = 1) %>% .[-1,]

# save SPOT and RNA ref file
write.table(as.data.frame(brca_st@assays$Spatial@counts), quote = F,row.names = T,sep = '\t',
            file = paste0(file_path,'BRCA_SPOT_Count.txt'))
# convert dgcMatrix to dataframe
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
tmp <- as_matrix(brca_sc@assays$RNA@counts)
write.table(tmp, quote = F,row.names = T,sep = '\t',file = paste0(file_path,'BRCA_SPOT_Ref_RNA.txt'))


#-----------------------------------------------------------
#    analyze outputs of uniPort
#-----------------------------------------------------------
source(paste0(file_path,'spatial_function.R'))

source('/data/yupines/kai/upload/spot/spatial_function.R')

ot <- read.table(paste0(file_path,'OT_BRCA_total_RNAref_KL05_OT05.txt'),sep = '\t',header = T, row.names = 1)
ot <- as.data.frame(t(ot))
ot_map <- mapCluster(ot,meta = brca_cluster, cluster = 'celltype_major')

# spatial scatter pie of cluster proportion
p <- stClusterPie(ot_map = ot_map, st = brca)

# single cluster proportion
stClusterExp(ot_map, brca, cluster = 'CAFs',cut = 0.25)
SpatialFeaturePlot(brca, features = c('CD3D','ERBB2'),alpha = c(0.8,0.9)) 












ot_map <- read.table('/data/yupines/kai/test/BRCA_OT_map_LT_T_i60k.txt',sep = '\t')
p <- stClusterPie(ot_map = ot_map, st = brca)

pdf('/data/yupines/kai/test/BRCA_celltype_LT_T_i60k.txt.pdf',width = 10,height = 7.5)
print(p)
dev.off()


