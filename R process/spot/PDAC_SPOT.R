library(scatterpie)
library(RColorBrewer)
library(grDevices)
library(Seurat)
library(tibble)
library(data.table)
set.seed(1234)
library(dplyr)
library(stringr)

file_path = '/Users/cao/Downloads/uniPort-main-2/R process/spot/'
#-----------------------------------------------------------
#    input data for uniPort
#-----------------------------------------------------------
# load st expression matrix
dataA = fread(paste0(file_path,"GSM3036911_PDAC-A-ST1-filtered.txt.gz"), header =T,check.names = F)
dataA = as.data.frame(dataA)
dataA = dataA %>% distinct(Genes,.keep_all = T) %>% column_to_rownames("Genes")

# load paired scRNA data
scdataA = fread(paste0(file_path,'GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz'),header = T) 
scdataA = as.data.frame(scdataA)
scdataA = scdataA[!duplicated(scdataA$Genes),]
rownames(scdataA) <- scdataA$Genes
scdataA <- scdataA[,-1]

names = colnames(scdataA)[1:ncol(scdataA)] %>% as.data.frame() %>% {colnames(.) <- 'raw_type';.}
names$cell = paste0('cell',1:ncol(scdataA))
names$cell_type = names$raw_type  
names$cell_type[str_detect(names$cell_type,'Ductal')] = 'Ductal'
names$cell_type[str_detect(names$cell_type,'Acinar cells')] = 'Acinar cells'
names$cell_type[str_detect(names$cell_type,'Cancer clone A')] = 'Cancer clone A'
names$cell_type[str_detect(names$cell_type,'Cancer clone B')] = 'Cancer clone B'
names$cell_type[str_detect(names$cell_type,'mDCs')] = 'mDCs'
names$cell_type[str_detect(names$cell_type,'Tuft cells')] = 'Tuft cells'
names$cell_type[str_detect(names$cell_type,'pDCs')] = 'pDCs'
names$cell_type[str_detect(names$cell_type,'Endocrine cells')] = 'Endocrine cells'
names$cell_type[str_detect(names$cell_type,'Endothelial cells')] = 'Endothelial cells'
names$cell_type[str_detect(names$cell_type,'Macrophages')] = 'Macrophages'
names$cell_type[str_detect(names$cell_type,'Mast cells')] = 'Mast cells'
names$cell_type[str_detect(names$cell_type,'T cells & NK cells')] = 'T & NK cells'
names$cell_type[str_detect(names$cell_type,'Monocytes')] = 'Monocytes'
names$cell_type[str_detect(names$cell_type,'RBCs')] = 'RBCs'
names$cell_type[str_detect(names$cell_type,'Fibroblasts')] = 'Fibroblasts'
colnames(scdataA) = paste0('cell',1:ncol(scdataA))

# save SPOT and RNA ref file for uniPort input
write.table(scdataA, quote = F,row.names = T, sep = '\t', file = paste0(file_path,'PDAC_scRNA.txt'))
write.table(dataA, quote = F,row.names = T, sep = '\t', file = paste0(file_path,'PDAC_ST.txt'))
write.table(names, quote = F,row.names = T, sep = '\t', file = paste0(file_path,'PDAC_scRNA_label.txt'))

# get coord from st data
ind <- as.data.frame(t(sapply(
  str_split(colnames(dataA), "x"), 
  function(x){
    x <- as.numeric(x)
    x <- as.vector(x)
  }))) %>% {
    names(.) <- c("row_ind", "col_ind")
    rownames(.) <- paste0(.$row_ind,"x",.$col_ind)
    rownames(.) <- paste0('X',rownames(.))
 ;.}


#-----------------------------------------------------------
#    analyze outputs of uniPort
#-----------------------------------------------------------
source(paste0(file_path,'spatial_function.R'))

ot <- read.table(paste0(file_path,'OT_PDAC.txt'),sep = '\t',header = T, row.names = 1)

rownames(names) <- names$cell
ot_map <- mapCluster(ot, meta = names, cluster = 'cell_type', min_cut = 0.15, balance = T)

# spatial scatter pie of cluster proportion
p <- stClusterPie(ot_map = ot_map, coord = ind, pie_scale = 0.8)
p

# single cluster proportion
p1 <- stClusterExp(ot_map, coord = ind, cluster = 'Cancer clone A',cut = 0.25)
p2 <- stClusterExp(ot_map, coord = ind, cluster = 'Ductal',cut = 0.25)

pdf('/data/yupines/kai/test/PDAC_CA_Ductal.pdf',width = 15,height = 6.5)
p1+p2
dev.off()

# gene expression
dt <- stGeneNorm(dataA)
stGeneExp(exp = dt, coord = ind, gene = c('CRISP3','TM4SF1'))
stGeneExp(exp = dt, coord = ind, gene = c('MUC5B'))





p1 <- stClusterExp(ot_map, brca, cluster = 'CAFs',cut = 0.15, point_size = 0.8)
p2 <- stClusterExp(ot_map, brca, cluster = 'Cancer.Epithelial',cut = 0.35, point_size = 0.8)

pdf('/data/yupines/kai/test/BRCA_CAFs_Cancer.Epithelial.pdf',width = 12,height = 5)
p1+p2
dev.off()








rownames(names) <- names$cell
ot <- read.table('/data/yupines/kai/test/OT_PDAC_LF.txt',sep = '\t',header = T, row.names = 1)
ot <- as.data.frame(t(ot))
rownames(ot) <- sapply(strsplit(rownames(ot),'\\.'),function(x)x[[1]])

ot_map <- mapCluster(ot, meta = names, cluster = 'cell_type', min_cut = 0.25,
                     balance = T)


# spatial scatter pie of cluster proportion
p <- stClusterPie(ot_map = ot_map, coord = ind, pie_scale = 0.8)

pdf('/data/yupines/kai/test/OT_PDAC_LF_T.pdf',height = 7.5,width = 10)
print(p)
dev.off()

