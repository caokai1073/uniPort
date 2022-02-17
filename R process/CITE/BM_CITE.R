library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)

file_path <- '/data/work/uniPort/BM_CITE/'

InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

write.table(as.data.frame(bm@assays$RNA@data), sep = '\t',quote = F,
            file = paste0(file_path,'bmcite_rna_processed.txt'))
            
write.table(as.data.frame(bm@assays$ADT@counts), sep = '\t',quote = F,
            file = paste0(file_path,'bmcite_ADT_ProName.txt'))
            
gene <- read.csv(paste0(file_path,'protein_gene_name.csv'))
rownames(bm@assays$ADT@counts) <- sapply(
  rownames(bm@assays$ADT@counts),
  function(x){gene$gene[gene$protein==x]}) %>% 
  {names(.) <-NULL;.}

write.table(as.data.frame(bm@assays$ADT@data),sep = '\t',quote = F,
            file = paste0(file_path,'bmcite_cite_processed.txt'))
            
meta <- bm@meta.data[,c("celltype.l1","celltype.l2","RNA.weight","ADT.weight","wsnn_res.2")]
write.table(meta,file = paste0(file_path,"meta_cite_combine.txt", sep='\t',quote = F))


