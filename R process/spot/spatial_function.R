
# make sure that rownames(ot) is the barcodes of single-cell reference data
# ref is the seurat Object
# st can be a spatial seurat object or a spatial expression matrix
# provided meta must have cell names as rownames
# default balance is TRUE
# `stGeneExp` function is only for non seuratObject data 
# `mapCluster` used for mapping single-cell cluster to corresponding spatial data
# `stClusterPie` used for ploting spatial scatter pie plot
# `stClusterExp` used for ploting single cluster proportion

theme_pie <- function(){
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20,hjust = 0.5),
        legend.position = 'right',
        legend.title=element_text(size=16),
        legend.text = element_text(size = 16)) 
}
theme_gene <- function(){
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = .5, face = "bold", size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
}

theme_cluster <- function(){
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = .5, face = "bold", size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank())
}

mapCluster <- function(ot, ref = NULL, cluster = NULL, meta = NULL, balance = F, min_cut = NULL){
  if(is(ref,"Seurat")){
    meta <- ref@meta.data
  } else {
    if(!is.null(meta)){
      meta <- meta
    } else {
      print('Please provide meta file of scRNA data')
    }
  }
  if(!is.null(cluster)){
    ref_cluster <- meta[cluster]
  } else {
    print('Please provide reference cluster name')
  }
  ref_cluster$cell <- rownames(ref_cluster)
  if(rownames(ot) %in% ref_cluster$cell){
    ot[] <- lapply(ot, as.numeric)
  } else {
    print('Please make sure that rownames of ot is the barcodes of scRNA reference data or rownames of meta file is cell name!')
  }
  if(isTRUE(balance)){
    ot_map <- lapply(unique(ref_cluster[,cluster]),function(x){
      cell = ref_cluster[ref_cluster[,cluster] == x,]$cell
      dt = t(as.data.frame(apply(ot[cell,],2,sum)))
      rownames(dt) <- x
      dt <- dt/length(cell)
      return(dt)
    }) %>% 
      do.call('rbind',.) %>% as.data.frame() %>% t()
  } else {
    ot_map <- lapply(unique(ref_cluster[,cluster]),function(x){
      cell = ref_cluster[ref_cluster[,cluster] == x,]$cell
      dt = t(as.data.frame(apply(ot[cell,],2,sum)))
      rownames(dt) <- x
      dt <- dt*(length(cell)/nrow(ref_cluster))
      return(dt)
    }) %>% 
      do.call('rbind',.) %>% as.data.frame() %>% t()
  }
  
  if(is.null(min_cut)){
    return(ot_map)
  } else {
    ot_map <- apply(ot_map,1,function(x){
      x[x < max(x)*min_cut] = 0
      return(x)
    }) %>% as.data.frame() %>% t()
    return(ot_map)
  }
}


stClusterPie <- function(ot_map, st = NULL, slice = NULL, coord = NULL, 
                         img_alpha = 0.8, pie_alpha = 0.9, pie_scale = 0.38,
                         color = colorRampPalette(brewer.pal(9,"Set1"))(ncol(cord_ot)-4)){
  if(is(st,"Seurat")){
    if (is.null(slice) && !is.null(names(st@images))){
      slice <- names(st@images)[1]
      warning(sprintf("Using slice %s for ploting", slice))
      img <- st@images[[slice]]@image
    } else {
      img <- st@images[[slice]]@image
    }
    img_grob <- grid::rasterGrob(
      matrix(rgb(img[,,1],img[,,2],img[,,3],alpha = img_alpha),nrow=dim(img)[1]),
      width = unit(1,"npc"))
    cord <- st@images[[slice]]@coordinates[,c('imagerow','imagecol')]
    cord$imagerow_scaled <- cord$imagerow*st@images[[slice]]@scale.factors$lowres
    cord$imagecol_scaled <- cord$imagecol*st@images[[slice]]@scale.factors$lowres
    cord_ot <- cbind(cord,ot_map)
    
    spatial_scatterPie_plot <- ggplot2::ggplot() + 
      ggplot2::annotation_custom(
        grob = img_grob, 
        xmin = 0, 
        xmax = ncol(img), 
        ymin = 0, 
        ymax = -nrow(img)) +
      scatterpie::geom_scatterpie(
        data = cord_ot, 
        ggplot2::aes(x = imagecol_scaled, 
                     y = imagerow_scaled), 
        cols = colnames(cord_ot)[5:ncol(cord_ot)], 
        color = NA,  
        alpha = pie_alpha, 
        pie_scale = pie_scale) +  
      ggplot2::ylim(nrow(img),0) + 
      ggplot2::xlim(0, ncol(img)) + 
      ggplot2::coord_fixed(ratio = 1,xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      scale_fill_manual('cluster',values = color) +
      theme_pie()
    return(spatial_scatterPie_plot)
    
  } else {
    if(!is.null(coord)){
      colnames(coord) <- c('row_ind','col_ind')
    } else {
      print('Please provide coordinates of st data!')
    }
    cord_ot <- cbind(coord,ot_map)
    
    spatial_scatterPie_plot <- ggplot() + 
      scatterpie::geom_scatterpie(
        aes(x= row_ind, y= col_ind),
        data = cord_ot,
        cols= colnames(cord_ot)[3:ncol(cord_ot)],
        color = NA,
        pie_scale = pie_scale) + 
      scale_fill_manual('cluster', values = colorRampPalette(brewer.pal(9,"Set1"))(ncol(cord_ot)-2)) +
      theme_pie()
    return(spatial_scatterPie_plot)
  }
}


stClusterExp <- function(ot_map, st = NULL, slice = NULL, coord = NULL, cluster, cut = NULL,
                         img_alpha = 0.8, point_size = 0.8){
  dt.use = t(apply(ot_map, 1, function(x) x/sum(x, na.rm = TRUE)))
  if(!is.null(cut)){
    dt.use[,cluster][dt.use[,cluster] < max(dt.use[,cluster])*cut] <- NA
  }
  if(is(st,"Seurat")){
    if (is.null(slice) && !is.null(names(st@images))){
      slice <- names(st@images)[1]
      warning(sprintf("Using slice %s for ploting", slice))
      img <- st@images[[slice]]@image
    } else {
      img <- st@images[[slice]]@image
    }
    img_grob <- rasterGrob(
      matrix(rgb(img[,,1],img[,,2],img[,,3],alpha = img_alpha),nrow=dim(img)[1]),
      width = unit(1,"npc"))
    cord <- st@images[[slice]]@coordinates[,c('imagerow','imagecol')]
    cord$imagerow_scaled <- cord$imagerow*st@images[[slice]]@scale.factors$lowres
    cord$imagecol_scaled <- cord$imagecol*st@images[[slice]]@scale.factors$lowres
    
    cord_ot <- cbind(cord,dt.use)
    
    spatial_cluster_plot <- ggplot2::ggplot(
      data = cord_ot,aes(x = imagecol_scaled, y = imagerow_scaled, color = get(cluster))) + 
      ggplot2::annotation_custom(grob = img_grob, xmin = 0,xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      ggplot2::geom_point(size = point_size,na.rm = TRUE) +
      ggplot2::scale_color_gradientn(na.value = "transparent",colors = rev(c('#d7191c','#fdae61','#ffffbf'))) + 
      ylim(nrow(img),0) + 
      xlim(0, ncol(img)) + 
      coord_fixed(ratio = 1,xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      labs(color = "Proportion") +
      ggtitle(cluster) +
      theme_pie()
    return(spatial_cluster_plot)
    
  } else {
    if(!is.null(coord)){
      colnames(coord) <- c('row_ind','col_ind')
    } else {
      print('Please provide coordinates of st data!')
    }
    dt.use = t(apply(ot_map, 1, function(x) x/sum(x, na.rm = TRUE)))
    dt.use <- cbind(coord,dt.use[,cluster,drop=F])
    max.val <- quantile(dt.use[, cluster],0.95)
    dt.use[, cluster] <- ifelse(dt.use[, cluster] > max.val, max.val, dt.use[, cluster])
    spatial_cluster_plot <- ggplot(
      dt.use, aes(row_ind, col_ind, color = get(cluster), size = get(cluster))) + 
      geom_point() + 
      scale_size(range = c(2, 6.5)) +
      ggtitle(cluster) + 
      scale_color_viridis(alpha = 0.82,breaks = c(min(dt.use[, cluster]),max.val*1.15), 
                          labels = c("low","high"),option = "inferno",
                          limits = c(min(dt.use[,cluster])/1.1, max.val*1.2)) + 
      theme_gene()
    return(spatial_cluster_plot)
  }
}

stGeneNorm <- function(st, assay_name = NULL, scale.factor = 10000, norm.method = "LogNormalize"){
  if(is(st,"Seurat") && !is.null(st@assays)){
    if(is.null(assay_name)){
      assay <- DefaultAssay(st)
      warning(sprintf("Using default assay %s", assay))
    } else {
      assay <- assay_name
    }
    stRNA <- NormalizeData(st, normalization.method = norm.method, scale.factor = scale.factor)
    stRNA <- ScaleData(stRNA, features = rownames(stRNA),do.center = F)
  } 
  else {
    stRNA <- CreateSeuratObject(counts = st)
    stRNA <- NormalizeData(stRNA, normalization.method = norm.method, scale.factor = scale.factor)
    stRNA <- ScaleData(stRNA, features = rownames(stRNA),do.center = F)
    assay <- DefaultAssay(stRNA)
  }
  dt <- as.data.frame(stRNA@assays$RNA@scale.data)
  return(dt)
}

stGeneExp <- function(exp, gene, coord = NULL, ncol = NULL, size_range = c(2, 6.5), color_alpha = 0.85){
  colnames(coord) <- c("row_ind","col_ind")
  cord_gene <- cbind(coord,t(exp[gene,]))
  
  if(isTRUE(length(gene) == 1)){
    fig <- ggplot(cord_gene, aes(row_ind, col_ind, color = get(gene), size = get(gene))) + 
      geom_point() + 
      scale_size(range = size_range) +
      ggtitle(gene) + 
      scale_color_viridis(alpha = color_alpha,
                          breaks = c(min(cord_gene[,gene]), max(cord_gene[,gene])), 
                          labels = c("low","high")) + 
      theme_gene()
    return(fig)
  } 
  else {
    if(!is.null(ncol)){
      ncol = ncol
    } else {
      if(length(gene) > 2){
        ncol = 3
      } else {
        ncol = 2
      }
    }
    fig.list <- lapply(colnames(cord_gene)[3:ncol(cord_gene)], function(gene){
      dt.use <- cord_gene
      dt.use <- dt.use[,c("row_ind","col_ind",gene)]
      ggplot(dt.use, aes(row_ind, col_ind, color = get(gene), size = get(gene))) + 
        geom_point() + 
        scale_size(range = size_range) +
        ggtitle(gene) + 
        scale_color_viridis(alpha = color_alpha,
                            breaks = c(min(dt.use[,gene]), max(dt.use[,gene])), 
                            labels = c("low","high")) + 
        theme_gene()
      })
    fig <- ggpubr::ggarrange(plotlist = fig.list, ncol = ncol)
    return(fig)
  }
}

