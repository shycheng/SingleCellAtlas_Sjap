library(tidyverse)
library(clusterProfiler)
#library(org.Sjaponicum.eg.db)
library(RColorBrewer)
library(ggsci)
#install.packages('/home/shycheng/genome/org.Sjaponicum.eg.db/',repos = NULL,type = 'source')

qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

col_vectors <- col_vector[-c(6,7,8)]

cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")

#' helper function used in stacked_vln_plot
#' plot.margin to adjust the white space between each plot. pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0, 0, -0, 0), "cm"), ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, cols = col_vectors[1:length(unique(obj@meta.data$Cell_Type))],... ) + 
    xlab("") + ylab(feature) + ggtitle("") + theme(legend.position = "none",axis.title.y  = element_text(size=8,face = 'bold'),plot.margin = plot.margin)
  return(p)
}
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

stacked_vln_plot<- function(obj, features,pt.size = 0, plot.margin = unit(c(-0, 0, -0, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  for (i in 1:(length(features)-1) ) {
    plot_list[[i]] <- plot_list[[i]] + theme( axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = plot.margin)
  }
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1,guides = 'collect')
  
  return (p)
}


barplot_per <- function(Seob,colValues = brewer.pal(7,'RdBu'),Cell_Type = "Cell_Type"){
  df = data.frame(cluster=Seob@meta.data[[Cell_Type]], day=Seob@meta.data[["orig.ident"]], Seob@reductions$umap@cell.embeddings,Cell_type = Seob@meta.data[[Cell_Type]])
  # g1 = ggplot_build(Umap_Female_seob_neoblast)
  # color_seurat = g1$data[[1]] %>% group_by(group) %>% dplyr::summarise(cols = dplyr::first(colour))
  bar_plot = function(df, position='stack') {
    df_bar = with(df, table(cluster, day))
    df_bar = as.data.frame(prop.table(df_bar,1) * 100) #%>%
    #arrange(day,desc(Freq))
    #df_bar$day <- factor(df_bar$day,levels = as.character(unique(Seob@meta.data$orig.ident)))
    df_bar$day <- factor(df_bar$day,levels = c('26dpi_F','22dpi_F','18dpi_F',"14dpi", 
                                               '18dpi_M', '22dpi_M','26dpi_M'))
    
    
    #df_bar$cluster <- factor(as.character(df_bar$cluster),levels = as.character(df_bar$cluster[1:length(levels(df_bar$cluster))]))
    #df_bar = df[,c('cluster', 'day')] %>% group_by(cluster,day) %>% summarise(freq = n())
    p = ggplot(df_bar, aes(fill=day, y=Freq, x=cluster))
    p = p + geom_bar(stat="identity", position=position,width = 0.8) + labs(fill='Day')
    p = p + theme(panel.background=element_blank())+ scale_fill_manual(values = colValues) +
      xlab('Clusters') +
      ylab('Percentage(%)') +
      ggtitle("Composition(%) of Clusters by Time") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = .5,face = 'bold'))
    return (p)
  }
  bar_plot(df)
}




barplot_cluster <- function(seob,umap_plot){
  require(patchwork)
  g1 = ggplot_build(umap_plot)
  color_seurat = g1$data[[1]] %>% group_by(group) %>% dplyr::summarise(cols = dplyr::first(colour))
  
  df = data.frame(cluster=seob@meta.data[["seurat_clusters"]], day=seob@meta.data[["orig.ident"]], seob@reductions$umap@cell.embeddings,Cell_type = seob@meta.data$Cell_Type)
  df_test = with(df, table(cluster,day))
  df_test = as.data.frame(prop.table(df_test,1) * 100)
  
  cell_type <- df %>% arrange(cluster) %>% 
    distinct(cluster,.keep_all = T) %>% 
    dplyr::select(Cell_type)
  
  df_type = with(df, table(day, Cell_type))
  df_type = as.data.frame(prop.table(df_type,1) * 100)
  n = nrow(df_type)/length(unique(df_type$Cell_type))
  df_type$cols <- rep(color_seurat$cols,each = n)
  # df_type = df_type %>% 
  #   arrange(day,Freq)
  cols = unique(df_type$cols)
  df_test$Cell_type <- factor(as.character(rep(cell_type$Cell_type,n)),levels = 
                                rev(as.character(unique(df_type$Cell_type))))
  df_test <- df_test %>% 
    arrange(desc(Cell_type),day,desc(Freq))
  df_test$cluster <- factor(as.character(df_test$cluster),
                            levels = as.character(
                              unique(df_test$cluster)
                            )
  )
  df_test$day <- factor(df_test$day,levels = c('26dpi_F','22dpi_F','18dpi_F',"14dpi", 
                                               '18dpi_M', '22dpi_M','26dpi_M'))
  
  p = ggplot(df_test, aes(fill=day, y=Freq, x=cluster)) 
  p = p + geom_bar(stat="identity", position='stack') + labs(fill='Day')
  p = p + theme(panel.background=element_blank())+ scale_fill_manual(values = brewer.pal(7,'RdBu')) +
    xlab('Clusters') +
    ylab('Percentage(%)') +
    ggtitle("Composition(%) of Clusters by Time")
  
  tab <- as.data.frame(table(cell_type$Cell_type))
  tab$Var1 <- factor(tab$Var1,levels = as.character(unique(df_type$Cell_type)))
  tab$x <- '1'
  corbar = ggplot(tab, aes(fill = Var1,y=Freq, x= x)) 
  corbar = corbar + geom_bar(stat="identity", position='stack') 
  corbar = corbar + theme(panel.background=element_blank())+ scale_fill_manual(values = cols) +
    xlab('') +
    ylab('')
  ps <- p + corbar + plot_layout(widths = c(19,1))
  return(ps)
}



barplot_celltype <- function(Seob,umap_plot,n){
  df = data.frame(cluster=Seob@meta.data[["seurat_clusters"]], day=Seob@meta.data[["orig.ident"]], Seob@reductions$umap@cell.embeddings,Cell_type = Seob@meta.data$Cell_Type)
  g1 = ggplot_build(umap_plot)
  color_seurat = g1$data[[1]] %>% group_by(group) %>% dplyr::summarise(cols = dplyr::first(colour))
  df_type = with(df, table(day, Cell_type))
  df_type = as.data.frame(prop.table(df_type,1) * 100)
  df_type$cols <- rep(color_seurat$cols,each = n)
  # df_type = df_type %>% 
  #   arrange(day,Freq)
  # df_type$Cell_type <- factor(as.character(df_type$Cell_type),levels = as.character(unique(df_type$Cell_type)))
  cols = unique(df_type$cols)
  #df_bar = df[,c('cluster', 'day')] %>% group_by(cluster,day) %>% summarise(freq = n())
  p = ggplot(df_type, aes(fill=Cell_type, y=Freq, x=day)) 
  p = p + geom_bar(stat="identity", position='stack',width = 0.68) + labs(fill='Cell type')
  p = p + theme(panel.background=element_blank())+ scale_fill_manual(values = cols) +
    xlab('') +
    ylab('Percentage(%)') +
    theme(axis.text.x = element_text(angle = 90, vjust = .5,hjust = -0.5))+ 
    ggtitle("Composition(%) of cell types by Time")
  return(p)
}

getFeaturePlots <- function(Seobj,id){
  p <- FeaturePlot(Seobj, 
                   reduction = "umap", 
                   features = id, 
                   #split.by = 'Sex',label.size = 3,ncol = 2,
                   label = T,min.cutoff = "q10", max.cutoff = "q90",
                   cols = c("grey", "red")) +NoAxes() + 
    labs(title = paste0("Gene:",id),subtitle  = paste0('Annotation:  ',ID2NAME[id,'Note']))
  return(p)
}

########################################################
library(reshape2)
library(ggpubr)
library(ggalluvial)

getCellprop_tab <- function(Seob,umap_plot,Cell_Type = "Cell_Type",UMAP="UMAP"){
  df = data.frame(cluster=Seob@meta.data[["seurat_clusters"]], day=Seob@meta.data[["orig.ident"]], Seob@reductions[[UMAP]]@cell.embeddings,Cell_type = Seob@meta.data[,Cell_Type])
  g1 = ggplot_build(umap_plot)
  color_seurat = g1$data[[1]] %>% group_by(group) %>% 
    dplyr::summarise(cols = dplyr::first(colour))
  df_type = with(df, table(day, Cell_type))
  df_type = as.data.frame(prop.table(df_type,1) * 100)
  colorValues <- c(rep(color_seurat$cols,each = nrow(df_type)/length(unique(df_type$Cell_type)))) # fill blank value
  df_type$cols <- c(colorValues,rep('grey',nrow(df_type) - length(colorValues)))
  return(df_type)
}

plot_alluvium <- function(seob ,umapplot,Cell_Type = 'NewCellType',filename,
                          width = 4.5,height = 6.5,UMAP='UMAP'){
  seob_CellProp <- getCellprop_tab(seob,umapplot,Cell_Type,UMAP)
  p_alluvium <- ggplot(seob_CellProp,
                       aes(x = day, stratum = Cell_type, alluvium = Cell_type,
                           y = Freq,
                           fill = Cell_type, label = Cell_type)) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_fill_manual(values = unique(seob_CellProp$cols)) + 
    geom_flow() +
    geom_stratum(stat = "alluvium",
                 color = "darkgray") +
    theme(legend.position = "right",   axis.line = element_blank(),panel.background=element_blank()) +
    ggtitle("Frequency of cells in different conditions")
  ggsave(filename = filename,width = width,height = height,plot = p_alluvium)
}



# ID2NAME <- id2name
# rownames(ID2NAME) <- gsub('_','-',rownames(id2name))
id2name <- read.delim("/home/shycheng/genome/Schitosoma_jp/Schistosoma_japonicum/SJ_GODB/sj_geneInfo_raw.tsv",sep = '\t',header = F)
colnames(id2name) <- c("id","Gene_Name","Note")
rownames(id2name) <- id2name$id
id2name$Gene_Name <- as.character(id2name$Gene_Name)
idx <- which(nchar(as.character(id2name$Gene_Name)) > 13)
id2name$Gene_Name[idx] <- as.character(rownames(id2name[idx,]))
id2name$Note <- gsub("\\(.*\\)","",id2name$Note)
id2name$id <- gsub('_','-',id2name$id)
rownames(id2name) <- id2name$id
ID2NAME <- id2name
rownames(ID2NAME) <- gsub('_','-',rownames(id2name))


simFP <- function(id,seob,min.cutoff = 'q1',max.cutoff = 'q99',label = FALSE,...){
  require(patchwork)
  id <- rownames(seob)[grepl(id,rownames(seob))]
  p <- Seurat::FeaturePlot(seob,features = id,label = label,min.cutoff = min.cutoff,max.cutoff = max.cutoff,raster=FALSE,...)
  p
}


multi_simFP <- function(ids,seob,min.cutoff = 'q1',max.cutoff = 'q99',label = FALSE,legend.position = "none",ncol = 3,...){
  purrr::map(ids,function(id){
    simFP(id = id,seob = seob,
          pt.size = 0.05,cols = c("grey", "red"),
          min.cutoff = 'q10',max.cutoff = 'q90',...) &
      ggtitle(label = paste0(id,"(",id2name[id,"Gene_Name"],")")) &
      theme(panel.border = element_rect(fill = NA,linewidth = 1,colour = 'black'),
            plot.title = element_text(face = "italic"),legend.position = legend.position)
  }) %>% patchwork::wrap_plots(ncol = ncol)
}




cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 








plot_umap <- function(seob, 
                      output_file,
                      group_by = "Cell_Type",
                      col_pals = c("#7FC97F", "#BEAED4", brewer.pal(4,'Oranges'), "#D95F02", "#7570B3",
                                   brewer.pal(5,'Blues'),"#E7298A","#1F78B4", "#B2DF8A",
                                   "#33A02C", "#FB9A99", "#E31A1C", "grey" ),
                      width=12, 
                      height=10, 
                      point_size=0.1,
                      show_labels=TRUE) {
  # 获取细胞类型和设置颜色
  
  cell_levels <- levels(seob@meta.data[,group_by])
  celltypes_pal <- col_pals
  names(celltypes_pal) <- cell_levels
  
  
  # 获取UMAP坐标
  umap <- seob@reductions$umap@cell.embeddings
  # 随机排序防止覆盖
  rand_ind <- sample.int(dim(umap)[1])
  
  col_pal <- celltypes_pal[as.character(seob@meta.data[rand_ind,][,group_by])]
  
  
  # 打开PDF设备
  pdf(output_file, width=width, height=height)
  
  # 设置绘图参数
  par(cex=2, las=1, mar=c(4,4,1,10), lwd=3)  # 增加右边距以容纳图例
  
  # 绘制UMAP散点图
  plot(umap[rand_ind,1], umap[rand_ind,2], pch=16, xlab='UMAP-1', ylab='UMAP-2', 
       frame=F, col=col_pal, cex=point_size)
  
  # 添加图例
  legend("topright", legend=cell_levels, fill=celltypes_pal, 
         border=NA, bty="n", cex=0.6, inset=c(-0.25, 0), xpd=TRUE)
  
  # 如果需要添加标签
  if(show_labels) {
    # 计算每个细胞类型的中心点
    centers <- lapply(cell_levels, function(ct) {
      cells <- which(seob@meta.data[,group_by] == ct)
      if(length(cells) > 0) {
        c(mean(umap[cells, 1]), mean(umap[cells, 2]))
      } else {
        c(NA, NA)
      }
    })
    
    # 添加标签
    for(i in 1:length(cell_levels)) {
      if(!any(is.na(centers[[i]]))) {
        text(centers[[i]][1], centers[[i]][2], labels=cell_levels[i], 
             col="black", font=1.5, cex=0.8)
      }
    }
  }
  
  # 关闭设备
  dev.off()
}









