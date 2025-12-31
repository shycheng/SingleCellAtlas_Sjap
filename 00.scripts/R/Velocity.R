# library(velocyto.R)
# 
# 
# # Loading and merge all samples' loom files ------------------------------------------------------
# rename_barcode <- function(x){
#   rownames(x) <- gsub('_','-',rownames(x))
#   colnames(x) <- gsub('-out:','_',colnames(x))
#   colnames(x) <- gsub('x','-1',colnames(x))
#   return(x)
# }
# 
# GeneNum <- nrow(seob@assays$RNA)
# GeneId <- rownames(seob@assays$RNA)
# 
# spliced.dat <- matrix(nrow = 9111)#行名为基因，列名为barcode
# unspliced.dat <- matrix(nrow = 9111)
# ambigous.dat <- matrix(nrow = 9111)
# 
# samples <- c('14dpi','18dpi_F','18dpi_M','22dpi_F','22dpi_M','26dpi_F','26dpi_M')
# for (s in samples){
#   dat <- read.loom.matrices(file=paste0("/public/shycheng/",s,"-out/","velocyto/",s,"-out.loom"))
#   sp <- dat$spliced
#   sp <- rename_barcode(sp)
#   sp <- sp[GeneId,]
#   spliced.dat <- cbind(spliced.dat, sp)
#   print(paste0(s," spliced finished!"))
#   usp <- dat$unspliced
#   usp <- rename_barcode(usp)
#   usp <- usp[GeneId,]
#   unspliced.dat <- cbind(unspliced.dat, usp)
#   print(paste0(s," unspliced finished!"))
#   amb <- dat$ambiguous
#   amb <- rename_barcode(amb)
#   amb <- amb[GeneId,]
#   ambigous.dat <- cbind(ambigous.dat, amb)
#   print(paste0(s," ambiguous finished!"))
# }
# 
# saveRDS(spliced.dat,file = "spliced.dat_loom.rds")
# saveRDS(ambigous.dat,file = "ambigous.dat_loom.rds")
# saveRDS(unspliced.dat,file = "unspliced.dat_loom.rds")



spliced.dat <- readRDS(file = "/home/shycheng/Projects/scRNA-seq/Trajectory_inference/spliced.dat_loom.rds")
ambigous.dat <- readRDS(file = "/home/shycheng/Projects/scRNA-seq/Trajectory_inference/ambigous.dat_loom.rds")
unspliced.dat <- readRDS(file = "/home/shycheng/Projects/scRNA-seq/Trajectory_inference/unspliced.dat_loom.rds")


# save files for scVelo
library(data.table)
save_filesFORscVelo <- function(seobj,filename,baseDIR){
  filesPATH = paste0(baseDIR,'/',filename,'_loom')
  dir.create(filesPATH)
  print(paste0('All files will be saved in ',filesPATH))
  #
  cells <- intersect(colnames(seobj),colnames(spliced.dat))
  print(paste0('Total cell number:',length(cells)))
  #
  ambigous.dat1 <- ambigous.dat[,c(cells)]
  spliced.dat1 <- spliced.dat[,c(cells)]
  unspliced.dat1 <- unspliced.dat[,c(cells)]
  #
  ambigousFile = paste0(filesPATH,'/',filename,'.ambigous.txt')
  splicedFile = paste0(filesPATH,'/',filename,'.spliced.txt')
  unsplicedFile = paste0(filesPATH,'/',filename,'.unspliced.txt')
  #
  data.table::fwrite(as.data.frame(ambigous.dat1),ambigousFile, sep="\t",row.names = T,col.names = T,verbose = T)
  print(paste0(ambigousFile,"... Finished!"))
  #
  data.table::fwrite(as.data.frame(spliced.dat1),splicedFile, sep="\t",row.names = T,col.names = T,verbose = T)
  print(paste0(splicedFile,"... Finished!"))
  #
  data.table::fwrite(as.data.frame(unspliced.dat1),unsplicedFile, sep="\t",row.names = T,col.names = T,verbose = T)
  print(paste0(unsplicedFile,"... Finished!"))
  
  ## 提取umap坐标和细胞类型注释
  seobj$barcode <- colnames(seobj)
  if ('Cell_Type' %in% colnames(seobj@meta.data)) {
    celltype <- 
      
      seobj@meta.data[,c('barcode','seurat_clusters',"Cell_Type")] 
  } else{
    celltype <- seobj@meta.data[,c('barcode','seurat_clusters')] 
  }
  write.csv(celltype, file=paste0(filesPATH,'/',filename,'.celltype.csv'),row.names = FALSE)
  
  genes <- rownames(spliced.dat1)
  write.csv(genes, file=paste0(filesPATH,'/',filename,'.genes.csv'))
  
  umaps <- seobj@reductions$umap@cell.embeddings
  write.csv(umaps,paste0(filesPATH,'/',filename,".umap.csv"))
  rm(ambigous.dat1,spliced.dat1,unspliced.dat1,
     ambigousFile,splicedFile,unsplicedFile
  )
  gc()
  print(paste0('All files has been saved in ',filesPATH))
}




