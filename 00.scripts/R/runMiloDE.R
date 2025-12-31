library(miloDE)
library(miloR)
library(scWGCNA)
library(scuttle)
suppressMessages(library(uwot))
library(scran)
suppressMessages(library(Seurat))
# plotting libraries
library(ggplot2)
library(viridis)
library(ggpubr)


sce = readRDS("01.rds_files/sce_mouseGastrulation.rds")
EmbryoCelltypeColours = readRDS("01.rds_files/EmbryoCelltypeColours.rds")


# downsample to few selected cell types
cts = c("Spinal cord" , "Haematoendothelial progenitors", "Endothelium" , "Blood progenitors 1" , "Blood progenitors 2")
sce = sce[, sce$celltype.mapped %in% cts]
# let's rename Haematoendothelial progenitors
sce$celltype.mapped[sce$celltype.mapped == "Haematoendothelial progenitors"] = "Haem. prog-s."


# delete row for tomato
sce = sce[!rownames(sce) == "tomato-td" , ]

sce = sce[,colData(sce)$sample %in% c(1,4)]

# add logcounts
sce = logNormCounts(sce)
sce@assays@data$logcounts <- as(sce@assays@data$logcounts, "CsparseMatrix")

# update tomato field to be more interpretable 
sce$tomato = sapply(sce$tomato , function(x) ifelse(isTRUE(x) , "Tal1_KO" , "WT"))

# for this exercise, we focus on 3000 highly variable genes (for computational efficiency)
dec.sce = modelGeneVar(sce)
hvg.genes = getTopHVGs(dec.sce, n = 3000)
sce = sce[hvg.genes , ]
# change rownames to symbol names
rowdata = as.data.frame(rowData(sce))
rownames(sce) = rowdata$SYMBOL



# add UMAPs
set.seed(32)
umaps = as.data.frame(uwot::umap(reducedDim(sce , "pca.corrected")))
# let's store UMAPs - we will use them for visualisation
reducedDim(sce , "UMAP") = umaps

# plot UMAPs, colored by cell types
umaps = cbind(as.data.frame(colData(sce)) , reducedDim(sce , "UMAP"))
names(EmbryoCelltypeColours)[names(EmbryoCelltypeColours) == "Haematoendothelial progenitors"] = "Haem. prog-s."
cols_ct = EmbryoCelltypeColours[names(EmbryoCelltypeColours) %in% unique(umaps$celltype.mapped)]



set.seed(32)

sce_milo = assign_neighbourhoods(sce , k = 20 , order = 2, 
                                 filtering = TRUE , reducedDim_name = "pca.corrected" , verbose = F)


nhoods_sce = nhoods(sce_milo)
# assign cell types for nhoods 
nhood_stat_ct = data.frame(Nhood = 1:ncol(nhoods_sce) , Nhood_center = colnames(nhoods_sce))
nhood_stat_ct = miloR::annotateNhoods(sce_milo , nhood_stat_ct , coldata_col = "celltype.mapped")



stat_auc = suppressWarnings(calc_AUC_per_neighbourhood(sce_milo ,  condition_id = "tomato"))


milo_res <- list(stat_auc = stat_auc,
                 sce = sce,
                 sce_milo = sce_milo
                 )

saveRDS(milo_res,"/home/shycheng/Projects/SingleCell_Sja2024/01.rds_files/milo_res.rds")

# de_stat = de_test_neighbourhoods(sce_milo ,
#                                  sample_id = "sample",
#                                  design = ~tomato,
#                                  covariates = c("tomato"),
#                                  subset_nhoods = stat_auc$Nhood[!is.na(stat_auc$auc)],
#                                  output_type = "SCE",
#                                  plot_summary_stat = TRUE,
#                                  layout = "UMAP", 
#                                  verbose = T)










