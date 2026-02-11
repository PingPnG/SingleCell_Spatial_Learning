#install.packages('Seurat')
library('Seurat')

###################name and load your dataset!! Copy file location and paste in-between ' '. Switch " \ " with " / ".
MDE_test <- readRDS(file = "PnG_face_seurat_subset_n10.RDS")
Idents(MDE_test) <- "celltype" ### adjust as original use subtype and hard to look




###########Display TSNE of data#######################
DimPlot(MDE_test)

DimPlot(MDE_test, label = TRUE, repel = TRUE) + NoLegend() 



###################Plot single gene #############################
FeaturePlot(MDE_test, "FADS2")

FeaturePlot(MDE_test, "KRT15", min.cutoff = 2, max.cutoff = NA)

#Plot Your Gene, random generate 10 genes and plot
sample(rownames(MDE_test), 10)

###################Plot Multiple gene show a group of genes as list#############################
FeaturePlot(MDE_test, features = c( "PTPRC", "KRT2", "KRT15", "PMEL"), min.cutoff = .08, ncol = 2)
##if gene is not expressed you will get an error and nothing will show out
FeaturePlot(MDE_test, features = c( "PTPRC", "KRT2", "UGCG", "PMEL", "FADS2", "FADS1"), min.cutoff = .08, ncol = 3)

#Plot 4 genes of your choice using three columns

######################Violin Plot################
VlnPlot(MDE_test, "LCE3E") + NoLegend() 

VlnPlot(MDE_test, features = c("COL6A1", "PECAM1")) + NoLegend()
####remove the points so it looks better 
VlnPlot(MDE_test, features = c("COL1A1", "COL17A1"),pt.size=0, ncol = 2, sort="decreasing")

#Plot 3 genes of your choice using three columns,high to low expression
VlnPlot(MDE_test, features = c("FADS2", "IL1A", "COL1A1"),pt.size=0, ncol = 3, sort="decreasing")

#####################Plot other MetaData besides Cell Type ########################
####list data for all the meta data and you can choose which type of meta data across different population
colnames(MDE_test@meta.data)

DimPlot(MDE_test, label = TRUE, repel = TRUE, group.by = "age_range")
###lovel this to show all the sub cell type of fibroblast 
DimPlot(MDE_test, label = TRUE, repel = TRUE, group.by = "Fb_sub", label.size = 2) + NoLegend()
#You can use label size to be 
DimPlot(MDE_test, label = TRUE, repel = TRUE, group.by = "age", label.size = 0)
#Use findmarker to identify all the cell biomarker
#####################Plot gene based on other MetaData ########################
VlnPlot(MDE_test, features = c("DCD", "DST","FADS2", "SLC52A1"), ncol = 2, pt.size=0, group.by = "age_range")
#sometimes look at if only one idividual have all the values. 
VlnPlot(MDE_test, "DCD", group.by = "age") ###look at the log transformed count , log based 10
VlnPlot(MDE_test, "NR4A1", group.by = "age")
################Xenium Data##############################
    #Clear Workspace#
#https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2
xenium<- readRDS(file = "scalp1_annRef.rds") 
#name of file everyone downloaded is scalp1_annRef.rds#
#rm("MDE_test")
VlnPlot(xenium, features = "TOP2A", group.by = "celltype3") + NoLegend()
colnames(xenium@meta.data)
dev.new() ####create a new window to make it bigger 
ImageDimPlot(xenium, cols = "polychrome", size = 0.75, group.by = "celltype3")
ImageFeaturePlot(xenium, features = c("FADS2"), size = 1, max.cutoff = .05)
####there are interactive DimPlot
##############VisiumHD Date#########################

visiumHD <- readRDS(file = "scalp_S2_final.rds")
VlnPlot(visiumHD, features = "KRT15", group.by = "celltype_final")+ NoLegend()
dev.new()
ImageDimPlot(visiumHD, cols = "glasbey", size = 0.75, group.by = "celltype_final", flip_xy = TRUE) ###use flip_xy to make it bigger so that so show better color
##when it is smoothout together, it will change to white and not look good, you can use perameter to pull specific region out
##subset function -- you can use subset out to only look at specific function
visiumHD <- UpdateSeuratObject(visiumHD)
SpatialFeaturePlot(visiumHD, features = c("KRT15"), alpha = c(0.1, 1), image.alpha=0)####image.alpha=0 get rid of the box
SpatialFeaturePlot(visiumHD, features = c("KRT1"), alpha = c(0.1, 1), image.alpha=0)+ DarkTheme()####image.alpha=0 get rid of the box
SpatialFeaturePlot(visiumHD, features = c("KRT79"), alpha = c(0.1, 1), image.alpha=0)####image.alpha=0 get rid of the box
SpatialFeaturePlot(visiumHD, features = c("KRT79"), alpha = c(0.1, 1), image.alpha=0)####image.alpha=0 get rid of the box

saveRDS(visiumHD, file = "scalp_S2_final_visiumHD.rds")
################here is how you can load other data into seurat
#https://satijalab.org/seurat/articles/pbmc3k_tutorial

##################need hair cycle data and change split.by################

hair <- readRDS(file = "Human_Anagen_HF.rds")
colnames(hair@meta.data)
DimPlot(hair)
DimPlot(hair, label = TRUE, repel = TRUE, group.by = "celltype", label.size = 3) + NoLegend()
FeaturePlot(hair, "TOP2A")
FeaturePlot(hair, "TOP2A", split.by = "Phase")
VlnPlot(hair, "TOP2A", split.by = "Phase")
###understand group by and split by
#VlnPlot(hair, "TOP2A", split.by = "Phase", group.by = "celltype")
VlnPlot(hair, "TOP2A", split.by = "G2M.Score") + NoLegend()

#####This way is too dirty
hair$phase.celltype <- paste(hair$Phase, hair$celltype, sep = ".")
VlnPlot(hair, features = "TOP2A", group.by = "phase.celltype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  NoLegend()
FeaturePlot(hair, "TOP2A", split.by = "phase.celltype")