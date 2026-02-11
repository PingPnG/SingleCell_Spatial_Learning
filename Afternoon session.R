library('Seurat')

###################name and load your dataset!! Copy file location and paste in-between ' '. Switch " \ " with " / ".
###RP treatment on skin
RP_explant <- readRDS(file = "Skin_Explant_RP_FlexSeq.rds")
DimPlot(RP_explant)

############## https://satijalab.org/seurat/articles/de_vignette ##########
###cell type identification is a manual process, 
##first look at the cluster and then do marker identification
##then based on the cell marker to define the cell type.
###usually start from cluster and then find marker , then manual annotation
### if you cluster differently the same cell can be named as different cell type
###cost for single cell from U michigan $900 per sample, sent 2mm biopsy , reembedded in 
###Paraphin and then cut 20X25nm sections for two used for the analysis, about half left 
###untouched

table(RP_explant@active.ident)
###pct1 percentage of the cell in this group, the next biggest group percentage is in pct2
Markers <- FindAllMarkers(object = RP_explant,
                       min.pct = 0.5,
                       logfc.threshold = 1,
                       only.pos=TRUE,
                       test.use = 'wilcox',
                       return.thresh= .01)

#You can also use more cell vs. the others by using , to seperate different cell type
Diff_Genes <- FindMarkers(object = RP_explant,
                             ident.1 = c("Keratinized Keratinocytes"), 
                             ident.2 = c("Basal Keratinocytes"),
                             logfc.threshold = 3,
                             min.pct = 0.1,
                             only.pos=FALSE,
                             test.use = 'wilcox')
#### if you are doing change to different test, using fold change as a large cut off to look at the differentiated genes

FeaturePlot(RP_explant, features = c( "IL37", "JAG2", "WNT10A","LCE1B"), ncol = 2)

VlnPlot(RP_explant, features = c( "IL37","LCE1B"), sort = "decreasing", ncol = 2)
###only show specific populations to only vidualize the specific cell types.
VlnPlot(RP_explant, features = "IL37", idents = c("Differentiated Keratinocytes",
                                                   "Keratinized Keratinocytes",
                                                   "Basal Keratinocytes"))

write.csv(Diff_Genes, "Diff_Genes.csv") ################saves inactive working directory###########


#########################################
colnames(RP_explant@meta.data)
Idents(RP_explant) <- "treatment"  ###set identity to specific metadata so that analysis is based on this
table(RP_explant$treatment)

DimPlot(RP_explant)
####Look at the treatment differences, change fold change cut off to smaller
Diff_Genes2 <- FindMarkers(object = RP_explant,
                          ident.1 = c("300 uM Retinyl propionate"), 
                          ident.2 = c("Vehicle"),
                          logfc.threshold = 1,
                          min.pct = 0.4,
                          only.pos=FALSE,
                          test.use = 'wilcox')

FeaturePlot(RP_explant, features = c( "KRT6A", "ZBTB16"), ncol = 2)


VlnPlot(RP_explant, features = "KRT6A", group.by = "celltype", split.by = "treatment")

VlnPlot(RP_explant, features = "KRT6A", group.by = "celltype", split.by = "treatment",
        idents = c("300 uM Retinyl propionate","Vehicle"))

##############################################
###how you create a new metadata type 
RP_explant$celltype_treatment<- paste(RP_explant$celltype, RP_explant$treatment, sep = "_")

Idents(RP_explant) <- "celltype_treatment"
table(RP_explant$celltype_treatment)
#####cell numbers can be artifact, might be just sample differences
Diff_Genes3 <- FindMarkers(object = RP_explant,
                          ident.1 = c("Fibroblasts_300 uM Retinyl propionate"), 
                          ident.2 = c("Fibroblasts_25 uM Clobetasol"),
                          logfc.threshold = .5,
                          min.pct = 0.4,
                          only.pos=FALSE,
                          test.use = 'wilcox')
VlnPlot(RP_explant, features = c("COL6A1", "COL6A2"),  group.by = "celltype", split.by = "treatment",
        idents = c("Fibroblasts_300 uM Retinyl propionate","Fibroblasts_25 uM Clobetasol"))

VlnPlot(RP_explant, features = c("COL6A1", "COL6A2"), split.by="orig.ident", group.by = "celltype", 
        idents = c("Fibroblasts_300 uM Retinyl propionate","Fibroblasts_25 uM Clobetasol"))

write.csv(Diff_Genes3, "Fibroblast_RP_vs_Clob.csv")

##############################Lets look at Diff gene file################################

              ###    https://toppgene.cchmc.org/enrichment.jsp    ####

                  ###    https://maayanlab.cloud/Enrichr/    ###




