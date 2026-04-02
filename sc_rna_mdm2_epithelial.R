library(Seurat)
library(scCustomize)
library(scater)
library(metap)
library(writexl)
setwd("/Users/tatalab/Documents/AKshara")
rds=readRDS("./mdm2KO rds file/mdm2_epithelial_epi_subset.Total_cells_BHT_and_Cont.integrated_reclus_supclus_06-11-23_final.rds")
DimPlot(rds)

rds@active.ident
rds<- PrepSCTFindMarkers(rds, assay = "SCT")
#find AT1 markers
at1_markers<- FindMarkers(rds, ident.1= "AT1", min.pct = .10, logfc.threshold = 1, only.pos = T)
at1_markers<- at1_markers[at1_markers$p_val_adj<0.01,]
at1_markers$Genes<- rownames(at1_markers)
write_xlsx(at1_markers, path = "./mdm2KO rds file/at1_marker_p_val_adj_0_01_logfc_1.xlsx")
FeaturePlot(rds, features = "Igfbp2", min.cutoff = "q10")

# find pats marker
pats_marker<- FindMarkers(rds, ident.1 = "PATS", min.pct = .10, logfc.threshold = 1, only.pos = T)
pats_marker<- pats_marker[pats_marker$p_val_adj < 0.01 ,]
pats_marker$Genes<- rownames(pats_marker)
write_xlsx(pats_marker, path = "./mdm2KO rds file/pats_marker_p_val_adj_0_01_logfc_1.xlsx")
FeaturePlot(rds, features = "S100a14", min.cutoff = "q10")

VlnPlot(rds, features = "Bdnf")

DefaultAssay(rds)<- "integrated"
#find AT1 markers int
at1_markers_int<- FindMarkers(rds, ident.1= "AT1", min.pct = .10, logfc.threshold = 1, only.pos = T)
at1_markers_int<- at1_markers_int[at1_markers_int$p_val_adj<0.01,]
at1_markers_int$Genes<- rownames(at1_markers_int)
write_xlsx(at1_markers_int, path = "./mdm2KO rds file/at1_marker_int_p_val_adj_0_01_logfc_1.xlsx")

# find pats marker
pats_marker_int<- FindMarkers(rds, ident.1 = "PATS", min.pct = .10, logfc.threshold = 1, only.pos = T)
pats_marker_int<- pats_marker[pats_marker_int$p_val_adj < 0.01 ,]
pats_marker_int$Genes<- rownames(pats_marker_int)
write_xlsx(pats_marker_int, path = "./mdm2KO rds file/pats_marker_int_p_val_adj_0_01_logfc_1.xlsx")
FeaturePlot(rds, features = "S100a14", min.cutoff = "q10")

