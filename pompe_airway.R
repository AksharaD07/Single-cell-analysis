library(Seurat)
library(scCustomize)
library(scater)
library(metap)

mxt_path<-"./Documents/AKshara/popmpe_airway/data/matrix.mtx.gz"
bar_path<-"./Documents/AKshara/popmpe_airway/data/barcodes.tsv.gz"
genes_path<-"./Documents/AKshara/popmpe_airway/data/genes.tsv.gz"

mat<- ReadMtx(mtx= mxt_path, features = genes_path, cells = bar_path)
pompe<- CreateSeuratObject(mat, min.cells = 3, min.features = 200)

## QC and normalization
pompe[["pct_mt"]]<- PercentageFeatureSet(pompe, pattern = "^MT-")
View(pompe@meta.data)
VlnPlot(pompe, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3)
pompe<- subset(pompe, subset=  nFeature_RNA > 200 & nFeature_RNA < 2500
                                  & pct_mt < 5 )
pompe<- NormalizeData(pompe, normalization.method = "LogNormalize", scale.factor = 10000)
pompe<- FindVariableFeatures(pompe, selection.method = "vst", nfeatures = 2000)

top10_gene<- head(VariableFeatures(pompe), 10)
plot1 <- VariableFeaturePlot(pompe)
plot2 <- LabelPoints(plot = plot1, points = top10_gene, repel = TRUE)
plot1 + plot2

# Dimensionality reduction
all_gene<- rownames(pompe)
pompe<- ScaleData(pompe, features = all_gene)
pompe<- RunPCA(pompe, features = VariableFeatures(pompe))
# Visualizations of PCA
print(pompe[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pompe, dims = 1:2, reduction = "pca")
DimPlot(pompe, reduction = "pca")
DimHeatmap(pompe, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pompe, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pompe)

pompe<- FindNeighbors(pompe, dims = 1:10)
pompe<- FindClusters(pompe, resolution = 0.5)
pompe<- RunUMAP(pompe, dims = 1:10)
DimPlot(pompe, reduction = "umap")
VlnPlot(pompe, features = "PDGRF")
FeaturePlot(pompe, features = "HAPLN3")
View(as.data.frame(VariableFeatures(pompe)))

X<-FindAllMarkers(pompe, min.pct = 0.25)
