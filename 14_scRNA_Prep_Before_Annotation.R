library(Seurat)
library(harmony)


## 01质控-------------------------------------

### 主要PercentageFeatureSet函数计算线粒体含量
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")### 人类使用pattern = "^MT-"，小鼠使用pattern = "^mt-"
metadata <- sce.all@meta.data### 该操作会在metadata数据里面增加一列叫做percent.mt

### 质控数据可视化，使用VlnPlot函数
### nFeature_RNA, number of Feature, 每个细胞中有多少个基因
### nCount_RNA, number of counts, 每个细胞中有多少个counts
### percent.mt, 我们自己增加的列,  每个细胞中线粒体基因的比例
VlnPlot(sce.all,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#过滤低质量细胞/基因
#过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
# 一般来说，在CreateSeuratObject的时候已经是进行了这个过滤操作
# 如果后期看到了自己的单细胞降维聚类分群结果很诡异，就可以回过头来看质量控制环节


### nFeature_RNA > 200
### percent.mt < 5
sce.all_filter <- subset(sce.all, subset = nFeature_RNA > 200 & percent.mt < 25)


## 02 常规操作--------------------------------------------------
sce.all_filter <- NormalizeData(sce.all_filter)
sce.all_filter <- FindVariableFeatures(sce.all_filter, selection.method = "vst", nfeatures = 2000)
sce.all_filter <- ScaleData(sce.all_filter, features = rownames(sce.all_filter))
sce.all_filter<- RunPCA(sce.all_filter, features = VariableFeatures(object = sce.all_filter),reduction.name = "pca")


## 03.批次矫正与否(根据情况，见印象笔记的单细胞笔记)-------------------------------------------------------------------------------
### 使用harmony进行批次矫正
### 默认reduction = "pca",group.by.vars参数输入的是批次信息,reduction.save 是结果保存的名称
sce.all_int <- RunHarmony(sce.all_filter,reduction = "pca",group.by.vars = "group",reduction.save = "harmony")


##进行RunUMAP
sce.all_int <- RunUMAP(sce.all_int, reduction = "harmony", dims = 1:30,reduction.name = "umap")


## 04 继续常规操作---------------------------------------------------------------------------------------------------
sce.all_int <- FindNeighbors(sce.all_int, reduction = "harmony", dims = 1:30)
sce.all_int <- FindClusters(sce.all_int, resolution = seq(0.2,1.2,0.1))##设置多个resolution选择合适的resolution
metadata <- sce.all_int@meta.data### 多个分辨率的分群信息会保存在metadata中

resolutions <- seq(0.2,1,0.1)##设置多个resolution选择合适的resolution

DimPlot(sce.all_int,reduction="umap",group.by=paste0("RNA_snn_res.",resolutions),ncol = 3,label = TRUE)

## 05选定一个resolution，准备进行细胞注释---------------------------------------
### 选定分辨率0.6作为分群
sce.all_int@meta.data$seurat_clusters <- sce.all_int@meta.data$RNA_snn_res.0.6
Idents(sce.all_int) <- "seurat_clusters"

### 作图展示
DimPlot(sce.all_int, reduction = "umap", label = T)

### 学习Dimplot的两个重要参数 group.by 和split.by
DimPlot(sce.all_int, reduction = "umap", group.by = "group")### group.by, 在一张图中展示多个信息
DimPlot(sce.all_int, reduction = "umap", split.by = "group")### split.by, 分成多个信息来展示


names(sce.all_int@reductions)### 当前的reduction信息
### FeaturePlot 如果不限定reduction，它是这样运作的：先找umap, 再找tsne, 再找pca

FeaturePlot(sce.all_int, features = "MS4A1", order = TRUE)


