
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
dir='resource/GSE245703/'
samples=list.files(dir)
samples

sceList=lapply(samples,function(pro){
  print(pro)
  sce=CreateSeuratObject(counts=Read10X_h5(file.path(dir,pro)),
                          project=gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro)),
                          min.cells=5,
                          min.features=300)#一个细胞至少有300个基因表达，一个基因至少在5个细胞中表达
  return(sce)
})
sceList
samples


sce.all=merge(x=sceList[[1]],
              y=sceList[-1],
              add.cell.ids= gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',samples)))

dim(sce.all)   #查看基因数和细胞总数

table(sce.all@meta.data$orig.ident)  #查看每个样本的细胞数


#添加分组
sce.all$group <- sce.all$orig.ident

#查看数据
table(sce.all@meta.data$group)
table(sce.all@meta.data$orig.ident)
table(sce.all_int@meta.data$seurat_clusters)


