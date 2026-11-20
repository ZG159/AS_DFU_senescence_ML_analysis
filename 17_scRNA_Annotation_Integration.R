#对人工注释和在线注释进行合并
ACT_annotation
print(ACT_annotation,n=16)

#根据两种注释方法，确定细胞大群的数目

#当然也可以用函数RenameIdents，参见"guozi/guozi03_multiple sample/scripts/03_celltype_annotaion.R"代码文件
#在这里，我们创建一个矩阵celltype
#参考（生信漫漫学：手动注释结合在线注释网站使用）https://mp.weixin.qq.com/s/TCDuIjKLlSjwOffWkla3Xg
celltype=data.frame(ClusterID=0:20,
                    celltype=0:20)
#定义细胞亚群
celltype[celltype$ClusterID%in%c(0,16),2]='Keratinocyte'
celltype[celltype$ClusterID%in%c(1),2]=' Basal cell'       
celltype[celltype$ClusterID%in%c(2),2]='Fibroblast cell'
celltype[celltype$ClusterID%in%c(3,13),2]='Pericyte'      #Pericyte = Pericyte_cell
celltype[celltype$ClusterID%in%c(4),2]='T cell'
celltype[celltype$ClusterID%in%c(5,20),2]='Monocyte'
celltype[celltype$ClusterID%in%c(6,14,17,18,19),2]='Endothelial cell'
celltype[celltype$ClusterID%in%c(7,13),2]='DC cell'
celltype[celltype$ClusterID%in%c(8),2]='Plasma cell'
celltype[celltype$ClusterID%in%c(9),2]='Macrophages cell'
celltype[celltype$ClusterID%in%c(10),2]='B cell'
celltype[celltype$ClusterID%in%c(11),2]='Neutrophil cell'
celltype[celltype$ClusterID%in%c(12),2]='SMC cell'
celltype[celltype$ClusterID%in%c(15),2]='Mast cell'


#现在把celltype的标签加到sce.all.int中去
sce.all_int@meta.data$celltype="NA"

for(i in 1:nrow(celltype)){
  print(i)
  sce.all_int@meta.data[which(sce.all_int@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
table(sce.all_int$celltype) #看一下注释的细胞类型



#完成绘图
Idents(sce.all_int)=sce.all_int$celltype #为了绘图
DimPlot(sce.all_int,group.by = "celltype",label = TRUE)

saveRDS(sce.all_int,file = "out/sce.all.int.Rdata")

