# 01 人工手动注释-------------------------------------------
#参考“生信漫漫学”的“单细胞数据集中常见marker基因汇总_2024版”（https://mp.weixin.qq.com/s/jzj8hLHdxHKqG5Re_KqVWg）


#以下为常见marker基因（生信技能树总结）：

T_cells <- c("CD3D", "CD3E", "CD8A", "PTPRC", "CD4", "CD2")
B_cells <- c("CD19", "CD79A", "MS4A1")
Plasma_cells <- c("IGHG1", "MZB1", "SDC1", "CD79A")
Monocytes_macrophages <- c("CD68", "CD163", "CD14")
NK_Cells <- c("FGFBP2", "FCGR3A", "CX3CR1")
Fibroblasts <- c("FGF7", "MME", "GSN", "LUM", "DCN")
Endothelial_cells <- c("PECAM1", "VWF")
Epi_tumor <- c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
Immune <- c("CD45", "PTPRC")
Epithelial_cancer <- c("EpCAM", "EPCAM")
Stromal <- c("CD10", "MME", "CD31", "PECAM1")
Mast <- c("CPA3", "CST3", "KIT", "TPSAB1", "TPSB2", "MS4A2")
#髓系
macrophages <- c("Adgre1", "Cd14", "Fcgr3", "LYZ", "CD68", "CD163")
cDCs <- c("Xcr1", "Flt3", "Ccr7", "CD1E")
pDCs <- c("Siglech", "Clec10a", "Clec12a", "CLEC4C")
monocytes <- c( "SPN", "CD300E")
neutrophils <- c("CSF3R", "S100A8", "CXCL3")

#下面是自己根据糖尿病足整理的marker基因：

B_cells <- c("CD19", "CD79A", "MS4A1")
T_cells <- c("CD3D", "CD3E", "CD8A", "PTPRC", "CD4", "CD2")
Monocytes <- c( "SPN", "CD300E")
Macrophages <- c("ADGRE1", "CD14","LYZ", "CD68", "CD163")
NK_cells <- c("FGFBP2", "FCGR3A", "CX3CR1")
Mast_cells <- c("CPA3", "CST3", "KIT", "TPSAB1", "TPSB2", "MS4A2")
Plasma_cells <- c("IGHG1", "MZB1", "SDC1")
DC_cells <- c("CLEC9A","CLEC10A","CST3","CADM1","COTL1","LYZ","DMXL2","FCER1A")
Neutrophils_cell <- c("CSF3R", "S100A8", "CXCL3")

Keratinocyte <- c("KRT1","KRT10","KRT5","KRT14")
Smooth_muscle_cell <- c("TAGLN","ACTA2")
Fibroblast_cell <- c("DCN","CFD","COL1A1","COL3A1")
vascular_endothelial_cell <- c("ACKR1","VEGFA","PECAM1","VWF")
Pericyte_cell <-  c("CSPG4", "PDGFRB", "TRPC6", "RGS5", "ABCC9", "KCNJ8", "NOTCH3", "HIGD1B", "COX4I2", "KCNK3")
Basal_cell <- c("KRT5","KRT14","CDH3","KRT15")

# 将所有 marker 基因整合到一个列表中
scRNA_marker_list1 <- list(
  SMCs_cell = Smooth_muscle_cell,
  KC_cell = Keratinocyte,
  Fibroblast_cell = Fibroblast_cell,
  Endothelial_cell = vascular_endothelial_cell,
  Pericyte_cell = Pericyte_cell
)

scRNA_marker_list2 <- list(
  T_cell = T_cells,
  B_cell = B_cells,
  Mono_cell = Monocytes,
  Macro_cell = Macrophages,
  NK_cell = NK_cells,
  Mast_cell = Mast_cells,
  Plasma_cell = Plasma_cells

)


scRNA_marker_list3 <- list(
  Neutrophils_cell = Neutrophils_cell,
  DC_cell = DC_cells
  
)
# 查看列表内容
print(scRNA_marker_list1)
print(scRNA_marker_list2)
print(scRNA_marker_list3)


p1 <- DotPlot(sce.all_int,assay="RNA",features=scRNA_marker_list1,
           group.by='RNA_snn_res.0.6'
)+theme(axis.text.x=element_text(angle=45,
                                     vjust=0.5,hjust=0.5))


p2 <- DotPlot(sce.all_int,assay="RNA",features=scRNA_marker_list2,
              group.by='RNA_snn_res.0.6'
)+theme(axis.text.x=element_text(angle=45,
                                 vjust=0.5,hjust=0.5))

p3 <- DotPlot(sce.all_int,assay="RNA",features=scRNA_marker_list3,
              group.by='RNA_snn_res.0.6'
)+theme(axis.text.x=element_text(angle=45,
                                 vjust=0.5,hjust=0.5))
p4 <- DotPlot(sce.all_int,assay="RNA",features=Basal_cell,
              group.by='RNA_snn_res.0.6'
)+theme(axis.text.x=element_text(angle=45,
                                 vjust=0.5,hjust=0.5))
p1
p2
p3
p4

