# 2025/07/10
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(patchwork)
library(ggplot2)
library(harmony)
library(pheatmap)
library(plyr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicFeatures)
library(DESeq2)
library(ggpubr)
library(DoubletFinder)
library(ReactomePA)
library(tidyverse)
library(biomaRt)
library(enrichplot)
library(ktplots)
library(viridis)
setwd("~/data/wuhao/scdata/Rrun/AllFigures")
options(Seurat.object.assay.version = "v3")
source("~/Code/LJ.function.R")

#===============================================================================
#================================= Figure 2 ====================================
#===============================================================================
setwd("~/data/wuhao/scdata/Rrun/Figure2")
options(Seurat.object.assay.version = "v3")
source("~/Code/LJ.function.R")

# TE_spheroid load
TEcounts<-Read10X(data.dir="~/data/wuhao/scdata/cellrangerout/TE_spheroid/outs/filtered_feature_bc_matrix")
dataTE<-CreateSeuratObject(counts = TEcounts, min.features = 200,project = "TE_spheroid")
dataTE@meta.data$tech<-"10xseq"
dataTE@meta.data$sampletype<-"TE_spheroid"
dataTE[["percent.mt"]] <- PercentageFeatureSet(dataTE, pattern = "^MT-")
dataTE <- subset(x = dataTE, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 10 & nCount_RNA < 250000)
plotQC(dataTE,"cut")
dataTE<-RMdoublets(dataTE)# rm doublets

# ICM load 
ICMcounts<-Read10X(data.dir="~/data/wuhao/scdata/cellrangerout/ICM/outs/filtered_feature_bc_matrix")
dataICM<-CreateSeuratObject(counts = ICMcounts, min.cells = 3, min.features = 200,project = "ICM")
dataICM@meta.data$tech<-"10xseq"
dataICM@meta.data$sampletype<-"ICM"
dataICM[["percent.mt"]] <- PercentageFeatureSet(dataICM, pattern = "^MT-")
dataICM <- subset(x = dataICM, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 10 & nCount_RNA < 250000)
plotQC(dataICM,"cut")
dataICM<-RMdoublets(dataICM)# rm doublets

# Blastoid load
BLcounts<-Read10X(data.dir="~/data/wuhao/scdata/cellrangerout/Blastoid/outs/filtered_feature_bc_matrix")
dataBD<-CreateSeuratObject(counts = BLcounts, min.cells = 3, min.features = 200,project = "Blastoid")
dataBD@meta.data$tech<-"10xseq"
dataBD@meta.data$sampletype<-"Blastoid"
dataBD[["percent.mt"]] <- PercentageFeatureSet(dataBD, pattern = "^MT-")
dataBD <- subset(x = dataBD, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 10 & nCount_RNA < 250000)
plotQC(dataBD,"cut")
dataBD<-RMdoublets(dataBD)# rm doublets

# use well known marker genes do celltype annotation
selfall<-merge(x=dataTE,y = c(dataICM,dataBD))
setwd("../selfrun2")
selfall.list <- SplitObject(selfall, split.by = "sampletype")
for (i in 1:length(x = selfall.list)) {
  selfall.list[[i]] <- NormalizeData(object = selfall.list[[i]], verbose = FALSE)
  selfall.list[[i]] <- FindVariableFeatures(object = selfall.list[[i]],
                                            selection.method = "vst",  nfeatures = 3000, verbose = FALSE)
  print(i)
}
selfall.anchors <- FindIntegrationAnchors(object.list = selfall.list, dims = 1:40,k.anchor = 5,k.filter = 20)
selfall <- IntegrateData(anchorset = selfall.anchors, dims = 1:100)
selfall <- ScaleData(selfall)
selfall <- RunPCA(object = selfall,npcs=200)
# rm outlaier cells
pcause=150
selfall <- FindNeighbors(object = selfall, dims = 1:pcause)
selfall <- FindClusters(object = selfall, resolution = 0.6)
selfall <- RunUMAP(object = selfall, dims = 1:pcause)
pdf(file = paste("./selfall_umap_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = selfall,reduction = "umap",group.by = "sampletype",pt.size = .3))
dev.off()
selfall@meta.data$filter<-"1"
selfall@meta.data$filter[selfall@reductions$umap@cell.embeddings[,1]< -8]<-"2"
selfall2<-subset(x = selfall, subset= filter== "1")
pdf(file = paste("selfall2_umap_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
DimPlot(object = selfall2,reduction = "umap",group.by = "sampletype",pt.size = 0.5)
dev.off()
selfall2.anchors <- FindIntegrationAnchors(object.list = selfall2.list, dims = 1:40,k.anchor = 5,k.filter = 20)
selfall2 <- IntegrateData(anchorset = selfall2.anchors, dims = 1:100)
selfall2 <- ScaleData(selfall2)
selfall2 <- RunPCA(object = selfall2,npcs=200)
selfall2 <- RunHarmony(selfall2,reduction="pca", group.by.vars = "sampletype", reduction.save="harmony")
for (i in 1:15) {
  pcause=i*10
  selfall2 <- FindNeighbors(object = selfall2,reduction="harmony", dims = 1:pcause)
  selfall2 <- FindClusters(object = selfall2, resolution = 1)
  selfall2 <- RunUMAP(object = selfall2, reduction="harmony", dims = 1:pcause)
  pdf(file = paste("./umap/selfall2_harmony_cluster",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = selfall2,reduction = "umap",label = T,pt.size = 0.3))
  dev.off()
  pdf(file = paste("./umap/selfall2_harmony_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = selfall2,reduction = "umap",group.by = "sampletype",pt.size = .3))
  dev.off()
  pdf(file = paste("./umap/selfall2_harmony_blastoidsct",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = selfall2,reduction = "umap",group.by = "blastoidsct",pt.size = .3))
  dev.off()
}
pcause=100
samplecolor<-c("TE_spheroid"="#4BC9CD","ICM"="#B799C7","Blastoid"="#94BF67")
selfall2 <- FindNeighbors(object = selfall2,reduction="harmony", dims = 1:pcause)
selfall2 <- FindClusters(object = selfall2, resolution = 1)
selfall2 <- RunUMAP(object = selfall2, reduction="harmony", dims = 1:pcause)
pdf(file = paste("./selfall2_harmony1_cluster",pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = selfall2,reduction = "umap",label = T,pt.size = 0.3))
dev.off()
pdf(file = paste("./selfall2_harmony_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = selfall2,reduction = "umap",group.by = "sampletype",pt.size = .3)+
        scale_color_manual(values = samplecolor))
dev.off()
pdf(file = paste("./selfall2_harmony_blastoidsct",pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = selfall2,reduction = "umap",group.by = "blastoidsct",pt.size = .3))
dev.off()
# dotplot
markerall<-read.table("blastocysts_markers.txt",header = T,sep = "\t")
for (i in 1:ncol(markerall)) {
  imarker<-cbind(markerall[,i],rep(colnames(markerall)[i],length(markerall[,i])))
  if (i==1) {
    markers<-imarker
  }else{
    markers<-rbind(markers,imarker)
  }
}
markers<-na.omit(markers)
markers<-as.matrix(markers)
markers<-unique(markers)
markers<-markers[markers[,1]!="",]
# dotplot show marker genes expression info
callmarkerdot(markers,selfall2)
# 1. Figure S2I Marker gene staining
markers<-markers[markers[,1]%in%rownames(selfall2@assays$RNA),]
for (i in 1:nrow(markers)) {
  pdf(file = paste("GeneUMAP/Alldata_FeaturePlot",markers[i,1],markers[i,2],".pdf",sep = "_"),width =5.5 ,height = 5)
  print(FeaturePlot(selfall2,slot = "data", features = markers[i,1],
                    ncol = 1,alpha=1,order = T)+scale_colour_gradientn(colours =viridis_pal()(100))
  )
  dev.off()
}
pdf(file = paste("./selfall2_harmony1_cluster",pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = selfall2,reduction = "umap",label = T,pt.size = 0.3))
dev.off()

# celltype annotation
selfall2@meta.data$cellcluster <- selfall2@active.ident
selfall2@meta.data$celltype    <- plyr::mapvalues(x=selfall2@meta.data$cellcluster, 
                                                  from=c(0:18),
                                                  to=c("EPI1","Mural.TE","EPI1","Polar.TE","Polar.TE","Mural.TE",
                                                  "EPI2","Mural.TE","EPI1","TE.fusion_competent","Mural.TE","TE",
                                                  "Polar.TE","IntermediateCell2","Hypoblast","IntermediateCell1",
                                                  "EPI2","IntermediateCell3","EPI2"))
# 2. Figure S2F and S2G
mycolors3<-c("EPI1"="#EB6363","EPI2"="#FF7B8D","Hypoblast"="#DA842A",
             "TE"="#008BF6","TE.fusion_competent"="#9490FF",
             "Polar.TE"="#BEA6F3","Mural.TE"="#00B6EB",
             "Human embryo"="#BDBDBD","Unknown"="#7E7E7E","IntermediateCell1"="#A7FF9F", #759F88
             "IntermediateCell2"="#65BD8D","IntermediateCell3"="#3EAC70"
)
pdf(file = "./selfall2_harmony_celltype.pdf",width = 10,height = 8)
print(DimPlot(object = selfall2,reduction = "umap",group.by = "celltype",label = T,pt.size = .5)+ 
        scale_color_manual(values = mycolors3)
)
dev.off()
pdf(file = paste("./selfall2_harmony_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = selfall2,reduction = "umap",group.by = "sampletype",pt.size = .3)+
        scale_color_manual(values = samplecolor))
dev.off()
Blastoids_b2c<-cbind(names(selfall2$celltype),as.matrix(selfall2$celltype),as.matrix(selfall2$sampletype))
Blastoids_b2c[,1] <- sub(".{2}$", "", Blastoids_b2c[,1])
colnames(Blastoids_b2c)<-c("Original_Barcodes","Predicted_Celltype","Sample")
write.table(Blastoids_b2c,file = "Figure S2F Barcodes to celltype.txt",sep = "\t",row.names = T,col.names = T)


# 3. Figure 2I heatmap
gl2<-read.table("GenelistforFigure2I.txt",header = F,sep = "\t")
selfall3<-subset(selfall2,subset=celltype%in%c("EPI1","Mural.TE","Polar.TE","EPI2",
                                               "TE.fusion_competent","TE","Hypoblast"))
selfall3$celltype2<-plyr::mapvalues(x=selfall3$celltype,
                                    from=c("EPI1","EPI2","Mural.TE","Polar.TE","TE.fusion_competent","Hypoblast"),
                                    to=c("EPI","EPI","Mural TE","Polar TE","TE fusion-competent","HYP"))
TopMatrix<-selfall3@assays$RNA@data[,]
celltype_specieslist<-unique(selfall3$celltype2)
for (i in 1:length(celltype_specieslist)) {
  icell<-TopMatrix[,selfall3$celltype2==celltype_specieslist[i]]
  if (i==1) {
    MeanMatrix<-rowMeans(icell)
  }else{
    MeanMatrix<-cbind(MeanMatrix,rowMeans(icell))
  }
}
colnames(MeanMatrix)<-celltype_specieslist
MarkerMeanMatrix<-MeanMatrix[rownames(MeanMatrix)%in%gl2[,1],]
MarkerMeanMatrix<-MarkerMeanMatrix[gl2[,1],unique(gl2[,2])]
pdf(file = "SelfallMarkermean_heatmap.pdf",width = 6,height = 10)
pheatmap(mat = MarkerMeanMatrix,scale = "row",cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(100),border_color = "white")
dev.off()

# 4. Figure S2 H cell type percent for each time
merge_bar_all<-cbind(as.matrix(selfall2$sampletype),as.matrix(selfall2$celltype))

for (i in unique(selfall2$sampletype)) {
  sign<-table(merge_bar_all[merge_bar_all[,1]==i,2])
  singletype_table<-cbind(rep(i,length(sign)),names(sign),as.matrix(sign))
  rownames(singletype_table)<-NULL
  if (i==unique(selfall2$sampletype)[1]) {
    table_type<-singletype_table
  }else{
    table_type<-rbind(table_type,singletype_table)
  }
}
alltable<-table_type
colnames(alltable)<-c("Stage","Celltype","Cellnumber")
table_sample_type<-as.data.frame(alltable)
table_sample_type$Stage<-factor(table_sample_type$Stage,levels = c("Blastoid","ICM","TE_spheroid"))
table_sample_type$Celltype<-as.factor(table_sample_type$Celltype)                                                             ##
table_sample_type$Cellnumber<-as.numeric(as.matrix(table_sample_type$Cellnumber))                                                 ##
alltable_percent<-ddply(table_sample_type,"Stage",transform,percent_weight=Cellnumber/sum(Cellnumber)*100)
write.csv(alltable_percent,file = "Thisstudy_allcells_CelltypePercent.csv",row.names = F,col.names = T)

mycolors3<-c("EPI1"="#EB6363","EPI2"="#FF7B8D","Hypoblast"="#DA842A",
             "TE"="#008BF6","TE.fusion_competent"="#9490FF",
             "Polar.TE"="#BEA6F3","Mural.TE"="#00B6EB",
             "Human embryo"="#BDBDBD","Unknown"="#7E7E7E","IntermediateCell1"="#A7FF9F", #759F88
             "IntermediateCell2"="#65BD8D","IntermediateCell3"="#3EAC70")

alltable_percent$Celltype<-factor(alltable_percent$Celltype,
                                  levels = rev(c("IntermediateCell1","IntermediateCell2","IntermediateCell3",
                                                 "EPI1","EPI2","Hypoblast","TE","Mural.TE","Polar.TE",
                                                 "TE.fusion_competent")))
p_bar2<-ggplot(alltable_percent,aes(x=Stage,y=percent_weight,fill=Celltype))+geom_bar(stat = "identity",width= 0.6)+
  theme_bw()+scale_fill_manual(values = mycolors3) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
  coord_cartesian(clip = 'off')
ggsave("This_study_allcells_celltype_percent.pdf", p_bar2, width=5 ,height=4)

alltable_percent$id<-alltable_percent$Celltype
alltable_percent$id<-as.factor(alltable_percent$id)
alltable_percent$id<-plyr::mapvalues(x=alltable_percent$id, 
                                     from=levels(alltable_percent$id),
                                     to=c(1:length(levels(alltable_percent$id))))
library(ggalluvial)
pdf(file = "This_study_allcells_celltype_alluvial.pdf",width = 6,height = 4)
ggplot(alltable_percent,
       aes(x = Stage, stratum = Celltype, alluvium = id, y = percent_weight,
           fill = Celltype, label = Celltype)) +
  scale_x_discrete(expand = c(.1, 0)) +
  geom_flow(stat = "alluvium", aes(flow = percent_weight), width = 1/3) +
  geom_stratum(alpha = .9, width = 1/3) +
  theme(legend.position = "none")+
  expand_limits(x = 1, y = 0) +theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = mycolors3)
dev.off()


# Figure 2H Downsampling and integration with human blastocysts data
setwd("./downsampling")
# load Petropoulos_2016_PMID_27062923 data
data1<-read.table("~/data/CellLines/humanE3_E7/Petropoulos_2016_PMID_27062923/counts.txt",sep = "\t",header = T,row.names = "X")
meta.data<-read.table("~/data/CellLines/humanE3_E7/Petropoulos_2016_PMID_27062923/Petropoulos_2016_PMID_27062923.meta.tsv",
                      sep = "\t",header = T)
data1<-count2tpm(data1,idType="SYMBOL")
meta.data_filter<-meta.data[!meta.data$devTime%in%c("E3","E4"),]
rownames(meta.data_filter)<-meta.data_filter$cell
data1<-data1[,rownames(meta.data_filter)]
data1_spa <- Matrix(as.matrix(data1), sparse = TRUE)
blastocyst1<-CreateSeuratObject(counts = data1_spa,meta.data = meta.data_filter, min.cells = 3, min.features = 200,project = "Petropoulos")
blastocyst1@meta.data$tech<-"10xseq"
blastocyst1@meta.data$sampletype<-"Petropoulos"
blastocyst1[["percent.mt"]] <- PercentageFeatureSet(blastocyst1, pattern = "^MT-")
blastocyst1 <- subset(x = blastocyst1, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 ) # & percent.mt < 15)
plotQC(blastocyst1,"cut")
blastocyst1<-RMdoublets(blastocyst1)
# load data Yanagida et al. 2021 PRJNA720968
# only use day 5 and day 6 cells and map embryo original celltype
data2<-read.table("~/data/CellLines/humanE3_E7/Yanagida_2021_PMID_33957081/Yanagida_2021_PMID_33957081.counts",
                  sep = ",",header = T,row.names = "Gene")
meta.data<-read.table("~/data/CellLines/humanE3_E7/Yanagida_2021_PMID_33957081/Yanagida_2021_PMID_33957081.meta.tsv",
                      sep = "\t",header = T)
data2<-count2tpm(data2,idType="SYMBOL")
rownames(meta.data)<-meta.data$cell
meta.data_filter<-meta.data[!meta.data$devTime%in%c("Gata3_D3","Gata3_D4"),]
rownames(meta.data_filter)<-meta.data_filter$cell
data2<-data2[,rownames(meta.data_filter)]
data2_spa <- Matrix(as.matrix(data2), sparse = TRUE)
blastocyst2<-CreateSeuratObject(counts = data2_spa,meta.data = meta.data, min.cells = 3, min.features = 200,project = "Yanagida")
blastocyst2@meta.data$tech<-"10xseq"
blastocyst2@meta.data$sampletype<-"Yanagida"
blastocyst2[["percent.mt"]] <- PercentageFeatureSet(blastocyst2, pattern = "^MT-")
blastocyst2 <- subset(x = blastocyst2, subset = nFeature_RNA > 200 & nFeature_RNA < 15000)# & percent.mt < 15)
blastocyst2 <- NormalizeData(object = blastocyst2, verbose = FALSE)
blastocyst2 <- FindVariableFeatures(object = blastocyst2,selection.method = "vst",nfeatures = 3000, verbose = FALSE)
plotQC(blastocyst2,"cut")
blastocyst2<-RMdoublets(blastocyst2)

refall<- merge(x = blastocyst1, y = blastocyst2)
refall$raw_annotation2 <- refall$raw_annotation
refall$raw_annotation2 <- plyr::mapvalues(x=refall$raw_annotation2, 
                                          from=c("EPI.PrE.INT"),
                                          to=c("Epiblast"))
set.seed(1)
downsampledcells<-sample(colnames(selfall2),ncol(refall),replace = F)
selfall_ds<-subset(selfall2,cells=downsampledcells)
selfall_ds2<-CreateSeuratObject(counts = selfall_ds@assays$RNA$counts,meta.data = selfall_ds@meta.data)
dsall<-merge(selfall_ds2,refall)
DefaultAssay(dsall)<-"RNA"
dsall.list <- SplitObject(dsall, split.by = "sampletype")
for (i in 1:length(x = dsall.list)) {
  dsall.list[[i]] <- NormalizeData(object = dsall.list[[i]], verbose = FALSE)
  dsall.list[[i]] <- FindVariableFeatures(object = dsall.list[[i]],
                                          selection.method = "vst",  nfeatures = 2000, verbose = FALSE)
  print(i)
}
dsall.anchors <- FindIntegrationAnchors(object.list = dsall.list, dims = 1:40,k.anchor = 5,k.filter = 20)
dsall <- IntegrateData(anchorset = dsall.anchors, dims = 1:100)
dsall <- ScaleData(dsall)
dsall <- RunPCA(object = dsall,npcs=200)
dsall <- RunHarmony(dsall,reduction="pca", group.by.vars = "sampletype", reduction.save="harmony")

for (i in 8:15) {
  pcause=i*10
  dsall <- FindNeighbors(object = dsall,reduction="harmony", dims = 1:pcause)
  dsall <- FindClusters(object = dsall, resolution = 0.6)
  dsall <- RunUMAP(object = dsall,reduction="harmony", dims = 1:pcause)
  pdf(file = paste("./umap/dsall_harmony_cluster",pcause,".pdf",sep = ""),width = 8,height = 6)
  print(DimPlot(object = dsall,reduction = "umap",label = T,pt.size = 1))
  dev.off()
  pdf(file = paste("./umap/dsall_harmony_sample",pcause,".pdf",sep = ""),width = 8,height = 6)
  print(DimPlot(object = dsall,reduction = "umap",group.by = "sampletype",pt.size = 1))
  dev.off()
  pdf(file = paste("./umap/dsall_harmony_blastoidsct",pcause,".pdf",sep = ""),width = 8,height = 6)
  print(DimPlot(object = dsall,reduction = "umap",group.by = "celltype",label = T,pt.size = 1))
  dev.off()
  pdf(file = paste("./umap/dsall_harmony_rawannotation",pcause,".pdf",sep = ""),width = 8,height = 6)
  print(DimPlot(object = dsall,reduction = "umap",group.by = "raw_annotation2",label = T,pt.size = 1))
  dev.off()
}
pcause = 80
dsall <- FindNeighbors(object = dsall,reduction="harmony", dims = 1:pcause)
dsall <- FindClusters(object = dsall, resolution = 8)
dsall <- RunUMAP(object = dsall, reduction="harmony", dims = 1:pcause)
dsall@meta.data$datatype<-"Reference"
dsall$datatype[dsall$sampletype%in%c(names(samplecolor))]="This Study"
pdf(file = paste("./dsall3_harmony_cluster",pcause,".pdf",sep = ""),width = 7,height = 4)
print(DimPlot(object = dsall,reduction = "umap",label = T,pt.size = 0.5,shape.by = "datatype")+
        scale_shape_manual(values = c(17,16)))
dev.off()
samplecolorref<-c(samplecolor,"Petropoulos"="#FF9288","Yanagida"="#FF81D6")
# 5. Figure 2H
dsall$raw_annotation2[is.na(dsall$raw_annotation2)]<-"This Study"
dsall$raw_annotation2<-factor(dsall$raw_annotation2,levels = rev(c("Epiblast","Hypoblast","ICM","Morula","Prelineage","TE","Unknown","This Study")))

mycolors<-c("Epiblast"="#EB6363","Hypoblast"="#DA842A","ICM"="#5D915E","Morula"="#3EAC70",
            "Prelineage"="lightgreen","TE"="#008BF6","Unknown"="#7E7E7E","This Study"="#BDBDBD")
dsall$datatype[dsall$datatype=="Reference"]<-"Human embryo"
pdf(file = paste("./dsall_harmony_reference_celltype",pcause,".pdf",sep = ""),width = 6.5,height = 5)
DimPlot(object = dsall,reduction = "umap",group.by = "raw_annotation2",label = T,
              order = c("Epiblast","Hypoblast","ICM","Morula","Prelineage","TE","Unknown","This Study"),pt.size = 2.5,shape.by = "datatype")+
        scale_shape_manual(values = c(17,16))+
  scale_color_manual(values = mycolors)
dev.off()
dsall$celltype[is.na(dsall$celltype)]<-"Human embryo"
mycolors2<-c("EPI1"="#EB6363","EPI2"="#FF7B8D","Hypoblast"="#DA842A","Intermediate.cell"="#3EAC70",
            "TE"="#008BF6","TE.fusion_competent"="#9490FF",
            "Polar.TE"="#BEA6F3","Mural.TE"="#00B6EB",
            "Human embryo"="#BDBDBD","Unknown"="#7E7E7E"
            )
pdf(file = paste("./dsall_harmony_this_study_celltype2",pcause,".pdf",sep = ""),width = 7,height = 5)
DimPlot(object = dsall,reduction = "umap",group.by = "celltype",label = T,
        order = c("EPI1","EPI2","Hypoblast","Intermediate.cell","TE","Mural.TE","Polar.TE","TE.fusion_competent","Unknown","Human embryo"),
        pt.size = 2.2,shape.by = "datatype")+
        scale_shape_manual(values = c(17,16))+
  scale_color_manual(values = mycolors2)
dev.off()
#===============================================================================
#================================ Figure 3 =====================================
#===============================================================================
setwd("~/data/wuhao/scdata/Rrun/AllFigures/Figure3")
# load Chimeric Blastoids and Mouse ICM (10x)
# load Chimeric Blastoids mouse
CBmcounts<-Read10X(data.dir="~/data/wuhao/scdata/cellrangerout/Chimeric_Blastoids_mouse/outs/filtered_feature_bc_matrix")
rownames(CBmcounts)<-paste0(rownames(CBmcounts),"-mouse")
CBm<-CreateSeuratObject(counts = CBmcounts, min.features = 200,project = "CB_mouse")
CBm@meta.data$tech<-"10xseq"
CBm@meta.data$sample2<-"CB_mouse"
# load Chimeric Blastoids human
CBhcounts<-Read10X(data.dir="~/data/wuhao/scdata/cellrangerout/Chimeric_Blastoids_human/outs/filtered_feature_bc_matrix")
rownames(CBhcounts)<-paste0(rownames(CBhcounts),"-human")
CBh_org<-CreateSeuratObject(counts = CBhcounts, min.features = 200,project = "CB_human")
CBh_org@meta.data$tech<-"10xseq"
CBh_org@meta.data$sample2<-"CB_human"
# species identify
CBh_org@meta.data$mouseCounts<-0
comcells<-intersect(names(CBh_org$mouseCounts),names(CBm$nCount_RNA))
CBh_org$mouseCounts[comcells]<-CBm$nCount_RNA[comcells]
CBh_org$species<-"unknown"
CBh_org$species[CBh_org$nCount_RNA*0.2 >CBh_org$mouseCounts]<-"Human"
CBh_org$species[CBh_org$mouseCounts*0.2>CBh_org$nCount_RNA]<-"Mouse"
CBm@meta.data$species<-"unknown"
CBm$species[names(CBh_org$species)[CBh_org$species=="Mouse"]]<-"Mouse"
CBm$species[setdiff(names(CBm$species),comcells)]<-"Mouse"
CBm2<-subset(CBm,subset=species=="Mouse")
CBm2[["percent.mt"]] <- PercentageFeatureSet(CBm2, pattern = "^mt-")
plotQC(CBm2,"raw")
CBm2 <- subset(x = CBm2, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
plotQC(CBm2,"cut")
CBm2<-RMdoublets(CBm2)
CBh2<-subset(CBh_org,subset=species=="Human")
CBh2[["percent.mt"]] <- PercentageFeatureSet(CBh2, pattern = "^MT-")
CBh2 <- subset(x = CBh2, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
plotQC(CBh2,"cut")
CBh2<-RMdoublets(CBh2)
# load mouse ICM
ICMmcounts<-Read10X(data.dir="~/data/wuhao/scdata/cellrangerout/MouseICM/outs/filtered_feature_bc_matrix")
rownames(ICMmcounts)<-paste0(rownames(ICMmcounts),"-mouse")
ICMm<-CreateSeuratObject(counts = ICMmcounts, min.features = 200,project = "ICM_mouse")
ICMm@meta.data$tech<-"10xseq"
ICMm@meta.data$sample2<-"ICM_mouse"
ICMm[["percent.mt"]] <- PercentageFeatureSet(ICMm, pattern = "^mt-")
ICMm <- subset(x = ICMm, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
plotQC(ICMm,"cut")
# rm doublets
ICMm<-RMdoublets(ICMm)
# merge alldata
commongenem<-intersect(rownames(CBm2),rownames(ICMm))
CBm2 <- subset(CBm2, features = commongenem)
ICMm <- subset(ICMm, features = commongenem)
CBall<-merge(x=CBm2,y=c(CBh2,ICMm))
plotQC(CBall,"cut")
CBall <- NormalizeData(object = CBall, verbose = FALSE)
CBall <- FindVariableFeatures(object = CBall,
                              selection.method = "vst",  nfeatures = 2000, verbose = FALSE)
CBall <- ScaleData(CBall)
CBall <- RunPCA(object = CBall,npcs=200)
pdf(file = "CBall_pca_use.pdf")
ElbowPlot(CBall, ndims = 100)
dev.off()
CBall <- RunHarmony(CBall,reduction="pca", group.by.vars = "sample2", reduction.save="harmony")
for (i in 1:10) {
  pcause=i*10
  CBall <- FindNeighbors(object = CBall,reduction="harmony", dims = 1:pcause)
  CBall <- FindClusters(object = CBall, resolution = 0.6)
  CBall <- RunUMAP(object = CBall, reduction="harmony", dims = 1:pcause)
  pdf(file = paste("./umap/CBall_rmintermediate_harmony_cluster",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = CBall,reduction = "umap",label = T,pt.size = 0.3))
  dev.off()
  pdf(file = paste("./umap/CBall_rmintermediate_harmony_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = CBall,reduction = "umap",group.by = "sample2",pt.size = .3))
  dev.off()
}
pcause=30
CBall <- FindNeighbors(object = CBall,reduction="harmony", dims = 1:pcause)
CBall <- FindClusters(object = CBall, resolution = 0.6)
CBall <- RunUMAP(object = CBall,reduction="harmony", dims = 1:pcause)
pdf(file = paste("CBall_rmintermediate_harmony_cluster",pcause,"_0.6.pdf",sep = ""),width = 6,height = 5)
DimPlot(object = CBall,reduction = "umap",label = T,pt.size = 0.3,order =names(CBall$sample2)[CBall$sample2=="CB_mouse"] )
dev.off()
pdf(file = paste("CBall_rmintermediate_harmony_sample",pcause,".pdf",sep = ""),width = 6,height = 5)
DimPlot(object = CBall,reduction = "umap",group.by = "sample2",pt.size = .5,order =names(CBall$sample2)[CBall$sample2=="CB_mouse"])
dev.off()
CBall$seuratcluster<-CBall@active.ident
CBall2<-subset(CBall,subset=seuratcluster!=15)
# Figure 3H and Figure S3A
pdf(file = paste("CBall2_rmintermediate_harmony_cluster",pcause,"_0.6.pdf",sep = ""),width = 6,height = 5)
DimPlot(object = CBall2,reduction = "umap",label = T,pt.size = 0.3,order =names(CBall$sample2)[CBall$sample2=="CB_mouse"] )
dev.off()
pdf(file = paste("CBall2_rmintermediate_harmony_sample",pcause,".pdf",sep = ""),width = 6,height = 5)
DimPlot(object = CBall2,reduction = "umap",group.by = "sample2",pt.size = .5,order =names(CBall$sample2)[CBall$sample2=="CB_mouse"])
dev.off()
CBall2$celltype_show<- plyr::mapvalues(x=CBall2@meta.data$seuratcluster,
                                       from=c(0:14),
                                       to=c("Mouse cells","TE 1","Polar TE","Mouse cells",
                                       "Mouse cells","TE 2","Mural TE","Mouse cells","Mouse cells","Mouse cells",
                                       "TE fusion-competent","TE 2","Mouse cells","hESC","Mouse cells"))
colorshow<-c("hESC"="#F76C6C","TE 1"="#0085ED","TE 2"="#96A841","Mural TE"="#43B158","Polar TE"="#6CC8CD",
             "TE fusion-competent"="#9587BE","Mouse cells"="#BDBDBD")
pdf(file = "CBall2_celltypeshow.pdf",width = 7,height = 5)
DimPlot(object = CBall2,reduction = "umap",group.by = "celltype_show",pt.size = .3,label = T)+
  scale_color_manual(values = colorshow)
dev.off()

CBall_b2c<-cbind(names(CBall2$celltype_show),as.matrix(CBall2$celltype_show),as.matrix(CBall2$sample2))
CBall_b2c[,1] <- sub(".{2}$", "", CBall_b2c[,1])
colnames(CBall_b2c)<-c("Original_Barcodes","Predicted_Celltype","Sample")
write.table(CBall_b2c,file = "Figure 3H Barcodes to celltype.txt",sep = "\t",row.names = T,col.names = T)
# Figure S3B identify human and mouse cells
Count_mat<-CBh_org@meta.data[,c(2,6,7)]
Count_mat<-Count_mat[intersect(rownames(Count_mat),substr(colnames(CBall), 1, nchar(colnames(CBall)) - 2)),]
Count_mat1<-as.matrix(cbind(Count_mat[,c(1,3)],"HumanReference"))
Count_mat2<-as.matrix(cbind(Count_mat[,c(2,3)],"MouseReference"))
colnames(Count_mat1)<-NULL
colnames(Count_mat2)<-NULL
Count_mat3<-rbind(Count_mat1,Count_mat2)
rownames(Count_mat3)<-NULL
Count_mat3<-as.data.frame(Count_mat3)
colnames(Count_mat3)<-c("AllCounts","Specie","Reference")
Count_mat3$AllCounts<-as.numeric(Count_mat3$AllCounts)
summary_df <- Count_mat3 %>%
  group_by(Reference, Specie) %>%
  summarise(mean_count = mean(AllCounts),
            sd_count = sd(AllCounts),
            .groups = "drop")
pdf(file = "celltype_identification_barplot.pdf",width = 7,height = 6)
ggplot(summary_df, aes(x = Reference, y = mean_count, fill = Specie)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
                position = position_dodge(width = 0.9),
                width = 0.3) +
  theme_classic() +
  labs(x = "Reference", y = "Mean Count", title = "Counts by Reference and Specie")
dev.off()
# human cell annotation 
markerall<-read.table("../mergeblastoids013_markers.txt",header = T,sep = "\t")
for (i in 1:ncol(markerall)) {
  imarker<-cbind(markerall[,i],rep(colnames(markerall)[i],length(markerall[,i])))
  if (i==1) {
    markers<-imarker
  }else{
    markers<-rbind(markers,imarker)
  }
}
markersh<-na.omit(markers)
markersh<-as.matrix(markersh)
markersh<-unique(markersh)
markersh2<-markersh[markersh[,1]!="",]
markersh2[,1]<-paste0(markersh2[,1],"-human")
callmarkerdot(markersh2,CBall,group.by="seuratcluster",fname="CBall-human_cluster")
markersh3<-markersh2[markersh2[,1]%in%rownames(CBall2),]
for (i in 1:nrow(markersh3)) {
  pdf(file = paste("./Featureplot/CBall_human_markergenes_staining",markersh3[i,1],".pdf",sep = ""),width =5 ,height = 4.5)
  print(FeaturePlot(CBall2,slot = "data", features = markersh3[i,1],cols = viridis_pal()(100)[21:100],order = T,pt.size=0.5,alpha=0.8))
  dev.off()
}
CBall2$celltype<- plyr::mapvalues(x=CBall2@meta.data$seuratcluster,
                                                  from=c(13,1,5,11,6,2,10),
                                                  to=c("EPI_h","TE1_h","TE2_h","TE2_h","Mural_h","PolarTE_h","TE_fusion_h"))

pdf(file = paste("CBall2_celltype",pcause,"_test.pdf",sep = ""),width = 8,height = 6)
DimPlot(object = CBall2,reduction = "umap",group.by = "celltype",pt.size = .3,label = T)
dev.off()
# Figure 3I human celltype dotplot
m3i<-read.table("GenelistforFigure3I.txt",header = F)
data.anno<-unique(m3i)
data.anno[,2]<-paste0(unique(data.anno[,2]),"-human")
data.usage <- DotPlot(CBall2,features = data.anno[,2], assay = 'RNA',group.by = "celltype")$data
data.usage<-data.usage[!data.usage$id%in%c(0,3,4,7,8,9,12,14),]
colnames(data.anno)<-c("label","features.plot")
data.anno<-as.data.frame(data.anno)
df.plot <- plyr::join(data.usage,data.anno)
df.plot<-na.omit(df.plot)
df.plot$features.plot<-sub("-human$","",df.plot$features.plot)
df.plot$id<-factor(df.plot$id)
df.plot$id<- plyr::mapvalues(x=df.plot$id,
                                  from=c("TE1_h","PolarTE_h","TE2_h","Mural_h",
                                         "TE_fusion_h","EPI_h"),
                                  to=c("TE 1","Polar TE","TE 2","Mural TE","TE fusion-competent","hESC"))
df.plot$id<-factor(df.plot$id,levels = c("hESC","TE 1","TE 2","Mural TE","Polar TE","TE fusion-competent"))
df.plot$features.plot<-factor(df.plot$features.plot,levels = m3i[,2] )

p <- ggplot(df.plot,aes(x=features.plot,y = as.numeric(id),size = pct.exp/100, color = avg.exp.scaled))+
  geom_point() + scale_size("% detected", range = c(0,6)) +  
  scale_color_gradientn(colours = c("white","red"),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") + cowplot::theme_cowplot() +
  ylab("") + xlab("Markers") + theme_bw() + 
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+ #复制 y轴 代替边框效果 
  facet_grid(~label, scales="free_x",space = "free")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Fig3Ihuman_dotplot.pdf",p,width =9,height = 4)
# mouse cell annotation 
markerall<-read.table("mergeblastoids013_markers.txt",header = T,sep = "\t")
for (i in 1:ncol(markerall)) {
  imarker<-cbind(markerall[,i],rep(colnames(markerall)[i],length(markerall[,i])))
  if (i==1) {
    markers<-imarker
  }else{
    markers<-rbind(markers,imarker)
  }
}
markers<-na.omit(markers)
markers<-as.matrix(markers)
markers<-unique(markers)
markersm<-markers[markers[,1]!="",]
library(stringr)
markersm[, 1] <- str_to_title(markersm[, 1])
markersm[,1]<-paste0(markersm[,1],"-mouse")
callmarkerdot(markersm,CBall2,group.by="seuratcluster",fname="CBall2-mouse_cluster")
CBall2$celltype<- plyr::mapvalues(x=CBall2@meta.data$celltype,
                                  from=c(0,3,4,7,8,12,14,9),
                                  to=c(rep("EPI_m",7),"PrE_m"))

pdf(file = paste("CBall2_celltype",pcause,".pdf",sep = ""),width = 8,height = 6)
DimPlot(object = CBall2,reduction = "umap",group.by = "celltype",pt.size = .3,label = T)
dev.off()

# Figure 3J CBm, ICMm, merge with mouse blastocysts 3.5 and 4.5 (PMID:30959515)
setwd("./CBmouse_withvivo")
datam3.53 <- read.table("./vivo3.54.5/Lib1-3_E3.5_counts.csv",
                       sep = ",",header=T,row.names = "X")
rownames(datam3.53)<-paste0("Lib1-3_E3.5_",sub("^b","",rownames(datam3.53)))
colnames(datam3.53)<-gsub("^MT\\.","MT-",colnames(datam3.53))
datam3.53_sp <- t(as.matrix(datam3.53))
datam3.53_sp <- as(datam3.53_sp, "dgCMatrix")
metadataall<-read.table("./vivo3.54.5/all_cells_metadata.csv",sep = ",",header = T,row.names = "index")
meta3.5<-metadataall[metadataall$Timepoint=="E3.5",]
meta3.53<-meta3.5[rownames(meta3.5)%in%rownames(datam3.53),]

m3.53<-CreateSeuratObject(counts = datam3.53_sp, meta.data = meta3.53, project = "Mousevivo")
m3.53@meta.data$tech<-"10xseq"
m3.53@meta.data$sample2<-"Mousevivo3.5"
m3.53[["percent.mt"]] <- PercentageFeatureSet(m3.53, pattern = "^MT-")
m3.53 <- subset(x = m3.53, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
m3.53 <- subset(m3.53,subset=Timepoint=="E3.5")
# load datam3.54
datam3.54 <- read.table("./vivo3.54.5/Lib1-4_E3.5_counts.csv",
                        sep = ",",header=T,row.names = "X")
rownames(datam3.54)<-paste0("Lib1-4_E3.5_",sub("^b","",rownames(datam3.54)))
colnames(datam3.54)<-gsub("^MT\\.","MT-",colnames(datam3.54))
datam3.54_sp <- t(as.matrix(datam3.54))
datam3.54_sp <- as(datam3.54_sp, "dgCMatrix")
meta3.54<-meta3.5[rownames(meta3.5)%in%rownames(datam3.54),]
m3.54<-CreateSeuratObject(counts = datam3.54_sp,meta.data = meta3.54, project = "Mousevivo")
m3.54@meta.data$tech<-"10xseq"
m3.54@meta.data$sample2<-"Mousevivo3.5"
m3.54[["percent.mt"]] <- PercentageFeatureSet(m3.54, pattern = "^MT-")
m3.54 <- subset(x = m3.54, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
m3.54 <- subset(m3.54,subset=Timepoint=="E3.5")
# load datam4.51
datam4.51 <- read.table("./vivo3.54.5/Lib1-1_E4.5_counts.csv",
                        sep = ",",header=T,row.names = "X")
rownames(datam4.51)<-paste0("Lib1-1_E4.5_",sub("^b","",rownames(datam4.51)))
colnames(datam4.51)<-gsub("^MT\\.","MT-",colnames(datam4.51))
datam4.51_sp <- t(as.matrix(datam4.51))
datam4.51_sp <- as(datam4.51_sp, "dgCMatrix")
meta4.5<-metadataall[metadataall$Timepoint=="E4.5",]
meta4.51<-meta4.5[rownames(meta4.5)%in%rownames(datam4.51),]
m4.51<-CreateSeuratObject(counts = datam4.51_sp, meta.data = meta4.51, project = "Mousevivo")
m4.51@meta.data$tech<-"10xseq"
m4.51@meta.data$sample2<-"Mousevivo4.5"
m4.51[["percent.mt"]] <- PercentageFeatureSet(m4.51, pattern = "^MT-")
m4.51 <- subset(x = m4.51, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
m4.51 <- subset(m4.51,subset=Timepoint=="E4.5")
# load datam4.5.2
datam4.52 <- read.table("./vivo3.54.5/Lib1-2_E4.5_counts.csv",
                        sep = ",",header=T,row.names = "X")
rownames(datam4.52)<-paste0("Lib1-2_E4.5_",sub("^b","",rownames(datam4.52)))
colnames(datam4.52)<-gsub("^MT\\.","MT-",colnames(datam4.52))
datam4.52_sp <- t(as.matrix(datam4.52))
datam4.52_sp <- as(datam4.52_sp, "dgCMatrix")
meta4.52<-meta4.5[rownames(meta4.5)%in%rownames(datam4.52),]
m4.52<-CreateSeuratObject(counts = datam4.52_sp,meta.data = meta4.52, project = "Mousevivo")
m4.52@meta.data$tech<-"10xseq"
m4.52@meta.data$sample2<-"Mousevivo4.5"
m4.52[["percent.mt"]] <- PercentageFeatureSet(m4.52, pattern = "^MT-")
m4.52 <- subset(x = m4.52, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
m4.52 <- subset(m4.52,subset=Timepoint=="E4.5")
ICMmcounts2<-ICMm@assays$RNA@counts
rownames(ICMmcounts2)<-toupper(gsub("-mouse$","",rownames(ICMmcounts2)))
ICMm3<-CreateSeuratObject(counts = ICMmcounts2, meta.data = ICMm@meta.data, project = "ICMmouse")
ICMm3$rmlist<-F
ICMm3$rmlist[sample(c(1:ncol(ICMm3)),ncol(ICMm3)*0.2)]<-T
ICMm4<-subset(ICMm3,subset=rmlist)
CBmcounts2<-CBm2@assays$RNA@counts
rownames(CBmcounts2)<-toupper(gsub("-mouse$","",rownames(CBmcounts2)))
CBm3<-CreateSeuratObject(counts = CBmcounts2, meta.data = CBm2@meta.data, project = "ChimericBmouse")
Mall<-merge(x=CBm3,y=c(ICMm4,m3.53,m3.54,m4.51,m4.52))
plotQC(Mall,"cut")
Mall.list <- SplitObject(Mall, split.by = "sample2")
for (i in 1:length(x = Mall.list)) {
  Mall.list[[i]] <- NormalizeData(object = Mall.list[[i]], verbose = FALSE)
  Mall.list[[i]] <- FindVariableFeatures(object = Mall.list[[i]],
                                        selection.method = "vst",  nfeatures = 2000, verbose = FALSE)
  print(i)
}
Mall.anchors <- FindIntegrationAnchors(object.list = Mall.list, dims = 1:60,k.anchor = 10,k.filter = 20)
Mall <- IntegrateData(anchorset = Mall.anchors, dims = 1:50,k.weight = 20)
Mall <- ScaleData(Mall)
Mall <- RunPCA(object = Mall,npcs=200)
pdf(file = "Mall_pca_use.pdf")
ElbowPlot(Mall, ndims = 100)
dev.off()
Mall <- RunHarmony(Mall,reduction="pca", group.by.vars = "sample2", reduction.save="harmony")
for (i in 1:10) {
  pcause=i*10
  Mall <- FindNeighbors(object = Mall,reduction="harmony", dims = 1:pcause)
  Mall <- FindClusters(object = Mall, resolution = 0.6)
  Mall <- RunUMAP(object = Mall, reduction="harmony", dims = 1:pcause)
  pdf(file = paste("./umap/Mall_rmintermediate_harmony_cluster",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = Mall,reduction = "umap",label = T,pt.size = 0.3))
  dev.off()
  pdf(file = paste("./umap/Mall_rmintermediate_harmony_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = Mall,reduction = "umap",group.by = "sample2",pt.size = .3))
  dev.off()
  pdf(file = paste("./umap/Mall_rmintermediate_harmony_celltype",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = Mall,reduction = "umap",group.by = "CellType",pt.size = .3))
  dev.off()
}
pcause=60
Mall <- FindNeighbors(object = Mall,reduction="harmony", dims = 1:pcause)
Mall <- FindClusters(object = Mall, resolution = 1)
Mall <- RunUMAP(object = Mall,reduction="harmony", dims = 1:pcause)
pdf(file = paste("Mall_rmintermediate_harmony_cluster",pcause,"_1.pdf",sep = ""),width = 6,height = 5)
DimPlot(object = Mall,reduction = "umap",label = T,pt.size = 0.1)
dev.off()
Mall$sample2<-as.factor(Mall$sample2)
pdf(file = paste("Mall_rmintermediate_harmony_sample",pcause,".pdf",sep = ""),width = 6,height = 5)
DimPlot(object = Mall,reduction = "umap",group.by = "sample2",pt.size = .5,order =c("CB_mouse","Mousevivo3.5","Mousevivo4.5","ICM_mouse"))
dev.off()
pdf(file = paste("./Mall_rmintermediate_harmony_celltype",pcause,".pdf",sep = ""),width = 6,height = 5)
print(DimPlot(object = Mall,reduction = "umap",group.by = "CellType",pt.size = .1))
dev.off()
Mall$seuratcluster<-Mall@active.ident
# Figure 3J
pdf(file = "Figure3J_sample_type.pdf",width = 6,height = 5)
DimPlot(object = Mall,reduction = "umap",group.by = "sample2",pt.size = .5,order =c("CB_mouse","Mousevivo3.5","Mousevivo4.5","ICM_mouse"))
dev.off()
# Fig S3D
Mall$CellType[Mall$CellType==""]="Unlabeled"
Mall$CellType[is.na(Mall$CellType)]="This Study"
mycolors<-c( "EPI" = "#F76C6C",
             "PrE" = "#00C853",
             "TE" = "#64B5F6",
             "Unlabeled"="#7E7E7E",
             "This Study"="#BDBDBD")
pdf(file = "./Figure S3D_reference_celltype.pdf",width = 6,height = 5)
DimPlot(object = Mall,reduction = "umap",group.by = "CellType",pt.size = 0.5,
              order = c("EPI","PrE","TE","Unlabeled","This Study"))+
  scale_color_manual(values = mycolors)
dev.off()
# Figure S3E all celltype
Mall$celltype2<- plyr::mapvalues(x=Mall@meta.data$seuratcluster,
                                 from=c(0:4,6,9,10,11,13,
                                        5,7,12,
                                        8),
                                 to=c(rep("EPI",10),"PrE","PrE","PrE","TE"))

pdf(file = "./Figure S3E_all_celltype.pdf",width = 6,height = 5)
DimPlot(object = Mall,reduction = "umap",group.by = "celltype2",pt.size = 0.5,
        order = c("EPI","PrE","TE","Unlabeled","This Study"))+
  scale_color_manual(values = mycolors)
dev.off()
# Figure 3K
m3k<-read.table("GenelistforFigure3k.txt",header = F,sep = "\t")
data.anno<-unique(m3k)
data.usage <- DotPlot(Mall,features = data.anno[,2], assay = 'RNA',group.by = "celltype2")$data
colnames(data.anno)<-c("label","features.plot")
data.anno<-as.data.frame(data.anno)
df.plot <- plyr::join(data.usage,data.anno)
df.plot<-na.omit(df.plot)
df.plot$id<-factor(df.plot$id,levels = c("EPI","PrE","TE"))
df.plot$features.plot<-factor(df.plot$features.plot,levels = m3k[,2] )
p <- ggplot(df.plot,aes(x=features.plot,y = as.numeric(id),size = pct.exp/100, color = avg.exp.scaled))+
  geom_point() + scale_size("% detected", range = c(0,6)) +  
  scale_color_gradientn(colours = c("white","red"),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") + cowplot::theme_cowplot() +
  ylab("") + xlab("Markers") + theme_bw() + 
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Fig3kmouse_dotplot.pdf",p,width =6,height = 4)
# Figure 3L
markers2<-m3k[,2]
for (i in 1:length(markers2)) {
  pdf(file = paste("./Featureplot/Mall_markergenes_staining",markers2[i],".pdf",sep = ""),width =5 ,height = 4.5)
  print(FeaturePlot(Mall,slot = "data", features = markers2[i],cols = viridis_pal()(100)[21:100],order = T,pt.size=0.5,alpha=0.8))
  dev.off()
}

# Figure S3E.2
merge_bar_all<-cbind(as.matrix(Mall$sample2),as.matrix(Mall$celltype2))

for (i in unique(Mall$sample2)) {
  sign<-table(merge_bar_all[merge_bar_all[,1]==i,2])
  singletype_table<-cbind(rep(i,length(sign)),names(sign),as.matrix(sign))
  rownames(singletype_table)<-NULL
  if (i==unique(Mall$sample2)[1]) {
    table_type<-singletype_table
  }else{
    table_type<-rbind(table_type,singletype_table)
  }
}
alltable<-table_type
colnames(alltable)<-c("Stage","Celltype","Cellnumber")
table_sample_type<-as.data.frame(alltable)
table_sample_type$Stage<-factor(table_sample_type$Stage,levels = c("CB_mouse","ICM_mouse","Mousevivo3.5","Mousevivo4.5"))
table_sample_type$Celltype<-as.factor(table_sample_type$Celltype)                                                             ##
table_sample_type$Cellnumber<-as.numeric(as.matrix(table_sample_type$Cellnumber))                                                 ##
alltable_percent<-ddply(table_sample_type,"Stage",transform,percent_weight=Cellnumber/sum(Cellnumber)*100)
write.csv(alltable_percent,file = "Mallsc_CelltypePercent.csv",row.names = F,col.names = T)
p_bar2<-ggplot(alltable_percent,aes(x=Stage,y=percent_weight,fill=Celltype))+geom_bar(stat = "identity",width= 0.6)+
  theme_bw()+scale_fill_manual(values = mycolors) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
  coord_cartesian(clip = 'off')
ggsave("FigureS3E2.Mall_celltype_percent.pdf", p_bar2, width=6 ,height=4)
# ==============================================================================
# ================================= Figure 4 ===================================
# ==============================================================================
# run CellphoneDB
Blastoids_data<-selfall2.list[[3]]
cpdbmeta<-as.matrix(Blastoids_data$celltype2)
cpdbmeta_df <- as.data.frame(cpdbmeta)
cpdbmeta_df$Cell <- rownames(cpdbmeta_df)
cpdbmeta_df<-cpdbmeta_df[,c(2,1)]
colnames(cpdbmeta_df)<-c("Cell","cell_type")
write.table(cpdbmeta_df,file=paste0(names(selfall2.list)[3],"_meta.txt"),sep = "\t",quote = F,row.names = F)
cpdbcounts<-round(Blastoids_data@assays$RNA@data,3)
cpdbcounts_df <- as.data.frame(cpdbcounts)
cpdbcounts_df$Gene <- rownames(cpdbcounts_df)
cpdbcounts_df<-cpdbcounts_df[,c("Gene",colnames(cpdbcounts))]
write.table(cpdbcounts_df,file = paste0(names(selfall2.list)[3],"_counts.txt"),sep = "\t",quote=F,row.names = F)
# shell run:   python HumanBlastoids_cellphonedb.py
setwd("~/data/wuhao/scdata/Rrun/AllFigures/Figure4 celltypecommunication")
mean   <- read.table("statistical_analysis_means_06_10_2025_174723_Blastoids.txt", sep = "\t",header = T)
pvalue <- read.table("statistical_analysis_pvalues_06_10_2025_174723_Blastoids.txt",sep = "\t")
colnames(pvalue)<-pvalue[1,]
pvalue<-pvalue[-1,]
colnames(pvalue)<-gsub("[ |()]",".",colnames(pvalue))
mean2<-mean[,c(1:13,which(colnames(mean)%in%c("Epiblast.Polar.TE","Epiblast.Mural.TE")))]
pl_FGF<-mean2$interacting_pair[grepl("Fibroblast growth factor",mean2$classification)]
pl_TGF<-mean2$interacting_pair[grepl("Transforming growth factor",mean2$classification)]
pl_WNT<-mean2$interacting_pair[grepl("Signaling by WNT",mean2$classification)]
pl_all<-c(pl_FGF,pl_TGF,pl_WNT)
mean3<-mean2[mean2$interacting_pair%in%pl_all,]
drawmatrix<-mean3[,c(2,13,14:ncol(mean3))]
pvalue2<-pvalue[,colnames(pvalue)%in%colnames(drawmatrix)]
pvalue2<-pvalue2[pvalue2[,1]%in%drawmatrix[,1],]
pvalue2<-pvalue2[,order(match(colnames(pvalue2),colnames(drawmatrix)))]
pvalue2<-pvalue2[order(match(pvalue2[,1],drawmatrix[,1])),]
pvalue2_melt<-melt(pvalue2,id.vars=c("interacting_pair","classification"))
drawmatrix_melt<-melt(drawmatrix,id.vars=c("interacting_pair","classification"))
if (identical(pvalue2_melt[,1],drawmatrix_melt[,1])&
    identical(pvalue2_melt[,3],drawmatrix_melt[,3])) {
  drawmatrix_melt<-cbind(drawmatrix_melt,as.numeric(pvalue2_melt[,4]))
}else{
  print("matrix not eactly match!")
}
colnames(drawmatrix_melt)<-c("gene_pair","classification","celltypes","score","pvalue")
drawmatrix_melt[,4]<-scale(drawmatrix_melt[,4])
drawmatrix_melt[,4]<-(drawmatrix_melt[,4] - min(drawmatrix_melt[,4])) / (max(drawmatrix_melt[,4]) - min(drawmatrix_melt[,4]))
drawmatrix_melt<-drawmatrix_melt[drawmatrix_melt[,5]<0.05,]
drawmatrix_melt<-drawmatrix_melt[!is.na(drawmatrix_melt[,4]),]
drawmatrix_melt<-drawmatrix_melt[drawmatrix_melt$classification!="",]
drawmatrix_melt<-drawmatrix_melt[drawmatrix_melt$classification!="Signaling by WNT inhibition",]
drawmatrix_melt<-drawmatrix_melt[!(grepl("^WNT9",drawmatrix_melt$gene_pair)|
                                     grepl("^WNT5",drawmatrix_melt$gene_pair)),]
drawmatrix_melt<-as.data.frame(drawmatrix_melt)
drawmatrix_melt$score<-as.numeric(drawmatrix_melt$score)
drawmatrix_melt$pvalue<-as.numeric(drawmatrix_melt$pvalue)
drawmatrix_melt$celltypes<-factor(drawmatrix_melt$celltypes,levels = c("Epiblast.Polar.TE","Epiblast.Mural.TE"))
# Figure 4A
p <- ggplot(drawmatrix_melt, aes(x=celltypes, y=fct_rev(gene_pair), size=-log10(pvalue+0.000001),colour=score )) + # fill=sore,size=pvalue
  geom_point()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1),strip.text.x = element_text(angle = 90, hjust = 1))+
  scale_colour_gradientn(colours = viridis_pal()(100)[51:100])
ggsave(p,filename = "HumanBlastoids_cellinteraction_selected2.pdf",width = 6,height = 8)
# Chimeric blastoids cellphonedb
imdata<-subset(Mall,subset=sample2=="CB_mouse")
cpdbcounts_m<-round(imdata@assays$RNA@data,3)
cpdbmeta_m  <-as.matrix(imdata$celltype2)
cpdbmeta_m <- as.data.frame(cpdbmeta_m)
cpdbmeta_m$Cell <- rownames(cpdbmeta_m)
cpdbmeta_m<-cpdbmeta_m[,c(2,1)]
colnames(cpdbmeta_m)<-c("Cell","cell_type")
cpdbmeta_m<-cpdbmeta_m[cpdbmeta_m[,2]!="TE",]
# extract CBh polar.TE and mural.TE
ihdata<-subset(CBall2,subset=sample2=="CB_human")
cpdbcounts_h<-round(ihdata@assays$RNA@data,3)
cpdbcounts_h<-cpdbcounts_h[grepl("-human$",rownames(cpdbcounts_h)),]
rownames(cpdbcounts_h)<-gsub("-human$","", rownames(cpdbcounts_h))
cpdbmeta_h  <-as.matrix(ihdata$celltype)
cpdbmeta_h <- as.data.frame(cpdbmeta_h)
cpdbmeta_h$Cell <- rownames(cpdbmeta_h)
cpdbmeta_h<-cpdbmeta_h[,c(2,1)]
colnames(cpdbmeta_h)<-c("Cell","cell_type")
cpdbmeta_h<-cpdbmeta_h[cpdbmeta_h[,2]!="EPI_h",]
cpdbmeta_h[cpdbmeta_h[,2]%in%c("TE1_h","TE2_h"),2]="TE_h"
cpdbmeta_h$cell_type[cpdbmeta_h$cell_type=="Mural_h"]<-"MuralTE_h"
# data merge
cpdbmeta_m[,2]<-paste0(cpdbmeta_m[,2],"_m")
cpdbmeta_all<-rbind(cpdbmeta_m,cpdbmeta_h)
Homogene <-intersect(rownames(cpdbcounts_h),rownames(cpdbcounts_m))
cpdbcounts_all<-cbind(cpdbcounts_h[Homogene,],cpdbcounts_m[Homogene,])
cpdbcounts_all<-cpdbcounts_all[,cpdbmeta_all$Cell]
cpdbcounts_all <- as.data.frame(cpdbcounts_all)
cpdbcounts_all$Gene <- rownames(cpdbcounts_all)
cpdbcounts_all<-cpdbcounts_all[,c("Gene",colnames(cpdbcounts_all)[1:(ncol(cpdbcounts_all)-1)])]
write.table(cpdbcounts_all,file = "ChimericB_counts.txt",sep = "\t",quote=F,row.names = F)
write.table(cpdbmeta_all, file="ChimericB_meta.txt", sep = "\t", quote = F, row.names = F)
# run Shell code: python ChimericBlastoids_cellphonedb.py
meanCB   <- read.table("statistical_analysis_means_07_17_2025_133857.txt", sep = "\t",header = T)
pvalueCB <- read.table("statistical_analysis_pvalues_07_17_2025_133857.txt",sep = "\t")
colnames(pvalueCB)<-pvalueCB[1,]
pvalueCB<-pvalueCB[-1,]
colnames(pvalueCB)<-gsub("[ |()]",".",colnames(pvalueCB))
meanCB2<-meanCB[,c(1:13,which(colnames(meanCB)%in%c("EPI.PolarTE_h","EPI.MuralTE_h")))]
drawmatrixCB<-meanCB2[meanCB2$interacting_pair%in%unique(drawmatrix_melt$gene_pair),c(2,13,14:ncol(meanCB2))]
pvalueCB2<-pvalueCB[,colnames(pvalueCB)%in%colnames(drawmatrixCB)]
pvalueCB2<-pvalueCB2[pvalueCB2[,1]%in%drawmatrixCB[,1],]
pvalueCB2<-pvalueCB2[,order(match(colnames(pvalueCB2),colnames(drawmatrixCB)))]
pvalueCB2<-pvalueCB2[order(match(pvalueCB2[,1],drawmatrixCB[,1])),]
pvalueCB2_melt<-melt(pvalueCB2,id.vars=c("interacting_pair","classification"))
drawmatrixCB_melt<-melt(drawmatrixCB,id.vars=c("interacting_pair","classification"))
if (identical(pvalueCB2_melt[,1],drawmatrixCB_melt[,1])&
    identical(pvalueCB2_melt[,3],drawmatrixCB_melt[,3])) {
  drawmatrixCB_melt<-cbind(drawmatrixCB_melt,as.numeric(pvalueCB2_melt[,4]))
}else{
  print("matrix not eactly match!")
}
colnames(drawmatrixCB_melt)<-c("gene_pair","classification","celltypes","score","pvalue")
drawmatrixCB_melt[,4]<-scale(drawmatrixCB_melt[,4])
drawmatrixCB_melt[,4]<-(drawmatrixCB_melt[,4] - min(drawmatrixCB_melt[,4])) / (max(drawmatrixCB_melt[,4]) - min(drawmatrixCB_melt[,4]))
drawmatrixCB_melt<-drawmatrixCB_melt[drawmatrixCB_melt[,5]<0.05,]
drawmatrixCB_melt<-drawmatrixCB_melt[!is.na(drawmatrixCB_melt[,4]),]
drawmatrixCB_melt<-drawmatrixCB_melt[drawmatrixCB_melt$classification!="",]
drawmatrixCB_melt<-as.data.frame(drawmatrixCB_melt)
drawmatrixCB_melt$score<-as.numeric(drawmatrixCB_melt$score)
drawmatrixCB_melt$pvalue<-as.numeric(drawmatrixCB_melt$pvalue)
drawmatrixCB_melt$classification<-factor(drawmatrixCB_melt$classification,levels = c("Signaling by WNT",
                                                                                 "Signaling by Transforming growth factor","Signaling by Fibroblast growth factor"))
merge_melt<-rbind(drawmatrix_melt,drawmatrixCB_melt)
merge_melt$score<-as.numeric(merge_melt$score)
merge_melt$pvalue<-as.numeric(merge_melt$pvalue)
merge_melt$classification<-factor(merge_melt$classification,levels = c("Signaling by WNT",
                                                                        "Signaling by Transforming growth factor","Signaling by Fibroblast growth factor"))
merge_melt$celltypes<-factor(merge_melt$celltypes,levels = c("Epiblast.Polar.TE","Epiblast.Mural.TE","EPI.PolarTE_h","EPI.MuralTE_h"))
# Figure 4A
p <- ggplot(merge_melt, aes(x=celltypes, y=fct_rev(gene_pair), size=-log10(pvalue+0.000001),colour=score )) + 
  geom_point()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1),strip.text.x = element_text(angle = 90, hjust = 1))+
  scale_colour_gradientn(colours = viridis_pal()(100)[51:100])
ggsave(p,filename = "Blastoids_cellinteraction_merged.pdf",width = 6,height = 8)
# TE cell integration
setwd("~/data/wuhao/scdata/Rrun/AllFigures/Fig4 celltypecommunication/AllTEintegration")
CBall2<-readRDS("../CBall2_backup_20250722.Rdata")
load("~/data/wuhao/scdata/Rrun/selfrun2/changehyporerun0206/Selfall2_2025.02.20.rdata")
# Chimeric TE
ChimericaTE<-subset(CBall2,subset=celltype%in%c("TE1_h","TE2_h","Mural_h","PolarTE_h","TE_fusion_h"))
ChimericaTE$sampletype=ChimericaTE$sample2
ChimericaTE<-subset(ChimericaTE,subset= percent.mt < 10)
DefaultAssay(ChimericaTE)<-"RNA"
ChimericaTEdata<-ChimericaTE@assays$RNA@counts
ChimericaTEdata<-ChimericaTEdata[grepl("-human$",rownames(ChimericaTEdata)),]
rownames(ChimericaTEdata)<-gsub("-human$","",rownames(ChimericaTEdata))
ChimericaTE2<-CreateSeuratObject(counts = ChimericaTEdata, min.features = 200,
                                 meta.data = ChimericaTE@meta.data,project = "CB_human")
# Human Blastoids TE
TE1<-subset(selfall2,subset=celltype%in%c("TE","TE.fusion_competent","Mural.TE","Polar.TE"))
DefaultAssay(TE1)<-"RNA"
hTE<-subset(TE1,subset=sampletype=="Blastoid")
# TE Spheroid 
TEspheroid<-subset(TE1,subset=sampletype=="TE_spheroid")
# Merge all TE
TEall<-merge(x=hTE,y=c(ChimericaTE2,TEspheroid))
TEall@active.ident<-as.factor(TEall$sampletype)
commongenes<-intersect(rownames(hTE),rownames(ChimericaTE2))
commongenes<-intersect(commongenes,rownames(TEspheroid))
TEall<-subset(TEall, features = commongenes)
TEall.list <- SplitObject(TEall, split.by = "sampletype")
for (i in 1:length(x = TEall.list)) {
  TEall.list[[i]] <- NormalizeData(object = TEall.list[[i]], verbose = FALSE)
  TEall.list[[i]] <- FindVariableFeatures(object = TEall.list[[i]],
                                          selection.method = "vst",  nfeatures = 2000, verbose = FALSE)
  print(i)
}
TEall.anchors <- FindIntegrationAnchors(object.list = TEall.list, dims = 1:20,k.anchor = 5,k.filter = 10)
TEall <- IntegrateData(anchorset = TEall.anchors, dims = 1:40,k.weight = 20)
TEall <- ScaleData(TEall)
TEall <- RunPCA(object = TEall,npcs=200)
pdf(file = "TEall_pca_use.pdf")
ElbowPlot(TEall, ndims = 100)
dev.off()
TEall$celltype<-plyr::mapvalues(x=TEall$celltype,
                  from=c("Mural_h","PolarTE_h","TE_fusion_h","TE1_h","TE2_h"),
                  to = c("Mural.TE","Polar.TE","TE.fusion_competent","TE","TE"))
for (i in 1:10) {
  pcause=i*10
  TEall <- FindNeighbors(object = TEall,dims = 1:pcause)
  TEall <- FindClusters(object = TEall, resolution = 0.6)
  TEall <- RunUMAP(object = TEall, dims = 1:pcause)
  pdf(file = paste("./umap/TEall_rmintermediate_cluster",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = TEall,reduction = "umap",label = T,pt.size = 0.3))
  dev.off()
  pdf(file = paste("./umap/TEall_rmintermediate_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = TEall,reduction = "umap",group.by = "sampletype",pt.size = .3))
  dev.off()
  pdf(file = paste("./umap/TEall_rmintermediate_celltypeold",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = TEall,reduction = "umap",group.by = "celltype",pt.size = .3))
  dev.off()
}
pcause=40
TEall <- FindNeighbors(object = TEall,dims = 1:pcause)
TEall <- FindClusters(object = TEall, resolution = 0.3)
TEall <- RunUMAP(object = TEall, dims = 1:pcause)
pdf(file = paste("TEall_rmintermediate_cluster",pcause,".pdf",sep = ""),width = 6,height = 5)
print(DimPlot(object = TEall,reduction = "umap",label = T,pt.size = 0.3))
dev.off()
pdf(file = paste("TEall_rmintermediate_sample",pcause,".pdf",sep = ""),width = 7,height = 5)
print(DimPlot(object = TEall,reduction = "umap",group.by = "sampletype",pt.size = .3))
dev.off()
pdf(file = paste("Eall_rmintermediate_celltypeold",pcause,".pdf",sep = ""),width = 6,height = 5)
print(DimPlot(object = TEall,reduction = "umap",group.by = "celltype",pt.size = .3))
dev.off()
TEall$seuratcluster<-TEall@active.ident
TEall$samplecelltype<-paste0(TEall$sampletype,TEall$celltype)
TEall@active.ident<-as.factor(TEall$sampletype)
TEall.ctlist<-SplitObject(TEall, split.by = "celltype")
hcg<-read.table("~/reference/genesets/human_protein-coding_gene_list.txt",header = T)
TopMatrix <- TEall@assays$RNA@data[,]
celltype_specieslist <- unique(TEall$samplecelltype)
for (i in 1:length(celltype_specieslist)) {
  icell<-TopMatrix[,TEall$samplecelltype==celltype_specieslist[i]]
  if (i==1) {
    MeanMatrix<-rowMeans(icell)
  }else{
    MeanMatrix<-cbind(MeanMatrix,rowMeans(icell))
  }
}
colnames(MeanMatrix)<-celltype_specieslist
muralTEDEGs<-FindAllMarkers(TEall.ctlist$Mural.TE,only.pos=T)
muralTEDEGs<-muralTEDEGs[muralTEDEGs$p_val_adj<0.05&abs(muralTEDEGs$avg_log2FC)>0.5,]
muralTEDEGs$state<-"Down"
muralTEDEGs$state[muralTEDEGs$avg_log2FC>0]<-"Up"
muralTEDEGs<-muralTEDEGs[rownames(muralTEDEGs)%in%hcg$symbol,]
write.csv(muralTEDEGs,file = "MuralTEDEGs.csv")
# Figure 4B TE DEGs heatmap
pdf(file = "muralTEup_eachcell_heatmap.pdf",width = 10,height = 10)
DoHeatmap(object = TEall.ctlist$Mural.TE,slot = "scale.data", features = muralTEDEGs$gene,group.by = "sampletype",label = TRUE,group.bar= T) +
  scale_fill_gradientn(colors = c("blue","white","red")) 
dev.off()
MMM<-MeanMatrix[muralTEDEGs$gene,unique(TEall.ctlist$Mural.TE$samplecelltype)]
annotation_col<-data.frame(Groups=colnames(MMM))
rownames(annotation_col)<-colnames(MMM)
MMM<-MMM[,c(2,3,1)]
pdf(file = "MuralTEMarkermean_heatmap.pdf",width = 6,height = 10)
pheatmap(mat = MMM,scale = "row",cluster_cols = F,cluster_rows = T,
         annotation_col = annotation_col,fontsize_row = 3,
         color = colorRampPalette(c("blue","white","red"))(100),border_color = "white")
dev.off()
# Figure 4 E function annotation
muralBlastoid_GOup<-callGO(muralTEDEGs$gene[muralTEDEGs$cluster=="Blastoid"&muralTEDEGs$state=="Up"],"human")
drawGOKEGG(muralBlastoid_GOup,icolor="Up")
muralBlastoid_KEGGup<-callKEGG(muralTEDEGs$gene[muralTEDEGs$cluster=="Blastoid"&muralTEDEGs$state=="Up"],"human")
drawGOKEGG(muralBlastoid_KEGGup,icolor="Up")
muralCB_GOup<-callGO(muralTEDEGs$gene[muralTEDEGs$cluster=="CB_human"&muralTEDEGs$state=="Up"],"human")
drawGOKEGG(muralCB_GOup,icolor="Up")
muralCB_KEGGup<-callKEGG(muralTEDEGs$gene[muralTEDEGs$cluster=="CB_human"&muralTEDEGs$state=="Up"],"human")
drawGOKEGG(muralCB_KEGGup,icolor="Up")
muralTES_GOup<-callGO(muralTEDEGs$gene[muralTEDEGs$cluster=="TE_spheroid"&muralTEDEGs$state=="Up"],"human")
drawGOKEGG(muralTES_GOup,icolor="Up")
muralTES_KEGGup<-callKEGG(muralTEDEGs$gene[muralTEDEGs$cluster=="TE_spheroid"&muralTEDEGs$state=="Up"],"human")
drawGOKEGG(muralTES_KEGGup,icolor="Up")
# PolarTE
polarTEDEGs<-FindAllMarkers(TEall.ctlist$Polar.TE,only.pos=T)
polarTEDEGs<-polarTEDEGs[polarTEDEGs$p_val_adj<0.05&abs(polarTEDEGs$avg_log2FC)>0.5,]
polarTEDEGs$state<-"Down"
polarTEDEGs$state[polarTEDEGs$avg_log2FC>0]<-"Up"
polarTEDEGs<-polarTEDEGs[rownames(polarTEDEGs)%in%hcg$symbol,]
write.csv(polarTEDEGs,file = "polarTEDEGs.csv")
MMM<-MeanMatrix[polarTEDEGs$gene,unique(TEall.ctlist$Polar.TE$samplecelltype)]
annotation_col<-data.frame(Groups=colnames(MMM))
rownames(annotation_col)<-colnames(MMM)
MMM<-MMM[rowSums(MMM)>0,]
pdf(file = "PolarTEMarkermean_heatmap.pdf",width = 6,height = 10)
pheatmap(mat = MMM,scale = "row",cluster_cols = T,cluster_rows = T,
         annotation_col = annotation_col,fontsize_row = 3,
         color = colorRampPalette(c("blue","white","red"))(100),border_color = "white")
dev.off()
polarBlastoid_GOup<-callGO(polarTEDEGs$gene[polarTEDEGs$cluster=="Blastoid"&polarTEDEGs$state=="Up"],"human")
drawGOKEGG(polarBlastoid_GOup,icolor="Up")
polarBlastoid_KEGGup<-callKEGG(polarTEDEGs$gene[polarTEDEGs$cluster=="Blastoid"&polarTEDEGs$state=="Up"],"human")
drawGOKEGG(polarBlastoid_KEGGup,icolor="Up")
polarCB_GOup<-callGO(polarTEDEGs$gene[polarTEDEGs$cluster=="CB_human"&polarTEDEGs$state=="Up"],"human")
drawGOKEGG(polarCB_GOup,icolor="Up")
polarCB_KEGGup<-callKEGG(polarTEDEGs$gene[polarTEDEGs$cluster=="CB_human"&polarTEDEGs$state=="Up"],"human")
drawGOKEGG(polarCB_KEGGup,icolor="Up")
polarTES_GOup<-callGO(polarTEDEGs$gene[polarTEDEGs$cluster=="TE_spheroid"&polarTEDEGs$state=="Up"],"human")
drawGOKEGG(polarTES_GOup,icolor="Up")
polarTES_KEGGup<-callKEGG(polarTEDEGs$gene[polarTEDEGs$cluster=="TE_spheroid"&polarTEDEGs$state=="Up"],"human")
drawGOKEGG(polarTES_KEGGup,icolor="Up")
# Figure S4A
meanCB   <- read.table("statistical_analysis_means_07_23_2025_124650.txt", sep = "\t",header = T)
pvalueCB <- read.table("statistical_analysis_pvalues_07_23_2025_124650.txt",sep = "\t")

colnames(pvalueCB)<-pvalueCB[1,]
pvalueCB<-pvalueCB[-1,]
colnames(pvalueCB)<-gsub("[ |()]",".",colnames(pvalueCB))
meanCB2<-meanCB[,c(1:13,which(colnames(meanCB)%in%c("MuralTE_h.EPI_m","MuralTE_h.PrE_m",
                                                    "PolarTE_h.EPI_m","PolarTE_h.PrE_m",
                                                    "TE_h.EPI_m","TE_h.PrE_m",
                                                    "TE_fusion_h.EPI_m","TE_fusion_h.PrE_m"
                                                    )))]
pl_FGF<-meanCB2$interacting_pair[grepl("Fibroblast growth factor",meanCB2$classification)]
pl_TGF<-meanCB2$interacting_pair[grepl("Transforming growth factor",meanCB2$classification)]
pl_WNT<-meanCB2$interacting_pair[grepl("Signaling by WNT",meanCB2$classification)]
pl_all<-c(pl_FGF,pl_TGF,pl_WNT)
meanCB3<-meanCB2[meanCB2$interacting_pair%in%pl_all,]
drawmatrixCB<-meanCB2[,c(2,13,14:ncol(meanCB3))]
pvalue2<-pvalueCB[,colnames(pvalueCB)%in%colnames(drawmatrixCB)]
pvalue2<-pvalue2[pvalue2[,1]%in%drawmatrixCB[,1],]
pvalue2<-pvalue2[,order(match(colnames(pvalue2),colnames(drawmatrixCB)))]
pvalue2<-pvalue2[order(match(pvalue2[,1],drawmatrixCB[,1])),]
pvalue2_melt<-melt(pvalue2,id.vars=c("interacting_pair","classification"))
drawmatrixCB_melt<-melt(drawmatrixCB,id.vars=c("interacting_pair","classification"))
if (identical(pvalue2_melt[,1],drawmatrixCB_melt[,1])&
    identical(pvalue2_melt[,3],drawmatrixCB_melt[,3])) {
  drawmatrixCB_melt<-cbind(drawmatrixCB_melt,as.numeric(pvalue2_melt[,4]))
}else{
  print("matrix not eactly match!")
}
colnames(drawmatrixCB_melt)<-c("gene_pair","classification","celltypes","score","pvalue")
drawmatrixCB_melt[,4]<-scale(drawmatrixCB_melt[,4])
drawmatrixCB_melt[,4]<-(drawmatrixCB_melt[,4] - min(drawmatrixCB_melt[,4])) / (max(drawmatrixCB_melt[,4]) - min(drawmatrixCB_melt[,4]))
drawmatrixCB_melt<-drawmatrixCB_melt[drawmatrixCB_melt[,5]<0.05,]
drawmatrixCB_melt<-drawmatrixCB_melt[!is.na(drawmatrixCB_melt[,4]),]
drawmatrixCB_melt<-drawmatrixCB_melt[drawmatrixCB_melt$classification!="",]
drawmatrixCB_melt<-as.data.frame(drawmatrixCB_melt)
drawmatrixCB_melt$score<-as.numeric(drawmatrixCB_melt$score)
drawmatrixCB_melt$pvalue<-as.numeric(drawmatrixCB_melt$pvalue)
MEcom<-intersect(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Mural.TE.Epiblast"]
  ,drawmatrixCB_melt$gene_pair[drawmatrixCB_melt$celltypes=="MuralTE_h.EPI_m"])
MHcom<-intersect(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Mural.TE.Hypoblast"]
                 ,drawmatrixCB_melt$gene_pair[drawmatrixCB_melt$celltypes=="MuralTE_h.PrE_m"])
PEcom<-intersect(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Polar.TE.Epiblast"]
                 ,drawmatrixCB_melt$gene_pair[drawmatrixCB_melt$celltypes=="PolarTE_h.EPI_m"])
PHcom<-intersect(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Polar.TE.Hypoblast"]
                 ,drawmatrixCB_melt$gene_pair[drawmatrixCB_melt$celltypes=="PolarTE_h.PrE_m"])
FEcom<-intersect(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="TE.fusion_competent.Epiblast"]
                 ,drawmatrixCB_melt$gene_pair[drawmatrixCB_melt$celltypes=="TE_fusion_h.EPI_m"])
FHcom<-intersect(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="TE.fusion_competent.Hypoblast"]
                 ,drawmatrixCB_melt$gene_pair[drawmatrixCB_melt$celltypes=="TE_fusion_h.PrE_m"])
comall<-c(MEcom,MHcom,PEcom,PHcom,FEcom,FHcom)
merge_melt<-rbind(drawmatrix_melt,drawmatrixCB_melt)
merge_melt$score<-as.numeric(merge_melt$score)
merge_melt$pvalue<-as.numeric(merge_melt$pvalue)
merge_melt<-merge_melt[merge_melt$gene_pair%in%comall,]
merge_melt$gene_pair<-factor(merge_melt$gene_pair,levels=unique(comall))
merge_melt$celltypes<-factor(merge_melt$celltypes,levels=c("Mural.TE.Epiblast","MuralTE_h.EPI_m",
"Polar.TE.Epiblast","PolarTE_h.EPI_m",
"Polar.TE.Hypoblast","PolarTE_h.PrE_m",
"TE.fusion_competent.Epiblast","TE_fusion_h.EPI_m",
"TE.fusion_competent.Hypoblast","TE_fusion_h.PrE_m"))
merge_melt<-merge_melt[!is.na(merge_melt$celltypes),]
merge_melt<-merge_melt[-2,]
merge_melt$classification<-as.factor(merge_melt$classification)
merge_melt <- merge_melt %>%
  arrange(classification, gene_pair) %>%
  mutate(gene_pair = fct_inorder(gene_pair))
write.csv(merge_melt,file = "HumanBlastoids and Chimeric Blastoids TE cellcommunication.csv")
annotation_df <- merge_melt %>%
  distinct(gene_pair, classification)
max_score <- max(merge_melt$score, na.rm = TRUE)
p <- ggplot(merge_melt, aes(x = celltypes, y = gene_pair)) + 
  geom_point(aes(size = -log10(pvalue + 1e-6), colour = score)) +
  geom_text(data = annotation_df, 
            aes(x = Inf, y = gene_pair, label = classification), 
            hjust = 0, size = 3, inherit.aes = FALSE) +
  scale_colour_gradientn(colours = viridis::viridis_pal()(100)[51:100]) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.margin = margin(5, 60, 5, 5),
    axis.text.y = element_text(size = 9)
  ) +
  coord_cartesian(clip = "off") 
ggsave(p,filename = "Common_Blastoids_cellinteraction_mergedAll2.pdf",width = 8,height = 10)
# ==============================================================================
# ================================ Figure 6 ====================================
# ==============================================================================
setwd("~/data/wuhao/scdata/Rrun/AllFigures/Figure6")
# load selfall2
load("Selfall2.rdata")
BD0 <-subset(x=selfall2,subset=sampletype=="Blastoid")
DefaultAssay(BD0)<-"RNA"
BD0@active.ident<-as.factor(BD0$sampletype)
BD0 <- subset(x=BD0, subset=celltype%in%c("EPI1","Mural.TE","Polar.TE","EPI2","TE.fusion_competent",
                                          "TE","Hypoblast"))
BD0 <- subset(x = BD0, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 10 & nCount_RNA < 150000)
BD0@assays$integrated=NULL
# load blastoids day1
library(Matrix)
mtx      <- readMM("./MergeBlastoidsday1/data/count_matrix.mtx")
genes    <- read.csv("./MergeBlastoidsday1/data/all_genes.csv")
barcodes <- read.csv("./MergeBlastoidsday1/data/cell_metadata.csv")
rownames(mtx)<-barcodes$bc_wells
colnames(mtx)<-genes$gene_name
mtx<-t(mtx)
mtx <- as(mtx, "dgCMatrix")
mtx <- mtx[!duplicated(rownames(mtx)), ]
rownames(barcodes)<-barcodes$bc_wells
BD13  <- CreateSeuratObject(counts=mtx,min.cells = 3,meta.data = barcodes,min.features=200,project="Blastoids_Day1")
BD13  <- subset(x=BD13, subset=sample%in%c("human_blastoid_implantation_day1",
                                           "human_blastoid_implantation_day3"))
BD13@meta.data$tech<-"dropseq"
BD13@meta.data$sampletype<-BD13@meta.data$sample
BD13[["percent.mt"]] <- PercentageFeatureSet(BD13, pattern = "^MT-")
BD13@active.ident<-as.factor(BD13$sampletype)
plotQC(BD13,"raw")
BD13 <- subset(x = BD13, subset = nFeature_RNA > 500 & nFeature_RNA < 15000 & percent.mt < 15 & nCount_RNA < 150000)
plotQC(BD13,"cut")
BD13<-RMdoublets(BD13)
BDall<-merge(x=BD0,y=BD13)
commongenes<-intersect(rownames(BD0),rownames(BD13))
BDall <- subset(BDall, features = commongenes)
plotQC(BDall,"cut")

BDall.list <- SplitObject(BDall, split.by = "sampletype")
for (i in 1:length(x = BDall.list)) {
  BDall.list[[i]] <- NormalizeData(object = BDall.list[[i]], verbose = FALSE)
  BDall.list[[i]] <- FindVariableFeatures(object = BDall.list[[i]],
                                          selection.method = "vst",  nfeatures = 3000, verbose = FALSE)
  print(i)
}
BDall.anchors <- FindIntegrationAnchors(object.list = BDall.list, dims = 1:60,k.anchor = 10,k.filter = 20)
BDall <- IntegrateData(anchorset = BDall.anchors, dims = 1:200)
BDall <- ScaleData(BDall)
BDall <- RunPCA(object = BDall,npcs=200)
pdf(file = "BDall_pca_use.pdf")
ElbowPlot(BDall, ndims = 100)
dev.off()
pcause=90
BDall <- FindNeighbors(object = BDall, dims = 1:pcause)
BDall <- FindClusters(object = BDall, resolution = 1.2)
BDall <- RunUMAP(object = BDall, dims = 1:pcause)
pdf(file = paste("BDall_rmintermediate_umap_cluster",pcause,"_1.2.pdf",sep = ""),width = 10,height = 8)
DimPlot(object = BDall,reduction = "umap",label = T,pt.size = 0.3)
dev.off()
pdf(file = paste("BDall_rmintermediate_umap_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
DimPlot(object = BDall,reduction = "umap",group.by = "sampletype",pt.size = .3)
dev.off()
pdf(file = paste("BDall_rmintermediate_umap_celltype",pcause,".pdf",sep = ""),width = 10,height = 8)
DimPlot(object = BDall,reduction = "umap",group.by = "celltype",label = T,pt.size = .3)
dev.off()
BDall@meta.data$mergecluster<-BDall@active.ident
c_20<-intersect(names(BDall@active.ident[BDall@active.ident==11]),
                rownames(BDall@reductions$umap@cell.embeddings)[BDall@reductions$umap@cell.embeddings[,1]>-7.5])
BDall@meta.data$mergecluster<-as.numeric(BDall@meta.data$mergecluster)
BDall$mergecluster[c_20]<-20
c_21<-intersect(names(BDall@active.ident[BDall@active.ident==15]),
                rownames(BDall@reductions$umap@cell.embeddings)[BDall@reductions$umap@cell.embeddings[,1]>7])
BDall@meta.data$mergecluster<-as.numeric(BDall@meta.data$mergecluster)
BDall$mergecluster[c_21]<-21
pdf(file = paste("BDall_rmintermediate_umap_cluster",pcause,"_1.2_changed.pdf",sep = ""),width = 10,height = 8)
DimPlot(object = BDall,reduction = "umap",group.by="mergecluster",label = T,pt.size = 0.3)
dev.off()
# cell type annotation
markerall<-read.table("mergeblastoids013_markers.txt",header = T,sep = "\t")
for (i in 1:ncol(markerall)) {
  imarker<-cbind(markerall[,i],rep(colnames(markerall)[i],length(markerall[,i])))
  if (i==1) {
    markers<-imarker
  }else{
    markers<-rbind(markers,imarker)
  }
}
markers<-na.omit(markers)
markers<-as.matrix(markers)
markers<-unique(markers)
markers<-markers[markers[,1]!="",]
callmarkerdot(markers,BDall,group.by="mergecluster",fname="BDall_cluster")
BDall$celltype2<-plyr::mapvalues(BDall$mergecluster,
                                 from = c(1:21),
                                 to   = c("Mural.TE","Mural.TE","Polar.TE","EPI1","EPI2",
                                 "Mural.TE","EEC1","EPI2","STB","EPI1","Polar.TE","EPI2",
                                 "TE_fusion","EEC2","TE","Polar.TE","EPI1","Hypoblast","EPI2","EPI2","Mural.TE"))
# Figure6A 
pdf(file = paste("BDall_rmintermediate_umap_celltype2_",pcause,"_1.2_changed_test.pdf",sep = ""),width = 10,height = 8)
DimPlot(object = BDall,reduction = "umap",group.by="celltype2",label = T,pt.size = 0.3)
dev.off()
# save cell type information
BDall_b2c<-cbind(names(BDall$celltype2),as.matrix(BDall$celltype2),as.matrix(BDall$sampletype))
BDall_b2c[BDall_b2c[,3]=="Blastoid",1] <- sub(".{2}$", "", BDall_b2c[BDall_b2c[,3]=="Blastoid",1])
colnames(BDall_b2c)<-c("Original_Barcodes","Predicted_Celltype","Sample")
write.table(BDall_b2c,file = "Figure 6A Barcodes to celltype.txt",sep = "\t",row.names = T,col.names = T)


# Figure 6B  MarkerMeanHeatmap
setwd("~/data/wuhao/scdata/Rrun/AllFigures/Fig6 mergeBlastoidsDay013/Figure6add")
m6b<-read.table("Figure6b.txt",header = F,sep = "\t")
TopMatrix<-BDall@assays$RNA@data[,]
celltype_specieslist<-unique(BDall$celltype2)
celltype_specieslist<-setdiff(celltype_specieslist,c("EEC1","EEC2"))
for (i in 1:length(celltype_specieslist)) {
  icell<-TopMatrix[,BDall$celltype2==celltype_specieslist[i]]
  if (i==1) {
    MeanMatrix<-rowMeans(icell)
  }else{
    MeanMatrix<-cbind(MeanMatrix,rowMeans(icell))
  }
}
colnames(MeanMatrix)<-celltype_specieslist
MarkerMeanMatrix<-MeanMatrix[rownames(MeanMatrix)%in%m6b_heatmapgenelist,]
colnames(MarkerMeanMatrix)<-c("Mural TE","EPI 1","EPI 2","Polar TE","Hypoblast","TE","STB","TE fusion-competent")
MarkerMeanMatrix<-MarkerMeanMatrix[m6b_heatmapgenelist,c(2,3,5,6,1,4,8,7)]

pdf(file = "Figure6Bheatmap.pdf",width = 10,height = 4)
pheatmap(mat = t(MarkerMeanMatrix),scale = "column",cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(100),border_color = "white")
dev.off()
# Figure S6B
pdf(file = "FigureS6B_BDall_sampletype.pdf",width = 11,height = 8)
DimPlot(object = BDall,reduction = "umap",group.by="sampletype",label = F,pt.size = 0.3)
dev.off()
# Figure S6C and Figure S6D
for (i in 1:nrow(m6b)) {
  pdf(file = paste("./markergenestaining/FeaturePlot",m6b[i,1],m6b[i,2],".pdf",sep = "_"),width =5.5 ,height = 5)
  print(FeaturePlot(BDall,slot = "data", features = m6b[i,1],cols = viridis_pal()(100)[21:100],pt.size = 0.5,
                    alpha=0.8,order = T))
  dev.off()
}
CGBgenes<-paste0("CGB",c(1:8))
CGBgenes<-intersect(rownames(BDall@assays$RNA),CGBgenes)
for (i in 1:length(CGBgenes)) {
  pdf(file = paste("./markergenestaining/FeaturePlot",CGBgenes[i],".pdf",sep = "_"),width =5.5 ,height = 5)
  print(FeaturePlot(BDall,slot = "data", features = CGBgenes[i],cols = viridis_pal()(100)[21:100],pt.size = 0.5,
                    alpha=0.8,order = T))
  dev.off()
}

# Run NM model no uteral cells
setwd("./Figure6/rmUterine_NM")
para.list <- list()
para.list$cellGroup <- TRUE ## whether providing group information for cells
para.list$runMiloR <- TRUE  ## whether run miloR aggregation 
para.list$cor.cutoff <- 0.3 ## correlation cutoff
ict<-setdiff(unique(BDall$celltype2),c("EEC1","EEC2"))
BDall_rm<-subset(BDall,subset=celltype2%in%ict)
BDall_rm.counts<-BDall_rm@assays$RNA@counts
BDall_rm.counts.meta <-data.frame(BDall_rm@meta.data)
BDall_rm.counts.meta$cell<-rownames(BDall_rm.counts.meta)
BDall_rm.counts.meta$EML <-"query"
BDall_rm.counts.meta$pj  <-"Blastoids_merge"
BDall_rm.counts.meta$group<-BDall_rm.counts.meta$sampletype
BDall_rm.counts.meta$group <- plyr::mapvalues(BDall_rm.counts.meta$group,
                                           from = c("human_blastoid_implantation_day1",
                                                    "human_blastoid_implantation_day3"),
                                           to   = c("Blastoids_iday1","Blastoids_iday3"))
rownames(BDall_rm.counts.meta)<-NULL
BDpredict_out_rm<-runNMModel(BDall_rm.counts,BDall_rm.counts.meta,filename="BDmerge13_rmUterine")
# Figure 6C and Figure S6E
drawNMModel(BDpredict_out_rm,filename="Blastoid",group.by="devTime")
drawNMModel(BDpredict_out_rm,filename="Blastoids_iday1",group.by="devTime")
drawNMModel(BDpredict_out_rm,filename="Blastoids_iday3",group.by="devTime")
drawNMModel(BDpredict_out_rm,filename="Blastoids_merge",group.by="pj")
ctcolor<-c("Zygote"="#999999","2–4 cell"="#F39B7E","8 cell"="#F8766C","AdvMes"="#E9842C",
           "Amnion"="#E71F18","Axial_Mes"="#084334","CTB"="#9CA700","DE"="#A3C8DD",
           "Epiblast"="#00B813","Erythroblasts"="#393A79","EVT"="#53B885",
           "ExE_Mes"="#0085ED","HEP"="#00EBEC","Hypoblast"="#7B95FF","ICM"="#BB81FF",
           "Morula"="#FF7E0E","Mesoderm"="#00C1A7","Prelineage"="#C49C93","PriS"="#F862DF",
           "STB"="#F7B6D2","TE"="#00B5ED","YSE"="#EA618E","Ambiguous"="#666666","low_cor"="#666666")
pdf(file = "BDall_RMuteral_nmmodel3_celltype.pdf",width = 4.8,height = 4)
print(ggplot()+geom_point(BDpredict_out_rm$umap2 %>% filter(!pj%in%"Blastoids_merge"),
                          mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",size=0.1,alpha=1)+
        geom_point(BDpredict_out_rm$umap2 %>% filter(pj%in%"Blastoids_merge"),mapping=aes(x=UMAP_1,y=UMAP_2,color=pred_EML),
                   size=0.2,alpha=0.8)+theme_bw()+scale_color_manual(values = ctcolor))
dev.off()

pdf(file = "BDall_RMuteral_nmmodel3_sample.pdf",width = 5,height = 4)
print(ggplot()+geom_point(BDpredict_out_rm$umap2 %>% filter(!pj%in%"Blastoids_merge"),
                          mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",size=0.1,alpha=1)+
        geom_point(BDpredict_out_rm$umap2 %>% filter(pj%in%"Blastoids_merge"),mapping=aes(x=UMAP_1,y=UMAP_2,color=devTime),
                   size=0.2,alpha=0.8)+theme_bw())
dev.off()

# Figure 6D cellphonedb: EEC1 ↔ mural TE，EEC2 ↔ polar TE
# on all pathway coincident with paper's figures
setwd("../Figure6/cellcommunication")
idata<-subset(BDall,subset=sampletype%in%c("human_blastoid_implantation_day1",
                                           "human_blastoid_implantation_day3"))
idata$sampletype <- plyr::mapvalues(idata$sampletype,
                                           from = c("human_blastoid_implantation_day1",
                                                    "human_blastoid_implantation_day3"),
                                           to   = c("Iday1","Iday3"))
idata<-subset(idata,subset=celltype2%in%c("EEC1","EEC2","Mural.TE","Polar.TE"))
idata$celltype3<-paste0(idata$celltype2,"_",idata$sampletype)

cpdbmeta<-as.matrix(idata$celltype3)
cpdbmeta_df <- as.data.frame(cpdbmeta)
cpdbmeta_df$Cell <- rownames(cpdbmeta_df)
cpdbmeta_df<-cpdbmeta_df[,c(2,1)]
colnames(cpdbmeta_df)<-c("Cell","cell_type")
write.table(cpdbmeta_df,file="ChimericB_meta.txt",sep = "\t",quote = F,row.names = F)
cpdbcounts<-round(idata@assays$RNA@data,3)
cpdbcounts_df <- as.data.frame(cpdbcounts)
cpdbcounts_df$Gene <- rownames(cpdbcounts_df)
cpdbcounts_df<-cpdbcounts_df[,c("Gene",colnames(cpdbcounts))]
write.table(cpdbcounts_df,file = "ChimericB_counts.txt",sep = "\t",quote=F,row.names = F)
# shell code: python ChimericBlastoids_cellphonedb.py
# Figure 6D
mean   <- read.table("../cellcommunication/statistical_analysis_means_07_25_2025_032835.txt", sep = "\t",header = T)
pvalue <- read.table("../cellcommunication/statistical_analysis_pvalues_07_25_2025_032835.txt",sep = "\t")
colnames(pvalue)<-pvalue[1,]
pvalue<-pvalue[-1,]
colnames(pvalue)<-gsub("[ |()]",".",colnames(pvalue))
mean2<-mean[,c(1:13,which(colnames(mean)%in%c("EEC1_Iday1.Mural.TE_Iday1","EEC1_Iday1.Polar.TE_Iday1","EEC2_Iday1.Mural.TE_Iday1","EEC2_Iday1.Polar.TE_Iday1",
                                              "Mural.TE_Iday1.EEC1_Iday1","Polar.TE_Iday1.EEC1_Iday1","Mural.TE_Iday1.EEC2_Iday1","Polar.TE_Iday1.EEC2_Iday1"
)))]
drawmatrix<-mean2[,c(2,13,14:ncol(mean2))]
pvalue2<-pvalue[,colnames(pvalue)%in%colnames(drawmatrix)]
pvalue2<-pvalue2[pvalue2[,1]%in%drawmatrix[,1],]
pvalue2<-pvalue2[,order(match(colnames(pvalue2),colnames(drawmatrix)))]
pvalue2<-pvalue2[order(match(pvalue2[,1],drawmatrix[,1])),]
pvalue2_melt<-melt(pvalue2,id.vars=c("interacting_pair","classification"))
drawmatrix_melt<-melt(drawmatrix,id.vars=c("interacting_pair","classification"))
if (identical(pvalue2_melt[,1],drawmatrix_melt[,1])&
    identical(pvalue2_melt[,3],drawmatrix_melt[,3])) {
  drawmatrix_melt<-cbind(drawmatrix_melt,as.numeric(pvalue2_melt[,4]))
}else{
  print("matrix not eactly match!")
}

colnames(drawmatrix_melt)<-c("gene_pair","classification","celltypes","score","pvalue")
drawmatrix_melt[,4]<-scale(drawmatrix_melt[,4])
drawmatrix_melt[,4]<-(drawmatrix_melt[,4] - min(drawmatrix_melt[,4])) / (max(drawmatrix_melt[,4]) - min(drawmatrix_melt[,4]))
drawmatrix_melt<-drawmatrix_melt[drawmatrix_melt[,5]<0.05,]
drawmatrix_melt<-drawmatrix_melt[!is.na(drawmatrix_melt[,4]),]
drawmatrix_melt<-drawmatrix_melt[drawmatrix_melt$classification!="",]
EEC1_gp<-setdiff(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="EEC1_Iday1.Polar.TE_Iday1"],
                 drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="EEC1_Iday1.Mural.TE_Iday1"]
)
EEC2_gp<-setdiff(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="EEC2_Iday1.Polar.TE_Iday1"],
                 drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="EEC2_Iday1.Mural.TE_Iday1"]
)
EEC1_gp2<-setdiff(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Polar.TE_Iday1.EEC1_Iday1"],
                  drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Mural.TE_Iday1.EEC1_Iday1"]
)
EEC2_gp2<-setdiff(drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Polar.TE_Iday1.EEC2_Iday1"],
                  drawmatrix_melt$gene_pair[drawmatrix_melt$celltypes=="Mural.TE_Iday1.EEC2_Iday1"]
)
comall<-c(intersect(EEC1_gp,EEC2_gp),intersect(EEC1_gp2,EEC2_gp2))

drawmatrix_melt<-drawmatrix_melt[drawmatrix_melt$gene_pair%in%comall,]
drawmatrix_melt<-as.data.frame(drawmatrix_melt)
drawmatrix_melt$score<-as.numeric(drawmatrix_melt$score)
drawmatrix_melt$pvalue<-as.numeric(drawmatrix_melt$pvalue)
drawmatrix_melt<-rbind(drawmatrix_melt,c("BMP7_BMR1A_ACR2A","Signaling by BMP","EEC2_Iday1.Mural.TE_Iday1",0,1))
drawmatrix_melt<-rbind(drawmatrix_melt,c("BMP7_BMR1A_ACR2A","Signaling by BMP","Mural.TE_Iday1.EEC2_Iday1",0,1))
drawmatrix_melt<-rbind(drawmatrix_melt,c("BMP7_BMR1A_ACR2A","Signaling by BMP","Mural.TE_Iday1.EEC1_Iday1",0,1))
drawmatrix_melt<-as.data.frame(drawmatrix_melt)
drawmatrix_melt$score<-as.numeric(drawmatrix_melt$score)
drawmatrix_melt$pvalue<-as.numeric(drawmatrix_melt$pvalue)
drawmatrix_melt$celltypes<-factor(drawmatrix_melt$celltypes,levels = c(
  "EEC2_Iday1.Mural.TE_Iday1","EEC1_Iday1.Mural.TE_Iday1",
  "EEC2_Iday1.Polar.TE_Iday1","EEC1_Iday1.Polar.TE_Iday1",
  "Mural.TE_Iday1.EEC2_Iday1","Mural.TE_Iday1.EEC1_Iday1",
 "Polar.TE_Iday1.EEC2_Iday1","Polar.TE_Iday1.EEC1_Iday1"))
p <- ggplot(drawmatrix_melt, aes(x=celltypes, y=fct_rev(gene_pair), size=-log10(pvalue+0.000001),colour=score )) + 
  geom_point()+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1),strip.text.x = element_text(angle = 90, hjust = 1))+
  scale_colour_gradientn(colours = viridis_pal()(100)[51:100])
ggsave(p,filename = "Figure6D_HumanBlastoidsimplantationday1_selected2.pdf",width = 8,height = 11)

# Figure S6J Calculated pseudotime trajectory for all TE cells by using monocle3
setwd("../Figure6/EMT")
BDTE<-subset(BDall_rm,subset=celltype2%in%c("TE","Mural.TE","Polar.TE","TE_fusion","STB"))
BDTE@meta.data$rm_noise=F
BDTE$rm_noise[BDTE@reductions$umap@cell.embeddings[,1]<0]=T
BDTE<-subset(BDTE,subset=rm_noise==F)
pdf(file = "BDTE_umap_celltype2_test.pdf",width = 10,height = 8)
DimPlot(object = BDTE,reduction = "umap",group.by="celltype2",label = T,pt.size = 0.3)
dev.off()
graycolor<-rep("lightgray",length(unique(BDall$celltype2)))
names(graycolor)<-unique(BDall$celltype2)
# Figure S6J.1
pdf(file = "BDTE_umap_gray.pdf",width = 6,height = 5)
DimPlot(object = BDall,reduction = "umap",group.by="celltype2",pt.size = 0.5)+scale_color_manual(values = graycolor)
dev.off()
BDTE_metadata <- BDTE@meta.data
BDTEmonocledata<-BDTE@assays$RNA$counts
saveRDS(BDTE_metadata,file = "./monocle3_result/BDTE_metadata.Rds")
saveRDS(BDTEmonocledata,file = "./monocle3_result/BDTEmonocledata.Rds")
int.embed <- Embeddings(BDTE, reduction = "umap")
saveRDS(int.embed,file = "./monocle3_result/BDTEint.embed.Rds")
library(monocle3)
library(tidyverse)
library(patchwork)
library(spdep)
library(cowplot)
all_metadata<-readRDS("BDTE_metadata.Rds")
monocledata<-readRDS("BDTEmonocledata.Rds")
gene_annotation <- data.frame(gene_short_name = rownames(monocledata))
rownames(gene_annotation) <- rownames(monocledata)
cds <- new_cell_data_set(monocledata,
                        cell_metadata = all_metadata,
                        gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, preprocess_method="PCA")
cds_BDTE <- cluster_cells(cds)
cds_BDall.embed <- cds_BDTE@int_colData$reducedDims$UMAP
int.embed <- readRDS("BDTEint.embed.Rds")
int_BDall.embed <- int.embed[rownames(cds_BDall.embed),]
cds_BDTE@int_colData$reducedDims$UMAP <- int_BDall.embed
cds_BDTE <- learn_graph(cds_BDTE, verbose = F, use_partition = F, close_loop = F)

pdf(file = "BDTEpseudotime_BDTEseuratUMAP_cluster.pdf",height = 8,width = 9)
p<-plot_cells(cds_BDTE,
             color_cells_by = "mergecluster",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=T,
             group_label_size=4,
             label_principal_points = TRUE,
             label_cell_groups = F,
             cell_size=0.2)
p$layers <- p$layers[-1]
print(p)
dev.off()

cds_BDTE <- order_cells(cds_BDTE, root_pr_nodes='Y_108')
pdf(file = "pseudotime_seuratUMAP_time.pdf",height = 4,width = 4)
p<-plot_cells(cds_BDTE,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
             label_leaves=F,
             label_branch_points=F,
             graph_label_size=1.5,
             group_label_size=4,cell_size=0.5)
p$layers <- p$layers[-1]
print(p)
dev.off()
# Figure S6J.2
emtg<-c("EPCAM","FN1","OCLN","SNAI1")
for (i in 1:length(emtg)) {
  igene<-emtg[i]
  p<-plot_genes_in_pseudotime(cds_BDTE[igene,],color_cells_by="celltype2")
  q<-ggplot(data = p$data,aes(pseudotime, log2(expression+1)))+
    geom_point(aes_string(color = "celltype2"),size=0.3)+geom_line(aes(x = pseudotime, y = log2(expectation+1)),data = p$data)+
    scale_color_manual(values = mycolors)+theme_bw()+ ggtitle(igene)
  eval(parse(text = paste("p",i,"=q",sep = "")))
}
pdf(file = "seuratUMAPgenetrand_EMTgenes.pdf",height = 5,width = 8)
plot_grid(p1,p2,p3,p4,ncol = 2)
dev.off()