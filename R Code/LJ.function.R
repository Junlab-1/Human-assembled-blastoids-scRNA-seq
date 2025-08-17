suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(tibble)
  library(miloR)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(igraph)
  library(batchelor)
  library(uwot)
  library(Seurat)
  library(DoubletFinder)
})

# QC function
plotQC<-function(seurat_obj,cutinfo="0"){
  if (!dir.exists("./QC")) dir.create("./QC", recursive = TRUE)
  obj_name <- deparse(substitute(seurat_obj))
  pdf(file = paste0("./QC/",obj_name,"_count_detail",cutinfo,".pdf"),width = 10,height = 8)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  feature_count.data<-cbind(VlnPlot(seurat_obj, features = "nFeature_RNA")$data,
                            VlnPlot(seurat_obj, features = "nCount_RNA")$data)[,-4]
  fc.plot<-ggplot(feature_count.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point()+theme_bw()+
    xlab("All counts per cell")+ylab("Gene numbers per cell")
  ggsave(fc.plot,filename = paste0("./QC/",obj_name,"_CountGenepercent",cutinfo,".pdf"),width = 8,height = 8)
}

# rm doublelets function
RMdoublets<-function(seurat_obj){
  datai_doublet <- seurat_obj
  obj_name <- deparse(substitute(seurat_obj))
  datai_doublet <- NormalizeData(object = datai_doublet, verbose = FALSE)
  datai_doublet <- FindVariableFeatures(object = datai_doublet,selection.method = "vst",nfeatures = 3000, verbose = FALSE)
  datai_doublet <- ScaleData(datai_doublet)
  datai_doublet <- RunPCA(object = datai_doublet,npcs=100)
  datai_doublet <- FindNeighbors(object = datai_doublet, dims = 1:50)
  datai_doublet <- FindClusters(object = datai_doublet, resolution = 0.5)
  datai_doublet <- RunUMAP(object = datai_doublet, dims = 1:50)
  sweep.res.list <- paramSweep_v3(datai_doublet, PCs = 1:10, sct = FALSE)
  for(j in 1:length(sweep.res.list)){
    if(length(sweep.res.list[[j]]$pANN[is.nan(sweep.res.list[[j]]$pANN)]) != 0){
      if(j != 1){
        sweep.res.list[[j]] <- sweep.res.list[[j - 1]]
      }else{
        sweep.res.list[[j]] <- sweep.res.list[[j + 1]]
      }
    }
  }
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk_v <- as.numeric(as.character(bcmvn$pK))
  pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
  nExp_poi <- round(0.03*ncol(datai_doublet))
  #doublet table
  datai_doublet <- doubletFinder_v3(datai_doublet, PCs = 1:10, pN = 0.25, 
                                    pK = pk_good, nExp = nExp_poi, 
                                    reuse.pANN = FALSE, sct = FALSE)
  colnames(datai_doublet@meta.data)[ncol(datai_doublet@meta.data)]="Doublet_type"
  # show doublets cells
  pdf(file = paste("./QC/",obj_name,"_UMAP_doublets.pdf",sep = ""),width = 7,height = 6)
  print(DimPlot(datai_doublet,reduction = "umap",group.by = "Doublet_type"))
  dev.off()
  if (identical(colnames(seurat_obj),colnames(datai_doublet))) {
    seurat_obj@meta.data$Doublet_type=datai_doublet@meta.data$Doublet_type
  }else(
    print("Cell Name is not same!")
  )
  seurat_obj<-subset(x = seurat_obj,subset= Doublet_type == "Singlet")
  return(seurat_obj)
}

library(patchwork)
library(viridis)
library(rlang)
# 3 run NM Lanner's model (PMID:39543283)
runNMModel<-function(alldata.counts,alldata.counts.meta,filename){
  source("~/reference/NM_model/HumanEarlyEmbryoRefProjection/HuEm_stable_reference_projection_tool/main.function2.R")
  if (ncol(alldata.counts)<200) {
    para.list$runMiloR=FALSE
  }
  if (para.list$runMiloR) {
    milo_out <- FunMiloCal(alldata.counts,alldata.counts.meta, temp.cal=TRUE)
  }else{
    milo_out <- FunMiloCal(alldata.counts,alldata.counts.meta, temp.cal=FALSE)
  }
  
  sf_out <- FunCalSF(milo_out)
  sf_out$query.sce.cor.out <- FunCalCor(sf_out) 
  sf_out$query.sce.cor.out.mean <- sf_out$query.sce.cor.out%>% gather(ref_cell,cor,-query_cell)%>%
    group_by(query_cell) %>% top_n(20,cor)%>% dplyr::summarise(cor_top_mean=mean(cor))
  sf_out$NWIN <- FunNWIN(sf_out)
  mnn.pairs.list <- list()
  for (ref_name in c("SPH2016","D3post","Meistermann_2021","CS7","nBGuo","Yan2013")) {
    print(ref_name)
    mnn.pairs.list[[ref_name]] <-FunCalMNN_each(sf_out,ref_name)
  }
  
  predict_out <-  FunProjCal(mnn.pairs.list,
                             query.sce.ob=sf_out$query.sce.ob,
                             query.sce.cor.out=sf_out$query.sce.cor.out,
                             temp.max=sf_out$NWIN$temp.max,
                             D2_umap_model=ref.umap,
                             Dmulti_umap_model=ref_umap_nDim_model,
                             cor.cutoff=para.list$cor.cutoff)
  
  predict_out$HS <- milo_out$HS
  predict_out$raw.meta <- milo_out$raw.meta
  predict_out$query.sce.cor.out.mean <- sf_out$query.sce.cor.out.mean
  predict_out <- FunPredAnno(predict_out,cor.cutoff=para.list$cor.cutoff)
  predict_out$full.anno %>% group_by(pred_EML) %>% dplyr::summarise(nCell=n_distinct(query_cell)) %>% arrange(desc(nCell))
  # add celltype information to predict_out$umap
  predict_out$umap2<-predict_out$umap
  predict_out$umap2$pred_EML<-NA
  predict_out$umap2$sub_pred_EML<-NA
  predict_out$umap2$pred_EML<-predict_out$anno$pred_EML[match(predict_out$umap2$cell,predict_out$anno$query_cell)]
  predict_out$umap2$sub_pred_EML<-predict_out$anno$sub_pred_EML[match(predict_out$umap2$cell,predict_out$anno$query_cell)]
  # save data
  saveRDS(predict_out,file = paste0(filename,"_predict_out.rdata"))
  return(predict_out)
}

# draw marker dotplot
callmarkerdot<-function(markers,obj,group.by=NULL,fname=NULL){
  data.anno<-unique(markers)
  data.usage <- DotPlot(obj,features = unique(data.anno[,1]), assay = 'RNA',group.by = group.by)$data
  if (is.null(fname)) {
    if (is.null(group.by)) {
      fname<- paste(deparse(substitute(obj)),"active.ident",sep = "_")
    }else{
      fname<- paste(deparse(substitute(obj)),deparse(substitute(group.by)),sep = "_")
    }
  }
  colnames(data.anno)<-c("features.plot","label")
  data.anno<-as.data.frame(data.anno)
  df.plot <- plyr::join(data.usage,data.anno)
  df.plot<-na.omit(df.plot)
  data.cluster<-"a"
  for (i in levels(df.plot$id)){
    i_list<-c(data.usage$avg.exp.scaled[data.usage$id == i],data.usage$pct.exp[data.usage$id == i])
    if (sum(data.cluster=="a")) {
      data.cluster=i_list
    }else{
      data.cluster=rbind(data.cluster,i_list)
    }
  }
  row.names(data.cluster)<-levels(df.plot$id)
  d <- dist(data.cluster, method = "euclidean")
  fit2 <- hclust(d, method="ward.D")
  pdf(file = paste0(fname,"_ClusterTree.pdf"),width = 8,height =6 )
  plot(fit2, col = "black", col.main = "#45ADA8", col.lab = "#7C8071",
       col.axis = "#F38630", lwd = 3, lty = 1, 
       axes = F, hang = -1)
  dev.off()
  cluster_order<-row.names(data.cluster)[fit2$order]
  df.plot$id <- factor(df.plot$id,levels = cluster_order)
  p <- ggplot(df.plot,aes(x=features.plot,y = as.numeric(id),size = (pct.exp^2)/100, color = avg.exp.scaled))+
    geom_point() + scale_size("% detected", range = c(0,6)) +  
    scale_color_gradientn(colours = c("white","red"),
                          guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                          name = "Average\nexpression") + cowplot::theme_cowplot() +
    ylab("") + xlab("Markers") + theme_bw() + 
    scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+ #复制 y轴 代替边框效果 
    facet_grid(~label, scales="free_x",space = "free")+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(fname,"_dotplot.pdf"),p,width = nrow(data.anno)*0.2,height = 6)
}

# draw function annotation result
drawGOKEGG<-function(ego_mat,icolor="Up"){
  obj_name <- deparse(substitute(ego_mat))
  if (nrow(ego_mat) == 0) {
    message(paste0(obj_name, " is empty. Skipping analysis."))
    return(NULL)
  }
  ego_mat<-ego_mat[1:min(10,nrow(ego_mat)),]
  p<-ggplot(ego_mat,aes(x=reorder(Description,Count),y=Count,fill=Count))+
    geom_bar(stat = "identity",colour="NA",width = 0.9,position = position_dodge(0.7))
  p<-p+xlab("")+scale_y_continuous(name = "Gene Count")+coord_flip()
  p<-p+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  if (icolor=="Up") {
    p<-p+scale_fill_gradient(low="#FF9289",high="#FF9289",name="Genes")
  }else if (icolor=="Down") {
    p<-p+scale_fill_gradient(low="#00DAE0",high="#00DAE0",name="Genes")
  }
  p<-p+scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(ego_mat$Count)/10)*10+1))
  p<-p+theme(text=element_text(size=16))+labs(title = obj_name)
  filename <- paste0(obj_name, ".pdf")
  if (file.exists(filename)) {
    filename <- paste0(obj_name, "_2.pdf")
  }
  ggsave(p,filename = filename,height = 5,width = 10)
}

