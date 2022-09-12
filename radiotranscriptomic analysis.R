library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(dplyr)
library(Scillus)
library(magrittr)
library(R.utils)
library(gplots)
library(ggplot2)
library(plyr)
library(tibble)
library(ggheatmap)

hcc <- readRDS("~hcc.rds")
Idents(hcc)="orig.ident"
hcc@meta.data["group.ident_new"] = "NA"
hcc@meta.data["orig.ident_new"] = "NA"
group_new <- read_csv("~group_new.csv")
hcc@meta.data["group.ident_new"]<- group_new$group[match(Idents(hcc), group_new$orig.ident)]
hcc@meta.data["orig.ident_new"]<- group_new$orig[match(Idents(hcc), group_new$orig.ident)]

Idents(hcc)="celltype_new"
plot1 <- DimPlot(hcc, reduction = "umap", label = TRUE,pt.size = 0.5,raster=FALSE)+ ggsci::scale_color_lancet()
plot2 <- DimPlot(hcc, reduction = "umap", group.by = "orig.ident_new",label = F,pt.size = 0.1,raster=FALSE)+ ggsci::scale_color_lancet()
ggsave("UMAP.png", plot = plot1, width = 7, height = 6)
ggsave("UMAP2.png", plot = plot2, width = 6, height = 6)

library(ggsci)
table1=table(hcc@meta.data$orig.ident_new,hcc@meta.data$celltype_new) 
table2=table(hcc@meta.data$group.ident_new,hcc@meta.data$celltype_new) 
balloonplot(table1)
balloonplot(table2)
stackData <- data.frame(sample=hcc$orig.ident_new,celltype=hcc$celltype_new)
immCelltype <- read.table("~cell_types.txt", header = T, row.names = NULL, sep = "\t")
stackData$celltype2 <- immCelltype$celltype2[match(stackData$celltype,immCelltype$celltype)]
stackData <-stackData[,-2]
stackData <- stackData%>%group_by(sample)%>%table()%>%as.data.frame.matrix()%>%t()
data <- as.data.frame(stackData)%>%rownames_to_column(var = "celltype2")%>%
  gather(key = "sample",value = "count",-celltype2)
plot1 <- ggplot(data, aes(fill=celltype2, y=count, x=sample)) + geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = pal_simpsons()(12))+ theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 0.9))+ xlab("")+ylab("percentage")+ ggsci::scale_color_lancet()
ggsave("celltype_percentage.png", plot = plot1, width = 8, height = 8)

Idents(hcc)="celltype_new"
hcc_b <-subset(hcc, idents = c("B cell", "Plasma cell"))
DefaultAssay(hcc_b) <- "integrated" 
hcc_b <- FindVariableFeatures(hcc_b, selection.method = "vst", nfeatures = 2000)
scale.genes <- rownames(hcc_b)
hcc_b <- ScaleData(hcc_b, features = scale.genes)

hcc_b <- RunPCA(hcc_b, features = VariableFeatures(hcc_b))
ElbowPlot(hcc_b, ndims=10, reduction="pca")
pc.num=1:5

hcc_b <- FindNeighbors(hcc_b, dims = pc.num)
hcc_b <- FindClusters(hcc_b, resolution = 0.3)
table(hcc_b@meta.data$seurat_clusters)
meta <- hcc_b@meta.data

hcc_b <- RunUMAP(hcc_b, dims = pc.num)
embed_umap <- Embeddings(hcc_b, 'umap') 
plot1 <- DimPlot(hcc_b, reduction = "umap", label = TRUE,pt.size = 2,raster=FALSE)+ ggsci::scale_color_lancet()

markers <- FindAllMarkers(hcc_b, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(markers, "hcc_b_markers.csv", row.names = F)
write.csv(markers.top10, "hcc_b_markers.top10.csv", row.names = F)

VlnPlot(hcc_b, features = c("CD79A","BANK1","MS4A1","TNFRSF13C","BCL11A"), pt.size=0,ncol = 2) #B细胞
VlnPlot(hcc_b, features = c("MZB1","SSR4"), pt.size=0,ncol = 2) #Plasma细胞

celltype_b <- read_csv("~celltype_b.csv")
hcc_b@meta.data$celltype_b = "NA"
for(i in 1:nrow(celltype_b)){hcc_b@meta.data[which(hcc_b@meta.data$seurat_clusters == celltype_b$ClusterID[i]),
                                              'celltype_b'] <- celltype_b$celltype[i]}
table(hcc_b@meta.data$celltype_b)
Idents(hcc_b)="celltype_b"
plot1 <- DimPlot(hcc_b, reduction = "umap", label = F,pt.size = 2,raster=FALSE)+ ggsci::scale_color_lancet()
plot2 <- DimPlot(hcc_b, reduction = "umap", group.by='group.ident_new', label = F,pt.size = 2,raster=FALSE)+ ggsci::scale_color_lancet()
ggsave("UMAP_b2.png", plot = plot1, width = 8, height = 7)
ggsave("UMAP_b3.png", plot = plot2, width = 7, height = 7)

plot1 = FeaturePlot(hcc_b, features = "MS4A1",reduction = "umap",pt.size = 3, ncol=1,raster=FALSE)+
  scale_color_gradient2(limits=c(1,3),low='White', mid = "LightGrey",high = 'Firebrick')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),axis.text = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=14)
  )
plot2 = FeaturePlot(hcc_b, features = "MZB1",reduction = "umap",pt.size = 3, ncol=1,raster=FALSE)+
  scale_color_gradient2(limits=c(1,4),low='White', mid = "Snow1",high = 'Firebrick')+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),axis.text = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=14)
  )
ggsave("UMAP_MS4A1.png", plot = plot1, width = 6, height = 7)
ggsave("UMAP_MZB1.png", plot = plot2, width = 6, height = 7)

library(clusterProfiler)
library(org.Hs.eg.db)
hcc_b <- readRDS("~hcc_b.rds")
Idents(hcc_b)= "group.ident_new"
DefaultAssay(hcc_b) <-"RNA"
dge.group <- FindMarkers(hcc_b, ident.1 = 'L', ident.2 = 'H')
sig_dge.group <- subset(dge.group, p_val_adj<0.01&abs(avg_log2FC)>1)

library(ggrepel)
dge.group$threshold = factor(ifelse(dge.group$p_val_adj < 0.01 & abs(dge.group$avg_log2FC) >= 1, 
                                    ifelse(dge.group$avg_log2FC >= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
dge.group$gene <- rownames(dge.group)
plot1 <- ggplot(dge.group,aes(x=avg_log2FC,y= -log10(p_val_adj),color=threshold))+
  geom_point(data = dge.group[dge.group$p_val_adj<0.01&abs(dge.group$avg_log2FC)>1,],size = 4)+ 
  geom_point(data = dge.group[dge.group$p_val_adj>0.01|abs(dge.group$avg_log2FC)<1,],size = 4)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+
  geom_text_repel(data = dge.group[dge.group$p_val_adj<0.1&abs(dge.group$avg_log2FC)>1,],
                  aes(label = gene),size = 3, color = "black",segment.color = "black", show.legend = FALSE )+
  ylab('-log10(P adj value)')+ xlab('log2(FoldChange)')+ 
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) + 
  geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.5) +
  theme_classic(base_line_size = 1)+
  theme(axis.title.x = element_text(size = 15,color = "black", face = "bold"),
        axis.title.y = element_text(size = 15,color = "black",face = "bold", vjust = 1.9, hjust = 0.5,angle = 90),
        legend.title = element_blank(),legend.text = element_text(color="black",size = 10, face = "bold"), 
        axis.text.x = element_text(size = 13, color = "black", face = "bold",vjust = 0.5,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 13,color = "black",face = "bold",vjust = 0.5,hjust = 0.5, angle = 0) 
  )
ggsave("volcano_b.png", plot = plot1, width = 6, height = 6)

library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
library(cogena)
library(BaseSet)
library(GSEABase)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

Idents(hcc_b) <- "celltype_b"
expr <- AverageExpression(hcc_b, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
head(expr)

gsva.res <- gsva(expr, genesets, method="gsva") 
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
pheatmap::pheatmap(gsva.res, show_colnames = T, 
                   scale = "row",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

expr1=as.matrix(hcc_b@assays$RNA@data)
geneset = getGmt("~h.all.v7.5.1.symbols.gmt")
hall <- gsva(expr1, gset.idx.list = geneset, kcdf="Gaussian",method = "zscore",parallel.sz=5)

de_gsva <- function(exprSet,meta,compare = NULL){
  allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  fit <- lmFit(exprSet,design)
  if(length(unique(meta))==2){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  return(allDiff)
}
meta <- hcc_b@meta.data[,c("group.ident_new")]
diff =de_gsva(exprSet = hall ,meta = meta,compare = "L-H")

idiff <-diff[["L-H"]]
df <- data.frame(ID = rownames(idiff), score = idiff$t)
df$group =sapply(1:nrow(idiff),function(x){
  if(idiff[x,"logFC"]>0 & idiff[x,"adj.P.Val"]<0.05){return("up")}
  else if(idiff[x,"logFC"]<0 & idiff[x,"adj.P.Val"]<0.05){return("down")}
  else{return("noSig")}
})

df = filter(df, abs(df$score)>1.5)
df$hjust = ifelse(df$score>0,1,0)
df$nudge_y = ifelse(df$score>0,-0.1,0.1)
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
limt = max(abs(df$score))
plot1 = ggplot(sortdf, aes(ID, score,fill=group))+geom_bar(stat = 'identity',alpha = 0.7)+
  scale_fill_manual(breaks=c("down","noSig","up"),values = c("#008020","grey","#08519C")) +
  geom_text(data = df, aes(label = ID, y = nudge_y),nudge_x =0,nudge_y =0,hjust =df$hjust,size = 5)+
  scale_y_continuous(limits=c(-limt,limt))+  coord_flip() +   theme_bw() +   theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+labs(x = "Hallmark pathways",y="t value of GSVA score",title = "Low versus High")+
  theme(plot.title = element_text(hjust = 0.5,size = 18),axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),axis.line = element_blank(),
        axis.ticks.y = element_blank(),legend.position = limt)
plot1
ggsave("gsva_b.png", plot = plot1, width = 11, height = 9)