library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(patchwork)
library(limma)
library(cowplot)
library(AnnotationHub)
library(ensembldb)
library(scales)
library(Repitools)
library(plyr)
library(gridExtra)

AS9232.data <- Read10X(data.dir = "32/")
AS9234.data <- Read10X(data.dir = "34/")

AS9232 <- CreateSeuratObject(counts = AS9232.data,project="ctl")
AS9234 <- CreateSeuratObject(counts = AS9234.data,project="vacc")

pigs <- merge(AS9232, y = AS9234, add.cell.ids = c("ctl","vacc"), project = "ARGILAGUETJOR_02")
View(pigs@meta.data)

### Refine gene annotation (replace some ensembl IDs with known gene names)

gtf <- rtracklayer::import('Sus_scrofa.Sscrofa11.1.103_KP055815.gtf')
gtf_df=annoGR2DF(gtf[,c("type","gene_id","gene_biotype","gene_name")])
gtf_gene = gtf_df %>%
  dplyr::filter(str_detect(type, "gene")) %>%
  mutate(feature = gene_name %>% 
           is.na %>%
           ifelse(gene_id, gene_name))
dict=read.table("ensembl_to_replace",h=F)
rownames(gtf_gene)=rownames(pigs[["RNA"]])
pigs[["RNA"]]@meta.features=gtf_gene[,c("chr","strand","gene_biotype")]
features_metadata=pigs[["RNA"]]@meta.features
rownames(features_metadata)=plyr::mapvalues(rownames(features_metadata),from=as.factor(dict$V1), to=as.factor(dict$V2))

### Calculate percentage of MT and virus RNAs

mito_genes=rownames(subset(features_metadata,gene_biotype=="Mt_rRNA"|gene_biotype=="Mt_tRNA"))
virus_genes=rownames(subset(features_metadata,chr=="KP055815.1"))
pigs$percent.mt <- PercentageFeatureSet(pigs, features=mito_genes)
pigs$percent.virus <- PercentageFeatureSet(pigs, features=virus_genes)

### Quality control plots

VlnPlot(pigs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.virus"),ncol = 2, pt.size=0)
FeatureScatter(pigs, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pigs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(pigs,  feature1 = "nFeature_RNA",feature2 = "percent.virus")
FeatureScatter(pigs, feature1 = "percent.mt", feature2 = "percent.virus")

### Filter low quality cells

pigs <- subset(pigs, subset = nFeature_RNA > 200 &
                 nFeature_RNA < 2000 & percent.mt < 5)

### SCT, integration, UMAP and cluster detection

pigs.split <- SplitObject(pigs, split.by = "orig.ident")
pigs.split <- lapply(X = pigs.split, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = pigs.split, nfeatures = 3000)
pigs.split <- PrepSCTIntegration(object.list = pigs.split, anchor.features = features)
pigs.anchors <- FindIntegrationAnchors(object.list = pigs.split, normalization.method = "SCT",anchor.features = features)
pigs.integrated <- IntegrateData(anchorset = pigs.anchors, normalization.method = "SCT") %>% 
                   RunPCA(verbose = FALSE) %>% 
                   RunUMAP(reduction = "pca", dims = 1:30) %>% 
                   FindNeighbors(reduction = "pca", dims = 1:30) %>% 
                   FindClusters(resolution = c(0.5,0.8),random.seed=123,save.SNN=TRUE)

#keep_feature <- Matrix::rowSums(counts(pigs.integrated) > 0) > 0
#pigs.integrated <- pigs.integrated[keep_feature, ]

DimPlot(pigs.integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(pigs.integrated, reduction = "umap", split.by = "orig.ident")
DimPlot(pigs.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

### Add cycle score and dissociation score

pigs.integrated <- CellCycleScoring(pigs.integrated,g2m.features=cc.genes$g2m.genes,s.features=cc.genes$s.genes, assay="SCT")
VlnPlot(pigs.integrated, features = c("S.Score","G2M.Score"), pt.size=0)

genes.dissoc <- c("ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL3", "CXCL2", "CXCL1", "DNAJA1", "DNAJB1", "DUSP1", 
                  "EGR1", "FOS", "FOSB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1A", "HSPA1B", 
                  "HSPA8", "HSPB1", "HSPE1", "HSPH1", "ID3", "IER2", "JUN", "JUNB", "JUND", "MT1X", "NFKBIA", 
                  "NR4A1", "PPP1R15A", "SOCS3", "ZFP36")

pigs.integrated<- AddModuleScore(pigs.integrated, features = list(genes.dissoc), ctrl.size =20, name = "genes_dissoc")

FeaturePlot(pigs.integrated, features ="percent.virus")
FeaturePlot(pigs.integrated, features ="genes_dissoc1")
FeaturePlot(pigs.integrated, features ="S.Score")
FeaturePlot(pigs.integrated, features ="G2M.Score")

### Conserved markers for each cluster

DefaultAssay(pigs.integrated) <- "RNA"
get_conserved <- function(cluster){
  FindConservedMarkers(pigs.integrated,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

conserved_markers <- map_dfr(0:14, get_conserved)

### Subclustering clusters 7,8,9, and 14

pigs.integrated <- FindSubCluster(pigs.integrated, "7", graph.name="integrated_nn",subcluster.name = "subcluster7",  resolution = 0.2, algorithm = 1)
pigs.integrated <- FindSubCluster(pigs.integrated, "8", graph.name="integrated_nn",subcluster.name = "subcluster8",  resolution = 0.3, algorithm = 1)
pigs.integrated <- FindSubCluster(pigs.integrated, "9", graph.name="integrated_nn",subcluster.name = "subcluster9",  resolution = 0.5, algorithm = 1)
pigs.integrated <- FindSubCluster(pigs.integrated, "14", graph.name="integrated_nn",subcluster.name = "subcluster14",  resolution = 0.3, algorithm = 1)

pigs.integrated$sub_clusters=as.character(Idents(pigs.integrated))

pigs.integrated$sub_clusters[Cells(subset(pigs.integrated, idents = 7))] <- pigs.integrated$subcluster7[Cells(subset(pigs.integrated, idents = 7))]
pigs.integrated$sub_clusters[Cells(subset(pigs.integrated, idents = 8))] <- pigs.integrated$subcluster8[Cells(subset(pigs.integrated, idents = 8))]
pigs.integrated$sub_clusters[Cells(subset(pigs.integrated, idents = 9))] <- pigs.integrated$subcluster9[Cells(subset(pigs.integrated, idents = 9))]
pigs.integrated$sub_clusters[Cells(subset(pigs.integrated, idents = 14))] <- pigs.integrated$subcluster14[Cells(subset(pigs.integrated, idents = 14))]

Idents(pigs.integrated)=pigs.integrated$sub_clusters

### UMAP by subclusters

DimPlot(pigs.integrated, reduction = "umap", group.by = "sub_clusters", label = TRUE, repel = TRUE,label.size = 3)+theme(legend.key.size = unit(0.2, 'cm'))

### Rename subcluster 9_1 as 9_1

pigs.integrated <- RenameIdents(object = pigs.integrated, '9_2' = "9_1")

### Conserved markers for each subcluster

Idents(pigs.integrated)=pigs.integrated@meta.data$sub_clusters
DefaultAssay(pigs.integrated) <- "RNA"
get_conserved <- function(cluster){
  FindConservedMarkers(pigs.integrated,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

conserved_markers <- map_dfr(levels(Idents(pigs.integrated)), get_conserved)

### Cell type annotation

pigs.integrated <- RenameIdents(pigs.integrated,'0'='CD4_Tcells_1','1'='Bcells_1','2' = 'CD8_Tcells','3'='Bcells_2',
                                       '4'='Bcells_3','5'='CD4CD8_Tcells','6'='Plasmablasts','7_0'='CTLs','7_1'='NKs',
                                       '8_0'='apop_Bcells','8_1'='apop_Tcells','8_2'='apop_cells','9_0'='Proliferating_Tcells','9_1'='GC_Bcells',
                                       '9_3'='Proliferating_CTLs','10'='CD4_Tcells_2',
                                       '11'='Plasma_cells','12'='gd_Tcells','13'='CD4_Tcells_3','14_0'='cDC_1','14_1'='macrophages','14_2'='cDC_2')
 
### UMAP by cell types (Figure 5A and Figure 7A)

DimPlot(pigs.integrated,label=F)
DimPlot(pigs.integrated,label=F,split.by = "orig.ident")

### Differential expression for each cell type

pigs.integrated$celltype.stim <- paste(Idents(pigs.integrated),pigs.integrated$orig.ident,sep="_")
pigs.integrated$celltype <- Idents(pigs.integrated)
Idents(pigs.integrated) <- "celltype.stim"
 
cell_types=levels(x = pigs.integrated$celltype)
 
for (i in cell_types){ 
   try({
     ident1 <- paste0(i,"_vacc")
     ident2 <- paste0(i,"_ctl")
     
     condition.diffgenes <- FindMarkers(pigs.integrated, ident.1 = ident1, ident.2= ident2, min.pct=0.25, logfc.threshold=0.25, slot="data")
     write.csv(condition.diffgenes, file=paste0(i,"_DE_vacc_vs_ctl_2reps.csv"))
   })
}

### UMAPs interesting markers (Figure 5C)

umap_genes=c("TYROBP","CST3","CTSL","C1QC","CD68","ENSSSCG00000036618","FCER1G","FSCN1","KLRB1","KLRK1","NKG7","PRF1","CD6","CD5","CD8A","CD8B","CD4","TRDC","CD2","CD3E","PRDM1","IRF4","PAX5","CD19","CD79A","CD79B","GCSAM","DUT","IRF8","HLA-DRA","ENSSSCG00000001455","XBP1","FOXP3","CXCL10","CD14","NCR1","TBX21","GATA3","CXCR5","CCR7","CD27")
umap_genes_q10=c("CD4","TRDC","CD14","CD6","CD27","CXCR5","PRDM1","PRF1")

for (i in umap_genes){
	FeaturePlot(pigs.integrated, features = i, min.cutoff = "q10",max.cutoff = "q90")
}
for (i in umap_genes_q10){
	FeaturePlot(pigs.integrated, features = i,max.cutoff = "q90")
}

### Violin plot percentage virus (Figure S12)

VlnPlot(pigs.integrated, features = "percent.virus", split.by = "orig.ident")+theme(axis.text = element_text(size = 6))

### Violin plots interferon genes (Figure 6 and Figure S11)

interferon_genes=c("MX2","IFIT1","IFITM3","ISG15","OAS2","ISG20","BST2","MX1","EIF2AK2","GBP2","IFI44L","IFIT2")
for (i in interferon_genes){
VlnPlot(pigs.integrated, features = i, split.by = "orig.ident")+theme(axis.text = element_text(size = 6))
}

### Dot plot subcluster markers (Figure 5B)

genes_dotplot=c("TYROBP","CST3","CTSL","C1QC","CD68","ENSSSCG00000036618","FCER1G","FSCN1","CD300C","KLRB1","KLRK1","NKG7","GZMA.1","FOXP3","CD6","CD5","CD8A","CD8B","CD4","TRDC","GATA3","CD2","CD3E","LEF1","TCF7","CCR7","CD27","PAX5","CD19","CD79A","CD79B","HLA-DRA","ENSSSCG00000001455","IRF4","XBP1","MZB1","JCHAIN","GCSAM")
Idents(pigs.integrated)=pigs.integrated@meta.data$celltype
levels(pigs.integrated)=c("cDC_1","cDC_2","macrophages","CD4_Tcells_1","CD4_Tcells_2","CD4_Tcells_3","CD4CD8_Tcells","CD8_Tcells","CTLs","NKs","Proliferating_CTLs","Proliferating_Tcells","gd_Tcells","Bcells_1","Bcells_2","Bcells_3","Plasmablasts","GC_Bcells","Plasma_cells","apop_cells","apop_Tcells","apop_Bcells")
dp=DotPlot(pigs.integrated,assay="RNA", features = genes_dotplot,dot.scale = 3,cluster.idents = FALSE)+
  theme(axis.text.x = element_text(angle=90))+theme(axis.text = element_text(size = 6,hjust=0.95,vjust=0.2))+
  coord_flip()+
  xlab('markers')+ ylab('cell_types')

### Violin plots GZMA.1,CXCL10,CCR7,CXCR4,PFN1 (Figure 8 and Figure S13)

VlnPlot(object = pigs_integrated, features = 'GZMA.1', idents =c('NKs','CTLs','Proliferating_CTLs'), split.by='orig.ident')+xlab('cell types')
VlnPlot(object = pigs_integrated, features = 'CXCL10', idents =c('Plasmablasts','CD4_Tcells_2','cDC_1'), split.by='orig.ident')+xlab('cell types')
VlnPlot(object = pigs_integrated, features = 'CCR7', idents =c('CTLs'), split.by='orig.ident')+xlab('cell type')+ theme(legend.position="none")
VlnPlot(object = pigs_integrated, features = 'CXCR4', idents =c('CTLs'), split.by='orig.ident')+xlab('cell type')+ theme(legend.position="none",axis.title.y = element_blank())
VlnPlot(object = pigs_integrated, features = 'PFN1', idents =c('CTLs'), split.by='orig.ident')+xlab('cell type')+ theme(legend.position="none",axis.title.y = element_blank())





