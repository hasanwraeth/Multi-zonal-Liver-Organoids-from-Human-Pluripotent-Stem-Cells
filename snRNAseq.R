#Library---------------
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(metap)
library(presto)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(velocyto.R)
library(pagoda2)
library(SeuratWrappers)
library(monocle3)
library(CellChat)

#10X------------------
HCLO <- Read10X(data.dir = "./outs/filtered_feature_bc_matrix/")
HCLO <- CreateSeuratObject(counts = HCLO, project = "HCLO", min.cells = 3, min.features = 200)
HCLO

HCLO[["percent.mt"]] <- PercentageFeatureSet(HCLO, pattern = "^MT-", assay='RNA')
head(HCLO@meta.data, 5)
VlnPlot(HCLO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(HCLO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HCLO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

HCLO <- subset(HCLO, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
HCLO <- NormalizeData(HCLO)#, normalization.method = "LogNormalize", scale.factor = 10000)
HCLO <- FindVariableFeatures(HCLO, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(HCLO), 100)

plot1 <- VariableFeaturePlot(HCLO)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(HCLO)
HCLO <- ScaleData(HCLO, features = all.genes)
#HCLO <- ScaleData(HCLO)
HCLO <- RunPCA(HCLO, features = VariableFeatures(object = HCLO))
VizDimLoadings(HCLO, dims = 1:4, reduction = "pca")
DimPlot(HCLO, reduction = "pca")
DimHeatmap(HCLO, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(HCLO, dims = 1:15, cells = 500, balanced = TRUE)

HCLO <- JackStraw(HCLO, num.replicate = 100)
HCLO <- ScoreJackStraw(HCLO, dims = 1:20)
JackStrawPlot(HCLO, dims = 1:15)
ElbowPlot(HCLO)

HCLO <- FindNeighbors(HCLO, dims = 1:15)
HCLO <- FindClusters(HCLO, resolution = 0.5)
head(Idents(HCLO), 5)
HCLO <- RunUMAP(HCLO, dims = 1:15)

                     
DimPlot(HCLO, reduction = "umap")

HCLO[["sub_annotations"]] <- Idents(HCLO)

saveRDS(HCLO, file = "./HCLO_cr.rds")


HCLO.markers <- FindAllMarkers(HCLO, only.pos = TRUE, min.pct = 0.25,
                               logfc.threshold = 0.25)
VlnPlot(HCLO, features = c("AFP"))
RidgePlot(HCLO, features = c("AFP"))
DotPlot(HCLO, features = c("ALB","FABP1","SOD1","DBI","GSS","ALDH1A2","GHR","HK1","WNT5A","HIF1A","GLS","CDH1","CPS1","HGD",'NFE2L2'))
FeaturePlot(HCLO, features = c("GSTA2"), order=T) #PAXBP1, TUT4 #ASGR1 #SULT1C2

HCLO.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(HCLO, features = top10$gene) + NoLegend()

DoHeatmap(HCLO, features = c("ALB","FABP1","SOD1","DBI","GSS",
                             "ALDH1A2","GHR","HK1","WNT5A","HIF1A",
                             "GLS","CDH1","CPS1","HGD",'NFE2L2')) + NoLegend()

ma5=HCLO.markers[HCLO.markers$cluster==8,]

HCLO <- SCTransform(HCLO, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = o,
  query = HCLO,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50)

pbmc3k <- MapQuery(
  anchorset = anchors,
  query = pbmc3k,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap")

#Ele----------------
library(ELeFHAnt)
human_tissues
EHCLO = CelltypeAnnotation(reference = o,
                                            query = HCLO, downsample = TRUE, downsample_to = 1000,
                                            validatePredictions = TRUE, annotationCol = "All_Integrated_Manual",
                                            species = "human", tissue = "Liver")
DimPlot(EHCLO, reduction = "umap")
DimPlot(EHCLO, group.by = "ELeFHAnt_Ensemble_CelltypePrediction", label = T, reduction = "umap", label.size = 6, repel = T) + NoLegend()
DimPlot(EHCLO, group.by = "ELeFHAnt_RF_CelltypePrediction", label = T, reduction = "umap", label.size = 6, repel = T) + NoLegend()
DimPlot(EHCLO, group.by = "ELeFHAnt_SVM_CelltypePrediction", label = T, reduction = "umap", label.size = 6, repel = T) + NoLegend()
DimPlot(EHCLO, group.by = "ELeFHAnt_LR_CelltypePrediction", label = T, reduction = "umap", label.size = 6, repel = T) + NoLegend()

#Integration-----------------------
hep.combined <- merge(o, y = HCLO, add.cell.ids = c("Liver","HCLO"),
                      project = "hep")
hep.list <- SplitObject(hep.combined, split.by = "orig.ident")
hep.list <- lapply(X = hep.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = hep.list)
hep.anchors <- FindIntegrationAnchors(object.list = hep.list, anchor.features = features)
hep.combined <- IntegrateData(anchorset = hep.anchors)
DefaultAssay(hep.combined) <- "integrated"

hep.combined <- ScaleData(hep.combined, verbose = FALSE)
hep.combined <- RunPCA(hep.combined, npcs = 30, verbose = FALSE)
#hep.combined <- RunHarmony(hep.combined, "sub_annotations")
hep.combined <- RunUMAP(hep.combined, reduction = "pca", dims = 1:15)
hep.combined <- FindNeighbors(hep.combined, reduction = "pca", dims = 1:15)
hep.combined <- FindClusters(hep.combined, resolution = 1) 
#resolution 0.4 to 1.2 for 3k cells controls cluster no., higher for more cells

saveRDS(hep.combined, file = "./hepcr.rds")

DimPlot(hep.combined, reduction = "umap", group.by = "sub_annotations")
DimPlot(hep.combined, reduction = "umap", label = TRUE, repel = TRUE)


FeaturePlot(hep.combined, features = c("APOC4")) #grey94, red, order =T

DotPlot(hep.combined, features = c("ALDH1A2","GHR","HIF3A","GLS","CDH1","ALB"), group.by="sub_annotations")



#Monocle3-----------------
#devtools::install_github("satijalab/seurat-wrappers")

cds <- SeuratWrappers::as.cell_data_set(HCLO)
head(colData(cds))
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
head(counts(cds))


recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- HCLO@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- HCLO@reductions$umap@cell.embeddings

cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj

cds <- learn_graph(cds, use_partition = F)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = colnames(cds[, clusters(cds) == "Unidentified"]))

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

head(pseudotime(cds), 10)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, ident, fill = ident)) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime), fill = ident)) + geom_boxplot()

#trace('calculateLW', edit = T, where = asNamespace("monocle3"))
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
FeaturePlot(HCLO, features = c("RPL22", "RERE", "ENO1", "KAZN"))
FeaturePlot(HCLO, features = c("AFP", "GATA4", "CDX2","DPP4"))

HCLO$pseudotime <- pseudotime(cds)
FeaturePlot(HCLO, features = "pseudotime")

RidgePlot(HCLO, features = c("RPL22", "RERE", "ENO1"),
          sort = T)

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("AFP"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = monocle3_pseudotime)

#saveRDS(cds, 'HCLO_cds.rds')
                     
#GSEA------------------
HCLO.genes = wilcoxauc(HCLO, 'seurat_clusters')  
dplyr::count(HCLO.genes, group)

listMarts()
ensembl = useEnsembl(biomart = 'ensembl')
View(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
View(listAttributes(ensembl))
View(listFilters(ensembl))

annoM = getBM(attributes = c('ensembl_gene_id','external_gene_name'),
              filters = c('external_gene_name'),
              values = unique(HCLO.genes$feature), mart = ensembl)

annodfM = left_join(HCLO.genes, annoM, 
                    by=c('feature' = 'external_gene_name'))
annodfM1 = annodfM[annodfM$padj < 0.05, ]
annodfM1 = annodfM1[abs(annodfM1$logFC) > 0.25, ]

ent_gene = getBM(attributes = c('entrezgene_id','external_gene_name'), 
                 filters = c('ensembl_gene_id'),
                 values = annodfM1$ensembl_gene_id, mart = ensembl)
ent_gene=ent_gene[order( ent_gene[,2] ),]
ent_gene = na.omit(ent_gene)
ent_gene1 =ent_gene$entrezgene_id
ent_gene1 = as.character(ent_gene1)
ent_uni = getBM(attributes = c('entrezgene_id','external_gene_name'), 
                filters = c('external_gene_name'),
                values = HCLO.genes$feature, mart = ensembl)
ent_uni=ent_uni[order(ent_uni[,2] ),]
ent_uni = na.omit(ent_uni)
ent_uni1 =ent_uni$entrezgene_id
ent_uni1 = as.character(ent_uni1)
ego = enrichGO(gene = ent_gene1, OrgDb = org.Hs.eg.db, 
               ont = "BP", universe = ent_uni1, readable = TRUE)

barplot(ego, showCategory = 10, font.size = 20)
dotplot(ego, showCategory = 20)


#annodfM2=annodfM1[annodfM1$group==0,]
#annodfM2=annodfM2[(annodfM2$feature %in% ent_gene$external_gene_name),]
annodfM1=annodfM1[order( annodfM1[,1] ),]

gl = (annodfM1$logFC)
names(gl)=ent_uni1
gl=sort(gl, decreasing = T)
gl = na.omit(gl)


kggl=gseGO(geneList=gl, 
           ont ="ALL", 
           keyType = "ENTREZID", 
           nPerm = 10000, 
           minGSSize = 3, 
           maxGSSize = 800, 
           pvalueCutoff = 0.05, 
           verbose = TRUE, 
           OrgDb = org.Hs.eg.db, 
           pAdjustMethod = "BY")

gseaplot2(kggl, geneSetID = 5720, 
          title=kggl$Description[5720], base_size=20)

View(as.data.frame(kggl))
kggl2=as.data.frame(kggl)
kggl2$sid=1:nrow(kggl)
View(kggl2)


#saveRDS(kggl, 'kggl.rds')

#loom and RNA velocity-------------------
#HCLO.loom <- as.loom(HCLO, "HCLO.loom", verbose = FALSE)
HCLO.loom

#emat <- HCLO.loom$spliced
# this dataset has already been pre-filtered, but this is where one woudl do some filtering
#emat <- emat[,colSums(emat)>=1e3]

HCLO <- RunVelocity(HCLO, deltaT = 1, kCells = 25, fit.quantile = 0.02)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

HCLO$barcode <- colnames(HCLO)
HCLO$UMAP_1 <- HCLO@reductions$umap@cell.embeddings[,1]
HCLO$UMAP_2 <- HCLO@reductions$umap@cell.embeddings[,2]
write.csv(HCLO@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(HCLO, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(HCLO@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)

HCLO <- ReadVelocity(file = "./velocyto/HCLO.loom")
HCLO <- as.Seurat(x = HCLO)
HCLO <- SCTransform(object = HCLO, assay = "spliced")
HCLO <- RunPCA(object = HCLO, verbose = FALSE)
HCLO <- FindNeighbors(object = HCLO, dims = 1:5)
HCLO <- FindClusters(object = HCLO)
HCLO <- RunUMAP(object = HCLO, dims = 1:5)
HCLO <- RunVelocity(object = HCLO, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = HCLO)))
names(x = ident.colors) <- levels(x = HCLO)
cell.colors <- ident.colors[Idents(object = HCLO)]
names(x = cell.colors) <- colnames(x = HCLO)
show.velocity.on.embedding.cor(emb = Embeddings(object = HCLO, reduction = "umap"), vel = Tool(object = HCLO, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
#saveRDS(HCLO, 'HCLO_veloR.rds')




















