#library------------------
library(dplyr)
library(Seurat)
library(patchwork)
#library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
#library(metap)
library(SingleCellExperiment)
library(miloR)
library(scrabbitr)
library(scater)
library(scran)
library(viridis)
library(anndata)

#sce conversion------------------
data <- read_h5ad("./Downloads/Hepatocytes_primary_ad.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
HCLOsc=as.SingleCellExperiment(data, assay='RNA')

plotReducedDim(HCLOsc, colour_by="source", dimred = "UMAP") 

#milo conversion------------------
HCLO_milo=Milo(HCLOsc)
HCLO_milo
HCLO_milo=buildGraph(HCLO_milo, k=10, d=30)
HCLO_milo <- makeNhoods(HCLO_milo, prop = 0.1, k = 10, d=30, refined = TRUE)
HCLO_milo <- buildNhoodGraph(HCLO_milo)

plotNhoodSizeHist(HCLO_milo)

p4 <- plotNhoodGraph(ad, size_range=c(0.1,3), node_stroke=0.1) + 
  scale_fill_viridis(name = "Nhood size", option = "viridis", direction = 1) 

#scrabbitr------------------

#save------------------
saveRDS(HCLO_milo, file = "./atlas_ad_milo.rds")

#create ortholog file--------
a=rownames(HCLO)
b=rownames(atlas)
c=as.data.frame(intersect(a,b))
names(a)[names(a) == "rownames(hep)"] <- "ref"
names(a)[names(a) == "rownames(atlas)"] <- "query"
rownames(c) = c[,1]

#Run pipeline---------------
out <- scrabbitr::calcNhoodSim(HCLO, atlas, c, sim_preprocessing="gene_spec", sim_measure="pearson",
                               hvg_join_type="intersection", max_hvgs=2000,
                               export_dir = "./milo/1/", 
                               verbose = TRUE)

# Load exported results------
nhood_sim <- as.matrix(fread("./milo/1/nhood_sim.tsv", sep="\t"), rownames=1)
r_vals <- fread("./milo/1/r_vals.tsv", sep="\t")
m_vals <- fread("./milo/1/m_vals.tsv", sep="\t")
out <- list(r_vals = r_vals, m_vals = m_vals, nhood_sim = nhood_sim)

nhood_sim[1:5,1:5]

# Calculate maximum correlations-----  
r_maxNhoods <- getMaxMappings(out$nhood_sim, 1, long_format=FALSE) # rabbit-mouse
m_maxNhoods <- getMaxMappings(out$nhood_sim, 2, long_format=FALSE) # mouse-rabbit
df_simFilt <- rbind(r_maxNhoods, m_maxNhoods)
options(repr.plot.width = 18, repr.plot.height = 8, repr.plot.res = 300)
p1 <- plotNhoodMaxSim(HCLO, r_maxNhoods) + scale_colour_viridis(option="magma")
p2 <- plotNhoodMaxSim(atlas, m_maxNhoods)
grid.arrange(p1,p2,nrow=1)

