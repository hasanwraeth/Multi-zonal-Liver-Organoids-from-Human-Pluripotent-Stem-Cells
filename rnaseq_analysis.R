#Day 28 -------------
#Load Packages #Enhanced Volcano

#library(readr)
#library(dplyr) has conflict with tidyverse
#library(magrittr)
library(BiocManager)
library(tximport)
library(DESeq2)
library(plotly)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

##NOT data.frame BUT tibble
#sample_table = read_csv("SraRunTable2.txt") %>% #View() #pipe symbol= %>%
#  select(`Sample Name`, source_name, Treatment,
#         Cell_Line, Cell_type, Time_point) %>%
#  slice(seq(1, 48, by=4)) #%>% View()
#, df[!duplicated(df$column),]
#install('tximport')

sample_table = read_csv("SraRunTable2.txt") %>% #select is dplyr function
  dplyr::select(`Sample Name`,
         source_name,
         Treatment,
         Cell_Line,
         Cell_type,
         Time_point) %>%
  unique()

'GSM4432378/quant.sf'
pull(sample_table, `Sample Name`)

loc = paste0(pull(sample_table, `Sample Name`), '/quant.sf')
names(loc) = pull(sample_table, `Sample Name`)
gene_map = read_csv("gene_map.csv", col_names = c('esntid', 'ensgid'))

count_data = tximport(
  files = loc,
  type = "salmon",
  tx2gene = gene_map,
  ignoreTxVersion = TRUE
)

#count_data['counts']
count_data$counts[1:6,]

#install('DESeq2')

#Day 30 -------------

sample_table2 = as.data.frame(sample_table)
colnames(sample_table2)[1] = "Sample"
conditions = c('control', 'z1', 'z3', 'mz')
conditions = rep(conditions, each = 3)
conditions = factor(conditions)
sample_table2$conditions = conditions

#y~x (y gene expression given x expt. design)
deseq_dataset = DESeqDataSetFromTximport(txi = count_data,
                                         colData = sample_table2,
                                         design =  ~ conditions)

counts(deseq_dataset)[1:6, 1:3]
count_data$counts[1:6, 1:3]

#normalization = geometric mean rows -> counts/gm of rows ->
#median of columns (size factor) -> raw counts/size factos
#normalization Factors (per gene) > size Factors (per column)
deseq_dataset = estimateSizeFactors(deseq_dataset)
normalizationFactors(deseq_dataset)
counts(deseq_dataset, normalized = TRUE) [1:6, 1:3]
boxplot(counts(deseq_dataset, normalized = TRUE))

#vst = ln transformed
vst = varianceStabilizingTransformation(deseq_dataset)
boxplot(assay(vst))

#observe cell line effect
plotPCA(vst, intgroup = 'conditions') +
  theme_bw()

#split expt. in two
dds1 = deseq_dataset[, 1:6]
dds2 = deseq_dataset[, 7:12]

#PCA for 2 grps
dds1 = estimateSizeFactors(dds1)
normalizationFactors(dds1)
vst1 = varianceStabilizingTransformation(dds1)
plotPCA(vst1, intgroup = 'conditions') +
  theme_bw()

dds2 = estimateSizeFactors(dds2)
normalizationFactors(dds2)
vst2 = varianceStabilizingTransformation(dds2)
plotPCA(vst2, intgroup = 'conditions') +
  theme_bw()

#Hierarchical clustering, euclidean distance
d = assay(vst1)
d = t(d)
d = dist(d)
h = hclust(d)
plot(h)

k = kmeans(t(assay(vst1)), centers = 2)
k$cluster


#Day 32 -----------------
#3 steps in DESeq2 analysis
#1) Estimate size factors (normalization)
#2) Estimate dispersion
#3) Apply statistics (Wald Test)

#dds1 = estimateDispersions(dds1)
#Doesn't work, so split experiment before importing data

# loc = paste0(pull(sample_table, `Sample Name`), '/quant.sf')
# names(loc) = pull(sample_table, `Sample Name`)
loc_nhbe = loc [1:6]
count_data_nhbe = tximport(
  files = loc_nhbe,
  type = "salmon",
  tx2gene = gene_map,
  ignoreTxVersion = TRUE
)
sample_table_nhbe = sample_table2[1:6, ]
sample_table_nhbe$conditions = factor(rep(c('control', 'treated'), each =
                                            3),
                                      levels = c('control', 'treated'))
dds_nhbe = DESeqDataSetFromTximport(txi = count_data_nhbe,
                                    colData = sample_table_nhbe,
                                    design =  ~ conditions)
dds_nhbe = estimateSizeFactors(dds_nhbe)
dds_nhbe = estimateDispersions(dds_nhbe)
plotDispEsts(dds_nhbe)

#Step 3
dds_nhbe = nbinomWaldTest(dds_nhbe)

##DESeq2 shortcut -----------------
#dds_nhbe = DESeq(dds_nhbe) does all 3 steps

result_nhbe = results(dds_nhbe)
#summary(1:100)
summary(result_nhbe)
View(as.data.frame(result_nhbe))

#result_table os a DataFrame not a data.frame!
result_nhbedf = as.data.frame(result_nhbe)

plotCounts(dds_nhbe, gene = 'ENSG00000000005', intgroup = 'conditions')
sum(complete.cases(result_nhbedf))

filter_df1 = result_nhbedf[complete.cases(result_nhbedf),]

#Filter results
#padj < 0.05
#log2FoldChange >1 or <-1

filter_df1$padj < 0.05
filter_df1.2 = filter_df1[filter_df1$padj < 0.05, ]
filter_df1.3 = filter_df1.2[abs(filter_df1.2$log2FoldChange) > 1, ]

plotMA(result_nhbe) #fix if padj not sorted plotMA(res[order(res$padj),])
#& = and | = or
filter_df1$test = filter_df1$padj < 0.05 &
  abs(filter_df1$log2FoldChange) > 1
filter_df1 = rownames_to_column(filter_df1, var = 'ensgene')

g = ggplot(filter_df1, aes(x = log2FoldChange, y = -log10(padj), name = ensgene )) +
  geom_point(aes(color=test), size = 1, alpha = 0.3) + 
  geom_vline(xintercept = 1, color = 'green', linetype=3) +
  geom_vline(xintercept = -1, color = 'green', linetype=3) +
  geom_hline(yintercept = -log10(0.05), color = 'blue', linetype=3)+
  #xlim(-3,3) + ylim(0,10) + #drops elements as scale is small
  scale_color_manual(values = c('black', 'red')) + theme_bw() +
  theme(legend.position = 'top')
# drop column filter_df1 = subset(filter_df1, select = -c(ID, BMI)) ---------
#install('plotly') #for interactive plots 
#ggplotly(g)

listMarts()
ensembl = useEnsembl(biomart = 'ensembl')
View(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
View(listAttributes(ensembl))
View(listFilters(ensembl))

getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version',
                     'ensembl_transcript_id', 'ensembl_transcript_id_version', 
                     'external_gene_name'), filters = c('ensembl_gene_id'),
      values = filter_df1$ensgene, mart = ensembl)

annotation = getBM(attributes = c('ensembl_gene_id', 'chromosome_name',
                                  'start_position', 'end_position', 'strand',
                                  'gene_biotype', 'description',
                                  'external_gene_name'), 
                   filters = c('ensembl_gene_id'),
                   values = filter_df1$ensgene, mart = ensembl)

annotated_df = left_join(filter_df1, annotation, 
                         by=c('ensgene' = 'ensembl_gene_id'))

g = ggplot(annotated_df, aes(x = log2FoldChange, y = -log10(padj), name = external_gene_name )) +
  geom_point(aes(color=test), size = 1, alpha = 0.3) + 
  geom_vline(xintercept = 1, color = 'green', linetype=3) +
  geom_vline(xintercept = -1, color = 'green', linetype=3) +
  geom_hline(yintercept = -log10(0.05), color = 'blue', linetype=3)+
  #xlim(-3,3) + ylim(0,10) + #drops elements as scale is small
  scale_color_manual(values = c('black', 'red')) + theme_bw() +
  theme(legend.position = 'top')
#ggplotly(g)

annotated_df1.2 = annotated_df[annotated_df$padj < 0.05, ]
annotated_df1.2 = annotated_df1.2[abs(annotated_df1.2$log2FoldChange) > 1, ]
degs = annotated_df1.2$ensgene

vst_nhbe = varianceStabilizingTransformation(dds_nhbe)
vst_nhbe_mat = assay(vst_nhbe)

hm_data = vst_nhbe_mat[degs,]
rownames(hm_data) = annotated_df1.2$external_gene_name
heatmap(hm_data)
pheatmap(hm_data, fontsize_row = 8, scale = 'row', cutree_cols = 2, 
         cutree_rows = 2) #, color = blue)
#display.brewer.all()
#blue = colorRampPalette(rev(brewer.pal(n = 12, name = "Blues")))(100)

#Hypergeometric test/ Fisher's exact test
#install('clusterProfiler')
#install('org.Hs.eg.db')
#GO enrichment, Kegg ------------------
ent_gene = getBM(attributes = c('entrezgene_id'), 
                   filters = c('ensembl_gene_id'),
                   values = annotated_df1.2$ensgene, mart = ensembl)
ent_gene =ent_gene$entrezgene_id
ent_gene = as.character(ent_gene)
ent_uni = getBM(attributes = c('entrezgene_id'), 
                 filters = c('ensembl_gene_id'),
                 values = annotated_df$ensgene, mart = ensembl)
ent_uni =ent_uni$entrezgene_id
ent_uni = as.character(ent_uni)
ego = enrichGO(gene = ent_gene, OrgDb = org.Hs.eg.db, 
               ont = "BP", universe = ent_uni)
View(summary(ego))
barplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)
cnetplot (ego, showCategory = 20)
goplot(ego, showCategory = 20)
keg = enrichKEGG(gene = ent_gene, universe = ent_uni)
View(summary(keg))

#Subsetting table based on key
#list = read_csv("Book2.csv")%>% #select is dplyr function
#  dplyr::select(`Sample Name`)
#tmp.subset = sample_table[sample_table$`Sample Name` %in% list$`Sample Name`,]

write_tsv(annotated_df1.2, "filtered_nhbe_results.txt")
