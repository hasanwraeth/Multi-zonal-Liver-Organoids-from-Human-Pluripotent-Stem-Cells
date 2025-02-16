setwd("C:/Users/Hasan/Desktop/Smith/bw")
setwd("C:/Users/Hasan/Desktop/Smith/macs2/SRR1642058/")

library("tidyverse")
library("DiffBind")
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v105)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(goseq)
library(rGREAT)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)


(bams <- dir("bam", "markDup.bam$", full.names = TRUE))
(Peaks <- dir("macs2", "narrowPeak$", full.names = TRUE, recursive = TRUE))
(SampleID <- basename(dirname(Peaks)))
(Condition <- sub(".rep.", "", SampleID))
(Replicate <- sub("^.*rep", "", SampleID))
Peakcaller <- "macs2"
PeakFormat <- "narrowPeak"
samples <- data.frame(SampleID=SampleID,
                      Condition=Condition,
                      Replicate=Replicate,
                      bamReads=bams,
                      Peaks=Peaks,
                      Peakcaller=Peakcaller,
                      PeakFormat=PeakFormat,
                      ScoreCol=5)

pf <- "DiffBind"
out <- "sample.csv"
dir.create(pf)
write.csv(samples, file.path(pf, out))
chip <- dba(sampleSheet = file.path(pf, out))
pdf(file.path(pf, "DiffBind.sample.correlation.pdf"), width = 9, height = 9)
plot(chip)
dev.off()
pdf(file.path(pf, "DiffBind.PCA.plot.pdf"))
dba.plotPCA(chip, DBA_CONDITION, label=DBA_ID)
dev.off()

chip <- dba.count(chip)
chip <- dba.normalize(chip) ## add for DiffBind 3.0
BLACKLIST <- FALSE ## fish data, no blacklist yet.
chip <- dba.blacklist(chip, blacklist=BLACKLIST, greylist=FALSE)
chip <- dba.contrast(chip, minMembers=1, categories = DBA_CONDITION)
chip <- dba.analyze(chip, bBlacklist = FALSE, bGreylist = FALSE)
chip.DB <- dba.report(chip, th=4)
head(chip.DB)
write.csv(chip.DB, file.path(pf, "DiffBind.results.csv"), row.names=FALSE)

#/ A minimal example on how to use EnrichedHeatmap together with rtracklayer.
#/ Required data are at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111902
#/ Go for the files :
#/ "GSM3045250_N_ATAC_Kaech_1_peaks.broadPeak.gz" and (genomic regions in BED-like format)
#/ "GSM3045250_N_ATAC_Kaech_1.bw"                     (the bigwig files with read counts)

require(EnrichedHeatmap)
require(rtracklayer)
require(circlize)
require(data.table)

#/ Load a BED file into R with data.table::fread(),
#/ then convert to GRanges format with GenomicRanges::GRanges():
targets <- makeGRangesFromDataFrame(
  df = fread("SRR1642058_peaks.narrowPeak", header = FALSE, data.table = FALSE),
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

#/ For this tutorial we take the first 10000 regions rather than the full file:
targets <- head(targets, 10000)

#/ We take the center of each region/peak and want to extend by 5kb each direction:
ExtendSize <- 2500
targets.extended  <- resize(targets, fix = "center", width = ExtendSize*2)

#/ We load the relevant parts of the bigwig file into R. 
#/ This is more efficient than reading the entire file. 
BigWig <- rtracklayer::import("SRR1642058.bw", 
                              format = "BigWig", 
                              selection = BigWigSelection(targets.extended))

#/ Create the normalizedMatrix that EnrichedHeatmap accepts as input.
#/ We use the targets center (width=1) because from what I understand normalizeMatrix
#/ does not allow to turn off its "extend" option. Therefore we trick it by simply
#/ providing the peak centers and then let the function extend it by our predefined window size.
normMatrix <- normalizeToMatrix(signal = BigWig, 
                                target = resize(targets, fix = "center", width = 1), 
                                background = 0, 
                                keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                target_ratio = 0,
                                mean_mode = "w0",       #/ see ?EnrichedHeatmap on other options
                                value_column = "score", #/ = the name of the 4th column of the bigwig
                                extend = ExtendSize)

#/ Make a color gradient that covers the range of normMatrix from 0 to the 99th percentile.
#/ The percentile avoids outliers to skew the heatmap:
col_fun = circlize::colorRamp2(quantile(normMatrix, c(0, .99)), c("white", "#005cff"))

#/ heatmap function:
EH <- EnrichedHeatmap( mat = normMatrix,
                       pos_line = FALSE, #/ no dashed lines around the start
                       border = FALSE,   #/ no box around heatmap
                       col = col_fun,    #/ color gradients from above
                       axis_name_gp = gpar(fontsize = 20),
                       column_title = "Smith", #/ column title 
                       column_title_gp = gpar(fontsize = 30, fontfamily = "sans"),
                       use_raster = TRUE, raster_quality = 10, raster_device = "tiff",
                       #/ turn off background colors
                       rect_gp = gpar(col = "transparent"),
                       #/ legend options:
                       heatmap_legend_param = list(labels_gp = gpar(fontsize = 20),
                                                   title_gp= gpar(fontsize = 20),title_position = "topcenter",
                                                   legend_direction = "horizontal",
                                                   title = "normalized counts"),
                       #/ options for the profile plot on top of the heatmap:
                       top_annotation = HeatmapAnnotation(
                         enriched = anno_enriched(ylim=c(0,1),
                                                  gp = gpar(col = "#005cff", lty = 1, lwd=4, ylab="occupancy score"),
                                                  axis_param=list(gp=gpar(size=40)), axis_title="occupancy score",
                                                  col="blue") )
) #/ end of EnrichedHeatmap function

#/ Save as pdf to disk:
#pdf("EnrichedHeatmap.pdf")

draw(EH,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "bottom",     #/ we want the legend below the heatmap
     annotation_legend_side = "bottom",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm") #/ some padding to avoid labels beyond plot borders
)

#dev.off()

macsPeaks = "SRR1642059_peaks.xls"
macsPeaks_DF <- read.delim(macsPeaks)
macsPeaks_DF[1:8, ]

macsPeaks_DF <- read.delim(macsPeaks, comment.char = "#")
macsPeaks_DF[1:2, ]

macsPeaks_GR <- GRanges(seqnames = macsPeaks_DF[, "chr"],
                        IRanges(macsPeaks_DF[,"start"],
                                macsPeaks_DF[, "end"]))
macsPeaks_GR
seqnames(macsPeaks_GR)
ranges(macsPeaks_GR)

mcols(macsPeaks_GR) <- macsPeaks_DF[, c("abs_summit", "fold_enrichment")]
macsPeaks_GR

macsPeaks_GR_np <- import("SRR1642059_peaks.narrowPeak", format = "narrowPeak")
macsPeaks_GR_np

seqlevelsStyle(macsPeaks_GR) <- "UCSC" #changes 1 to chr1

peakAnno <- annotatePeak(macsPeaks_GR, tssRegion = c(-500, 500), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                         annoDb = "org.Hs.eg.db")

class(peakAnno)
peakAnno_GR <- as.GRanges(peakAnno)
peakAnno_DF <- as.data.frame(peakAnno)
peakAnno_GR[1:2, ]
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
upsetplot(peakAnno, vennpie = F)

annotatedPeaksGR_TSS <- peakAnno_GR[peakAnno_GR$annotation == "Promoter",]
genesWithPeakInTSS <- unique(annotatedPeaksGR_TSS$geneId)
genesWithPeakInTSS[1:2]
allGeneGR <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
allGeneGR[1:2, ]
allGeneIDs <- allGeneGR$gene_id

GO_result <- enrichGO(gene = genesWithPeakInTSS, universe = allGeneIDs, 
                      OrgDb = org.Hs.eg.db, ont = "BP")
GO_result_df <- data.frame(GO_result)
GO_result_df[1:5, ]

GO_result_plot <- pairwise_termsim(GO_result)
emapplot(GO_result_plot, showCategory = 20)

msigdbr_collections()
msig_t2g <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)
msig_t2g <- msig_t2g[, colnames(msig_t2g) %in% c("gs_name", "entrez_gene")]
msig_t2g[1:3, ]

hallmark <- enricher(gene = genesWithPeakInTSS, universe = allGeneIDs, TERM2GENE = msig_t2g)
hallmark_df <- data.frame(hallmark)
hallmark_df[1:3, ]

allGenesForGOseq <- as.integer(allGeneIDs %in% genesWithPeakInTSS)
names(allGenesForGOseq) <- allGeneIDs
allGenesForGOseq[1:3]

pwf = nullp(allGenesForGOseq, "hg38", "knownGene", plot.fit = FALSE)
Myc_hallMarks <- goseq(pwf, "hg38", "knownGene", gene2cat = data.frame(msig_t2g))
Myc_hallMarks[1:3, ]

great_Job <- submitGreatJob(macsPeaks_GR, species = "hg38", version = "4.0.4", request_interval = 1)
availableCategories(great_Job)

macsSummits_GR <- GRanges(seqnames(macsPeaks_GR), IRanges(macsPeaks_GR$abs_summit,
                                                          macsPeaks_GR$abs_summit), score = macsPeaks_GR$fold_enrichment)
macsSummits_GR <- resize(macsSummits_GR, 100, fix = "center")
macsSummits_GR

macsSummits_GR2=macsSummits_GR
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(macsSummits_GR2)) %in% main.chroms)
macsSummits_GR2<- macsSummits_GR2[keep.peaks, ]
peaksSequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, macsSummits_GR2)
names(peaksSequences) <- paste0(seqnames(macsSummits_GR2), ":", start(macsSummits_GR2),
                                "-", end(macsSummits_GR2))

peaksSequences[1:2, ]
writeXStringSet(peaksSequences, file = "SRR1642059.fa")

pfm <- getMatrixByName(JASPAR2020, name = " KAT3B")
pfm

#cat myhg19Genes.bed sed 's/^/chr/' > myhg19Genes_new.bed

bed = "myhg19Genes_new.bed"
bed=read.delim(bed, header = F)
seqlevelsStyle(bed) 
