# SIGNAC INSTALL SUCKS. USE A DOCKER IMAGE OR .SIF file
# How to initiate an interactive R session in docker:
#docker run --rm -it -v /home/bizon/Desktop/projects/darpa_project/data/atac_asap_test/atac1/outs/:/my_folder  timoast/signac:latest 
# How to get it going in singulatiry as a standalone script:
#singularity exec --bind /hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac signac_latest.sif Rscript process_ATAC.R > proc_asap303.out

# NOTE, installing these biocmanager things will ask for updates...
# DO NOT update any of these things. Takes forever

library(Signac)
library(Seurat)
#BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
BiocManager::install(c('EnsDb.Mmusculus.v79'))
BiocManager::install("biovizBase")
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(SeuratDisk)
# Budha had this in here, I'm not quite sure why...
# might end up breaking the script who knows.
library(biovizBase)

# Currently all of our data is in this directory
#setwd('/home/bizon/Desktop/projects/darpa_project/data/atac_asap_test/atac1/outs')

# Steal the code line-by-line from Signac
# When running in docker, we can move stuff just to "my_folder"
counts <- Read10X_h5(filename = "out_y303/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "out_y303/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
# Need to check to see what minCells and minFeatures is doing
# Not sure what "min cells" and "min features" actually does. But apparently it is too stringent for
# whatever data I am dealing with....
# They suggest 200 minimum features, I had to set it to 1...
# And we STILL only get 3k cells, well short of the ASAPseq
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = 'out_y303/outs/fragments.tsv.gz',
  min.cells = 10,
  # Budha set this to 200 so lets keep that.
  min.features = 200
)

data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# extract gene annotations from EnsDb
# Moderately slow, especially when we end up getting larger files...
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# Budha commented this out, lets do what he does. He knows R
#genome(annotations) <- "mm10"

Annotation(data) <- annotations

## Evidently these are QC metrics
# compute nucleosome signal score per cell
data <- NucleosomeSignal(object = data)

# compute TSS enrichment score per cell
data <- TSSEnrichment(object = data, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100
data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

# Apparently in R you save figures by pre-opening the pdf...
pdf(file = "out_y303/outs/qc1.pdf", width = 8,height = 4)
DensityScatter(data, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
# dev.off does the printing, I guess
dev.off()

data$high.tss <- ifelse(data$TSS.enrichment > 3, 'High', 'Low')
pdf(file = "out_y303/outs/qc2.pdf", width = 8,height = 4)
TSSPlot(data, group.by = 'high.tss') + NoLegend()
# dev.off does the printing, I guess
dev.off()

# FOR SOME REASON, this is failing. I think it might have something to do with the nucleosome signal being 0 across all data?
# Not sure exactly what that means for this...
#data$nucleosome_group <- ifelse(data$nucleosome_signal > 1, 'NS > 1', 'NS < 1')
#pdf(file = "out_y303/outs/qc3.pdf", width = 8,height = 4)
#FragmentHistogram(object = data, group.by = 'nucleosome_group')
# dev.off does the printing, I guess
#dev.off()

pdf(file = "out_y303/outs/qc3.pdf", width = 16,height = 4)
VlnPlot(
  object = data,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

# So this is how we're supposed to QC, but it seems like with this pilot data (non DARPA), we won't be getting ANY of these
# to work (nowhere near the expected number of peaks)
data <- subset(
  x = data,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
#data

## This final normalization step might be where we then output our files.
data <- RunTFIDF(data)

##########################################################################################
# So this is a function Budha grabbed, but apparently it pulls out what we need.

#Get the findRegion function from Signac repo:
#Function taken from here: https://github.com/stuart-lab/signac/blob/master/R/utilities.R
FindRegion <- function(
  object,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}

###########################################################################################
# All Budha's code for running things...

mat <- GetAssayData(object = data, assay = "peaks", slot = "data")
mat <- as.matrix(mat)

mode <- 'Multi' # Set to Single or Multi, depending on the requirement of analysis

regions <- row.names(mat)
gene_info <- rep("NA", length(regions))
gene_meta <- data.frame(regions, gene_info)

# Certain patches of the code in the loop is taken from the code of the Signac function AnnotatioPlot
# I will just trust that this works... sure Budha...

for (k in 1:length(regions)) {
  
  # Might want to look at RNA instead of "peaks" for this script.
  # run both to see what works... Compare output differences.
  region <- FindRegion(object = data, region = regions[k], sep = c("-", "-"), 
        assay = "peaks", extend.upstream = 0, extend.downstream = 0)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
    
  annotation <- Annotation(object = data)
    
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  genes.keep <- unique(x = annotation.subset$gene_name)
  
  if (length(genes.keep)>0) {
    if (length(genes.keep)==1) {
      gene_meta$gene_info[k] <- genes.keep
    }
    if (length(genes.keep)>1) {
      if (mode == 'Single') {
      gene_meta$gene_info[k] <- genes.keep[1]
    }
    
    if (mode == "Multi") {
      gene_meta$gene_info[k] <- paste(genes.keep, collapse = ", ")
      print(k)
    }
    }
    
  }
}

write.csv(mat, "mat.csv")
write.csv(gene_meta, "gene_meta.csv")

SaveH5Seurat(data, 'out_y303/outs/atac_process303.h5', overwrite = FALSE, verbose = TRUE)

# This conversion step *might* be where I'm running into an issue. Basically for whatever reason,
# the H5Seurat file above cannot be read in to python. So I think the convert step is necessary...
# Maybe the annotations could be saved as a CSV?

Convert("out_y303/outs/atac_process303.h5.h5seurat",dest="h5ad")

###########################################################################################
# EVERYTHING BELOW HERE RUN JUST THIS FIRST TIME. COMMENT ALL OF IT OUT LATER
# WE WANT THIS DONE IN PYTHON, NOT R...
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunSVD(data)

# This is a plot
#DepthCor(pbmc)

data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3)

pdf(file = "out_y303/outs/umap.pdf", width = 8,height = 8)
DimPlot(object = data, label = TRUE) + NoLegend()
dev.off()

gene.activities <- GeneActivity(data)

data[['RNA']] <- CreateAssayObject(counts = gene.activities)
data <- NormalizeData(
  object = data,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)
DefaultAssay(data) <- 'RNA'

pdf(file = "out_y303/outs/expression.pdf", width = 12,height = 8)
FeaturePlot(
  object = data,
  features = c('Ms4a1', 'Cd3d', 'Lef1', 'Nkg7', 'Trem1', 'Cd4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()

DefaultAssay(data) <- 'peaks'
# The identifiers used to find differentially expressed peaks
# are from some anchor-finding analysis they do with matching RNAseq
#da_peaks <- FindMarkers(
#  object = data,
#  ident.1 = "CD4 Naive",
#  ident.2 = "CD14+ Monocytes",
#  test.use = 'LR',
#  latent.vars = 'nCount_peaks'
#)

pdf(file = "out_y303/outs/coverage.pdf", width = 16,height = 8)
CoveragePlot(
  object = data,
  region = "chr2-87011729-87035519",
  annotation = FALSE,
  peaks = FALSE
)
dev.off()

pdf(file = "out_y303/outs/coverage2.pdf", width = 16,height = 8)
CoveragePlot(
  object = data,
  region = "Cd3d",
  annotation = FALSE,
  peaks = FALSE
)
dev.off()

SaveH5Seurat(data, 'out_y303/outs/atac_process303.h5', overwrite = FALSE, verbose = TRUE)

# This conversion step *might* be where I'm running into an issue. Basically for whatever reason,
# the H5Seurat file above cannot be read in to python. So I think the convert step is necessary...
# Maybe the annotations could be saved as a CSV?

Convert("out_y303/outs/atac_process303.h5.h5seurat",dest="h5ad")