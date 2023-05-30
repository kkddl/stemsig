

library('hdf5r')
library(Seurat)
sce <- Read10X_h5(filename = '~/Glioma_GSE103224_expression.h5')
sce = CreateSeuratObject(sce,project = "GSE103224")
sce

library(readr)
metadata <- read_delim("~/Glioma_GSE103224_CellMetainfo_table.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

row.names(metadata) <- metadata$Cell


###scCancer mRNAsi
stem.sig.file <- system.file("txt", "pcbc-stemsig.tsv", package = "scCancer")
stem.sig <- read.delim(stem.sig.file, header = FALSE, row.names = 1)

X = GetAssayData(sce,'counts')
common.genes <- intersect(rownames(stem.sig), rownames(X))
X <- X[common.genes, ]
stem.sig <- stem.sig[common.genes, ]

s <- apply(X, 2, function(z) {cor(z, stem.sig, method = "sp", use = "complete.obs")})
names(s) <- colnames(X)

s <- s - min(s)
s <- s / max(s)
s = as.data.frame(s)

sce <- AddMetaData(object = sce, metadata = s)

###
meta2 = metadata[,c(1,4,5,6)]
rownames(meta2) = meta2$Cell
sce <- AddMetaData(object = sce, metadata = meta2)

#
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
D = GetAssayData(object = sce, slot = "counts")
sce$log10GenesPerUMI <- log10(sce$nFeature_RNA)/log10(sce$nCount_RNA)
##"nFeature_RNA", "nCount_RNA", "percent.mt"
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.01)
sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

##normalizing the data 
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection) 
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sce), 10)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)

# 
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
#
print(sce[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sce, dims = 1:2, reduction = "pca")
DimPlot(sce, reduction = "pca",label = T)
DimHeatmap(sce, dims = 1:12, cells = 500, balanced = TRUE)

sce <- JackStraw(sce, num.replicate = 100)
sce <- ScoreJackStraw(sce, dims = 1:10)
JackStrawPlot(sce, dims = 1:10)
ElbowPlot(sce)



sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.5)
head(Idents(sce), 5)
#
set.seed(123)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
DimPlot(sce,reduction = "tsne",label=T)
sce <- RunUMAP(object = sce, dims = 1:15, do.fast = TRUE)
DimPlot(sce, reduction = "umap",label = T)

#plot metadata
DimPlot(sce, group.by="Celltype..major.lineage.", label=T, label.size=5, reduction='tsne')
DimPlot(sce, group.by="orig.ident", label=T, label.size=5, reduction='tsne')
DimPlot(sce, group.by="Celltype..malignancy.", label=T, label.size=5, reduction='tsne')

library(ggplot2)
library(dplyr)
mat1=s
colnames(mat1)="score"
mat2=Embeddings(sce,"tsne")
mat3=merge(mat2,mat1,by="row.names")
head(mat3)
mat3%>%ggplot(aes(tSNE_1,tSNE_2))+geom_point(aes(color=score))+
  scale_color_gradient(low = "grey",high = "purple")+theme_classic()

ggsave("CCR7.2.pdf",device = "pdf",width = 13.5,height = 12,units = "cm")


sce@meta.data$malig_stem = ifelse(sce@meta.data$Celltype..malignancy. == 'Malignant cells',
                                  ifelse(sce@meta.data$s > 0.5, 'malig-high', 'malig-low'), 'non-malig')

kids <- purrr::set_names(levels(factor(sce@meta.data$malig_stem)))
length(kids)

sids <- purrr::set_names(levels(factor(sce@meta.data$orig.ident)))
length(sids)




table(sce@meta.data$orig.ident)

n_cells <- as.numeric(table(sce@meta.data$orig.ident))
m <- match(sids, sce@meta.data$orig.ident)
ei <- data.frame(sce@meta.data[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select(-"malig_stem")


groups <- sce@meta.data[, c("malig_stem","orig.ident")]
library(Matrix.utils)
pb <- aggregate.Matrix(t(sce@assays$RNA@counts), groupings = groups, fun = "sum")
pb[1:6, 1:6]
library(Matrix.utils)
library(magrittr)
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)
str(pb)
options(width = 100)
table(sce@meta.data$malig_stem, sce@meta.data$orig.ident)
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}
de_samples <- purrr::map(1:length(kids), get_sample_ids) %>%
  unlist()

samples_list <- purrr::map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}
library(purrr)
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

gg_df <- data.frame(malig_stem = de_cluster_ids,
                    orig.ident = de_samples)

metadata <- gg_df %>%
  dplyr::select(malig_stem, orig.ident) 
metadata    
metadata$malignancy = ifelse(substr(metadata$malig_stem,1,2) == 'ma', 'malig','non-malig')
metadata$score = ifelse(substr(metadata$malig_stem,7,8) == 'hi','high',
                        ifelse(substr(metadata$malig_stem,7,8) == 'lo','low','not'))
clusters <- levels(factor(metadata$malignancy))
clusters
clusters[1]

cluster_metadata <- metadata[which(metadata$malignancy == clusters[1]), ]
head(cluster_metadata)

rownames(cluster_metadata) <- cluster_metadata$orig.ident
head(cluster_metadata)

counts <- data.frame(pb[[1]],pb[[2]])


library(DESeq2)
metadata$grp = paste(metadata$orig.ident,rep(c(0,1,2),each = 8),sep = '.')
cluster_metadata = metadata[1:16,]
dds <- DESeqDataSetFromMatrix(round(counts), 
                              colData = cluster_metadata, 
                              design = ~ score)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA

DESeq2::plotPCA(rld, intgroup = "orig.ident")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
library(pheatmap)
pheatmap(rld_cor, annotation = cluster_metadata[,'grp'])

pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(cluster_metadata, pca$x)
ggplot(df) + geom_point(aes(x=PC6, y=PC7, color = orig.ident))

##
dds <- DESeq(dds)

plotDispEsts(dds)
##Results
# Output results of Wald test for contrast for stim vs ctrl
levels(factor(cluster_metadata$score))[2]
levels(factor(cluster_metadata$score))[1]

contrast <- c("score", levels(factor(cluster_metadata$score))[1], levels(factor(cluster_metadata$score))[2]) 

# resultsNames(dds)
relevel(dds$score,ref = 'low')
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)
resultsNames(dds)

res <- lfcShrink(dds, 
                 contrast =  contrast,
                 res=res,type = 'ashr')
library(tibble)
# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  as_tibble()
# Check results output
res_tbl
padj_cutoff <- 0.05

sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

sig_res

write.csv(sig_res,
          paste0("results", clusters[1], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
library(readr)
sig_res <- read_csv("resultsmalig__vs__sig_genes.csv")

gene = modProbes[modProbes %in% sig_res$gene]


####220913
library(EnhancedVolcano)
EnhancedVolcano(sig_res,
          lab = sig_res$gene,
             x = 'log2FoldChange',
          y = 'padj',pCutoff = 5e-02,)

genelist = c('PECAM1','VWF', 'ENG','AIF1', 'C5AR1','CD68','GFAP', 'AQP4', 'ALDOC','OLIG1', 'OLIG2', 'PDGFRA', 'DLL3','CD24','STMN2','DCX'
)
# astrocytes (GFAP, AQP4, ALDOC) and OPCs (OLIG1, OLIG2, PDGFRA, DLL3),  neuroblasts like CD24 and STMN2,DCX. 
DotPlot(sce, features = genelist)


##cellchat
{
Sys.setenv(RETICULATE_PYTHON=) 
sce@commands$FindClusters


sce@meta.data$cellchat = ifelse(sce@meta.data$malig_stem == 'malig-high', 'malig-high',
                                  ifelse(sce@meta.data$malig_stem == 'malig-low', 'malig-low', 
                                         ifelse(sce@meta.data$Celltype..major.lineage. == 'Neuron', 'neuron',
                                                ifelse(sce@meta.data$Celltype..major.lineage. == 'Endothelial', 'endothelial','mono/macro'))))
                                
data.input = sce@assays$RNA@data
identity = data.frame(group = sce$cellchat,row.names = names(sce$Cell))
unique(identity$group) # check the cell labels

#
cellchat = createCellChat(object = data.input, meta = sce@meta.data,group.by = 'cellchat')
summary(cellchat)
levels(cellchat@idents)
str(cellchat)
groupSize = as.numeric(table(cellchat@idents))

CellChatDB = CellChatDB.human
str(CellChatDB)

colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation) 
CellChatDB.use = CellChatDB 
cellchat@DB = CellChatDB.use 


cellchat = subsetData(cellchat)
future::plan('multiprocess',workers=4)

cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)

cellchat = projectData(cellchat, PPI.human)

cellchat = computeCommunProb(cellchat, raw.use = F, population.size = T) 
#filter out communication with few cells in certain cell groups
cellchat = filterCommunication(cellchat, min.cells = 10)
df.net = subsetCommunication(cellchat)
write.csv(df.net,'net_lr.csv')


## @netP
cellchat = computeCommunProbPathway(cellchat)
df.netp = subsetCommunication(cellchat, slot.name = 'netP')
write.csv(df.netp,'net_pathway.csv')


##
cellchat = aggregateNet(cellchat)
#
groupSize = as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = T)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = 'Number of interactions')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = 'Interaction weights/strength')

#number of interaction
mat = cellchat@net$count
par(mfrow = c(3,3), xpd = T)
for(i in 1:nrow(mat)) {
  # i=1
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2, 
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i]
                   )
}
##
mat = cellchat@net$weight
par(mfrow = c(3,3), xpd = T)
for(i in 1:nrow(mat)){
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i,]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2, arrow.size = 0.1,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}



#
cellchat@netP$pathways #
pathways.show <- c("TGFb")  

#hierarchy plot 
levels(cellchat@idents) #show all celltype
vertex.receiver = c(4,5) # a numeric vector giving the index of the celltype as targets
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = 'hierarchy')

#circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = 'circle')

#chord plot 

#heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = 'Reds')


#
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TGFb = extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F) 

#hierarchy plot
LR.show = pairLR.TGFb[1,]
vertex.receiver = c(4,5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver, layout = 'hierarchy')

#circle plot
netVisual_individual(cellchat,  signaling = pathways.show, pairLR.use = LR.show, layout = 'circle')
netVisual_individual(cellchat,  signaling = pathways.show, pairLR.use = LR.show, layout = 'chord')



pathways.show.all = cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(4,5) 
dir.create('all_pathways_com_circle')
setwd('all_pathways_com_circle')

for( i in 1:length(pathways.show.all)) {
  # visualize communication network associated with both signaling pathway and individual LR pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c('pdf'),
            vertex.receiver = vertex.receiver, layout = 'circle')
  #compute and visualize the contribution of each ligand receptor pair to the overall signaling pathway
  gg = netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename = paste0(pathways.show.all[i], '_LR_contri.pdf'),
         plot = gg, width = 5, height = 2.5, units = 'in', dpi = 300)}
setwd('../')
  
  
  
##
  #
  levels(cellchat@idents)
  #show all the significant interactions (LR pairs)
  
  p = netVisual_bubble(cellchat, sources.use = c(2,3), targets.use = c(1,4), remove.isolate = F)
  ggsave('malig-high-low.pdf', p, width = 8, height = 12)

  #
  netVisual_bubble(cellchat, sources.use = c(2,3), targets.use = c(1,4),
                   signaling = c("CCL","CXCL"), remove.isolate = F)
  
  #
  pairLR.use = extractEnrichedLR(cellchat, signaling = c('CCL','CXCL','TGFb'))
  netVisual_bubble(cellchat, sources.use = c(2,3), targets.use = c(1,4),
                   pairLR.use = pairLR.use, remove.isolate = T)

}