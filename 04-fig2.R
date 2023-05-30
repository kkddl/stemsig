# Fig3 
KM
library(survival)
library(survminer)
mrnasi_dat2$sigroup = ifelse(mrnasi_dat2$score>median(mrnasi_dat2$score),'H','L')
mrnasi_dat2$OS = as.numeric(mrnasi_dat2$OS)
fit<-survfit(Surv(OS.time,OS) ~ sigroup, data=mrnasi_dat2)
KMsurvival_plot<-ggsurvplot(fit,pval = TRUE, 
                            risk.table = "abs_pct", 
                            risk.table.y.text.col = T,
                            risk.table.y.text = FALSE,
                            xlab = "Time",   
                            surv.median.line = "hv", 
                            ncensor.plot = FALSE, 
                            palette = c('#e41a1c','#377eb8'), 
                            ggtheme = theme_classic()
)
#F2a
KMsurvival_plot
lgg_pfi = na.omit(lgg_pfi)
colnames(lgg_pfi)[1] = 'Tumor_Sample_Barcode'
pfi = merge(mrnasi_dat2,lgg_pfi,by='Tumor_Sample_Barcode',all = F)
pfi$PFI = as.factor(pfi$PFI)
library(cmprsk)
library(splines)
attach(pfi)
crmod <- cuminc(PFI.time,PFI,pfi$sigroup)
crmod

#F2b
f2b = plot(crmod,xlab = 'Time', ylab = 'Survival',lwd=2,lty=1,
     col = c('#e41a1c','#377eb8','#4daf4a','#984ea3'))
{
  pfi$PFI[pfi$PFI == '2'] ='1'
  pfi$PFI = as.numeric(pfi$PFI)
  fit<-survfit(Surv(PFI.time,PFI) ~ sigroup, data=pfi)
  f2b<-ggsurvplot(fit,pval = TRUE,
                              risk.table = "abs_pct",  
                              risk.table.y.text.col = T,
                              risk.table.y.text = FALSE,
                              xlab = "Time",  
                              surv.median.line = "hv", 
                              ncensor.plot = FALSE, 
                              palette = c('#e41a1c','#377eb8'), 
                              ggtheme = theme_classic()
)
                              }
#F2c to do

colnames(tpms) = substr(colnames(tpms),1,12)
f2cinput = tpms[,mrnasi_dat2$Tumor_Sample_Barcode]
f2cinput = as.data.frame(f2cinput)
f2cinput <- f2cinput[rowSums(f2cinput==0)< dim(f2cinput)[2],]
f2cinput <- f2cinput[rowSums(f2cinput>=1)>=200,]

f2c_input = t(scale(t(log2(f2cinput+1))))
f2c_input[f2c_input>2]=2
f2c_input[f2c_input < -2]= -2
f2c_input = as.data.frame(f2c_input)

sigroup<-mrnasi_dat2[,'sigroup'] 
sigroup<-factor(sigroup,levels=c('L','H'))

pvalues <- sapply(1:nrow(f2cinput),function(i){
  data<-cbind.data.frame(gene=as.numeric(f2cinput[i,]),sigroup)
  p=wilcox.test(gene~sigroup, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")

# Calculate fold-change for each gene

conditionsLevel<-levels(sigroup)
dataCon1=f2cinput[,c(which(sigroup==conditionsLevel[1]))]
dataCon2=f2cinput[,c(which(sigroup==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results based on FDR threshold

outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(f2cinput)
outRst=na.omit(outRst)
fdrThres=0.01
input = outRst[outRst$FDR<fdrThres & abs(outRst$log2foldChange)>1.5,]

sigroup = as.data.frame(sigroup)
rownames(sigroup) = mrnasi_dat2$Tumor_Sample_Barcode
library(pheatmap)
f2c<-pheatmap(f2c_input[rownames(input),],
                            color = colorRampPalette(c('#377eb8','#f7f7f7','#e41a1c'))(100), 
                            annotation_col = sigroup, 
                            show_colnames =FALSE, 
                            show_rownames = F, 
                            cluster_cols=FALSE,  
                            cluster_rows=TRUE)  
f2c


#F2d to do 

##KEGG
nrDEG = input
library( "clusterProfiler" )
library( "org.Hs.eg.db" )
nrDEG$SYMBOL <- rownames(nrDEG)
df <- bitr( rownames( nrDEG ), fromType = "SYMBOL", toType = c( "ENTREZID" ), 
            OrgDb = org.Hs.eg.db )
head( df )
nrDEG = merge( nrDEG, df, by = 'SYMBOL' )
head( nrDEG )

gene_up = nrDEG[ nrDEG$log2foldChange > 0 , 'ENTREZID' ] 
gene_down = nrDEG[ nrDEG$log2foldChange < 0 , 'ENTREZID' ]
gene_diff = c( gene_up, gene_down )
gene_all = as.character(nrDEG[ ,'ENTREZID'] )
g_list = list( gene_up = gene_up, gene_down = gene_down, gene_diff = gene_diff)

kk.up <- enrichKEGG(gene = gene_up,
                    organism = 'hsa',
                    universe = gene_all,
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01)
kk.dowm <- enrichKEGG(gene = gene_down,
                      organism = 'hsa',
                      universe = gene_all,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01)

kegg_down_dt <- as.data.frame(kk.dowm)
kegg_up_dt <- as.data.frame( kk.up )
down_kegg <- kegg_down_dt[ kegg_down_dt$pvalue < 0.01, ]
down_kegg$group <- 'down_pathway'
up_kegg <- kegg_up_dt[ kegg_up_dt$pvalue < 0.01, ]
up_kegg$group <- 'up_pathway'
dat = rbind(up_kegg,down_kegg)
dat$pvalue = -log10(dat$pvalue)
dat$group =  factor(dat$group)

library(ggpubr)
ggbarplot(dat,x = 'Description',y = 'pvalue',
          fill = 'group',
          color = 'white',
          palette = 'jco',
          sort.val = 'asc',
          xlab = 'Pathway names',
          ylab = '-log10 P-value',
          title = 'Pathway enrichment') +
  rotate() +
  theme_minimal()

##GO
{

  go <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all") 
  barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
  go1 <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="bp") 
  go1 <- simplify(go1)
  go2 <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="cc") 
  go2 <- simplify(go2)
  go3 <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="mf") 
  go3 <- simplify(go3)
  library(patchwork);library(tidyverse)
  dotplot(go1) / dotplot(go2)/dotplot(go3)
  
  p = dotplot(go1) + scale_y_discrete(labels=function(x) str_wrap(x, width=20)) 
  p
  p / dotplot(go2)/dotplot(go3)
  }
BP <- enrichGO( gene          =  gene_diff,
                universe      =  gene_all,
                OrgDb         =  org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           =  'BP',
                pAdjustMethod = "BH",
                pvalueCutoff  =  0.05,
                qvalueCutoff  =  0.05,
                readable      =  TRUE)
barplot(BP,showCategory=20)
dotplot(BP,showCategory=20)

CC <- enrichGO( gene          =  gene_diff,
                universe      =  gene_all,
                OrgDb         =  org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           =  'CC',
                pAdjustMethod = "BH",
                pvalueCutoff  =  0.05,
                qvalueCutoff  =  0.05,
                readable      =  TRUE)
barplot(CC,showCategory=20)
dotplot(CC,showCategory=20)

MF <- enrichGO( gene          =  gene_diff,
                universe      =  gene_all,
                OrgDb         =  org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           =  'MF',
                pAdjustMethod = "BH",
                pvalueCutoff  =  0.05,
                qvalueCutoff  =  0.05,
                readable      =  TRUE)
barplot(MF,showCategory=20) 
dotplot(MF,showCategory=20) +
  scale_x_continuous(limits = c(0,0.08), breaks = c(0.00,0.04,0.08))
#

#F2e 
##oncoplot
maf = data.table::as.data.table(read.csv(file = './rawdata/TCGA.LGG.mutect.somatic.maf',
                                         header = T,sep = '\t',stringsAsFactors = F,comment.char = '#'))
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 12)
maf = maf[maf$Hugo_Symbol %in% nrDEG$SYMBOL]
cliformaf = sigroup
cliformaf$Tumor_Sample_Barcode = rownames(cliformaf)

maf = maf[-which(maf$Tumor_Sample_Barcode %in% c('TCGA-DU-6392')),]
mafh = maf[maf$Tumor_Sample_Barcode %in% rownames(sigroup)[sigroup$sigroup == 'H'],]
mafl = maf[maf$Tumor_Sample_Barcode %in% rownames(sigroup)[sigroup$sigroup == 'L'],]

library(maftools)
maf_lgg=read.maf(maf = maf,clinicalData = cliformaf)
maf_lgg@clinical.data$sigroup = as.factor(maf_lgg@clinical.data$sigroup)
subtypecolors = c("#e41a1c", "#377eb8",'#4daf4a')
names(subtypecolors) = levels(maf_lgg@clinical.data$sigroup)
phecolors = list(sigroup = subtypecolors)


getSampleSummary(maf_lgg)
getClinicalData(maf_lgg)
plotmafSummary(maf=maf_lgg, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
f2e2 = oncoplot(maf=maf_lgg, clinicalFeatures = 'sigroup', bgCol = "#CCCCCC",borderCol = 'white',annotationColor = phecolors,
         annoBorderCol = "#ebebeb", sortByAnnotation = T,removeNonMutated=T)

clin.low <- subset(cliformaf, sigroup=="L")$Tumor_Sample_Barcode
clin.high <- subset(cliformaf, sigroup=="H")$Tumor_Sample_Barcode

lgg.low <- subsetMaf(maf=maf_lgg, tsb=clin.low, isTCGA=TRUE)
lgg.high <- subsetMaf(maf=maf_lgg, tsb=clin.high, isTCGA=TRUE)

comp <- mafCompare(m1=lgg.low, m2=lgg.high, m1Name="low_risk", m2Name="high_risk", minMut=5)
forestPlot(mafCompareRes=comp, pVal=0.05, color=c("royalblue","maroon" ), geneFontSize=0.8)

genes <- comp$results$Hugo_Symbol[1:10]
#F2e
f2e = coOncoplot(m1=lgg.low, m2=lgg.high, m1Name="low_risk", m2Name="high_risk")









#F2f 

library(TCGAbiolinks)
query.lgg.nocnv <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Copy number variation",
  legacy = TRUE,
  file.type = "nocnv_hg19.seg",
  sample.type = c("Primary Tumor")
)

GDCdownload(query.lgg.nocnv)

cnvMatrix <- GDCprepare(query.lgg.nocnv)
head(cnvMatrix)


gdac.root <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/" 
file <- paste0(gdac.root, "genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt")
if(!file.exists(basename(file))) 
  downloader::download(file, file.path("~/Downloads", basename(file))) 

markersMatrix <-  readr::read_tsv(
  file.path("~/Downloads", basename(file)), 
  col_names = FALSE, col_types = "ccn", 
  progress = FALSE)
head(markersMatrix)

library(dplyr)
cnvMatrix %<>%
  filter(abs(Segment_Mean) > 0.3) %>%
  mutate(label = if_else(Segment_Mean < -0.3, 0, 1)) %>%
  dplyr::select(-Segment_Mean) %>%
  `names<-`(c("Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")) %>%
  mutate(Chromosome = ifelse(Chromosome == "X", 23, Chromosome),
         Chromosome = ifelse(Chromosome == "Y", 24, Chromosome),
         Chromosome = as.integer(Chromosome)) %>%
  as.data.frame()
cnvMatrixH = cnvMatrix[substr(cnvMatrix$Sample.Name,1,12) %in% cliformaf[which(cliformaf$sigroup == 'H'),2],]
cnvMatrixL = cnvMatrix[substr(cnvMatrix$Sample.Name,1,12) %in% cliformaf[which(cliformaf$sigroup == 'L'),2],]





commonCNV <- readr::read_tsv(
  file.path("~/Downloads", basename(file)), 
  progress = FALSE
) %>%
  mutate(markerID = paste(Chromosome, Start, sep = ":"))


markersMatrix_filtered <- markersMatrix %>%
  `names<-`(c("Probe.Name", "Chromosome", "Start")) %>%
  mutate(Chromosome = ifelse(Chromosome == "X", 23, Chromosome),
         Chromosome = ifelse(Chromosome == "Y", 24, Chromosome),
         Chromosome = as.integer(Chromosome)) %>%
  mutate(markerID = paste(Chromosome, Start, sep = ":")) %>%
  filter(!duplicated(markerID)) %>%

  filter(!markerID %in% commonCNV$markerID) %>%
  dplyr::select(-markerID)

#recurrent CNV
library(gaia)
library(TCGAbiolinks)

##### H
set.seed(200)
markers_obj <- load_markers(as.data.frame(markersMatrix_filtered))
nbsamplesH <- length(unique(cnvMatrixH$Sample.Name))
cnv_objH <- load_cnv(cnvMatrixH, markers_obj, nbsamplesH)
suppressWarnings({
  cancer <- "LGG"
  resultsH <- runGAIA(
    cnv_objH, markers_obj,
    output_file_name = paste0("~/Downloads/GAIA_", cancer, "_fltH.txt"),
    aberrations = -1,
    chromosomes = -1,
    approximation = TRUE,
    num_iterations = 5000,
    threshold = 0.25
  )
})

#
RecCNVH <- as_tibble(resultsH) %>%
  mutate_at(vars("q-value"), as.numeric) %>%
  mutate_at(vars(-"q-value"), as.integer) %>% {
    minval = min(.$`q-value`[which(.$`q-value` != 0)])
    mutate(., `q-value` = ifelse(`q-value` == 0, minval, `q-value`))
  } %>%
  mutate(score = -log10(`q-value`)) %>%
  as.data.frame()

threshold <- 0.3
gaiaCNVplot(RecCNVH,threshold)
#####

##### L

set.seed(200)
markers_obj <- load_markers(as.data.frame(markersMatrix_filtered))
nbsamplesL <- length(unique(cnvMatrixL$Sample.Name))
cnv_objL <- load_cnv(cnvMatrixL, markers_obj, nbsamplesL)
suppressWarnings({
  cancer <- "LGG"
  resultsL <- runGAIA(
    cnv_objL, markers_obj,
    output_file_name = paste0("~/Downloads/GAIA_", cancer, "_fltL.txt"),
    aberrations = -1,
    chromosomes = -1,
    approximation = TRUE,
    num_iterations = 5000,
    threshold = 0.25
  )
})

#
RecCNVL <- as_tibble(resultsL) %>%
  mutate_at(vars("q-value"), as.numeric) %>%
  mutate_at(vars(-"q-value"), as.integer) %>% {

    minval = min(.$`q-value`[which(.$`q-value` != 0)])
    mutate(., `q-value` = ifelse(`q-value` == 0, minval, `q-value`))
  } %>%
  mutate(score = -log10(`q-value`)) %>%
  as.data.frame()

threshold <- 0.3
gaiaCNVplot(RecCNVL,threshold)
#####L

save(resultsH, resultsL, RecCNVH, RecCNVL, threshold, file = paste0("~/prj01/", cancer,"_CNV_results.rda"))

#recurrent CNV
library(GenomicRanges)


genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg19") %>%
  filter(external_gene_name != "" & chromosome_name %in% c(1:22, "X", "Y")) %>%
  mutate(
    chromosome_name = ifelse(
      chromosome_name == "X", 23, ifelse(chromosome_name == "Y", 24, chromosome_name)
    ),
    chromosome_name = as.integer(chromosome_name)
  ) %>%
  arrange(chromosome_name, start_position) %>%
  dplyr::select(c("external_gene_name", "chromosome_name", "start_position","end_position")) %>%
  `names<-`(c("GeneSymbol","Chr","Start","End"))

genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
save(genes_GR, genes, file = "~/genes_GR.rda", compress = "xz")


#CNV 
#####H
# 
sCNVH <- dplyr::select(RecCNVH, c(1:4, 6)) %>%
  filter(`q-value` <= threshold) %>%
  arrange(1, 3) %>%
  `names<-`(c("Chr","Aberration","Start","End","q-value"))
#  GenomicRanges 
sCNVH_GR <- makeGRangesFromDataFrame(sCNVH, keep.extra.columns = TRUE)
# 
hits <- findOverlaps(genes_GR, sCNVH_GR, type = "within")
sCNVH_ann <- cbind(sCNVH[subjectHits(hits),], genes[queryHits(hits),])
library(glue)
# 
AmpDel_genesH <- as_tibble(sCNVH_ann, .name_repair = "universal") %>%
  mutate(
    AberrantRegion = glue("{Chr...1}:{Start...3}-{End...4}"),
    GeneRegion = glue("{Chr...7}:{Start...8}-{End...9}")
  ) %>%
  dplyr::select(c(6, 2, 5, 10, 11)) %>%
  mutate(Aberration = if_else( Aberration == 0, "Del", "Amp"))
#####H

#####L
# 
sCNVL <- dplyr::select(RecCNVL, c(1:4, 6)) %>%
  filter(`q-value` <= threshold) %>%
  arrange(1, 3) %>%
  `names<-`(c("Chr","Aberration","Start","End","q-value"))
# 
sCNVL_GR <- makeGRangesFromDataFrame(sCNVL, keep.extra.columns = TRUE)
# 
hits <- findOverlaps(genes_GR, sCNVL_GR, type = "within")
sCNVL_ann <- cbind(sCNVL[subjectHits(hits),], genes[queryHits(hits),])
library(glue)
# 
AmpDel_genesL <- as_tibble(sCNVL_ann, .name_repair = "universal") %>%
  mutate(
    AberrantRegion = glue("{Chr...1}:{Start...3}-{End...4}"),
    GeneRegion = glue("{Chr...7}:{Start...8}-{End...9}")
  ) %>%
  dplyr::select(c(6, 2, 5, 10, 11)) %>%
  mutate(Aberration = if_else( Aberration == 0, "Del", "Amp"))
#####L

AmpDel_genesH = AmpDel_genesH[AmpDel_genesH$GeneSymbol %in% rownames(input),]
AmpDel_genesL = AmpDel_genesL[AmpDel_genesL$GeneSymbol %in% rownames(input),]

save(RecCNVH, AmpDel_genesH, RecCNVL, AmpDel_genesL, file = paste0("~/prj01/", cancer, "_CNV_results2.rda"))

#
#cnv
load("~/LGG_CNV_results2.rda")
library(dplyr)
cnvH <- as_tibble(RecCNVH) %>%
  # threshold = 0.3
  filter(`q-value` <= threshold) %>%
  dplyr::select(c(1, 3, 4, 2)) %>%
  # 
  mutate(
    Chromosome = ifelse(
      Chromosome == 23, "X", ifelse(Chromosome == 24, "Y", Chromosome)
    ),
    Chromosome = paste0("H_chr", Chromosome)
  ) %>%
  mutate(CNV = 1) %>%
  `names<-`(c("Chromosome","Start_position","End_position","Aberration_Kind","CNV"))

cnvL <- as_tibble(RecCNVL) %>%
  # threshold = 0.3
  filter(`q-value` <= threshold) %>%
  dplyr::select(c(1, 3, 4, 2)) %>%
  # 
  mutate(
    Chromosome = ifelse(
      Chromosome == 23, "X", ifelse(Chromosome == 24, "Y", Chromosome)
    ),
    Chromosome = paste0("L_chr", Chromosome)
  ) %>%
  mutate(CNV = 1) %>%
  `names<-`(c("Chromosome","Start_position","End_position","Aberration_Kind","CNV"))

cnv = rbind(cnvH,cnvL)

{
library(circlize)
  cnvH <- as_tibble(RecCNVH) %>%

    filter(`q-value` <= threshold) %>%
    dplyr::select(c(1, 3, 4, 2)) %>%
    mutate(
      Chromosome = ifelse(
        Chromosome == 23, "X", ifelse(Chromosome == 24, "Y", Chromosome)
      ),
      Chromosome = paste0("chr", Chromosome)
    ) %>%
    mutate(CNV = 1) %>%
    `names<-`(c("Chromosome","Start_position","End_position","Aberration_Kind","CNV"))
  
  cnvL <- as_tibble(RecCNVL) %>%

    filter(`q-value` <= threshold) %>%
    dplyr::select(c(1, 3, 4, 2)) %>%
    mutate(
      Chromosome = ifelse(
        Chromosome == 23, "X", ifelse(Chromosome == 24, "Y", Chromosome)
      ),
      Chromosome = paste0("chr", Chromosome)
    ) %>%
    mutate(CNV = 1) %>%
    `names<-`(c("Chromosome","Start_position","End_position","Aberration_Kind","CNV"))
  
  #
  mut_type <- c(
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",'In_Frame_Del', 'In_Frame_Ins','Splice_Site'
  )
  mut <- filter(mafl, Variant_Classification %in% mut_type) %>%
    mutate(
      mut.id = glue::glue("{Chromosome}:{Start_Position}-{End_Position}|{Reference_Allele}"),
      .before = Hugo_Symbol
    ) %>%
    dplyr::select(c("Chromosome","Start_Position","End_Position","Variant_Classification","Hugo_Symbol")) %>%
    distinct_all() %>% {
      
      typeNames <<- unique(.$Variant_Classification)
      type = structure(1:length(typeNames), names = typeNames)
      mutate(., type = type[Variant_Classification])
    } %>%
    dplyr::select(c(1:3,6,4,5))
  
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
#  CNV 
colors <- c("forestgreen","firebrick")
circos.genomicTrackPlotRegion(
  cnvL, ylim = c(0, 1.2),
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region, value, ytop.column = 2,
      ybottom = 0, col = colors[value[[1]] + 1],
      border = "white"
    )
    cell.xlim = get.cell.meta.data("cell.xlim")
    circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
  }
)
# 
colors <- structure(c("blue","green","red", "gold"), names = typeNames[1:4])
circos.genomicTrackPlotRegion(
  mut, ylim = c(1.2, 4.2),
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(
      region, value, cex = 0.8, pch = 16, 
      col = colors[value[[2]]], ...)
  }
)
circos.clear()

#
legend(
  -0.2, 0.2, bty = "n", y.intersp = 1,
  c("Amp", "Del"), pch = 15,
  col = c("firebrick", "forestgreen"),
  title = "CNVs", text.font = 1,
  cex = 0.8, title.adj = 0
)
legend(
  -0.2, 0, bty = "n", y.intersp = 1,
  names(colors), pch = 20, col = colors,
  title = "Mutations", text.font = 1,
  cex = 0.8, title.adj = 0
)
}
