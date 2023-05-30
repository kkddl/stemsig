

#f5a boxplot+violinplot+dotplot
{
  dat_for_f5 = mrnasi_dat2
dat_for_f5$gc = paste0('t',dat_for_f5$gc)

compare_means(score ~ gc,  data = dat_for_f5)
my_comparisons_f5 = list(c('t1','t2'),c('t2','t3'),c('t1','t3'))
f5a1 <- ggplot(dat_for_f5, mapping = aes(x = gc,  y = score, fill=gc)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons_f5,label = 'p.signif')+ 
  stat_compare_means(label.y = 1.6)

f5a1

f5a2 <- ggplot(dat_for_f5, mapping = aes(x = gc,  y = log2(riskscore+1), fill=gc)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = log2(riskscore+1)),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'riskscore') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons_f5,label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.6) 
f5a2

f5a3 <- ggplot(dat_for_f5, mapping = aes(x = gc,  y = Age, fill=gc)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = Age) ,width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a')) +
  labs(x = '', y = 'Age') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons_f5,label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.6) 
f5a3

f5a4 = ggplot(data = dat_for_f5)+
  geom_bar(aes(x=gc,fill=as.factor(OS)),position = "fill")+
  scale_fill_manual(values=c('#F3DE35','#3A2348'))+theme_bw()+
  theme(panel.grid = element_blank())
f5a4

f5a5 = ggplot(data = dat_for_f5)+
  geom_bar(aes(x=gc,fill=Gender),position = "fill")+
  scale_fill_manual(values=c('#F3DE35','#3A2348'))

f5a6 = ggplot(data = dat_for_f5)+
  geom_bar(aes(x=gc,fill=Subtype),position = "fill")+
  scale_fill_manual(values=color2)+theme_bw()+
  theme(panel.grid = element_blank(),legend.key.size = unit(0.5, "inches")) 
f5a6

f5a7 = ggplot(data = dat_for_f5)+
  geom_bar(aes(x=gc,fill=Diagnosis),position = "fill")+
  scale_fill_manual(values=color2)+theme_bw()+
  theme(panel.grid = element_blank())
f5a7

f5a <- plot_grid(f5a1,f5a2,f5a3,
                 align = "v",axis = "l",
                 ncol = 3)

}

#
{
color1 = c('#2171A9','#A2C9D9')
color2 = c('#e41a1c','#316587','#42AC6A','#F3DE35','#3A2348')
f5e2 = ggplot(data = dat_for_f5)+
  geom_bar(aes(x=gc,fill=MGMT),position = "fill",width = 0.7)+
  scale_fill_manual(values=color1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()
f5e2

f5e3 = ggplot(data = dat_for_f5[dat_for_f5$IDH != 'NA',])+
  geom_bar(aes(x=gc,fill=IDH),position = "fill")+
  scale_fill_manual(values=color1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()
f5e3

f5e4 = ggplot(data = dat_for_f5)+
  geom_bar(aes(x=gc,fill=pq),position = "fill")+
  scale_fill_manual(values=color1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()
f5e4

f5e5 = ggplot(data = dat_for_f5[dat_for_f5$ATRX != 'NA',])+
  geom_bar(aes(x=gc,fill=ATRX),position = "fill")+
  scale_fill_manual(values=color1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()
f5e5

f5e6 = ggplot(data = dat_for_f5[dat_for_f5$TERT != 'NA',])+
  geom_bar(aes(x=gc,fill=TERT),position = "fill")+
  scale_fill_manual(values=color1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()


f5e <- plot_grid(f5e1,f5e2,f5e3,f5e4,f5e5,f5e6,
                 align = "v",axis = "l",
                 ncol = 3)

#TMB
#tmb
x = tmb(maf = maf_lgg)
x = x[x$Tumor_Sample_Barcode %in% rownames(dat_for_f5),]
x$gc = dat_for_f5[match(x$Tumor_Sample_Barcode,rownames(dat_for_f5)),'gc']
x=na.omit(x) 


compare_means(total_perMB_log ~ gc,  data = x)
my_comparisons = list(c('t1','t2'),c('t2','t3'),c('t1','t3'))
f5e1 <- ggplot(x, mapping = aes(x = gc,  y = total_perMB_log, fill=gc)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = total_perMB_log),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'TMB') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.5) 
f5e1
}

#CNV
{
  
  cnv = load('./CNV_tcga_lgg.rda')
  lgg_seg = eval(parse(text=cnv))
  lgg_seg = lgg_seg[,-1]
  lgg_seg = lgg_seg[,c('Sample','Chromosome','Start','End','Num_Probes','Segment_Mean')]
  tumor_seg = lgg_seg[substr(lgg_seg$Sample,14,15) == '01',]
  write.table(tumor_seg,'seg.txt',sep = '\t',col.names = T,row.names = F)
  lgg_marker_file = read.delim('./snp6.na35.remap.hg38.subset.txt.gz')
  lgg_marker_file = lgg_marker_file[lgg_marker_file$freqcnv =='FALSE',]
  lgg_marker_file = lgg_marker_file[,c(1,2,3)]
  write.table(lgg_marker_file,'lgg_marker_file.txt',sep='\t',col.names = T,row.names = F)
  
  library(copynumber)
  
  tumor_seg$Chromosome = as.numeric(tumor_seg$Chromosome)
  tumor_seg = na.omit(tumor_seg)
  tumor_seg = as.data.frame(tumor_seg)
  tumor_seg$arm = rep('p',77214)
  tumor_seg = tumor_seg[,c(1,2,7,3,4,5,6)]
  colnames(tumor_seg) = c('sampleID','chrom','arm','start.pos','end.pos','n.probes',
                          'mean')
  
  
  t1 = rownames(mrnasi_dat2)[mrnasi_dat2$gc == 't1']
  t2 = rownames(mrnasi_dat2)[mrnasi_dat2$gc == 't2']
  t3 = rownames(mrnasi_dat2)[mrnasi_dat2$gc == 't3']
  t1_tumorseg = tumor_seg[substr(tumor_seg$sampleID,1,12) %in% t1,]
  t2_tumorseg = tumor_seg[substr(tumor_seg$sampleID,1,12) %in% t2,]
  t3_tumorseg = tumor_seg[substr(tumor_seg$sampleID,1,12) %in% t3,]
  
 plotFreq(segments=t1_tumorseg, thres.gain=0.2,thres.loss=-0.1)
  plotFreq(segments=t2_tumorseg, thres.gain=0.2,thres.loss=-0.1)
plotFreq(segments=t3_tumorseg, thres.gain=0.2,thres.loss=-0.1)

}


{
#F5c
##oncoplot

maf = data.table::as.data.table(read.csv(file = './TCGA.LGG.mutect.somatic.maf',
                                         header = T,sep = '\t',stringsAsFactors = F,comment.char = '#'))
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 12)
gc = as.data.frame(mrnasi_dat2$gc)
rownames(gc) = rownames(mrnasi_dat2)
cliformaf =gc
cliformaf$Tumor_Sample_Barcode = rownames(cliformaf)


#
maf = maf[-which(maf$Tumor_Sample_Barcode %in% c('TCGA-DU-6392')),]

library(maftools)
maf_lgg=read.maf(maf = maf,clinicalData = cliformaf)
names(maf_lgg@clinical.data)[1] = 'gc'
maf_lgg@clinical.data$gc = as.factor(maf_lgg@clinical.data$gc)
subtypecolors = c('#e41a1c','#316587','#42AC6A','#F3DE35')
names(subtypecolors) = levels(maf_lgg@clinical.data$gc)
phecolors = list(gc = subtypecolors)

getSampleSummary(maf_lgg)
getClinicalData(maf_lgg)
plotmafSummary(maf=maf_lgg, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
f5c = oncoplot(maf=maf_lgg, clinicalFeatures = 'gc', bgCol = "#CCCCCC",borderCol = 'white',annotationColor = phecolors,
                annoBorderCol = "#ebebeb", sortByAnnotation = T,removeNonMutated=T)

##high low risk group maf
clin1 <- subset(cliformaf, gc=="1")$Tumor_Sample_Barcode
clin2 <- subset(cliformaf, gc=="2")$Tumor_Sample_Barcode
clin3 = subset(cliformaf, gc=="3")$Tumor_Sample_Barcode
clin4 = subset(cliformaf, gc=="4")$Tumor_Sample_Barcode

lgg1 <- subsetMaf(maf=maf_lgg, tsb=clin1, isTCGA=TRUE)
lgg2 <- subsetMaf(maf=maf_lgg, tsb=clin2, isTCGA=TRUE)
lgg3 = subsetMaf(maf=maf_lgg, tsb=clin3, isTCGA=TRUE)
lgg4 = subsetMaf(maf=maf_lgg, tsb=clin4, isTCGA=TRUE)


comp <- mafCompare(m1=lgg1, m2=lgg2,  m1Name="subtype1", m2Name="subtype2",minMut=5)
forestPlot(mafCompareRes=comp, pVal=0.05, color=c("royalblue","maroon" ), geneFontSize=0.8)

genes <- comp$results$Hugo_Symbol[1:10]
#F5c2
f5c2 = coOncoplot(m1=lgg1, m2=lgg2,m1Name="subtype1", m2Name="subtype2") #只能画两个


}

#Fig6 ESTIMATE purity 
{
mrnasi_dat2$`ABSOLUTE purity` = as.numeric(as.factor(mrnasi_dat2$`ABSOLUTE purity`))
mrnasi_dat2$`ESTIMATE immune score` = as.numeric(as.factor(mrnasi_dat2$`ESTIMATE immune score`))
mrnasi_dat2$`ESTIMATE stromal score`=as.numeric(as.factor(mrnasi_dat2$`ESTIMATE stromal score`))
mrnasi_dat2$gc = paste0('t',mrnasi_dat2$gc)
names(mrnasi_dat2)[16:18] = c('purity','immunescore','stromalscore')
compare_means(immunescore ~ gc,  data = mrnasi_dat2)
my_comparisons = list(c('t1','t2'),c('t2','t3'),c('t3','t1'))


f6a1 <- ggplot(mrnasi_dat2, mapping = aes(x = gc,  y = purity, fill=gc)) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = purity),width = 0.25,outlier.colour=NA,alpha=0.9) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'purity') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.5) 
f6a1


f6a2 <- ggplot(mrnasi_dat2, mapping = aes(x = gc,  y = immunescore, fill=gc)) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = immunescore),width = 0.25,alpha=0.9,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a')) +
  labs(x = '', y = 'immune score') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.5)
f6a2

f6a3 <- ggplot(mrnasi_dat2, mapping = aes(x = gc,  y = stromalscore, fill=gc)) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = stromalscore),width = 0.25,alpha=0.9,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a')) +
  labs(x = '', y = 'stromal score') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.5)
f6a3
plot_grid(f6a1,f6a2,f6a3,ncol = 3,align = 'h')
}


#cibersort 
#import ciberst
{ library(readr)
  ciberst <- read_delim("CIBERSORTx_Job2_Results.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
  ciberst = tibble::column_to_rownames(ciberst,var = 'Mixture')
ciberst = ciberst[,-(23:25)]
ciberst = ciberst[rownames(mrnasi_dat2),] 
ciberst = na.omit(ciberst)
cforcb = mrnasi_dat2[rownames(ciberst),'gc']

##
library(tidyr);library(dplyr)
re = ciberst
dat <- re %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
a = dat %>% group_by(Cell_type)  %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)
dat$Group = mrnasi_dat2[dat$Sample, 'gc']
library(ggpubr)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
f6b = ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3'))+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
f6b
}

#
{##risk score & immune infiltration correlation
all(rownames(sort_clinic)==rownames(re))
cor1 = cbind(sort_clinic$RiskScore,re)
names(cor1) = makeNames(names(cor1))
colnames(cor1)[1] = 'risk_score'

c7 = ggplot(cor1, aes(x = risk_score, y = T.cells.CD4.naive)) +
  geom_point(colour='blue',size = 1) +
  geom_smooth(method = lm,colour='red',se=F)+
  stat_cor(method = "pearson")+theme_test()

library(patchwork)
c1+c6+c5+
  c7+c3+c4+plot_annotation(title = "", tag_levels = "A")
}

{#F6D 
my_plot = function(x){
exp_f6d = data.frame(tpms[x,rownames(mrnasi_dat2)])
exp_f6d$gc = mrnasi_dat2$gc
names(exp_f6d)[1]='exp'
exp_f6d$exp = log2(exp_f6d$exp+1)

compare_means(exp ~ gc,  data = exp_f6d)
my_comparisons = list(c('s1','s2'),c('s2','s3'),c('s3','s4'))
ggplot(exp_f6d, mapping = aes(x = gc,  y = exp, fill=gc)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = gc,  y = exp),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'Relative expression') +theme_classic()+ ggtitle(x)+
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif')+stat_compare_means()
}
my_list = lapply(c('PDCD1','CD274','PDCD1LG2','CTLA4','CD80','CD86'), my_plot)
library(cowplot)
plot_grid(plotlist = my_list, align = "h", 
          nrow = 2)
}

# TIDE 
tide = t(apply(log2(tpms+1),1,function(x)x-mean(x)))
write.table(tide,file = 'tide.txt',row.names = T,col.names = T,sep = '\t')
#import lgg-tide.csv
library(readr)
lgg_tide <- read_csv("lgg-tide.csv")
tide = lgg_tide[match(rownames(mrnasi_dat2),lgg_tide$Patient) ,]
all(tide$Patient == rownames(mrnasi_dat2))
all(rownames(gc) == tide$Patient)
tide$gc = mrnasi_dat2$gc

color2 = c('#316587','#42AC6A','#F3DE35','#3A2348')
f6e = ggplot(data = tide)+
  geom_bar(aes(x=gc,fill=Responder),position = "fill",width = 0.7)+
  scale_fill_manual(values=color2)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()
f6e

mrnasi_dat2$tide = tide$Responder
f6f = ggplot(mrnasi_dat2, mapping = aes(x = tide,  y = log2(riskscore+1), fill=tide)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = tide,  y = log2(riskscore+1)),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'riskscore') +theme_classic()+ 
  stat_compare_means(label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 3) 
f6f


# f6h drug
library(oncoPredict)
dir='~/DataFiles//DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
table(colSums(testExpr)) 

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = prot_expr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]
dim(testPtype)
names(testPtype)
rownames(testPtype) = testPtype$V1
testPtype = testPtype[rownames(mrnasi_dat2),]
mrnasi_dat2$risk_group = risk_group

for_f6h = cbind(testPtype$temozolomide,mrnasi_dat2$gc)
for_f6h = as.data.frame(for_f6h)
colnames(for_f6h) = c('IC50','gc')
for_f6h$gc = as.character(for_f6h$gc)

f6h = ggplot(for_f6h, mapping = aes(x = gc,  y = IC50, fill=gc)) +
  
  geom_boxplot(aes(x = gc,  y = IC50),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'IC50') +theme_classic()+ 
  stat_compare_means(label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 3) +
  scale_y_continuous(limits = c(0, 1e+07))
f6h


