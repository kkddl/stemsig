#fig8
#Volcano plot of Spearman’s coefficient between the expression of compounds targets and the GCScore
{  #F8a
  library(readxl)
  drugtgt <- read_excel
  tgt = drugtgt$`Target genes`[!duplicated(drugtgt$`Target genes`)]
tgt = f2cinput[rownames(f2cinput) %in% tgt,]
all(colnames(tgt) == rownames(mrnasi_dat2))
tgt = as.data.frame(t(log2(tgt+1)))
for_f7a = cbind(tgt,mrnasi_dat2$riskscore)

cor_gene_rs=data.frame()
for(i in names(tgt)){
    x=as.numeric(tgt[,i])
    y=as.numeric(mrnasi_dat2[,"riskscore"])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
      cor_gene_rs=rbind(cor_gene_rs,cbind(tgt=i,cor,pvalue))
}
cor_gene_rs$cor = as.numeric(cor_gene_rs$cor)
cor_gene_rs$pvalue = as.numeric(cor_gene_rs$pvalue)
cor_gene_rs$group = 'not-sig'
cor_gene_rs$group[which((cor_gene_rs$pvalue < 0.05) & (cor_gene_rs$cor > 0.4) )] = 'red'
cor_gene_rs$group[which((cor_gene_rs$pvalue < 0.05) & (cor_gene_rs$cor < -0.4))] = 'blue'
cor_gene_rs$pvalue = -log10(cor_gene_rs$pvalue)


f8a = ggscatter(cor_gene_rs,x='cor',y='pvalue',
          color = 'group',
          palette = c('#2f5688','#BBBBBB','#CC0000'),
          repel = T,
          size = 1,
          label = cor_gene_rs$tgt,
          font.label = 8,
          xlab = 'correlation',
          ylab = '-log10(p-value)')+theme_classic()+
  geom_hline(yintercept = 20, linetype='dashed',color = '#BBBBBB') +
  geom_vline(xintercept = c(-0.4,0.4),linetype='dashed',color = '#BBBBBB')+geom_rug(aes(color = group))#
}

#f8b Volcano plot of Spearman’s coefficient between the CERES score of compounds targets and the GCScore
# import sample_info and crispr_gene_effect
{ #F8b
  library(readr)
  sample_info <- read_csv("")
  glcl = sample_info[which(sample_info$"lineage_subtype" == 'glioma'),'DepMap_ID']

head = fread('./',sep=',',header = T)
head = tibble::column_to_rownames(head,var = 'DepMap_ID')
names(head) = sapply(strsplit(names(head),split = '\\('),function(x) x[1]) #
names(head) = gsub(' ','',names(head))

glclexp = na.omit(head[unlist(glcl),])
tgt = drugtgt$`Target genes`[!duplicated(drugtgt$`Target genes`)]
glclexp = glclexp[,colnames(glclexp) %in% tgt] #CERES score 

#CCLE_expression log10(TPM+1)
clexp = fread('./CCLE_expression.csv',sep=',',header = T)
clexp[1:5,1:5]
clexp = tibble::column_to_rownames(clexp,var = 'V1')
names(clexp) = sapply(strsplit(names(clexp),split = '\\('),function(x) x[1]) #
names(clexp) = gsub(' ','',names(clexp))

clexp = na.omit(clexp[rownames(glclexp),])
glclexp = glclexp[rownames(clexp),]
clrs = predict(step.multi_COX,type = "risk",newdata =clexp) #compute risk score of cell lines
clrs = as.data.frame(clrs)
all(rownames(clrs)==rownames(glclexp))

cor_ceres_rs=data.frame()
for(i in names(glclexp)){
  x=as.numeric(glclexp[,i])
  y=as.numeric(clrs[,"clrs"])
  corT=cor.test(x,y)
  cor=corT$estimate
  pvalue=corT$p.value
  cor_ceres_rs=rbind(cor_ceres_rs,cbind(glclexp=i,cor,pvalue))
}
cor_ceres_rs$cor = as.numeric(cor_ceres_rs$cor)
cor_ceres_rs$pvalue = as.numeric(cor_ceres_rs$pvalue)
cor_ceres_rs$group = 'not-sig'
cor_ceres_rs$group[which((cor_ceres_rs$pvalue < 0.05) & (cor_ceres_rs$cor > 0.3) )] = 'red'
cor_ceres_rs$group[which((cor_ceres_rs$pvalue < 0.05) & (cor_ceres_rs$cor < -0.3))] = 'blue'
cor_ceres_rs$pvalue = -log10(cor_ceres_rs$pvalue)


f8b = ggscatter(cor_ceres_rs,x='cor',y='pvalue',
          color = 'group',
          palette = c('#2f5688','#BBBBBB','#CC0000'),
          repel = T,
          size = 1,
          label = cor_ceres_rs$glclexp,
          font.label = 8,
          xlab = 'correlation',
          ylab = '-log10(p-value)')+theme_classic()+
  geom_hline(yintercept = -log10(0.05), linetype='dashed',color = '#BBBBBB') +
  geom_vline(xintercept = c(-0.3,0.3),linetype='dashed',color = '#BBBBBB')+geom_rug(aes(color = group))
}

{
#f8c Four GDSC-related compounds were identified by Spearman correlation analysis between the GCScore and AUC value
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
testPtype <- fread('./DrugPredictions-gdsc2.csv', data.table = F)
testPtype[1:4, 1:4]
dim(testPtype)
names(testPtype)
rownames(testPtype) = testPtype$V1
testPtype = testPtype[rownames(mrnasi_dat2),]
mrnasi_dat2$riskscore = RiskScore
all(rownames(mrnasi_dat2)==rownames(testPtype))
testPtype = testPtype[,-1]

cor_auc_rs=data.frame()
for(i in names(testPtype)){
  x=as.numeric(testPtype[,i])
  y=as.numeric(mrnasi_dat2[,"riskscore"])
  corT=cor.test(x,y)
  cor=corT$estimate
  pvalue=corT$p.value
  cor_auc_rs=rbind(cor_auc_rs,cbind(testPtype=i,cor,pvalue))
}
cor_auc_rs$cor = as.numeric(cor_auc_rs$cor)
cor_auc_rs$pvalue = as.numeric(cor_auc_rs$pvalue)
cor_auc_rs$group = 'not-sig'
cor_auc_rs$group[which((cor_auc_rs$pvalue < 0.05) & (cor_auc_rs$cor > 0.3) )] = 'red'
cor_auc_rs$group[which((cor_auc_rs$pvalue < 0.05) & (cor_auc_rs$cor < -0.3))] = 'blue'
cor_auc_rs$pvalue = -log10(cor_auc_rs$pvalue)
#boxplot

newinput = testPtype[,cor_auc_rs[which(cor_auc_rs$cor < -0.2),'testPtype']]

all(rownames(newinput)==rownames(mrnasi_dat2))
newinput$risk_group = mrnasi_dat2$risk_group

nn = melt(newinput, id.vars = 'risk_group') 
nn$value = as.numeric(nn$value)


f8c2 <- ggplot(data=nn,aes(x=variable,y=value,fill=risk_group))+ 
  geom_boxplot()+labs(x='',y='IC50')+
  scale_x_discrete(breaks=colnames(newinput)[1:(ncol(newinput)-1)],labels=colnames(newinput)[1:(ncol(newinput)-1)])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "top",axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  stat_compare_means(method = "t.test",aes(label = ..p.signif..))

f8c1<- ggplot(data = cor_auc_rs[which(cor_auc_rs$cor < -0.2),],aes(cor,forcats::fct_reorder(testPtype,cor,.desc = T))) +
  geom_segment(aes(xend=0,yend=testPtype),linetype = 2) +
  geom_point(aes(size=pvalue),col = '#e41a1c') +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.2, -0.4),
                  expand = expansion(mult = c(0.01,.1))) +
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())
}

{ save(prot_expr,file = 'prot_expr.Rdata')
  #f8d Four CTRP2-related compounds were identified by Spearman correlation analysis between the GCScore and AUC value
  library(oncoPredict)
  dir='~/DataFiles//DataFiles/Training Data/'
  CTRP2_Expr = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
  CTRP2_Res = readRDS(file = file.path(dir,"CTRP2_Res.rds"))
  CTRP2_Res <- exp(CTRP2_Res) 
  CTRP2_Expr = log10(CTRP2_Expr+1)
  
  testExpr<- CTRP2_Expr[,sample(1:ncol(CTRP2_Expr),10)]
  #
  table(colSums(testExpr)) 
  
  calcPhenotype(trainingExprData = CTRP2_Expr,
                trainingPtype = CTRP2_Res,
                testExprData = prot_expr,
                batchCorrect = 'eb',  #   "eb" for ComBat  
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 10, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData' )
  
  library(data.table)
  testPtype <- fread('', data.table = F)
  testPtype[1:4, 1:4]
  dim(testPtype)
  names(testPtype)
  rownames(testPtype) = testPtype$V1
  testPtype = testPtype[rownames(mrnasi_dat2),]
  mrnasi_dat2$riskscore = RiskScore
  all(rownames(mrnasi_dat2)==rownames(testPtype))
  testPtype = testPtype[,-1]
  
  cor_ctrp_rs=data.frame()
  for(i in names(testPtype)){
    x=as.numeric(testPtype[,i])
    y=as.numeric(mrnasi_dat2[,"riskscore"])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    cor_ctrp_rs=rbind(cor_ctrp_rs,cbind(testPtype=i,cor,pvalue))
  }
  cor_ctrp_rs$cor = as.numeric(cor_ctrp_rs$cor)
  cor_ctrp_rs$pvalue = as.numeric(cor_ctrp_rs$pvalue)
  cor_ctrp_rs$group = 'not-sig'
  cor_ctrp_rs$group[which((cor_ctrp_rs$pvalue < 0.05) & (cor_ctrp_rs$cor > 0.3) )] = 'red'
  cor_ctrp_rs$group[which((cor_ctrp_rs$pvalue < 0.05) & (cor_ctrp_rs$cor < -0.3))] = 'blue'
  cor_ctrp_rs$pvalue = -log10(cor_ctrp_rs$pvalue)
  
  library(reshape2)
  newinput = testPtype[,cor_ctrp_rs[cor_ctrp_rs$group == 'blue','testPtype']]
  
  all(rownames(newinput)==rownames(mrnasi_dat2))
  newinput$risk_group = mrnasi_dat2$risk_group
  
  nn = melt(newinput, id.vars = 'risk_group') 
  nn$value = as.numeric(nn$value)
  
  f8d2 <- ggplot(data=nn,aes(x=variable,y=value,fill=risk_group))+ 
    geom_boxplot()+labs(x='',y='IC50')+
    scale_x_discrete(breaks=colnames(newinput)[1:(ncol(newinput)-1)],labels=colnames(newinput)[1:(ncol(newinput)-1)])+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "top",axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
    scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
    stat_compare_means(method = "t.test",aes(label = ..p.signif..))
 f8d2 
 
 f8d1<- ggplot(data = cor_ctrp_rs[which(cor_ctrp_rs$cor < -0.3),],aes(cor,forcats::fct_reorder(testPtype,cor,.desc = T))) +
   geom_segment(aes(xend=0,yend=testPtype),linetype = 2) +
   geom_point(aes(size=pvalue),col = '#e41a1c') +
   scale_size_continuous(range =c(2,8)) +
   scale_x_reverse(breaks = c(0, -0.2, -0.4),
                   expand = expansion(mult = c(0.01,.1))) + #左右留空
   theme_classic() +
   labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
   theme(legend.position = "bottom", 
         axis.line.y = element_blank())
}

{# import PRISM data secondary-screen-dose-response...csv
  prism = fread('',header = T)
  prism2 = prism[prism$depmap_id %in% glcl$DepMap_ID,]

  prism2 = prism2[,c('depmap_id','auc','name')]
  prism2 = na.omit(prism2)
  prism3 = tidyr::spread(prism2, name, auc,fill = NA)
  
}
cowplot::plot_grid(f8c1, f8c2, f8d1, f8d2, labels=c("C", "", "D", ""), 
          ncol=2, 
          rel_widths = c(2, 2)) 

#f8f
{
  input_f8f = cor_auc_rs[cor_auc_rs$group != 'not-sig',]
  input_f8f$testPtype = unlist(strsplit(input_f8f$testPtype,split='_'))[c(1,3,5,7,9,11,13,15,17,19)]
  input_f8f$targ = c('PARP1/2','PDGFR/VEGFR','ABL1','TNKS1','NTRK1','IGF1R','HDAC8','PORCN','IDH1','IGF1R') 
 input_f8f$pathway = c('IGF signaling','RTK signaling','TCR signaling','Wnt signaling','RTK signaling','IGF signaling','Chromatin histone acetylation','Wnt signaling','Metabolism','RTK signaling')
  
 par(pin = c(5,2))
 f8f = ggplot(input_f8f, aes(x=testPtype,y=pathway)) +
   geom_point(aes(size=pvalue,color=cor)) +
   scale_color_gradientn('Correlation', 
                         colors=colorRampPalette(c('steelblue','white','darkred'))(100)) + 
   theme_bw() +
   theme(panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text=element_text(size=14, colour = "black"),
         axis.text.x = element_text(angle = 90, hjust = 1),
         axis.text.y = element_text(size=12, colour = "black"),
         axis.title = element_blank(),
         panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))+coord_cartesian(clip = 'off')

 f8fa = annoSegment(object = f8f,
             annoPos = 'top',
             xPosition = c(1:10),
             segWidth = 0.8,
             addText = T,
             textLabel = input_f8f$targ,
             textRot = 45,textHVjust = 1,vjust=1,textCol =rep('black',10),
             textSize = 10)+ theme(plot.margin=unit(rep(1.5,4),'cm')) 
 
 inputbar = data.frame(group = input_f8f$group, 
                       pathway = input_f8f$pathway,
                       num = rep(1,10))

 par(pin = c(2,3))
 bar <- ggplot(inputbar, aes(x=pathway, y=num, fill=group)) + geom_col() +
   labs(x='',y='Num of drugs')+
   coord_flip()+theme_bw()+scale_fill_manual(values = c('steelblue','darkred'))+
   theme(panel.grid =element_blank())+ 
   theme(axis.text = element_blank()) +   ## 
   theme(panel.border = element_blank()) +   ## 
   theme(axis.line = element_line(size=0.2, colour = "black"))

 bar
 
 cowplot::plot_grid(f8fa, bar, align = 'vh')
 }
