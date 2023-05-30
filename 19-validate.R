
cgga
cggacli 
library(limma)


forf7d1 = prot_exprc2[gfor_ml,] 
forf7d1 = as.data.frame(t(forf7d1))

test2_logistic<-predict(mult.model,newdata = forf7d1) 
forf7d1 $gc = test2_logistic
kmf7d1 = cggacli2[match(rownames(forf7d1),cggacli2$CGGA_ID),]
all(kmf7d1$CGGA_ID==rownames(forf7d1))
kmf7d1$gc = forf7d1$gc
names(kmf7d1)[7] = 'OS.time'
names(kmf7d1)[8] = 'OS'
{
  fit<-survfit(Surv(OS.time,OS) ~ gc, data=kmf7d1)
  f7g1<-ggsurvplot(fit,pval = TRUE, 
                              risk.table = "abs_pct",  
                              risk.table.y.text.col = T,
                              risk.table.y.text = FALSE,
                              xlab = "Time",  
                              surv.median.line = "hv", 
                              ncensor.plot = FALSE, 
                              # legend.labs =
                              #   c("low risk", "high risk"),    
                              palette = c('#e41a1c','#377eb8' ,'#4daf4a','#984ea3'), ###  
                              ggtheme = theme_classic()
  )
  
  f7_input = prot_exprc2[,kmf7d1$CGGA_ID]
  f7_input = t(scale(t(log2(f7_input+1))))
  f7_input[f7_input>2]=2
  f7_input[f7_input < -2]= -2
  f7_input = as.data.frame(f7_input)
  
  kmf7d1 = kmf7d1[order(kmf7d1$gc),]
  rownames(kmf7d1) =kmf7d1$CGGA_ID
  f7_input = f7_input[,rownames(kmf7d1)]
  gc = as.data.frame(kmf7d1$gc)
  rownames(gc) = rownames(kmf7d1)
  names(gc) = 'gc'
  f7f<-pheatmap(f7_input[gfor_ml,], #
                color = colorRampPalette(c('#377eb8','#f7f7f7','#e41a1c'))(100), 
                annotation_col = data.frame(kmf7d1[,c(3,4,5,11:13)],gc), #
                show_colnames =FALSE, #
                show_rownames = T, #
                cluster_cols=F,  
                cluster_rows=T) 
}

{#cgga325
  cgga2 <- read_delim("~CGGA.mRNAseq_325.RSEM-genes.20200506.txt.zip", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
  cggacli2 <- read_delim("~/CGGA.mRNAseq_325_clinical.20200506.txt.zip", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
  library(limma)
  load("")
  genecode=human_geneInfo_genecode_v25
  prot_genes = genecode[which(genecode$type == 'protein_coding'),'symbol']
  prot_exprc2=cgga2[cgga2$Gene_Name %in% prot_genes,]
  prot_exprc2 = avereps(prot_exprc2[,-1],ID = prot_exprc2$Gene_Name)
  #prot_expr
  names_tc 
  prot_exprc2 = prot_exprc2[rownames(prot_exprc2) %in% names_tc,]
  fpkmToTpm <- function(fpkm)
  {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  prot_exprc2 <- apply(prot_exprc2,2,fpkmToTpm)
  cggacli2 = cggacli2[cggacli2$Grade %in% c('WHO II','WHO III'), ]
  prot_exprc2 = prot_exprc2[,cggacli2$CGGA_ID]
  
}

