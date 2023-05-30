
library(ggplot2)

# F1a 
score_dat = as.data.frame(score)
score_dat$ID = substr(rownames(score_dat),1,12)
score_dat = score_dat[order(score_dat$score),]
score_dat$score = as.numeric(score_dat$score)
rownames(score_dat)= score_dat[,'ID']

for_anno = clinic_surv
for_anno$ID = rownames(for_anno)

mrnasi_dat = merge(score_dat,for_anno,by = 'ID')
mrnasi_dat2 = mrnasi_dat[,-1];rownames(mrnasi_dat2)=mrnasi_dat$ID
mrnasi_dat2 = mrnasi_dat2[order(mrnasi_dat2$score),]
mrnasi_dat2[,2] = as.numeric(mrnasi_dat2[,2])
mrnasi_dat2[,3] = as.factor(mrnasi_dat2[,3])
mrnasi_dat2[,7] = as.factor(mrnasi_dat2[,7])
mrnasi_dat2[,9] = as.factor(mrnasi_dat2[,9])
mrnasi_dat2[,11] = as.factor(mrnasi_dat2[,11])
mrnasi_dat2[,12] = as.factor(mrnasi_dat2[,12])
mrnasi_dat2[,14] = as.factor(mrnasi_dat2[,14])
mrnasi_dat2[,19] = as.factor(mrnasi_dat2[,19])

mrnasi_dat2$ID = 1:nrow(mrnasi_dat2)

#p1
p1=ggplot(mrnasi_dat2, aes(x=ID, y=score)) +
  geom_line(color="#3077B2", size=2) +geom_area( fill="#3077B2", alpha=0.4) +
  theme_classic()+labs(x="Patient ID",y="mRNAsi")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0)) #调整坐标轴起点
p1


minput = t(mrnasi_dat2[,1])
colnames(minput) = rownames(mrnasi_dat2);rownames(minput) = 'score'
library(pheatmap)

#p2
# p2 = pheatmap(minput,
#         annotation_col = data.frame(mrnasi_dat2[,c(2,3,7,9,11,12,14,19)]), 
#          show_colnames =FALSE, 
#          show_rownames = F, 
#          cluster_cols=FALSE,  
#          cluster_rows=FALSE,
#         legend = F,cellheight = 10)

library(cowplot)

#03-fig1

My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

{#实心p1
p1 <- ggplot(mrnasi_dat2, aes(x=ID, y=score))+
  geom_bar(stat="identity",col = "#3077B2")+
  My_Theme1 +
  labs(y = "mRNAsi") +
  scale_x_continuous(expand = c(0,0)) 

colnames(mrnasi_dat2)

colnames(mrnasi_dat2) = c('score','Age','Gender','OS','OS.time','Grade','Histology','Diagnosis','IDH',
                          'IDH/codel','pq','MGMT','TERT','ATRX','BRAF','ABSOLUTE purity','ESTIMATE immune score',
                          "ESTIMATE stromal score",'Subtype','ID')

p2 <- ggplot(mrnasi_dat2,aes(ID,1))+
  geom_tile(aes(fill = Histology))+
  My_Theme2+
  labs(y = "Histology")+
  scale_fill_manual(values = c("#74C065","#B096C2","#FAA96C","#EFF882","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) + #不留空
  theme(legend.direction = "vertical", legend.box = "vertical") 

p3 <- ggplot(mrnasi_dat2,aes(ID,1))+
  geom_tile(aes(fill = IDH))+
  My_Theme2+
  labs(y = "IDH status")+
  scale_fill_manual(values = c("#84CBA8","#F2FA9C","#B1A5C8","#FA5E5C")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.direction = "vertical", legend.box = "vertical") #不留空

p4 <- ggplot(mrnasi_dat2,aes(ID,1))+
  geom_tile(aes(fill = Gender ))+
  My_Theme2+
  labs(y = "Gender")+
  scale_fill_manual(values=c("#E00F0A","#3D6197","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.direction = "vertical", legend.box = "vertical") #不留空  

p5 <- ggplot(mrnasi_dat2,aes(ID,1))+
  geom_tile(aes(fill = pq))+
  My_Theme2 +
  labs(y = "1p/19q codeletion")+
  scale_fill_manual(values=c("#64B685","#FC6D4C","#6D6466")) +
  scale_x_continuous(expand = c(0,0))  +
  theme(legend.direction = "vertical", legend.box = "vertical")#不留空

p6 = ggplot(mrnasi_dat2,aes(ID,1))+
  geom_tile(aes(fill = MGMT))+
  My_Theme2 +scale_fill_manual(values=c("#64B685","#FC6D4C","#6D6466")) +
  labs(y = "MGMT promoter status")+
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.direction = "vertical", legend.box = "vertical")

p7 = ggplot(mrnasi_dat2,aes(ID,1))+
  geom_tile(aes(fill = Subtype))+
  My_Theme2 +scale_fill_manual(values=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')) +
  labs(y = "Subtype")+
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.direction = "vertical", legend.box = "vertical")


legend_a <- get_legend(p2+theme(legend.position = "bottom"))
legend_b <- get_legend(p3+theme(legend.position = "bottom"))
legend_c <- get_legend(p4+theme(legend.position = "bottom"))
legend_d <- get_legend(p5+theme(legend.position = "bottom"))
legend_e <- get_legend(p6+theme(legend.position = "bottom"))
legend_f <- get_legend(p7+theme(legend.position = "bottom"))

library(patchwork)
blank_p <- plot_spacer() + theme_void() # create a blank plot for legend alignment 
leg = plot_grid(legend_a,legend_b,legend_c,legend_d,legend_e,legend_f,nrow = 1, axis = 'lt',align = 'v')

p <- plot_grid(p1,p2,p3,p4,p5,p6,p7,leg,
               align = "v",axis = "l",
               ncol = 1, rel_heights = c(4,1,1,1,1,1,1,4)
               )
}

#F1a finished
p
#################################################################################

#F1b 
library(ggpubr)
{
  dat_for_subtype = mrnasi_dat2[mrnasi_dat2$Subtype != 'NA',]
compare_means(score ~ Subtype,  data = dat_for_subtype)
my_comparisons = list(c('CL','PN'),c('CL','ME'),c('ME','NE'),c('NE','PN'),c('ME','PN'))
f1b <- ggplot(dat_for_subtype, mapping = aes(x = Subtype,  y = score, fill=Subtype)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = Subtype,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif')+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.55) 
f1b

f1b2 = ggplot(mrnasi_dat2, mapping = aes(x = Gender,  y = score, fill=Gender)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = Gender,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(label = 'p.signif')
f1b2

dat_for = mrnasi_dat2[!is.na(mrnasi_dat2$Grade),]
f1b3 = ggplot(dat_for, mapping = aes(x = Grade,  y = score, fill=Grade)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = Grade,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')
f1b3

mrnasi_dat2$OS = as.factor(mrnasi_dat2$OS)
f1b4 = ggplot(mrnasi_dat2, mapping = aes(x = OS,  y = score, fill=OS)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = OS,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')
f1b4

dat_for = mrnasi_dat2[mrnasi_dat2$IDH != 'NA',]
f1b5 = ggplot(dat_for, mapping = aes(x = IDH,  y = score, fill=IDH)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = IDH,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')
f1b5

f1b6 <- ggplot(mrnasi_dat2, mapping = aes(x = pq,  y = score, fill=pq)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = pq,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')
f1b6

f1b7 <- ggplot(mrnasi_dat2, mapping = aes(x = MGMT,  y = score, fill=MGMT)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = MGMT,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')
f1b7

dat_for = mrnasi_dat2[mrnasi_dat2$TERT != 'NA',]
f1b8 <- ggplot(dat_for, mapping = aes(x = TERT,  y = score, fill=TERT)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = TERT,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')
f1b8

dat_for = mrnasi_dat2[mrnasi_dat2$ATRX != 'NA',]
f1b9 <- ggplot(dat_for, mapping = aes(x = ATRX,  y = score, fill=ATRX)) +
  geom_violin(trim = T,alpha = 0.8) +
  geom_jitter(position = position_jitter(0.2), alpha=0.3)+
  geom_boxplot(aes(x = ATRX,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')
f1b9

f1ba <- plot_grid(f1b,f1b2,f1b3,f1b4,f1b5,f1b6,f1b7,f1b8,f1b9,
                 align = "v",axis = "l",
                 ncol = 3)
}
########################################
#F1c 
#TMB
maf = data.table::as.data.table(read.csv(file = '.TCGA.LGG.mutect.somatic.maf',
                                         header = T,sep = '\t',stringsAsFactors = F,comment.char = '#'))
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 12)
mrnasi_dat2$Tumor_Sample_Barcode = rownames(mrnasi_dat2)
library(maftools)
maf_lgg=read.maf(maf = maf,clinicalData = mrnasi_dat2)

tmutbur = tmb(maf_lgg, captureSize = 50, logScale = TRUE)
mut = merge(tmutbur,mrnasi_dat2,by = 'Tumor_Sample_Barcode',all = F)
mut = mut[mut$Tumor_Sample_Barcode != 'TCGA-DU-6392',]
mut$logTMB = log2(mut$total)
{
f1c1 = ggplot(mut, aes(x=ID, y=score))+
  geom_bar(stat="identity",col = "#3077B2")+
  My_Theme1 +
  labs(y = "mRNAsi") +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.direction = "vertical", legend.box = "vertical")

f1c2 = ggplot(mut,aes(ID,1))+
  geom_tile(aes(fill = logTMB))+
  My_Theme2+
  labs(y = "TMB")+
  scale_x_continuous(expand = c(0,0))+
  theme(legend.direction = "vertical", legend.box = "vertical")

f1c3 = ggplot(mut,aes(ID,1))+
  geom_tile(aes(fill = IDH ))+
  My_Theme2+
  labs(y = "IDH")+
  scale_fill_manual(values=c("#E00F0A","#3D6197","#6D6466")) +
  scale_x_continuous(expand = c(0,0))+
  theme(legend.direction = "vertical", legend.box = "vertical")


f1c4 = ggplot(mut,aes(ID,1))+
  geom_tile(aes(fill = TERT ))+
  My_Theme2+
  labs(y = "TERT")+
  scale_fill_manual(values=c("#E00F0A","#3D6197","#6D6466")) +
  scale_x_continuous(expand = c(0,0))+
  theme(legend.direction = "vertical", legend.box = "vertical")

f1c5 = ggplot(mut,aes(ID,1))+
  geom_tile(aes(fill = ATRX ))+
  My_Theme2+
  labs(y = "ATRX")+
  scale_fill_manual(values=c("#E00F0A","#3D6197","#6D6466")) +
  scale_x_continuous(expand = c(0,0))+
  theme(legend.direction = "vertical", legend.box = "vertical")

legend_a <- get_legend(f1c1+theme(legend.position = "bottom"))
legend_b <- get_legend(f1c2+theme(legend.position = "bottom"))
legend_c <- get_legend(f1c3+theme(legend.position = "bottom"))
legend_d <- get_legend(f1c4+theme(legend.position = "bottom"))
legend_e <- get_legend(f1c5+theme(legend.position = "bottom"))

leg = plot_grid(legend_a,legend_b,legend_c,legend_d,legend_e,legend_f,nrow = 1, axis = 'lt',align = 'v')

f1c <- plot_grid(f1c1,f1c2,f1c3,f1c4,f1c5,leg,
               align = "v",axis = "l",
               ncol = 1, rel_heights = c(4,1,1,1,1,4)
               )
}
mut$tmbg = ifelse(mut$total > median(mut$total), 'h','l')
f1d = ggplot(mut, mapping = aes(x = tmbg,  y = score, fill=tmbg)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = tmbg,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#fdae61','#2c7bb6')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(method = 't.test')

mrnasi_dat2$agegroup = ifelse(mrnasi_dat2$Age>=60,'>=60',
                              ifelse(mrnasi_dat2$Age>50,'50-59',
                                     ifelse(mrnasi_dat2$Age>40,'40-49','<=40')))
mrnasi_dat2$agegroup = factor(mrnasi_dat2$agegroup, levels=c('<=40','40-49','50-59','>=60'),ordered = T)

f1b1 = ggplot(mrnasi_dat2, mapping = aes(x = agegroup,  y = score, fill=agegroup)) +
  geom_violin(trim = T,alpha = 0.5) +
  geom_jitter(position = position_jitter(0.3), alpha=0.3)+
  geom_boxplot(aes(x = agegroup,  y = score),width = 0.25,outlier.colour=NA) +
  scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3')) +
  labs(x = '', y = 'mRNAsi') +theme_classic()+ 
  stat_compare_means(label.y = 1.55)

