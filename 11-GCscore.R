#multivariate cox
{
multi_input.expr<-as.data.frame(t(tpms[gene,]))
multi_input.expr = multi_input.expr[rownames(mrnasi_dat2),]

library(survival)

mySurv<-with(mrnasi_dat2,Surv(OS.time,OS))

multi_COX<-coxph(mySurv ~ ., data=multi_input.expr)  

summary(multi_COX)
step.multi_COX=step(multi_COX,direction = "both")  
step.multi_COX 
}

{
  multi_input.expr = as.data.frame(t(tpms[names(step.multi_COX$coefficients),rownames(mrnasi_dat2)]))
}

#riskscore ROC
{
RiskScore<-predict(step.multi_COX,type = "risk",newdata =multi_input.expr)
risk_group<-ifelse(RiskScore>=median(RiskScore),'high','low')
risk_group
all(rownames(RiskScore) == rownames(risk_group))
all(rownames(mrnasi_dat2)== rownames(risk_group))
KM.input2<-cbind(mrnasi_dat2[,c("OS","OS.time")],RiskScore,risk_group)

library(survival)
library(survminer)

str(KM.input2)
fit2<-survfit(Surv(OS.time,OS) ~ risk_group, data=KM.input2)
print(fit2) #
summary(fit2) #
summary(fit2)$table #
KMsurvival_plot2<-ggsurvplot(fit2,pval = TRUE, 
                             risk.table = "abs_pct",  
                             risk.table.y.text.col = T,
                             risk.table.y.text = FALSE,
                             xlab = "Time in months",  
                             surv.median.line = "hv", 
                             ncensor.plot = FALSE, 
                             palette = c('#e41a1c','#377eb8'), 
                             ggtheme = theme_classic())


library(timeROC)

time_ROC_input2<-KM.input2
time_ROC2<-timeROC(T=time_ROC_input2$OS.time, 
                   delta=time_ROC_input2$OS, 
                   marker=time_ROC_input2$RiskScore, 
                   cause=1, 
                   weighting = "marginal", 
                   times = c(365,1095,1825), 
                   ROC=TRUE,
                   iid=TRUE) #计算AUC

time_ROC2 
#ROC
summary(time_ROC2) 
time_ROC2$AUC
library(ggplot2)
time_ROC2$TP
summary(time_ROC2)

time_ROC.res2<-data.frame(TP_1year=time_ROC2$TP[,1], 
                          FP_1year=time_ROC2$FP[,1],  
                          TP_3year=time_ROC2$TP[,2],  
                          FP_3year=time_ROC2$FP[,2], 
                          TP_5year=time_ROC2$TP[,3], 
                          FP_5year=time_ROC2$FP[,3]) 
time_ROC2$AUC 
TimeROC_plot2<-ggplot()+
  geom_line(data=time_ROC.res2,aes(x=FP_1year,y=TP_1year),size=1,color="red")+
  geom_line(data=time_ROC.res2,aes(x=FP_3year,y=TP_3year),size=1,color="blue")+
  geom_line(data=time_ROC.res2,aes(x=FP_5year,y=TP_5year),size=1,color="black")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2 
  )+
  theme_bw()+
  annotate("text",x=0.75,y=0.25,size=4.5,
           label=paste0("AUC at 1 years = ",round(time_ROC2$AUC[[1]],3)),color="red")+
  annotate("text",x=0.75,y=0.15,size=4.5,
           label=paste0("AUC at 3 years = ",round(time_ROC2$AUC[[2]],3)),color="blue")+
  annotate("text",x=0.75,y=0.05,size=4.5,
           label=paste0("AUC at 5 years = ",round(time_ROC2$AUC[[3]],3)),color="black")+
  labs(x="False positive rate",y="True positive rate")+
  theme(axis.text=element_text(face="bold", size=11,  color="black"),
        axis.title=element_text(face="bold", size=14, color="black")) 
TimeROC_plot2 
}




###risk scoreDEG
{
f2cinput = tpms[,mrnasi_dat2$Tumor_Sample_Barcode]

RiskScore<-predict(step.multi_COX,type = "risk",newdata =multi_input.expr)
risk_group<-ifelse(RiskScore>=median(RiskScore),'high','low')

risk_group<-factor(risk_group,levels=c('low','high'))

all(rownames(risk_group)==colnames(f2cinput))

p_rs <- sapply(1:nrow(f2cinput),function(i){
  data<-cbind.data.frame(gene=as.numeric(f2cinput[i,]),risk_group)
  p=wilcox.test(gene~risk_group, data)$p.value
  return(p)
})
fdr_rs=p.adjust(p_rs,method = "fdr")

# Calculate fold-change for each gene

conditionsLevel_rs<-levels(risk_group)
dataCon1=f2cinput[,c(which(risk_group==conditionsLevel_rs[1]))]
dataCon2=f2cinput[,c(which(risk_group==conditionsLevel_rs[2]))]
foldChange_rs=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results based on FDR threshold

outRst_rs<-data.frame(log2foldChange=foldChange_rs, pValues=p_rs, FDR=fdr_rs)
all(rownames(outRst_rs)==rownames(f2cinput))

outRst_rs=na.omit(outRst_rs)
fdrThres=0.01
input_rs = outRst_rs[outRst_rs$FDR<fdrThres & abs(outRst_rs$log2foldChange)>1.2,]
}