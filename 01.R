

###########
library(openxlsx)
library(tidyverse)

#data import
fpkm <- read.table("./rawdata/TCGA-LGG.htseq_fpkm.tsv.gz",sep = "\t" , header = T, 
                   #row.names = "Ensembl_ID",
                   stringsAsFactors = FALSE ,
                   check.names = FALSE)

probeMap[1:4,]
library(dplyr)
expr <- fpkm %>% 
  inner_join(probeMap, by = c("Ensembl_ID" = "id")) %>% 
  select(gene , starts_with("TCGA") )
cli <- read_delim("rawdata/TCGA-LGG.GDC_phenotype.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)#import dataset导入
cli[1:4,1:4]
#survival data
surv <- read.table("./rawdata/TCGA-LGG.survival.tsv.gz", header = T,
                   stringsAsFactors = FALSE)

cli_surv <- cli %>% 
  inner_join(surv,by = c("submitter_id.samples" = "sample")) %>% 
  select(submitter_id,age_at_index.demographic,gender.demographic,
         neoplasm_histologic_grade,primary_diagnosis.diagnoses,OS,OS.time)

cli_surv2 <- tcga_all %>% 
  inner_join(cli_surv,by = c("Case" = "submitter_id")) %>% 
  select(Case,age_at_index.demographic,gender.demographic,OS,OS.time,
         neoplasm_histologic_grade,Histology,primary_diagnosis.diagnoses,`IDH status`,`IDH/codel subtype`,`1p/19q codeletion`,`MGMT promoter status`,
         "TERT promoter status",'ATRX status',"BRAF V600E status","ABSOLUTE purity","ESTIMATE immune score","ESTIMATE stromal score","Transcriptome Subtype")
cli_surv3 = cli_surv2[!duplicated(cli_surv2$Case),]
clinic_surv = tibble::column_to_rownames(cli_surv3,var='Case')

load("./rawdata/human_geneInfo_genecode_v25.rda")
genecode=human_geneInfo_genecode_v25
prot_genes=genecode[genecode$type == 'protein_coding','symbol']
prot_expr=expr[expr$gene %in% prot_genes,]
prot_expr=na.omit(prot_expr)

prot_expr = prot_expr[,c(1,which(substr(colnames(prot_expr),14,16)=='01A'))]
library(limma)
prot_expr = avereps(prot_expr[,-1],ID = prot_expr$gene)




prot_expr = prot_expr[names_tc,]

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
fpkm = 2^prot_expr-1
tpms <- apply(fpkm,2,fpkmToTpm)
colSums(tpms)

## tpms TCGALGG tpm mat prot_exprc CGGA693 tpm matrix