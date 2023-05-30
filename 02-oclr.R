
#
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapser)

synLogin(email = "", password = "")
synRNA <- synGet( "syn2701943", downloadLocation = "~/Downloads/PCBC")

library(tidyverse)

exp <- read_delim(file = 'rnaseq_norm.tsv',delim='\t') %>%
  separate(col = "tracking_id", sep = "\\.", into = c("Ensembl_ID", "suffix")) %>%
  dplyr::select(-suffix) %>%
  column_to_rownames("Ensembl_ID") %>%
  as.matrix()

library(org.Hs.eg.db)
unimap <- mapIds(
  org.Hs.eg.db, keys = rownames(exp), keytype = "ENSEMBL", 
  column = "SYMBOL", multiVals = "filter"
)
data.exp <- exp[names(unimap),]
rownames(data.exp) <- unimap

synMeta <- synTableQuery("SELECT UID, Diffname_short FROM syn3156503")
library(dplyr)
library(tibble)
metaInfo <- synMeta$asDataFrame() %>%
  dplyr::select(UID, Diffname_short) %>%
  column_to_rownames("UID") %>%
  filter(!is.na(Diffname_short))

X <- data.exp  
y <- metaInfo[colnames(X), ]
names(y) <- colnames(X)
head(y)

#model construct
install.packages('gelnet')
library(gelnet)

# 
X <- data.exp
m <- apply(X, 1, mean)
X <- X - m
# 
sc <- which(y == "SC")
X.sc <- X[, sc]
X.or <- X[, -sc]

model.RNA <- gelnet(t(X.sc), NULL, 0, 1)
head(model.RNA$w)
save(X, y, model.RNA, file = "model.rda")

#prediction
load("~/prj01/model.rda")
common <- intersect(names(model.RNA$w), rownames(tpms))
X <- tpms[common, ]
w <- model.RNA$w[common]

#mRNAsi
score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)

#function
predict.mRNAsi <- function(exp, modelPath='model.rda') {
  load(modelPath)
  
  common <- intersect(names(model.RNA$w), rownames(exp))
  X <- exp[common, ]
  w <- model.RNA$w[common]
  
  score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
  score <- score - min(score)
  score <- score / max(score)
}
# score <- predict.mRNAsi(exp, '~/Downloads/PCBC/model.rda')
