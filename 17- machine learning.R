
library(readr)
consensus_cluster <- read_csv("", 
                              col_names = FALSE)
all(rownames(mrnasi_dat2) == consensus_cluster$X1)
mrnasi_dat2 = mrnasi_dat2[consensus_cluster$X1,]
mrnasi_dat2$gc = consensus_cluster$X2
# caret
library(caret)
n0 = apply( tpms<0.5,1,sum)
i0 = which(n0 > 200)
ftpms = tpms[-i0,]
data = ftpms[rownames(ftpms) %in% gforcluster,rownames(mrnasi_dat2)]
data = data.frame(t(data))
data$gc = as.factor(mrnasi_dat2$gc)
{
  
}
set.seed(3031)
intrain <- createDataPartition(y = data[,'gc'], p= 0.7)[[1]]

# seperate test and training set

training <- data[intrain,]
testing <- data[-intrain,]

{
# Boruta
library(Boruta)
set.seed(1)
boruta <- Boruta(x=training[,-ncol(training)], y=training[,'gc'], pValue=0.01, mcAdj=T, 
                 maxRuns=1000)
boruta
table(boruta$finalDecision)
Boruta::plotImpHistory(boruta)

library(dplyr)
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}
boruta.variable.imp <- boruta.imp(boruta)

library(YSX)#
sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
           xtics_angle = 90)
}
#
boruta.imp <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta")
g1 = boruta.imp$Item



#svm 
library(e1071)
set.seed(2016)
model = svm(formula = gc~., data = training)
# 
pred <- predict(model,training)
tab <- table(Predicted = pred,Actual = training$gc)
1-sum(diag(tab))/sum(tab) #


set.seed(2016)
model = svm(formula = gc~., data = training,kernel='linear')
pred <- predict(model,training)

tab <- table(Predicted = pred,Actual = training$gc)
1-sum(diag(tab))/sum(tab) #
sum(diag(tab))/sum(tab)#

#
results_test<-predict(model,testing,type="class")
res_test<-table(results_test,testing$gc)
1-sum(diag(res_test))/sum(res_test)


summary(model)


set.seed(124)
tmodel <- tune(svm,gc~.,data = training,kernel='linear',
               ranges = list(epsilon = seq(0,1,0.1,),
                             cost = 2^(2:4)))
plot(tmodel)

mymodel <- tmodel$best.model
summary(mymodel)

pred <- predict(mymodel,testing)
tab <- table(Predicted = pred,Actual = testing$gc)
tab
1-sum(diag(tab))/sum(tab)

#feature selection 
set.seed(150)
rfeCNTL <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
svm.features <- rfe(training[, -ncol(training)], training[, 'gc'], sizes = c(1:(ncol(training)-1)), 
                     rfeControl = rfeCNTL, method="svmLinear")

 g2 = svm.features[['optVariables']]
 
 
#xgboost 
library(xgboost)
?xgboost
xgb_training = training
xgb_training$gc = as.numeric(xgb_training$gc)-1
dtrain <- xgb.DMatrix(as.matrix(xgb_training[,-ncol(xgb_training) ]) , label= xgb_training$gc %>% as.matrix)

bst = xgboost(data = dtrain, max.depth = 6, eta = 0.2, nrounds = 20,num_class=3, nthread =2, objective = 'multi:softmax')
pred <- predict(bst, dtrain)

pre_xgb = round(pred)
table(xgb_training$gc,pre_xgb)

tab <- table(Predicted = pre_xgb, Actual = xgb_training$gc)
tab
1-sum(diag(tab))/sum(tab)

#test
xgb_testing = testing
xgb_testing$gc =  as.numeric(xgb_testing$gc)-1
dtest=xgb.DMatrix(as.matrix(xgb_testing[,-ncol(xgb_testing)]) , label= xgb_testing$gc %>% as.matrix)    
pred_test = predict(bst,dtest)
pre_test = round(pred_test)
table(xgb_testing$gc,pre_test)
#
tab <- table(Predicted = pre_test, Actual = xgb_testing$gc)
tab
1-sum(diag(tab))/sum(tab)#


#ROC
library(multiROC)
library(pROC)
xgboost_roc <- multiclass.roc(xgb_training$gc,as.numeric(pre_xgb))

{
  rf_pred = data.frame(pred_test)
  rf_pred = model.matrix(~as.factor(pred_test), rf_pred)
colnames(rf_pred) = paste(colnames(rf_pred),'_pred_RF')
rf_pred = rf_pred[,-1]
true_label = model.matrix(~as.factor(gc), xgb_testing)

colnames(true_label) = gsub('.*?\\.','',colnames(true_label))
colnames(true_label) = paste(colnames(true_label),'_true')
final_df = cbind(true_label,rf_pred)
final_df = as.data.frame(final_df[,-1])
final_df$`factor(gc)1 _true` = as.factor(final_df$`factor(gc)1 _true`)
final_df$`factor(gc)2 _true` = as.factor(final_df$`factor(gc)2 _true`)

names(final_df)[1:2] = c('gc1_true','gc2_true')
names(final_df)[3:4] =c('gc1_pred_RF','gc2_pred_RF')

roc_res <- multi_roc(final_df, force_diag=T)
plot_roc_df <- plot_roc_data(roc_res)

require(ggplot2)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))
}

# xgb.cv 
best_param = list()
best_seednumber = 1234
best_logloss = Inf
best_logloss_index = 0

# 
for (iter in 1:50) {
  param <- list(objective = "multi:softmax",     
                max_depth = sample(6:10, 1),               
                eta = runif(1, .01, .3),                   
                gamma = runif(1, 0.0, 0.2),                
                subsample = runif(1, .6, .9),             
                colsample_bytree = runif(1, .5, .8), 
                min_child_weight = sample(1:40, 1),
                max_delta_step = sample(1:10, 1),
                num_class =3
  )
  cv.nround = 50                                   
  cv.nfold = 5                                     
  seed.number = sample.int(10000, 1)[[1]]
  set.seed(seed.number)
  mdcv <- xgb.cv(data=dtrain, params = param, nthread=6, metrics=c("auc",'mlogloss',"merror"),
                 nfold=cv.nfold, nrounds=cv.nround, watchlist = list(),
                 verbose = F, early_stop_round=8, maximize=FALSE)
  
  min_logloss = min(mdcv$evaluation_log[,test_mlogloss_mean])
  min_logloss_index = which.min(mdcv$evaluation_log[,test_mlogloss_mean])
  
  if (min_logloss < best_logloss) {
    best_logloss = min_logloss
    best_logloss_index = min_logloss_index
    best_seednumber = seed.number
    best_param = param
  }
}

(nround = best_logloss_index)
set.seed(best_seednumber)
best_seednumber
(best_param)               

# auc | rmse | error 
xgb_plot=function(input,output){
  history=input
  train_history=history[,1:6]%>%mutate(id=row.names(history),class="train")
  test_history=history[,7:12]%>%mutate(id=row.names(history),class="test")
  colnames(train_history)=c("logloss.mean","logloss.std","auc.mean","auc.std","error.mean","error.std","id","class")
  colnames(test_history)=c("logloss.mean","logloss.std","auc.mean","auc.std","error.mean","error.std","id","class")
  
  his=rbind(train_history,test_history)
  his$id=his$id%>%as.numeric
  his$class=his$class%>%factor
  
  if(output=="auc"){ 
    auc=ggplot(data=his,aes(x=id, y=auc.mean,ymin=auc.mean-auc.std,ymax=auc.mean+auc.std,fill=class),linetype=class)+
      geom_line()+
      geom_ribbon(alpha=0.5)+
      labs(x="nround",y=NULL,title = "XGB Cross Validation AUC")+
      theme(title=element_text(size=15))+
      theme_bw()
    return(auc)
  }
  
  
  if(output=="mlogloss"){
    rmse=ggplot(data=his,aes(x=id, y=logloss.mean,ymin=logloss.mean-logloss.std,ymax=logloss.mean+logloss.std,fill=class),linetype=class)+
      geom_line()+
      geom_ribbon(alpha=0.5)+
      labs(x="nround",y=NULL,title = "XGB Cross Validation mlogloss")+
      theme(title=element_text(size=15))+
      theme_bw()
    return(rmse)
  }
  
  if(output=="merror"){
    error=ggplot(data=his,aes(x=id,y=error.mean,ymin=error.mean-error.std,ymax=error.mean+error.std,fill=class),linetype=class)+
      geom_line()+
      geom_ribbon(alpha=0.5)+
      labs(x="nround",y=NULL,title = "XGB Cross Validation ERROR")+
      theme(title=element_text(size=15))+
      theme_bw()
    return(error)
  }
  
}

xgb_plot(mdcv$evaluation_log[,-1]%>%data.frame,"auc")
xgb_plot(mdcv$evaluation_log[,-1]%>%data.frame,"mlogloss")
xgb_plot(mdcv$evaluation_log[,-1]%>%data.frame,"merror")

#
model <- xgb.train(data=dtrain, params=best_param, nrounds=nround, nthread=6, watchlist = list())
#Importance
importanceRaw <- xgb.importance(feature_names=colnames(dtrain), model = model)
library(Ckmeans.1d.dp)
xgb.ggplot.importance(importanceRaw)     # importance 

# #--------------------------------------------------------------------------------------
# feature selection    
cum_impt=data.frame(names=importanceRaw$Feature,impt=cumsum(importanceRaw$Importance))
cum_impt=filter(cum_impt,cum_impt$impt<0.9)
selected_feature<-cum_impt$names

train=dplyr::select(xgb_training,selected_feature)
dtrain<- xgb.DMatrix(data=train%>%as.matrix,label= xgb_training$gc%>%as.matrix)

model <- xgb.train(data=dtrain, params=best_param, nrounds=nround, nthread=2, watchlist = list())
impselect <- xgb.importance(feature_names=colnames(dtrain), model = model)
xgb.ggplot.importance(impselect)
# #--------------------------------------------------------------------------------------
model
g3 = model[['feature_names']]



##Lasso
library(lars) 
library(glmnet)
x = as.matrix(log2(training[,-ncol(training)]+1))
training$gc = as.numeric(training$gc)
y=training$gc
set.seed(123)
model_lasso <- glmnet(x, y, family="multinomial", nlambda=50, alpha=1)
print(model_lasso)

head(coef(model_lasso, s=c(model_lasso$lambda[29],0.009)))
#plot.glmnet(model_lasso, xvar = "norm", label = TRUE)
plot(model_lasso, xvar="lambda", label=TRUE)
set.seed(123)
cv_fit <- cv.glmnet(x=x, y=y, alpha = 1, nlambda = 1000)
plot(cv_fit)
#plot.cv.glmnet(cv_fit)
# 
c(cv_fit$lambda.min,cv_fit$lambda.1se) 
set.seed(123)
model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(y ,lasso.prob)
dat=as.data.frame(re[,1:2])
colnames(dat)=c('event','prob')
dat$event=as.factor(dat$event)
library(ggpubr) 
p <- ggboxplot(dat, x = "event", y = "prob",
               color = "event", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

library(ROCR)
library(glmnet)
library(caret)
# calculate probabilities for TPR/FPR for predictions
pred <- predict(re[,2], re[,1])
perf <- performance(pred,"tpr","fpr")
performance(pred,"auc") # shows calculated AUC for model
plot(perf,colorize=FALSE, col="black") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
head(fit$beta)

g4 = choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]

for_ml = intersect(g1,g2)
for_ml2 =  intersect(g3,g4) 
gfor_ml = intersect(for_ml,for_ml2)
save(g1,g2,g3,g4,gfor_ml,gforcluster,input_rs, file='17.Rdata')


## venn plot

data_ls = list(g1,g2,g3,g4)
str(data_ls)
names(data_ls) = c('G1','G2','G3','G4')

{library(magrittr)
  library(readr)
  library(purrr)
  library(VennDiagram)
  library(colorfulVennPlot)
  library(RColorBrewer)
  library(grDevices)
  library(Cairo)
  library(stringr)
  library(tibble)
  library(tidyr)
# 
number_area <- 2^length(data_ls) - 1 # 

# 
intToBin <- function(x){
  if (x == 1)
    1
  else if (x == 0)
    NULL
  else {
    mod <- x %% 2
    c(intToBin((x-mod) %/% 2), mod)
  }
}

x_area <- seq(number_area) %>% map(., ~intToBin(.x)) %>% 
  map_chr(., ~paste0(.x, collapse = "")) %>% 
  map_chr(., ~str_pad(.x, width = length(data_ls), side = "left", pad = "0"))

#

# 
G_union <- data_ls$G1 %>% union(data_ls$G2) %>% 
  union(data_ls$G3) %>% union(data_ls$G4)


area_calculate <- function(data_ls, character_area){
  character_num <- 1:4 %>% map_chr(., ~substr(character_area, .x, .x)) %>% 
    as.integer() %>% as.logical()
  
  element_alone <- G_union
  for (i in 1:4) {
    element_alone <- 
      if (character_num[i]) {
        intersect(element_alone, data_ls[[i]])
      } else {
        setdiff(element_alone, data_ls[[i]])
      }
  }
  return(element_alone)
}

# 
element_ls <- map(x_area, ~area_calculate(data_ls = data_ls, character_area = .x))

# 
quantity_area <- map_int(element_ls, length)

# 
percent_area <- (quantity_area / sum(quantity_area)) %>% round(3) # 
percent_area <- (percent_area * 100) %>% paste0("%")

#
length_pallete <- max(quantity_area) - min(quantity_area) + 1
color_area <- colorRampPalette(brewer.pal(n = 7, name = "YlGn"))(length_pallete) 
color_tb <- tibble(quantity = seq(min(quantity_area), max(quantity_area), by = 1),
                   color = color_area) 

nest1 <- tibble(quantity = quantity_area, percent = percent_area, area = x_area) %>% 
  group_by(quantity) %>% nest() %>% 
  left_join(color_tb, by = "quantity") %>% 
  arrange(quantity) %>% 
  unnest(cols = c(data))

# Venn
regions <- nest1$quantity
names(regions) <- nest1$area

CairoPDF(file = "venn_num.pdf", width = 8, height = 6)
plot.new()
f7b=plotVenn4d(regions, Colors = nest1$color, Title = "", 
           labels = paste0("G", 1:4)) 
}


#logistics回归
library(lattice)
library(ggplot2)
library(caret)
library(e1071)
library(nnet)
library(pROC)
gfor_ml 
train_data = training[,c(gfor_ml,'gc')]
train_data$gc = as.factor(train_data$gc)
test_data = testing[,c(gfor_ml,'gc')]
test_data$gc = as.factor(test_data$gc)

mult.model<-multinom(gc~.,data=train_data)
summary(mult.model)

train_logistic<-predict(mult.model,newdata = train_data) 
table(train_data$gc,train_logistic)
conMat3 = confusionMatrix(factor(train_logistic),factor(train_data$gc))
# 
z <- summary(mult.model)$coefficients/summary(mult.model)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1))*2
p

# 
exp(coef(mult.model))
# 
head(pp<-fitted(mult.model))
# 
pre_logistic<-predict(mult.model,newdata = test_data) 

# 
 table(test_data$gc,pre_logistic)
 # 
 conMat4<-confusionMatrix(factor(pre_logistic),factor(test_data$gc))
 
 
 
 
#ROC 
 mn_pred <- predict(mult.model, test_data, type = 'prob') #test_data / train_data
 mn_pred <- data.frame(mn_pred)
 colnames(mn_pred) <- paste(colnames(mn_pred), "_pred_MN")
 
 test_data$gc = as.factor(test_data$gc)
 true_label = model.matrix(~gc, test_data) 
 true_label = true_label[,-1]
 true_label <- data.frame(true_label)
 colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
 colnames(true_label) <- paste(colnames(true_label), "_true")
 final_df <- cbind(true_label, mn_pred)
 names(final_df)[1:2] = c('2_true','3_true')
 names(final_df)[3:5] = c('1_pred_RF','2_pred_RF','3_pred_RF')
 
 roc_res <- multi_roc(final_df, force_diag=T)
 pr_res <- multi_pr(final_df, force_diag=T)
 
 plot_roc_df <- plot_roc_data(roc_res)
 plot_pr_df <- plot_pr_data(pr_res)
 
 require(ggplot2)
 ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
   geom_path(aes(color = Group), size=1.5) +
   geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                colour='grey', linetype = 'dotdash') +
   theme_bw() + 
   theme(plot.title = element_text(hjust = 0.5), 
         legend.justification=c(1, 0), legend.position=c(.95, .05),
         legend.title=element_blank(), 
         legend.background = element_rect(fill=NULL, size=0.5, 
                                          linetype="solid", colour ="black"))
 
 ggplot(plot_pr_df, aes(x=Recall, y=Precision)) + 
   geom_path(aes(color = Group), size=1.5) + 
   theme_bw() + 
   theme(plot.title = element_text(hjust = 0.5), 
         legend.justification=c(1, 0), legend.position=c(.95, .05),
         legend.title=element_blank(), 
         legend.background = element_rect(fill=NULL, size=0.5, 
                                          linetype="solid", colour ="black"))
 