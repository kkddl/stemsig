##clustering #https://www.jianshu.com/p/a85faae7fea6 #https://www.jianshu.com/p/302e4e1fb9d9
# library(ConsensusClusterPlus)
# 
# res <- ConsensusClusterPlus(as.matrix(f2c_input[rownames(input),]), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "pam", corUse = "complete.obs", distance = "pearson",seed=123456, plot="png", writeTable=T)
# 
# # 5.用PAC的方法确定最佳聚类数
# #   面积最小值对应K为最佳K
# Kvec = 2:maxK
# x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
# PAC = rep(NA,length(Kvec)) 
# names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
# for(i in Kvec){
#   M = res[[i]]$consensusMatrix
#   Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
#   PAC[i-1] = Fn(x2) - Fn(x1)
# } 
# 
# optK = Kvec[which.min(PAC)]  # 理想的K值 = 3
# 
# icl = calcICL(res,
#               title="./",
#               plot="pdf")
