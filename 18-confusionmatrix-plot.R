{
tab_classify <- as.data.frame.array(conMat3 $table)
### 
blue   <- "#204F8D"
lblue  <- "#498EB9"
dwhite <- "#B6D1E8"
white  <- "#E6EAF7"

pdf("confusionmatrix.pdf",width = 5,height = 5)

par(bty="n", mgp = c(2,0.5,0), mar = c(5.1,6.1,4.1,2.1),tcl=-.25, font.main=3)
par(xpd=NA)
# y
plot(c(0,ncol(tab_classify)),c(0,nrow(tab_classify)), # 
     col = "white", # 
     xlab = "",xaxt = "n", # 
     ylab = "",yaxt = "n") # 


axis(2, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE) 
text(y = 0.5:(ncol(tab_classify)-0.5),
     par("usr")[1], 
     labels = rownames(tab_classify)[nrow(tab_classify):1], 
     srt = 0, pos = 2, xpd = TRUE)
mtext("Prediction", side=2, line = 4.5) 

# x
axis(1, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE) 
text(x = 0.5:(ncol(tab_classify)-0.5),
     par("usr")[1] - 0.2, 
     labels = colnames(tab_classify), 
     srt = 45, pos = 1, xpd = TRUE)
mtext("Actual", side=1, line = 3.5) 

# 
input_matrix <- as.matrix(tab_classify) 
mat.max = max(input_matrix) 
unq.value <- unique(sort(as.vector(input_matrix))) 
rbPal <- colorRampPalette(c(white,dwhite,lblue,blue)) 
col.vec <- rbPal(max(unq.value) + 1) 
col.mat <- matrix(NA,byrow = T,ncol = ncol(input_matrix),nrow = nrow(input_matrix)) 
# 
for (i in 1:nrow(col.mat)) {
  for (j in 1:ncol(col.mat)) {
    col.mat[i,j] <- col.vec[input_matrix[i,j] + 1]
  }
}
# 
x_size <- ncol(input_matrix)
y_size <- nrow(input_matrix)

my_xleft = rep(c(0:(x_size-1)),each = x_size) 
my_xright = my_xleft + 1 
my_ybottom = rep(c((y_size-1):0),y_size) 
my_ytop = my_ybottom + 1 
rect(xleft = my_xleft,
     ybottom = my_ybottom,
     xright = my_xright,
     ytop = my_ytop,
     col=col.mat, #
     border = F) # 

text(my_xleft + 0.5,my_ybottom + 0.5,input_matrix, cex = 1.3) #

invisible(dev.off())
}
