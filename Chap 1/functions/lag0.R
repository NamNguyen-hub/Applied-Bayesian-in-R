lag0 <- function(x,p){

R = nrow(x)
C = ncol(x)

x1 = as.matrix(x[1:(R-p),])
out = rbind(matrix(0,p,C), x1)
  return(out)
  }