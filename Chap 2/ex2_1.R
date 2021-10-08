# library <- function(x){
#   x = toString(substitute(x))
#   if(!require(x,character.only=TRUE)){
#     install.packages(x, repos="http://cran.us.r-project.org")
#     base::library(x,character.only=TRUE)
#   }}
rm(list=ls())
library('rstudioapi')
setwd(dirname(getActiveDocumentContext()$path))
getwd()


