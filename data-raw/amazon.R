## code to prepare `amazon` dataset goes here
library(R.matlab)
library(Matrix)
out <- readMat('../data/amazon.mat')
Xt <- Matrix(out$Xt)
amazon = order_cols_by_freq(Xt, 50)
# idx <- order(colMeans(Xt), decreasing = T)
# d <- min(50, nrow(Xt))
# amazon = Xt[, idx[1:d]]
usethis::use_data(amazon, overwrite = TRUE)
