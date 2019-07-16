devtools::install_github("ryananeff/superNOVA")

rm(list = ls(all.names = TRUE))
set.seed(12345)

library("DGCA")
library("superNOVA")

options(stringsAsFactors = FALSE)

message("loading data")
message(Sys.time())

data(darmanis) #loads daramanis
data(design_mat) #loads design_mat

rownames(design_mat) = colnames(darmanis)

chow_res = chowCor(matA = darmanis, design_mat=design_mat,compare = c("oligodendrocyte", "neuron"))
ddcor_res = flattenChow(chow_res)
