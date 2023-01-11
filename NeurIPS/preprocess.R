library(rhdf5)
library(Matrix)
setwd("~/")  #path is hidden in the demo

#------------- multiome train data ---------
#dna(ATAC) data 
h5ls('train_multi_inputs.h5')
dna_data <- hdf5r::H5File$new('train_multi_inputs.h5', mode = 'r')
rna_data <- hdf5r::H5File$new('train_multi_targets.h5', mode = 'r')
# # ------------------- create sparse matrix ------------------
dna_count_all <- as(dna_data[["train_multi_inputs/block0_values"]][1:228942,1: 105942],"dgCMatrix")
dna_count_all@Dimnames <- list(dna_data[["train_multi_inputs/axis0"]][],dna_data[["train_multi_inputs/axis1"]][])
# dim(dna_count_all) #228942 105942

rna_count_all <- as(rna_data[["train_multi_targets/block0_values"]][1:23418, 1: 105942],"dgCMatrix")
rna_count@Dimnames <- list(rna_data[["train_multi_targets/axis0"]][],rna_data[["train_multi_targets/axis1"]][])
# # dim(rna_count) #23418 105942
# save.image("./rda_dataset/multiome_counts_matrix.RData")

#------------- multiome test data ---------
h5ls('test_multi_inputs.h5')
dna_data <- hdf5r::H5File$new('test_multi_inputs.h5', mode = 'r')
# # ------------------- create sparse matrix ------------------
dna_count_test <- as(dna_data[["test_multi_inputs/block0_values"]][1:228942,1:55935],"dgCMatrix")
dna_count_test@Dimnames <- list(dna_data[["test_multi_inputs/axis0"]][],dna_data[["test_multi_inputs/axis1"]][])

format(object.size(dna_count_test),units="GB")
save(dna_count_test,file="./rda_dataset/multiome_test.RData")
#------------ gene activity matrix------------
h5ls("./Nov9/gene_activity/train_multi_inputs_gam_2000.h5ad")
# h5dump("./Nov9/gene_activity/train_multi_inputs_gam_2000.h5ad")
data <-hdf5r::H5File$new(filename = "./Nov9/gene_activity/train_multi_inputs_gam_2000.h5ad",mode = 'r')

# as.integer(length(data[["/obs/_index"]][])):as.integer(length(data[["/var/_index"]][]))
# gam <- new("dgCMatrix", i =data[["/X/indices"]][], p = data[["/X/indptr"]][],x = data[["/X/data"]][], Dim = c(as.integer(105942),as.integer(25406)))
gam <- new("dgCMatrix", i =data[["/X/indices"]][], p = data[["/X/indptr"]][],x = data[["/X/data"]][], Dim = c(as.integer(length(data[["/obs/_index"]][])),as.integer(length(data[["/var/_index"]][]))))
gam <- t(gam)
rownames(gam) <- data[["/var/_index"]][]
colnames(gam) <- data[["/obs/_index"]][]

format(object.size(gam),units = "GB")


# h5ls("./Nov9/gene_activity/train_multi_inputs_gam10000.h5ad")
