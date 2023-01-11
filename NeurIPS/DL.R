library(Seurat)
# remotes::install_github("satijalab/seurat", "feat/dictionary", quiet = TRUE)
library(SeuratDisk)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(rhdf5)
library(Matrix)
args = commandArgs(T)
setwd("~/") #path is hidden in this demo

load("./rda_dataset/multiome_counts_matrix.RData")
rna_count_all <- rna_count
load("./savedata/annotations.RData")

# #----------------------- metadata -----------------------
meta <- read.csv("metadata.csv",header = TRUE)
head(meta) 
#donor 13176 27678 31800 32606 
donors <- c(13176,27678, 31800, 32606)
#i = 1
i = args[1]
ind_donor <- meta$cell_id[meta$donor == donors[i]]

# ------------------- split data by days ------------------
# dim(dna_count_all)[2] #105942
Prediction <- list()
ind0 = which(colnames(dna_count_all) %in% ind_donor) #35396
Cor <- c()

for(j in 1:1){
print(c("j=",j))
ii <- seq(j,length(ind0),3)
ind = ind0[ii]
dna_count <- dna_count_all[,ind]
dna_count_day7 <- dna_count[,colnames(dna_count) %in% cellid_day7]
# print(c("dim-dna-day7",dim(dna_count_day7)) )
dna_count_day234 <- dna_count[,colnames(dna_count) %in% c(cellid_day2,cellid_day3,cellid_day4)]


rna_count <- rna_count_all[,ind]
rna_count_day7 <- rna_count[,colnames(rna_count) %in% cellid_day7]
# print(c("dim-rna-day7",dim(rna_count_day7)) )
rna_count_day234 <- rna_count[,colnames(rna_count) %in% c(cellid_day2,cellid_day3,cellid_day4)]
# print(c("dim-rna-day234",dim(rna_count_day234)) )
#------------
# ------------- create multiome Seurat -----------------
rna_count_multi <- rna_count_day234

obj.multi <- CreateSeuratObject(counts = rna_count_multi) 
#print(c("dim-obj.multi",dim(obj.multi)) )

dna_count_multi <- dna_count_day234
grange.counts <- StringToGRanges(rownames(dna_count_multi), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- dna_count_multi[as.vector(grange.use), ]  
# print(c("#features removed not in standardchrom",dim(dna_count_day234)[1] - dim(atac_counts)[1] ))

# Add in ATAC-seq data as ChromatinAssay object
#but without fragment file????
chrom_assay <- CreateChromatinAssay(
counts = atac_counts,
sep = c(":", "-"),
genome = 'hg38',
min.cells = 10,
annotation = annotations
)
obj.multi[["ATAC"]] <- chrom_assay

#------------- create ref -------------------
# load("./savedata/objrna.RData")
rna_count_ref <- rna_count_day234
obj.rna =  CreateSeuratObject(counts = rna_count_ref) 
obj.rna = SCTransform(object = obj.rna,conserve.memory = TRUE) %>% RunPCA() %>% RunUMAP(dims = 1:50, return.model = TRUE) 

# save(obj.rna,file="./savedata/objrna_001.RData")
#format(object.size(obj.rna),units="GB")
# -------------- preprocess ------------------------
DefaultAssay(obj.multi) <- "RNA"
obj.multi <- SCTransform(obj.multi,conserve.memory = TRUE) 

# multiome-ATAC & ATAC query is already TF-IDF
DefaultAssay(obj.multi) <- "ATAC"
#obj.multi <- RunTFIDF(obj.multi)
obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
#--------------- query -----------

#dna_count_query <- dna_count_day7
dna_count_query <- dna_count[,colnames(dna_count) %in% cellid_day7]
#  dna_count_query <- dna_count_all[,colnames(dna_count) %in% cellid_day7]


atac_assay <- CreateChromatinAssay(
  counts = dna_count_query,
  sep = c(":", "-"),
  annotation = annotations
)

obj.atac  <- CreateSeuratObject(counts = atac_assay,assay = 'ATAC')

#------- predict------
dims.atac <- 2:50
dims.rna <- 1:50
DefaultAssay(obj.multi) <-  "RNA"
DefaultAssay(obj.rna) <- "SCT"

obj.rna.ext <- PrepareBridgeReference(reference = obj.rna,
                                      bridge = obj.multi, 
                                      reference.reduction = "pca",
                                      reference.dims = dims.rna,
                                      normalization.method = "SCT"
) # normalization.method for reference


#query 38157 22082
bridge.anchor <- FindBridgeTransferAnchors(extended.reference = obj.rna.ext, 
                                           query = obj.atac,
                                           reduction = "lsiproject",
                                           dims = dims.atac
)

trans <- TransferData(anchorset =  bridge.anchor, 
                refdata = GetAssayData(obj.rna[['RNA']]),
                query = obj.atac,weight.reduction = "Bridge.reduc" )


# print(dim(trans[["id"]]) == dim(rna_count_day7))
Prediction[[j]] <- trans[["id"]][]
# save(Prediction,file=paste0("./savedata/trans_donor_pred_",donors[i],"_",j,".RData"))



#------------- calculate correlation matrix ---------------
print("calculate Cor")
test_count<- rna_count_day7
C <- c()
for(m in 1:(dim(test_count)[2])){
  C[m]<- cor(as.vector(trans[["id"]][,m]),test_count[,m])
}

Cor <- c(Cor,C)
save(Cor,file=paste0("./savedata/trans_donor_cor_",donors[i],"_",j,".RData"))
}

save(Prediction,Cor,file=paste0("./savedata/trans_donor_predcor_",donors[i],".RData"))





