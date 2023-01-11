setwd("~/")
library(doParallel)
library(qlcMatrix)
library(sparseMatrixStats)
registerDoParallel(cores=127)

print("start loading")
load("./rda_dataset/multiome_counts_matrix.RData")
rna_count_all <- rna_count # 23418
rm(rna_count)
load("./Nov9/rna_rank.RData") #load in B_rank full matrix for RNA
print("finish loading ")
#--- filter out peaks less than 5% cells ------
n <- dim(dna_count_all)[2]
# sum(dna_count_all[2,] != 0) /n
A <- t(dna_count_all)
nonzero <- diff(A@p)
sum(nonzero/n >= 0.05)

filtered_peak <- dna_count_all[which(nonzero/n >= 0.05),] #25514 filtered_peak

print("start")
frag <- 5000
dna_matrix <- filtered_peak
N <- c(0,seq(frag,dim(rna_count_all)[1],frag),(dim(rna_count_all)[1])) #for rna
M <- c(0,seq(frag,dim(dna_matrix)[1],frag), (dim(dna_matrix)[1]))#dna
# N <- c(1,seq(frag,dim(rna_count_all)[1],frag),dim(rna_count_all)[1]) #for rna
# M <- c(1,seq(frag,dim(dna_count_all)[1],frag), dim(dna_count_all)[1])#dna

COR <- foreach(j=2:length(N),.combine=cbind) %dopar% {
# COR <- foreach(j=2:3,.combine=cbind) %dopar% {
# Correlation <- list()
# j = 2
    print(paste0("j=",j))
    # n <- N[j] #rna
    # B <- rna_count_all[N[j-1]:N[j],]
    # B_rank <- rowRanks(B,ties.method = "average")
    B_rank_part <- B_rank[(N[j-1]+1):(N[j]),]
    
    ans <- foreach(i=2:length(M),.combine=rbind) %dopar% {
#   ans <- foreach(i=2:5,.combine=rbind) %dopar% {
# i =2
        print(paste0("i=",i))
        # m <- M[i] #dna
        A <- dna_matrix[(M[i-1]+1):M[i],]
        A_rank <- rowRanks(A,ties.method = "average")
        print("start compute cor")
        COR_1 <- corSparse(t(A_rank),t(B_rank_part))
        rownames(COR_1) <- rownames(A)
        colnames(COR_1) <- rownames(rna_count_all)[(N[j-1]+1):N[j]]
        COR_1
        # Correlation[[i-1]] <- COR_1
    }
    ans
    
}

# format(object.size(COR_1),units = "MB")

save(COR,file = "./Nov9/COR/cor_filtered_all.RData")
