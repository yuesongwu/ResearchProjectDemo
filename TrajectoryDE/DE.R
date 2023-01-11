library(edgeR)
library(dplyr)
library(ggplot2)
library(DESeq2)
setwd("./DE")

#this DE gene test demo contains 
#edgeR 
#DESeq2
#t-test 
#edgeR + Wilcox test(this one need large sample size, not suitable for this case)

# -------------------- edgeR -------------------------
data <- read.delim("gene_expected_count.annot.txt")
ind=55488:55492
WT_Batf2 =  data[-ind,5:12]
#remove all zero genes
rownames(WT_Batf2)=data$gene_id[-ind]
WT_Batf2 <- WT_Batf2[-which(apply(WT_Batf2,1,sum) == 0),]


Counts = WT_Batf2
# colnames(Counts) = c(paste0("WT",1:4),paste0("Baft2",1:4))
head(Counts)

group = factor(c(rep("WT",4), rep("Baft",4)))
y <- DGEList(counts=Counts, genes=rownames(Counts),group = group)
y
keep <- filterByExpr(y,group =group)
y <- y[keep,,keep.lib.sizes=FALSE]
dim(y)
y <- calcNormFactors(y)
# mice <- as.factor(1:8)
# design <- model.matrix(~group+mice)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
res=topTags(lrt,n= "Inf")

res = as.data.frame(res) 
res = res%>%arrange(PValue)%>% left_join(data%>%select(gene_id,external_gene_name),by = c("genes" = "gene_id") )
res = res[,c(1,7,2:6)]
head(res)
write.csv(res,file = "./edgeR_result.csv",quote = FALSE)

##add
top = 10
DE_genes_count = res[1:top,c(1,5,6) ] %>% left_join(data%>%select(gene_id,external_gene_name,paste0("X6972.WG.",1:8)),by = c("genes" = "gene_id") )

df = cbind.data.frame(genes = rep(DE_genes_count[,4],each =8),counts = log(matrix(t(DE_genes_count[,-c(1:4)])) +1 ),group = c(rep("WT",4) , rep("Baft",4)))
p2 <- ggplot(df, aes(x=genes,y = counts,fill = group )) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))
p2 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
summary(decideTestsDGE(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col=2)
# res = res%>%arrange(res$logFC)
# res =res%>%left_join(data%>%select(external_gene_name,gene_id,description),by = c("genes"="gene_id"))
# res = res %>% select(external_gene_name,genes,description,logFC ,logCPM,PValue,FDR)

# head(res)

deGenes <- decideTestsDGE(lrt,p.value=0.1)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes,cex = 0.8)
abline(h=c(-1, 1), col=2)
?plotSmear



# -------------- Deseq2-----------------
data <- read.delim("gene_expected_count.annot.txt")

ind=55488:55492

coldata =data.frame(condition = factor(c(rep("WT",4),rep("Batf2",4))))
WT_Batf2 =  data[-ind,5:12]
head(WT_Batf2)

#remove all zero genes
rownames(WT_Batf2)=data$gene_id[-ind]
WT_Batf2 <- WT_Batf2[-which(apply(WT_Batf2,1,sum) == 0),]


WT_Batf2 <- DESeqDataSetFromMatrix(countData = WT_Batf2,
                              colData = coldata,
                              design= ~ condition)
WT_Batf2
# dds@assays@data@listData$counts
WT_Batf2$condition <- relevel(WT_Batf2$condition, ref = "WT")
# dds$condition
WT_Batf2 <- DESeq(WT_Batf2)
res <- results(WT_Batf2)

# --------------------- t-test ---------------------------------------
TPM <- read.delim("gene_TPM.annot.txt")
ind=55488:55492
#remove all zero genes
WT_Batf2 =  TPM[-ind,5:12]
rownames(WT_Batf2)=TPM$gene_id[-ind]
WT_Batf2 <- WT_Batf2[-which(apply(WT_Batf2,1,sum) == 0),]
head(WT_Batf2)

errind <- c(22390,31589)
WT_Batf2[errind,]
# WT_Batf2 <- WT_Batf2[-errind,]
pvalues <- sapply(1:nrow(WT_Batf2), function(i){
if(i %in% errind ) return(NA)
p <- t.test(WT_Batf2[i,1:4],WT_Batf2[i,5:8])$p.value
return(p)
})

WT_Batf2 <- cbind(WT_Batf2,pvalues = pvalues)
WT_Batf2$fdr <- p.adjust(WT_Batf2$pvalues, method = "fdr")
WT_Batf2 = WT_Batf2%>%rownames_to_column()%>%arrange(pvalues)%>%left_join(TPM%>%select(external_gene_name,gene_id),by = c("rowname" = "gene_id") )


# DE_genes_list <- WT_Batf2$external_gene_name[1:10]
# DE_genes_count <- data.frame(id = DE_genes_list ) %>% left_join(TPM%>%select(external_gene_name,paste0("X6972.WG.",1:8)),by = c("id" = "external_gene_name") )
log(WT_Batf2[4,2:9]+1)
log(WT_Batf2[5,2:9]+1)

df = cbind.data.frame(genes = rep(WT_Batf2$external_gene_name[1:10],each =8),counts = log(matrix(t(WT_Batf2[1:10,2:9])) +1),group = c(rep("WT",4) , rep("Baft",4)))
p2 <- ggplot(df, aes(x=genes,y = counts,fill = group )) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))
p2 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


# Calculate the fold-change for each gene
dataCon1 <- WT_Batf2[,2:5]
dataCon2 <- WT_Batf2[,6:9]
foldChanges <-   log2(rowMeans(dataCon2)+1) - log2(rowMeans(dataCon1)+1) 
WT_Batf2$logFC <- foldChanges
WT_Batf2$pvalues_sign = ifelse(WT_Batf2$logFC>=0,WT_Batf2$pvalues,-WT_Batf2$pvalues)
WT_Batf2=WT_Batf2 %>% arrange(pvalues_sign)

write.csv(WT_Batf2%>%select(external_gene_name,pvalues_sign), file = "DE_t_symbol.csv",row.names = FALSE,quote = FALSE)


write.csv(WT_Batf2%>%select(rowname,external_gene_name,logFC,pvalues,fdr)%>%arrange(pvalues), file = "ttest_result.csv",row.names = FALSE,quote = FALSE)



# ----------edgeR + Wilcox ------------------------
library(edgeR)
library(dplyr)
setwd("./DE")
data <- read.delim("gene_expected_count.annot.txt")
ind=55488:55492
WT_Batf2 =  data[-ind,5:12]
#remove all zero genes
rownames(WT_Batf2)=data$gene_id[-ind]
WT_Batf2 <- WT_Batf2[-which(apply(WT_Batf2,1,sum) == 0),]


Counts = WT_Batf2
# colnames(Counts) = c(paste0("WT",1:4),paste0("Baft2",1:4))
head(Counts)

group = factor(c(rep("WT",4), rep("Baft",4)))
y <- DGEList(counts=Counts, genes=rownames(Counts),group = group)
y
keep <- filterByExpr(y,group =group)
y <- y[keep,,keep.lib.sizes=FALSE]
dim(y)
y <- calcNormFactors(y)
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)
# Run the Wilcoxon rank-sum test for each gene
conditions = group
pvalues <- sapply(1:nrow(count_norm), function(i){
data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
p <- wilcox.test(gene~conditions, data)$p.value
return(p)
})
fdr <- p.adjust(pvalues, method = "fdr")
# Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
# write.table(outRst[outRst$FDR<fdrThres,], file = "examples/examples.WilcoxonTest.rst.tsv", sep="\t", quote = F, row.names = T, col.names = T)
dim(outRst[outRst$FDR<fdrThres,])


range(pvalues )
range(fdr)

