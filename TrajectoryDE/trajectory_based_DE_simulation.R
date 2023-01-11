library(MASS)
library(scuttle)
library(RColorBrewer)
library(PseudotimeDE)
library(mgcv)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tibble)
library(scales)
library(slingshot)
library(parallel)
library(irlba)
library(mclust)
library(stringr)
library(abind)
setwd("~/") 


corenum = 15
options(mc.cores = corenum)

args = commandArgs(T)
seed = as.integer(args[1])

set.seed(seed)

#--------------------------- functions ------------------------
cubic_func = function(x,theta,pos = 0,factor=1){
    sin(theta*(x+pos))*factor
}

sub_index = function(x,rate = 0.8) {
  sample(x = c(1:dim(sce)[2]), size = rate*dim(sce)[2], replace = FALSE)
}
subsampling = function(x, sce) {
  sce <- sce[, x]
  
  rd <- irlba::prcomp_irlba(t(assays(sce)$logcounts), scale. = FALSE)$x[, 1:2]
  # rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
  reducedDims(sce) <- SimpleList(PCA = rd)

  # fit <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "GMM")
  fit <- slingshot(sce, reducedDim = 'PCA')
  # scale pseudotime to [0,1]fit <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "GMM")
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, single_ori_tbl, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime) #if start end reverse
  }
  tbl
}


print("set true model")
n = 1000 #n cells
p = 20 #p genes; 1:10 non-DE; 11:20 DE 
S = 1
#--------------------------------------------------------
#------------- for Y ~ NB(mean,disp)-----------------
#
### dispersion
disp_true = runif(p,0.1,1) # phi # var =  mean + mean^2/phi
###------- mean link -------
PT_true = matrix(sort(runif(n,0,1)),nrow = n, ncol = 1) # pseudotime_true; single lineage
beta0_true = runif(p,2,4)  #! for non-DE, mean_true = exp(beta0) 
######------for the true spline, assume 0 for non-DE; cubic polynomial for DE genes
cubic_list = list()
cubic_list[[1]] = cubic_func(PT_true,theta= pi,pos = -0.4,factor = -10 ) # decreasing (pi)
cubic_list[[2]] = cubic_func(PT_true,theta= pi,pos = 0.4,factor = -10 ) # increasing (pi)
cubic_list[[3]]= cubic_func(PT_true,theta=pi,factor = 10) # increasing and decreasing (pi)
cubic_list[[4]]= cubic_func(PT_true,theta=pi,pos = 0.2,factor = 10) # increasing and decreasing---(pi)
cubic_list[[5]]= cubic_func(PT_true,theta=pi,pos = -0.2,factor = 10) # increasing--- and decreasing (pi)
f = cbind(matrix(rep(0,n*p/2),nrow = n, ncol = p/2),do.call(cbind,cubic_list),do.call(cbind,cubic_list))
mean_true = exp(t(apply(f, 1,function(x){x+beta0_true })))
count = sapply(1:p,function(j,mu,theta){rnegbin(mu[,j],mu = mu[,j], theta = theta[j])},mean_true, disp_true)
rownames(count) = paste0("cl",1:n)
colnames(count) = paste0("G",1:p)
print(count[1:5,1:5])
#
#------------- for Y ~ NB(mean,disp) end -----------------
#---------------------------------------------------------------

print("start testing preparation")
#--------------------------------------------------------------------------------------------

genes_cal = c(paste0("G",c(1:2,5:6,8:9)),"G11","G12","G19","G20")
pv_array = array(dim = c(length(genes_cal),S,6),dimnames = list(genes_cal,1:S,c("emp_pv_0.8", "para_pv_0.8","emp_pv_all", "para_pv_all","fix_pv","theory_pv")))  # gene X Simulation X pvtype
pv_PT = array(dim = c(length(genes_cal),S,3),dimnames = list(genes_cal,1:S,c("pv", "intercept","true_beta")))  # gene X Simulation X pvtype
skip = c()
multi_lineage_flag = c()  #true for multi_lineage
for (itr in 1:S){
    print(itr)
    sce <- SingleCellExperiment(assays = List(counts = t(count)))
    #--------- preprocess  ---------------------
    #! didn't filter
    # geneFilter <- apply(assays(sce)$counts,1,function(x){
    #     sum(x >= 3) >= 10
    # })
    # sce <- sce[geneFilter, ]
    # assays(sce)$norm <- FQnorm(assays(sce)$counts) #742 300
    sce = logNormCounts(sce)

    # PCA
    pca <- prcomp(t(log1p(assays(sce)$logcounts)), scale. = FALSE)
    rd1 <- pca$x[,1:2]
    plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
    reducedDims(sce) <- SimpleList(PCA = rd1)
    ## mclust 
    cl1 <- Mclust(rd1)$classification
    plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
    colData(sce)$GMM <- cl1
    # PI by slingshot 
    print("get PI for the original data")
    sce2 <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "GMM")
    multi_lineage_flag[itr] = !is.null(colData(sce2)$slingPseudotime_2)
    sce <- slingshot(sce, reducedDim = 'PCA')  #! without clusterLabels, it forces one lineage
    single_ori_tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(sce)$slingPseudotime_1))
    # plot(single_ori_tbl$pseudotime,type = 'l')
    #--------- Tests  ---------------------
    # PsudoDE 80 percent
    print("get PI for each subsamples")
    single_index_0.8 <- mclapply(seq_len(n), sub_index,rate=0.8)
    # single_index_all <- mclapply(seq_len(n), sub_index,rate=1)
    # get PI for each subsamples 
    flag = tryCatch(single_sub_tbl_0.8  <- mclapply(single_index_0.8 , subsampling , sce = sce),warning = function(c){1},error = function(c){1})
    if(is.numeric(flag)) {
        skip = c(skip,itr)
        next
    }
    # print(system.time(single_sub_tbl_all <- mclapply(single_index_all , subsampling , sce = sce))) # 100s 15 cores
    single_sub_tbl_all = vector("list", n)
    single_sub_tbl_all = lapply(1:n,function(i){ single_ori_tbl })
    # ---------perform test
    print("pseudoDE subsamples pv")
    print(system.time(res_0.8 <- PseudotimeDE::runPseudotimeDE(gene.vec = genes_cal,
                                        ori.tbl = single_ori_tbl,
                                        sub.tbl = single_sub_tbl_0.8, 
                                        mat = sce, ## You can also use a matrix or SeuratObj as the input
                                        model = "nb",
                                        mc.cores = corenum)))                                      
    print(system.time(res_all <- PseudotimeDE::runPseudotimeDE(gene.vec = genes_cal,
                                        ori.tbl = single_ori_tbl,
                                        sub.tbl = single_sub_tbl_all, 
                                        mat = sce, 
                                        model = "nb",
                                        mc.cores = corenum)))    #    ori.tbl and    sub.tbl[[i]] are the same                                

    #get theoretical pv for sub_all
    print("theory")
    theory_pv = c()
    pv_pt = c()
    interc_pt = c()
    pseudotime =  single_ori_tbl$pseudotime
    indices = which(rownames(sce)%in% genes_cal )
    print(system.time(for(i in 1:length(indices) ){
        j = indices[i]
        expv <- assays(sce)$counts[rownames(sce)[j], single_ori_tbl$cell]
        dat <- cbind(pseudotime, expv) %>% as_tibble()
        fit <- gam(expv ~ s(pseudotime, k = 6, bs = "cr"), family = nb(link = "log"),
                            data = dat,
                            knots = list(pseudotime = c(0:5/5)), control = list(nthreads = 1))
        theory_pv[i] = summary(fit)$s.pv

        #--------------------- if pseudotime_true is known, theoretical performance --------------
        dat2 = data.frame(expv = count[,j],pseudotime  = PT_true)
        fit2 <- gam(expv ~ s(pseudotime, k = 6, bs = "cr"), family = nb(link = "log"),
                                data = dat2,
                                knots = list(pseudotime = c(0:5/5)), control = list(nthreads = 1))
        pv_pt[i] = summary(fit2)$s.pv              
        interc_pt[i]= summary(fit2)$p.coeff 

    }))

    pv_array[,itr,1] = res_0.8$emp.pv  
    pv_array[,itr,2] = res_0.8$para.pv
    pv_array[,itr,3] = res_all$emp.pv
    pv_array[,itr,4] = res_all$para.pv  
    pv_array[,itr,5] = res_all$fix.pv
    pv_array[,itr,6] = theory_pv

    pv_PT[,itr,1] = pv_pt
    pv_PT[,itr,2] = interc_pt
    pv_PT[,itr,3] = beta0_true[indices]

}

save(pv_array,count,skip,pv_PT,multi_lineage_flag,file = paste0("./data/pv_100_",seed,".rda"))



