#######################
# create null random matrices and fit spiec easi
###########################

# desired clustering at which to generate networks
# either ASV (set to "") or cOTU (set to "_cOTU")
cluster<-c("")
#cluster<-c("_cOTUs")

# packages and functions
    library(tidyverse)
    library(SpiecEasi)

# randomise matrices either "between" or "within" samples 
# preferred method is "within". As stated in https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176751 this 
# is the "most basic null model is the element-wise permutation of the OTU abundance matrix." However this perturbation is only 
# performed within samples to retain the intra-sample species abundance relationship (SARs; an important conserved property in real communities).
# shufflign all elements across samples would result in unrealistic samples where species abundances do not fit observed SARs.
# Nonetheless SARs across the entire dataset are lost in this process. This is equivalent to the "neutral approximation" 
# that taxa are ecologically indistiguishable and thus may be common in one sample though rare in another (see https://pubmed.ncbi.nlm.nih.gov/22341498/).
# The null networks thus generated allow the question to be posed: which properties of this network are generated by niche processes including 1)
# biotic interactions among nodes, and 2) shared components of nodes' niches that lead to non-interaction based cooccurance.
randomiseMat<-function(mat, method, rngseed=1){
    set.seed(rngseed)
    if (method=="between") {
        rand <- mat[sample(nrow(mat)),]
        rownames(rand)<-rownames(mat)  
    } else if (method=="within") {
        rand<-mat
        for (i in 1:nrow(mat)){
            rand[i,]<-mat[i,sample(ncol(mat))]
        }
    }
    return(rand)
}

fitSE<-function(mat_list, lmr=5e-04){
    se<-spiec.easi(
        mat_list, 
        method='mb', 
        sel.criterion = "bstars",
        verbose=TRUE,
        lambda.min.ratio=lmr, 
        nlambda=100, 
        pulsar.select=TRUE, 
        pulsar.params = list(rep.num=50, thresh=0.05, ncores=4)
        )
    print(paste0("Stability: ", getStability(se)))
    return(se)
}

# 1) load networking matrices - ready to input into spiec-easi
Matrices<-readRDS(paste0("Outputs/NetworkMatrices", cluster, ".RDS"))

all<-Matrices[["all"]]
coastal<-Matrices[["coastal"]]
oceanic<-Matrices[["oceanic"]]


# 7) create random matrices
all_rand<-lapply( all, randomiseMat, method="within")
coastal_rand<-lapply( coastal, randomiseMat, method="within")
oceanic_rand<-lapply( oceanic, randomiseMat, method="within")

# 8) run spiec easi

getMaxReps<- function(N, maxRAM=12) {(maxRAM * 1e+9) / (N^2*8) }
getN<-function(mat_list) {  sum(unlist(lapply(mat_list, ncol)))   }

n<-getN(all)
n
getMaxReps(n)

se_all<-fitSE(all)
saveRDS( se_all, paste0("Outputs/se_all", cluster, ".RDS"))

se_all_rand<-fitSE(all_rand)
saveRDS( se_all_rand, paste0("Outputs/se_all_rand", cluster, ".RDS"))



################################################################
# unused coastal and oceanic split (not enough data to fit)
# multiple cluster options therefore are not implemented
se_coastal<-fitSE(coastal)
saveRDS( se_coastal, "Outputs/se_coastal.RDS")

se_coastal_rand<-fitSE(coastal_rand)
saveRDS( se_coastal_rand, "Outputs/se_coastal_rand.RDS")

se_oceanic<-fitSE(oceanic)
saveRDS( se_oceanic, "Outputs/se_oceanic.RDS")

se_oceanic_rand<-fitSE(oceanic_rand, lmr=5e-05)
saveRDS( se_oceanic_rand, "Outputs/se_oceanic_rand.RDS")
