#######################
# create matrices for spiec-easi input
###########################

# packages and functions
    library(phyloseq)
    library(tidyverse)

applyfilter_ps<-function(ps, minProportion=0.0001){
   totalReads<-sum(sample_sums(ps))
    ps_filt<-ps %>%
            filter_taxa(., function(x) sum(x) > minProportion*totalReads, TRUE) %>%
            filter_taxa(., function(x) sum(x>0)>=8, TRUE)
    print(paste0("Filtering out OTUs with abundance less than ",  minProportion*totalReads  ))
    print(paste0( "Proportion of total reads retained: ", sum(sample_sums(ps_filt))/totalReads))
    return(ps_filt)
}

alignNames<-function(names){
    names %>%
        str_replace("OLO4", "OLO04") %>%
        str_replace("PAL1$", "PAL01") %>%
        str_replace("TUT1$", "TUT01") %>%
        str_replace("JAR8", "JAR08") %>%
        str_replace("KIN5", "KIN05") %>%
        str_replace("JAR4", "JAR04")  %>%
        str_replace("OFU3", "OFU03")  %>%
        str_replace("OFU4", "OFU04") %>%
        str_replace("ROS4", "ROS04")  %>%
        str_replace("BAA1", "BAA01") %>%
        str_replace("ART1", "ART01") %>%
        str_replace("KOA1", "KOA01") 
}

filterByRealm<-function(mat, realm){ 
    newmat<-mat[rownames(mat) %in% realmEncoding$ARMS[realmEncoding$realm==realm],] 
    newmat<-newmat[,colSums(newmat)!=0]
    print(paste0(dim(newmat)[2], " OTUs remaining in realm ", realm, " ; ", dim(mat)[2]-dim(newmat)[2], " removed"))
    return(newmat)
}

# 1) load data
ps_l<-readRDS(file.path("Outputs", "AaronsARMS.RDS"))
metab<-read.csv("Data/Hartmann_ARMS_100_cleaned_600_forJames_matchIDs.csv") %>%
    select(-Metabs)

# 2) rarefy data
ps16<-rarefy_even_depth(ps_l[["16s"]], rngseed=1)
psc<-rarefy_even_depth(ps_l[["COI"]], rngseed=1)
psm<-transform_sample_counts(otu_table(metab, taxa_are_rows=TRUE), function(x) ceiling(x/10) ) # reduce to same order of magnitude as metabarcoidng data
psm<-rarefy_even_depth(psm, rngseed=1)


# 3) apply filtering to reduce data to manageable quantity
ps16_filt<-applyfilter_ps(ps16)
psc_filt<-applyfilter_ps(psc)
psm_filt<-applyfilter_ps(psm)


# 4) get matrices for spieceasi
EsvMat<-ps16_filt %>%
    otu_table %>% 
    t %>%
    unclass
EsvMat<-EsvMat[order(rownames(EsvMat)),]

CoiMat<-psc_filt %>%
    otu_table %>% 
    t %>%
    unclass
CoiMat<-CoiMat[order(rownames(CoiMat)),]

metabMat<-psm_filt %>%
    otu_table %>% 
    t %>%
    unclass
metabMat<-metabMat[order(rownames(metabMat)),]

# 5) align matrices - indented lengthy data cleaning section

    # first correct and align vectors of rownames
    rn16<-rownames(EsvMat) %>%
        #corrections
        str_replace("DMSOA", "") %>%
        str_replace("DMSO", "") %>%
        str_replace("\\.", "") %>%
        str_replace("\\.", "") %>%
        str_replace("SES", "")  %>%
        str_replace("BS", "") %>%
        str_replace("_", "") %>% 
        alignNames

    rncoi<-rownames(CoiMat)%>%
        #corrections
        str_replace("DMSO", "") %>%
        str_replace("BS", "") %>%
        str_replace("_", "") %>%
        str_replace("_", "") %>% 
        alignNames


    rnmetab<-rownames(metabMat)%>%
        #corrections
        str_replace("VIP_", "") %>%
        str_replace("INDO_", "") %>%
        str_replace("LI_", "") %>%
        str_replace("SI_", "") %>%
        str_replace("_", "") %>%
        str_replace("NAP1B", "BNAP1B") %>% 
        alignNames

    #then check they are correct manually
    rn16[!rn16 %in% rncoi ]
    rn16[!rn16 %in% rnmetab ]

    rncoi[!rncoi %in% rn16 ]
    rncoi[!rncoi %in% rnmetab ]

    rnmetab[!rnmetab %in% rn16 ]
    rnmetab[!rnmetab %in% rncoi ]

    #then reassign them 
    rownames(metabMat)<-rnmetab
    rownames(CoiMat)<-rncoi
    rownames(EsvMat)<-rn16

    #reorder to alphabetical now we have final row names
        EsvMat<-EsvMat[order(rownames(EsvMat)),]
        CoiMat<-CoiMat[order(rownames(CoiMat)),]
        metabMat<-metabMat[order(rownames(metabMat)),]

    # Get a string of the common column names using intersect function (twice)
    common_column_names <- intersect(intersect(rn16, rncoi), rnmetab)

    # Filter out columns whose name isn't in common_column_names
    EsvMat <- EsvMat[rownames(EsvMat) %in% common_column_names,]
    CoiMat <- CoiMat[rownames(CoiMat) %in% common_column_names,]
    metabMat <- metabMat[rownames(metabMat) %in% common_column_names,]


# 6) subset to oceanic and continental
#categorise oceanic and coastal (coastal=1, oceanic =2)
realmEncoding<-data.frame(ARMS=common_column_names, realm=c(rep(1, 23), rep(2, 15), rep(1, 8), rep(2, 3)))

all<-list(EsvMat, CoiMat, metabMat)
coastal<-lapply( all, filterByRealm, 1)
oceanic<-lapply( all, filterByRealm, 2)


# 7) save networking matrices - ready to input into spiec-easi
saveRDS( list("all"=all, "coastal"=coastal, "oceanic"=oceanic), "Outputs/NetworkMatrices.RDS")
