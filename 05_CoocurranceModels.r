#######################
# cooccurance models & nulls
###########################

# packages and functions
    library(phyloseq)
    library(tidyverse)
    library(gridExtra)
    library(SpiecEasi)

applyfilter_ps<-function(ps, minProportion=0.001){
   totalReads<-sum(sample_sums(ps))
    ps_filt<-ps %>%
            filter_taxa(., function(x) sum(x) > minProportion*totalReads, TRUE) %>%
            filter_taxa(., function(x) sum(x>0)>=8, TRUE)
    print(paste0("Filtering out OTUs with abundance less than ",  minProportion*totalReads  ))
    print(paste0( "Proportion of total reads retained: ", sum(sample_sums(ps_filt))/totalReads))
    return(ps_filt)
}

applyfilter_mat<-function(metab, minProportion=0.001){
   totalReads<-sum(colSums(metab))
    metab_filt<-metab[rowSums(metab)>minProportion*totalReads & rowSums(metab>0)>=8,]
    print(paste0("Filtering out OTUs with abundance less than ",  minProportion*totalReads  ))
    print(paste0( "Proportion of total reads retained: ", sum(colSums(metab_filt))/totalReads))
    return(metab_filt)
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

randomiseMat<-function(mat){
    rand <- mat[sample(nrow(mat)),]
    rownames(rand)<-rownames(mat)  
    return(rand)
}

fitSE<-function(mat_list){
    se<-spiec.easi(
        mat_list, 
        method='mb', 
        sel.criterion = "bstars",
        verbose=TRUE,
        lambda.min.ratio=5e-04, 
        nlambda=100, 
        pulsar.select=TRUE, 
        pulsar.params = list(rep.num=10, thresh=0.05, ncores=4)
        )
    print(paste0("Stability: ", getStability(se)))
    return(se)
}

getMultiNetworkDataframe<-function(se){
        
    #### MB method ####
    # Symmetrize the MB assymetric graph by taking the maximum value of the edge pairs e_ij / e_ji
    sebeta <- symBeta(getOptBeta(se), mode='maxabs')
    elist <- Matrix::summary(sebeta)

    # Calculate the number of columns - i.e. number of OTUs/metabolites 
    n_bac <- ncol(EsvMat)
    n_euk <- ncol(CoiMat)
    n_met <- ncol(metabMat)

    # Create a new column that encodes the types of association between two nodes 
    # Nodes are numbered in order from bac to euk to mets 

    elist <- elist %>%
    mutate(interaction.type = case_when(
        i %in% 1:n_bac & 
        j %in% 1:n_bac ~ "Bacteria-Bacteria",

        i %in% (n_bac+1):(n_bac+n_euk) & 
        j %in% (n_bac+1):(n_bac+n_euk) ~ "Eukaryote-Eukaryote",
        
        i %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) & 
        j %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Metabolite-Metabolite",

        i %in% 1:n_bac & 
        j %in% (n_bac+1):(n_bac+n_euk) ~ "Bacteria-Eukaryote",
        j %in% 1:n_bac & 
        i %in% (n_bac+1):(n_bac+n_euk) ~ "Bacteria-Eukaryote",

        i %in% 1:n_bac & 
        j %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Bacteria-Metabolite",
        j %in% 1:n_bac & 
        i %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Bacteria-Metabolite",

        i %in% (n_bac+1):(n_bac+n_euk) & 
        j %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Eukaryote-Metabolite",
        j %in% (n_bac+1):(n_bac+n_euk) & 
        i %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Eukaryote-Metabolite",

        ))
    return(elist)
}

getEdgeDf<-function(se.test, se.rand){
    df1<-cbind(getMultiNetworkDataframe(se.test), group="Observed")
    df2<-cbind(getMultiNetworkDataframe(se.rand), group="Null")
    df<-rbind(df1, df2)
    df<-as.data.frame(df)
    df<-as_tibble(df) %>%
        mutate(value=as.numeric(x))
    return(df)
}

makeNetworkFigure<-function(edge_df){
    edge_df %>% 
    ggplot(data=., aes(x=value))+
        geom_histogram( aes(fill=group),alpha=0.5, position="identity",  binwidth=0.05)+
                scale_fill_manual(values=c("black", "blue"))+

        facet_wrap(~interaction.type)+
        geom_vline(xintercept=0)+
        theme_classic()+
        scale_x_continuous(limits=c(-0.5, 1), expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab ("Edge Weights")+ ylab("Frequency")+
        theme(legend.title=element_blank(), strip.background=element_blank(), 
            strip.text = element_text(size = 12))
}

# 1) load data
ps_l<-readRDS(file.path("Outputs", "AaronsARMS.RDS"))
metab<-read.csv("Data/Hartmann_ARMS_100_cleaned_600_forJames_matchIDs.csv") %>%
    select(-Metabs)

# 2) rarefy data
ps16<-rarefy_even_depth(ps_l[["16s"]], rngseed=1)
psc<-rarefy_even_depth(ps_l[["COI"]], rngseed=1)


# 3) apply filtering to reduce data to manageable quantity
ps16_filt<-applyfilter_ps(ps16)
psc_filt<-applyfilter_ps(psc)
metab_filt<-applyfilter_mat(metab)


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

metabMat<-metab_filt %>% t
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


# 7) create random matrices
all_rand<-lapply( all, randomiseMat)
coastal_rand<-lapply( coastal, randomiseMat)
oceanic_rand<-lapply( oceanic, randomiseMat)

# 8) run spiec easi

getMaxReps<- function(N, maxRAM=12) {(maxRAM * 1e+9) / (N^2*8) }
getMaxReps(N=dim(EsvMat)[2]+dim(CoiMat)[2]+dim(metabMat)[2])

se_all<-fitSE(all)
se_all_rand<-fitSE(all_rand)
se_coastal<-fitSE(coastal)
se_coastal_rand<-fitSE(coastal_rand)
se_oceanic<-fitSE(oceanic)
se_oceanic_rand<-fitSE(oceanic_rand)



# 9) export plotting dataframe 
edge_all_df<-getEdgeDf(se_all, se_all_rand)
edge__coastal_df<-getEdgeDf(se_coastal, se_coastal_rand)
edge_oceanic_df<-getEdgeDf(se_oceanic, se_oceanic_rand)



# make plots - move to next script
p_all<-makeNetworkFigure(edge_all_df)
ggsave("Figures/Figure_all.pdf", p_all,  width = 12, height = 8)
p_coastal<-makeNetworkFigure(edge__coastal_df)
ggsave("Figures/Figure_coastal.pdf", p_coastal,  width = 12, height = 8)
p_oceanic<-makeNetworkFigure(edge_all_df)
ggsave("Figures/Figure_oceanic.pdf", p_oceanic,  width = 12, height = 8)
