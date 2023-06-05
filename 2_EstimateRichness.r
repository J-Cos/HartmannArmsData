# 2_EstimateRichness
# This script applies a multiple rarefying approach to estimate the ESV richness on each ARMS unit. Outputs 
# both means, sd, se. and all replicate richness valeus as seperate csvs. 

#packages
    library(tidyverse)
    library(phyloseq)
    library(BioinformaticsPackage) # see this here (not yet on CRAN): https://github.com/J-Cos/BioinformaticsPackage

#functions
    GetAaronsRichnessDataframe<-function(psList=ps_l, Depth=Depths, Gene, NumberOfReplicates){

        smallestsamples<-ps_l[[Gene]] %>%sample_sums %>% sort %>% head


        if (as.vector(smallestsamples[1])==(Depth[[Gene]]+1)) {
            Richnesses<-BioinformaticsPackage::GetMultiRarefiedRichnessEstimates(psList[[Gene]], RarefyDepth=Depth[[Gene]], NumReplicates=NumberOfReplicates)
        } else{print(paste0("Warning: ", Gene, " data has changed"))}

        RichnessReplicates<-Richnesses %>% 
            select(ARMS, Observed, RandomSeed) %>%
            pivot_wider(names_from=RandomSeed, values_from=Observed) 
            
        RichnessStats<-RichnessReplicates %>%
            rowwise() %>%
            mutate(
                mean = mean(c_across(where(is.numeric))),
                SD = sd(c_across(where(is.numeric))),
                SE = SD/NumberOfReplicates) %>%
            select(ARMS, mean, SD, SE)

        return(list(
                "Stats"=RichnessStats,
                "Replicates"=RichnessReplicates
            )
        )
    }

GetRarefiedReplicatePsList<-function(psList=ps_l, Depth=Depths, Gene, NumberOfReplicates){

        smallestsamples<-ps_l[[Gene]] %>%sample_sums %>% sort %>% head


        if (as.vector(smallestsamples[1])==(Depth[[Gene]]+1)) {
            ps_rarefiedReplicate_l<-list()
            for (i in 1:NumberOfReplicates) {
                        ps_rarefiedReplicate_l[[i]]<-ps_l[[Gene]] %>%
                            rarefy_even_depth(  physeq=., sample.size = Depth[[Gene]], rngseed = i, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) 
            }
        } else{
            print(paste0("Warning: ", Gene, " data has changed"))
        }

        return(ps_rarefiedReplicate_l)
}
    

GetEndemicsDataframe<-function(psList){
    df_list<-list()
    for (i in 1:length(psList)){

        ps<-psList[[i]]
        df_list[[i]]<-data.frame(row.names=sample_names(ps),ARMS=sample_data(ps)$ARMS, NumberEndemics=NA, Replicate=i)


        for (sample in sample_names(ps)) {
            
            ESVsPresentInSample<- taxa_names(ps)    [(otu_table(ps)[,sample]>0)]
            
             df_list[[i]][sample,"NumberEndemics"]<-ps %>%
                prune_taxa(ESVsPresentInSample, .) %>%
                prune_samples(!sample_names(ps) %in% sample, .) %>%
                taxa_sums %>%
                `<`(1) %>%
                sum
        }
    }

    df<-bind_rows(df_list) %>%
        as_tibble %>%
        pivot_wider(names_from=Replicate, values_from=NumberEndemics)

    return(df)
}

#1) load data
    ps_l<-readRDS(file.path("Outputs", "AaronsARMS.RDS"))

# define parameters
    Depths<-list(   "COI"=39397,
                    "16s"=55296) #one fewer than minimum depth on each
    Replicates<-100

#2) multirarefy
    #get sample depths
        MatchedCOIARMS_indices<-sample_data(ps_l[["COI"]])$ARMS %in% sample_data(ps_l[["16s"]])$ARMS
        NumberMatchedARMS<-sum(MatchedCOIARMS_indices)
        NonmathcedCoiArms<-sample_data(ps_l[["COI"]])$ARMS[!MatchedCOIARMS_indices]

    #get multirarefied dataframes
        RichnessCoiDf_list<-GetAaronsRichnessDataframe(Gene="COI", NumberOfReplicates=Replicates)
        Richness16sDf_list<-GetAaronsRichnessDataframe(Gene="16s", NumberOfReplicates=Replicates)

#3) reformat and output to csv
    RichnessStats<-full_join(RichnessCoiDf_list[["Stats"]],Richness16sDf_list[["Stats"]], by="ARMS", suffix=c("_COI", "_16s")) %>%
        mutate(NumberRarefyReplicates=Replicates)

    write.csv(RichnessStats, file.path("Outputs", "RichnessStats.csv"))
    write.csv(RichnessCoiDf_list[["Replicates"]], file.path("Outputs", "COIRarefyReplicates.csv"))
    write.csv(Richness16sDf_list[["Replicates"]], file.path("Outputs", "16sRarefyReplicates.csv"))

# 4) get endemics and output to csv
    psCOI_rarefiedReplicate_l<-GetRarefiedReplicatePsList(psList=ps_l, Depth=Depths, Gene="COI", NumberOfReplicates=Replicates)
    EndemicCOIDf<-GetEndemicsDataframe(psCOI_rarefiedReplicate_l)
    ps16_rarefiedReplicate_l<-GetRarefiedReplicatePsList(psList=ps_l, Depth=Depths, Gene="16s", NumberOfReplicates=Replicates)
    Endemics16Df<-GetEndemicsDataframe(ps16_rarefiedReplicate_l)

    EndemicCOIStats<-EndemicCOIDf %>%
            rowwise() %>%
            mutate(
                mean = mean(c_across(where(is.numeric))),
                SD = sd(c_across(where(is.numeric))),
                SE = SD/Replicates) %>%
            select(ARMS, mean, SD, SE)

    Endemic16Stats<-Endemics16Df %>%
            rowwise() %>%
            mutate(
                mean = mean(c_across(where(is.numeric))),
                SD = sd(c_across(where(is.numeric))),
                SE = SD/Replicates) %>%
            select(ARMS, mean, SD, SE)

    EndemicityStats<-full_join(EndemicCOIStats,Endemic16Stats, by="ARMS", suffix=c("_COI", "_16s")) %>%
        mutate(NumberRarefyReplicates=Replicates)
    
    write.csv(EndemicityStats, file.path("Outputs", "EndemicityStats.csv"))
    write.csv(EndemicCOIDf, file.path("Outputs", "endemics_COIRarefyReplicates.csv"))
    write.csv(Endemics16Df, file.path("Outputs", "endemics_16sRarefyReplicates.csv"))
