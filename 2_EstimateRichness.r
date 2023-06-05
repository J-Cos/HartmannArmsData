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

GetSiteEndemicsDataframe<-function(psList){
    #gets the number of site endemic ESVs per ARMS. These are definied as ESVs occuring on at least 2 ARMS at a site but on no ARMS outside a site.

    df_list<-list()
    for (i in 1:length(psList)){
        print(paste0("Getting endemics for replicate ", i))
        ps<-psList[[i]]
        
        #get site ids matching metabolite data as used by Aaron
            LongSiteNames<-sample_data(ps)$Site %>% as.vector
            LongSiteNames<-substr(LongSiteNames,1,nchar(LongSiteNames)-1)

            EmmaSiteNames<-substr(LongSiteNames,1,nchar(LongSiteNames)-1)
            DitaSiteNames<-substr(LongSiteNames,2,nchar(LongSiteNames))

            ShortSiteNames<-LongSiteNames
            ShortSiteNames[grepl("\\d", LongSiteNames)]<-EmmaSiteNames[grepl("\\d", LongSiteNames)] # replace those that still contain numbers (emma's) with emma's site names
            ShortSiteNames[!grepl("\\d", LongSiteNames)]<-DitaSiteNames[!grepl("\\d", LongSiteNames)] # replace those that no longer contain numbers (Dita's) with dita's site names

            sample_data(ps)$Site <-ShortSiteNames
        
        #make dataframe list to populate
            df_list[[i]]<-data.frame( row.names=sample_names(ps), ARMS=sample_data(ps)$ARMS, Site=sample_data(ps)$Site,  NumberEndemics=NA, Replicate=i)

        #get a list of endemic ESVs for each site (those occuring on more than 1 ARMS in a site but nowhere else)
            SiteEndemics_list<-list()
            for (site in unique(sample_data(ps)$Site) ) {

                if (sum(sample_data(ps)$Site==site)>1){

                    samplesFromSite<-sample_names(ps)[sample_data(ps)$Site==site]

                    ESVsOnMoreThan1ArmsAtSite<- taxa_names(ps)    [     rowSums(otu_table(ps)[,samplesFromSite]>0)>1        ]
                    
                    ESVsOnMoreThan1ArmsAtSiteTrueIfEndemic<-ps %>%
                        prune_taxa(ESVsOnMoreThan1ArmsAtSite, .) %>%
                        prune_samples(!sample_names(ps) %in% samplesFromSite, .) %>%
                        taxa_sums %>%
                        `<`(1) 

                    SiteEndemics_list[[site]]<-ESVsOnMoreThan1ArmsAtSite[ESVsOnMoreThan1ArmsAtSiteTrueIfEndemic]
                }
            }

        # populate dataframe list with number of endemic esvs per site        
            for (sample in sample_names(ps)) {
                
                ESVsPresentInSample<- taxa_names(ps)    [(otu_table(ps)[,sample]>0)]
                Site<-sample_data(ps)[sample,]$Site
                if (!is.null(SiteEndemics_list[[Site]])) {
                    df_list[[i]][sample,"NumberEndemics"]<-sum(ESVsPresentInSample %in%SiteEndemics_list[[Site]])
                }
            }

    }
    #combine list of dataframes into single dataframe
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
    SiteEndemicCOIDf<-GetSiteEndemicsDataframe(psCOI_rarefiedReplicate_l)

    ps16_rarefiedReplicate_l<-GetRarefiedReplicatePsList(psList=ps_l, Depth=Depths, Gene="16s", NumberOfReplicates=Replicates)
    Endemics16Df<-GetEndemicsDataframe(ps16_rarefiedReplicate_l)
    SiteEndemic16Df<-GetSiteEndemicsDataframe(ps16_rarefiedReplicate_l)

    EndemicCOIStats<-EndemicCOIDf %>%
            rowwise() %>%
            mutate(
                mean = mean(c_across(where(is.numeric))),
                SD = sd(c_across(where(is.numeric))),
                SE = SD/Replicates) %>%
            select(ARMS, mean, SD, SE)

   SiteEndemicCOIStats<-SiteEndemicCOIDf %>%
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

   SiteEndemic16Stats<-SiteEndemic16Df %>%
            rowwise() %>%
            mutate(
                mean = mean(c_across(where(is.numeric))),
                SD = sd(c_across(where(is.numeric))),
                SE = SD/Replicates) %>%
            select(ARMS, mean, SD, SE)

    EndemicityStats<-full_join(EndemicCOIStats,Endemic16Stats, by="ARMS", suffix=c("_COI", "_16s")) %>%
        mutate(NumberRarefyReplicates=Replicates)
    
    SiteEndemicityStats<-full_join(SiteEndemicCOIStats,SiteEndemic16Stats, by="ARMS", suffix=c("_COI", "_16s")) %>%
        mutate(NumberRarefyReplicates=Replicates)

    #save arms endemics
        write.csv(EndemicityStats, file.path("Outputs", "EndemicityStats.csv"))
        write.csv(EndemicCOIDf, file.path("Outputs", "endemics_COIRarefyReplicates.csv"))
        write.csv(Endemics16Df, file.path("Outputs", "endemics_16sRarefyReplicates.csv"))

    #save site endemics
        write.csv(SiteEndemicityStats, file.path("Outputs", "SiteEndemicityStats.csv"))
        write.csv(SiteEndemicCOIDf, file.path("Outputs", "SiteEndemics_COIRarefyReplicates.csv"))
        write.csv(SiteEndemic16Df, file.path("Outputs", "SiteEndemics_16sRarefyReplicates.csv"))
