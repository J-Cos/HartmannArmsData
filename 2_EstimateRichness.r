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

    MatchAaronsSiteNames<-function(ps){
        #get site ids matching metabolite data as used by Aaron
            LongSiteNames<-sample_data(ps)$Site %>% as.vector
            LongSiteNames<-substr(LongSiteNames,1,nchar(LongSiteNames)-1)

            EmmaSiteNames<-substr(LongSiteNames,1,nchar(LongSiteNames)-1)
            DitaSiteNames<-substr(LongSiteNames,2,nchar(LongSiteNames))

            ShortSiteNames<-LongSiteNames
            ShortSiteNames[grepl("\\d", LongSiteNames)]<-EmmaSiteNames[grepl("\\d", LongSiteNames)] # replace those that still contain numbers (emma's) with emma's site names
            ShortSiteNames[!grepl("\\d", LongSiteNames)]<-DitaSiteNames[!grepl("\\d", LongSiteNames)] # replace those that no longer contain numbers (Dita's) with dita's site names

            sample_data(ps)$Site <-ShortSiteNames
        return(ps)
    }
    
    GetRarefiedReplicatePsList<-function(psList=ps_l, Depth=Depths, Gene, NumberOfReplicates){

            smallestsamples<-psList[[Gene]] %>%sample_sums %>% sort %>% head


            if (as.vector(smallestsamples[1])==(Depth[[Gene]]+1)) {
                ps_rarefiedReplicate_l<-list()
                for (i in 1:NumberOfReplicates) {
                            ps_rarefiedReplicate_l[[i]]<-psList[[Gene]] %>%
                                rarefy_even_depth(  physeq=., sample.size = Depth[[Gene]], rngseed = i, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) 
                }
            } else{
                print(paste0("Warning: ", Gene, " data has changed"))
            }

            return(ps_rarefiedReplicate_l)
    }

    PrunePhyloseqToARMsForEndemicAnalysisfunction<-function(ps_l_item, DesiredArmsForEndemicsAnalysis){ 
        ps_l_item<- ps_l_item %>%
                    MatchAaronsSiteNames %>%
                    prune_samples(sample_data(.)$ARMS %in% DesiredArmsForEndemicsAnalysis, .) 
        
        SiteWith3Arms<-names(which(table(sample_data(ps_l_item)$Site)==3))
        ps_l_item<-prune_samples(sample_data(ps_l_item)$Site %in% SiteWith3Arms, ps_l_item) %>%
            prune_taxa(taxa_sums(.)>0, .)                                  
    }    

    GetEndemicsDataframe<-function(psList){
        df_list<-list()
        for (i in 1:length(psList)){
            print(paste0("Getting ARMS endemics for replicate ", i))
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
        #requires that all sites have exactly three arms! these should be aligned with the arms used in the metabolite endemics calculation.
        #gets the number of site endemic ESVs per ARMS. These are definied as ESVs occuring on at least 2 ARMS at a site but on no ARMS outside a site.

        df_list<-list()
        for (i in 1:length(psList)){
            print(paste0("Getting site endemics for replicate ", i))
            ps<-psList[[i]]
                        
            #make dataframe list to populate
                df_list[[i]]<-data.frame( row.names=sample_names(ps), ARMS=sample_data(ps)$ARMS, Site=sample_data(ps)$Site,  NumberEndemics=NA, Replicate=i)

            #get a list of endemic ESVs for each site (those occuring on more than 1 ARMS in a site but nowhere else)
                SiteEndemics_list<-list()
                for (site in unique(sample_data(ps)$Site) ) {

                    if (sum(sample_data(ps)$Site==site)==3){

                        samplesFromSite<-sample_names(ps)[sample_data(ps)$Site==site]

                        ESVsOnMoreThan1ArmsAtSite<- taxa_names(ps)    [     rowSums(otu_table(ps)[,samplesFromSite]>0)>1        ]
                        
                        ESVsOnMoreThan1ArmsAtSiteTrueIfEndemic<-ps %>%
                            prune_taxa(ESVsOnMoreThan1ArmsAtSite, .) %>%
                            prune_samples(!sample_names(ps) %in% samplesFromSite, .) %>%
                            taxa_sums %>%
                            `<`(1) 

                        SiteEndemics_list[[site]]<-ESVsOnMoreThan1ArmsAtSite[ESVsOnMoreThan1ArmsAtSiteTrueIfEndemic]
                    } else {
                        print("ARMS not filtered appropriately for endemics analysis")
                        return(NULL)
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

    GetEndemicsStatisticalSummary<-function(EndemicDf){
        EndemicDf %>%
            rowwise() %>%
            mutate(
                mean = mean(c_across(where(is.numeric))),
                SD = sd(c_across(where(is.numeric))),
                SE = SD/Replicates) %>%
            select(ARMS, mean, SD, SE)
    }

#1) load data
    ps_l<-readRDS(file.path("Outputs", "AaronsARMS.RDS"))

# define parameters
    Depths<-list(   "COI"=39397,
                    "16s"=55296) #one fewer than minimum depth on each
    Replicates<-100

#2) save community matrices
    psc<-prune_samples(sample_data(ps_l[["COI"]])$ARMS %in% sample_data(ps_l[["16s"]])$ARMS, ps_l[["COI"]])
    ps16<-prune_samples(sample_data(ps_l[["16s"]])$ARMS %in% sample_data(ps_l[["COI"]])$ARMS, ps_l[["16s"]])

    sample_sums(psc) %>% sort
    sample_sums(ps16) %>% sort
    
    pscr<-rarefy_even_depth(psc)
    ps16r<-rarefy_even_depth(ps16)

    GetCommMatNamedByARMS<-function(ps) {
        cm<-otu_table(ps) 
        colnames(cm)<-    as.vector(sample_data(ps)[,"ARMS"])$ARMS
        return(cm)
    }

    write.csv(GetCommMatNamedByARMS(pscr), "Outputs/CommunityMatrix_COI.csv")
    write.csv(GetCommMatNamedByARMS(ps16r), "Outputs/CommunityMatrix_16s.csv")

#3) multirarefy

    #get sample depths
        MatchedCOIARMS_indices<-sample_data(ps_l[["COI"]])$ARMS %in% sample_data(ps_l[["16s"]])$ARMS
        NumberMatchedARMS<-sum(MatchedCOIARMS_indices)
        NonmathcedCoiArms<-sample_data(ps_l[["COI"]])$ARMS[!MatchedCOIARMS_indices]

    #get multirarefied dataframes
        RichnessCoiDf_list<-GetAaronsRichnessDataframe(Gene="COI", NumberOfReplicates=Replicates)
        Richness16sDf_list<-GetAaronsRichnessDataframe(Gene="16s", NumberOfReplicates=Replicates)

#4) reformat and output to csv
    RichnessStats<-full_join(RichnessCoiDf_list[["Stats"]],Richness16sDf_list[["Stats"]], by="ARMS", suffix=c("_COI", "_16s")) %>%
        mutate(NumberRarefyReplicates=Replicates)

    write.csv(RichnessStats, file.path("Outputs", "RichnessStats.csv"))
    write.csv(RichnessCoiDf_list[["Replicates"]], file.path("Outputs", "COIRarefyReplicates.csv"))
    write.csv(Richness16sDf_list[["Replicates"]], file.path("Outputs", "16sRarefyReplicates.csv"))

# 5) get endemics and output to csv
    #prune to endemics
        DesiredArmsForEndemicsAnalysis<-c(
            "ARSG1A", "ARSG1B", "ARSG1C", 
            "ASEU1A", "ASEU1B", "ASEU1C", 
            "BCEN1A", "BCEN1B", "BCEN1C",
            "BHOR1A", "BHOR1B", "BHOR1C",
            "BNAP1A", "BNAP1B", "BNAP1C",
            "JAR11", "JAR04", "JAR08",
            "KIN13A", "KIN13B", "KIN13C",
            "SKOT1A", "SKOT1B", "SKOT1C",
            "SPRM1A", "SPRM1B", "SPRM1C",
            "OFU04", "OFU03B", "OFU03C",
            "ROS04", "ROS19", "ROS25",
            "TUT01", "TUT11C", "TUT22",
            "ART01A", "ART01B", "ART02A",
            "BAA01A", "BAA01B", "BAA02A",
            "BAO01B", "BAO01C", "BAO02A",
            "TWI01A", "TWI01B", "TWI01C"
        )

        ps_endemic_l<-lapply( ps_l, PrunePhyloseqToARMsForEndemicAnalysisfunction, DesiredArmsForEndemicsAnalysis)

    #get rarefying depths
        lapply(lapply(ps_endemic_l, sample_sums), min)
        Depths<-list(   "COI"=39397,
                        "16s"=60642) #one fewer than minimum depth on each

    #muti-rarefy and calculate endemics
        psCOI_rarefiedReplicate_l<-GetRarefiedReplicatePsList(psList=ps_endemic_l, Depth=Depths, Gene="COI", NumberOfReplicates=Replicates)
        EndemicCOIDf<-GetEndemicsDataframe(psCOI_rarefiedReplicate_l)
        SiteEndemicCOIDf<-GetSiteEndemicsDataframe(psCOI_rarefiedReplicate_l)

        ps16_rarefiedReplicate_l<-GetRarefiedReplicatePsList(psList=ps_endemic_l, Depth=Depths, Gene="16s", NumberOfReplicates=Replicates)
        Endemics16Df<-GetEndemicsDataframe(ps16_rarefiedReplicate_l)
        SiteEndemic16Df<-GetSiteEndemicsDataframe(ps16_rarefiedReplicate_l)

        EndemicCOIStats<-GetEndemicsStatisticalSummary(EndemicCOIDf)
        SiteEndemicCOIStats<-GetEndemicsStatisticalSummary(SiteEndemicCOIDf) 
        Endemic16Stats<-GetEndemicsStatisticalSummary(Endemics16Df)
        SiteEndemic16Stats<-GetEndemicsStatisticalSummary(SiteEndemic16Df)

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

