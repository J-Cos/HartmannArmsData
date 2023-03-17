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

#1) load data
    ps_l<-readRDS(file.path("Outputs", "AaronsARMS.RDS"))

# define parameters
    Depths<-list(   "COI"=39397,
                    "16s"=55296) #one fewer than minimum depth on each
    Replicates<-1000

#2) multirarefy
    #get sample depths
        MatchedCOIARMS_indices<-sample_data(ps_l[["COI"]])$ARMS %in% sample_data(ps_l[["16s"]])$ARMS
        NumberMatchedARMS<-sum(MatchedCOIARMS_indices)
        NonmathcedCoiArms<-sample_data(ps_l[["COI"]])$ARMS[!MatchedCOIARMS_indices]

    #get multirarefied dataframes
        RichnessCoiDf_list<-GetAaronsRichnessDataframe(Gene="COI", NumberOfReplicates=Replicates)
        Richness16sDf_list<-GetAaronsRichnessDataframe(Gene="16s", NumberOfReplicates=Replicates)

#3) reformat and output to csv
    RichnessStats<-left_join(RichnessCoiDf_list[["Stats"]],Richness16sDf_list[["Stats"]], by="ARMS", suffix=c("_COI", "_16s")) %>%
        mutate(NumberRarefyReplicates=Replicates)

    write.csv(RichnessStats, file.path("Outputs", "RichnessStats.csv"))
    write.csv(RichnessCoiDf_list[["Replicates"]], file.path("Outputs", "COIRarefyReplicates.csv"))
    write.csv(Richness16sDf_list[["Replicates"]], file.path("Outputs", "16sRarefyReplicates.csv"))
