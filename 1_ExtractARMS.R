#packages
    library(tidyverse)
    library(phyloseq)

#functions
GetSmallerDuplicateSampleNames<-function(ps){   
    duplicateARMS<-sample_data(ps)$ARMS[sample_data(ps)$ARMS%>% duplicated]
    smallerDuplicates<-c()
    for (ARMS in duplicateARMS){
        duplicateSamples<-sample_names(ps)[sample_data(ps)$ARMS ==ARMS]
        largestDuplicate<-which.max(sample_sums(ps)[duplicateSamples])
        smallerDuplicates<-c(smallerDuplicates, duplicateSamples[-largestDuplicate])
    }
    return(smallerDuplicates)
}

PrunePhyloseqToAaronsData<-function(phyloseqRDS){
    ps<-readRDS(file=phyloseqRDS) %>% 
            prune_samples(sample_data(.)$Fraction=="BS" | sample_data(.)$Fraction=="Sessile",.) %>% # "BS" in 16s and "Sessile" in COI - Aarons data is all sessile
            prune_samples(sample_data(.)$ARMS %in% HartmannARMS,.) %>%
            prune_samples(!sample_names(.) %in% GetSmallerDuplicateSampleNames(.),.) %>%
            prune_taxa(taxa_sums(.)>0, .)
    return(ps)
}


#1) load data
    #get list of Aarons ARMS
        HartmannARMS<-read.csv(file.path("Data", "HartmannARMS.csv")) %>% select("Matching.ARMS") %>% unlist %>% as.vector

#2) extract those arms only from COI and 16s datasets
        ps16Aaron<-PrunePhyloseqToAaronsData(phyloseqRDS=file.path("/home/j/Dropbox/CrossPacific_Paper", "Outputs", "ps16.RDS"))
        psCoiAaron<-PrunePhyloseqToAaronsData(phyloseqRDS=file.path("/home/j/Dropbox/CrossPacific_Paper", "Outputs", "psCOI.RDS"))

#3) save output
    saveRDS(list("16s"=ps16Aaron, "COI"=psCoiAaron), file.path("Outputs", "AaronsARMS.RDS"))