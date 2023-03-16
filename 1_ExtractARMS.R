#packages
    library(tidyverse)
    library(phyloseq)

#1) load data
    #get list of Aarons ARMS
        HartmannARMS<-read.csv(file.path("Data", "HartmannARMS.csv")) %>% select("Matching.ARMS") %>% unlist %>% as.vector
    #extract those arms only from COI and 16s datasets
        ps16Aaron<-readRDS(file=file.path("/home/j/Dropbox/CrossPacific_Paper", "Outputs", "ps16.RDS")) %>% 
            prune_samples(sample_data(.)$Fraction=="BS",.) %>%
            prune_samples(sample_data(.)$ARMS %in% HartmannARMS,.) %>%
            prune_taxa(taxa_sums(.)>0, .)
        psCoiAaron<-readRDS(file=file.path("/home/j/Dropbox/CrossPacific_Paper", "Outputs", "psCOI.RDS")) %>% 
            prune_samples(sample_data(psCOI)$Fraction=="Sessile",.) %>%
            prune_samples(sample_data(.)$ARMS %in% HartmannARMS,.)  %>%
            prune_taxa(taxa_sums(.)>0, .)
            
#2) save output
    saveRDS(list("16s"=ps16Aaron, "COI"=psCoiAaron), file.path("Outputs", "AaronsARMS.RDS"))