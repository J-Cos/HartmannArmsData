
MakeModelSelectionTableForPrinting<-function(modelList ){

        b<-unlist(lapply(modelList, MuMIn::AICc))
        names(b)<-1:length(b)
        AICorder<-as.integer(names(sort(b)))
        R2s<-signif(unlist(lapply(modelList, MuMIn::r.squaredGLMM))[c(TRUE, FALSE)], 2)[AICorder]


        table<-MuMIn::model.sel(modelList, rank="AICc") %>% 
          as.data.frame %>%
          select(-`(Intercept)`, -logLik, -AICc) %>%
          mutate(R2= round(R2s, 2))%>%
          mutate(weight=round(weight*100, 2) ) %>%
          mutate(delta=round(delta, 2))
        
        if ("family" %in% names(table)) {
          table<-select(table, -family)
        }
        
        names(table)[4:length(names(table))]<-c( "K", "deltaAICc", "AICc weight", "marginal R2")

    return(table)
  }