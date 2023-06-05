# 3d_SiteEndemics vs Site_Endemics
# Frequentist model selection


#packages
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(rstan)
  library(brms)
  library(tidybayes)
  library(ggridges)
  #library(ggstance)
  #library(ggmcmc)
  library(RColorBrewer)
  library(cowplot)
  library(modelr)
  library(GGally)
  library(sjPlot)
  library(ggnewscale)

  library(lme4)
  library(lmerTest)

#functions
    source(file.path("Code", "Functions", "Functions.R"))

#1) load data
    df<-read.csv(file.path("Outputs", "SiteEndemicityStats.csv")) %>%
        left_join(read.csv(file.path("Data", "MetaboliteData.csv"))) %>%
        left_join(read.csv(file.path("Data", "nceas_scores.csv"))) %>%
        as_tibble

    metabEndemics<-read.csv(file.path("Data", "Hartmann_metabolite_endemics.csv")) %>%
        select(ARMS, site.endemic)

    df<-left_join(df, metabEndemics, by=c("Metab_Sample"= "ARMS"))

#2) Correlation matrices
  df %>% 
    select(mean_COI, mean_16s, site.endemic, nceas_score)%>%
    GGally::ggpairs()

  ggsave("Figures/SiteEndemics_v_SiteEndemics_correlations.png", width = 8, height = 7)


#3) metab-esv: frequentist linear models (as a sanity check)
  plotList<-c(FitLinearModel(df, "site.endemic", "mean_16s"),
    FitLinearModel(df, "site.endemic", "mean_16s", EcoregionRandomIntercept=TRUE),
    FitLinearModel(df, "site.endemic", "mean_COI"),
    FitLinearModel(df, "site.endemic", "mean_COI", EcoregionRandomIntercept=TRUE))

  png("Figures/SiteEndemics_v_SiteEndemics_LinearModels.png", height = 16.6, width = 11.7, units = 'in', res = 300)
    cowplot::plot_grid(plotlist=plotList, nrow=4)
  dev.off()

#4) stress-metab: frequentist linear model

  #both genes
    dat<-df %>%
      select(mean_COI, mean_16s, site.endemic, nceas_score) %>%
      filter(complete.cases(.))
    
    dat %>% PlotVIF("site.endemic")

    MuMIn::dredge(lm(site.endemic ~ ., data = dat, na.action = "na.fail"))

    # mixed model
    dat<-df %>%
      select(mean_COI, mean_16s, site.endemic, nceas_score, Ecoregion) %>%
      filter(complete.cases(.))  #%>%
      #mutate_at(c("mean_16s", "mean_COI", "site.endemic", "nceas_score"), ~(scale(.) %>% as.vector))

      fm1<-glmer.nb(site.endemic~ nceas_score + (1|Ecoregion), data = dat)
      fm2<-glmer.nb(site.endemic~ mean_16s + (1|Ecoregion), data = dat)
      fm3<-glmer.nb(site.endemic~ mean_COI + (1|Ecoregion), data = dat)
      fm4<-glmer.nb(site.endemic~ nceas_score + mean_16s + (1|Ecoregion), data = dat)
      fm5<-glmer.nb(site.endemic~ nceas_score + mean_COI + (1|Ecoregion), data = dat)
      fm6<-glmer.nb(site.endemic~ mean_COI + mean_16s + (1|Ecoregion), data = dat)
      fm7<-glmer.nb(site.endemic~ nceas_score + mean_COI + mean_16s + (1|Ecoregion), data = dat)
      fm8<-glmer.nb(site.endemic~ (1|Ecoregion), data = dat)

      MuMIn::model.sel(list(fm1, fm2, fm3, fm4, fm5, fm6, fm7, fm8)) %>% 
      as.data.frame %>%
      write.csv("Outputs/FrequentistModelSelection_SiteEndemics_BothGenes.csv")
      car::qqPlot(residuals(fm1))

  #COI only
    dat<-df %>%
      select(mean_COI, site.endemic, nceas_score) %>%
      filter(complete.cases(.))

    dat %>% PlotVIF(outcome="site.endemic")

    MuMIn::dredge(lm(site.endemic ~ ., data = dat, na.action = "na.fail"))

    # mixed model
    dat<-df %>%
      select(mean_COI, site.endemic , nceas_score, Ecoregion) %>%
      filter(complete.cases(.)) #%>%
      #mutate(site.endemic=log(site.endemic+0.00000001))
      #mutate_at(c("mean_COI", "site.endemic", "nceas_score"), ~(scale(.) %>% as.vector))

    fm1<-glmer.nb(site.endemic ~ nceas_score + (1|Ecoregion), data = dat, family="Gamma")
    fm2<-glmer.nb(site.endemic ~ mean_COI + (1|Ecoregion), data = dat)
    fm3<-glmer.nb(site.endemic ~ nceas_score + mean_COI + (1|Ecoregion), data = dat)
    fm4<-glmer.nb(site.endemic ~ (1|Ecoregion), data = dat)

    MuMIn::model.sel(list(fm1, fm2, fm3, fm4), rank="AICc") %>% 
      as.data.frame %>%
      write.csv("Outputs/FrequentistModelSelection_SiteEndemics_COI.csv")

    car::qqPlot(residuals(fm1))


  #16s only 
    dat<-df %>%
      select(mean_16s, site.endemic, nceas_score) %>%
      filter(complete.cases(.))
  
    dat %>% PlotVIF(outcome="site.endemic")

    MuMIn::dredge(lm(Metabolite_richness ~ ., data = dat, na.action = "na.fail"))

    # mixed model
    dat<-df %>%
      select(mean_16s, site.endemic, nceas_score, Ecoregion) %>%
      filter(complete.cases(.)) #%>%
      #mutate_at(c("mean_16s", "Metabolite_richness", "nceas_score"), ~(scale(.) %>% as.vector))

    fm1<-glmer.nb(site.endemic~ nceas_score + (1|Ecoregion), data = dat)
    fm2<-glmer.nb(site.endemic~ mean_16s + (1|Ecoregion), data = dat)
    fm3<-glmer.nb(site.endemic~ nceas_score + mean_16s + (1|Ecoregion), data = dat)
    fm4<-glmer.nb(site.endemic~ (1|Ecoregion), data = dat)

    MuMIn::model.sel(list(fm1, fm2, fm3, fm4)) %>% 
      as.data.frame %>%
      write.csv("Outputs/FrequentistModelSelection_SiteEndemics_16s.csv")
    car::qqPlot(residuals(fm2))


  #fit best model to all data
    mod<-fm2
      summary(mod)

        OutputSubdirectory<-file.path("Outputs", "FrequentistMixedModel_SiteEndemics")
        dir.create(OutputSubdirectory)


        sjPlot::tab_model(mod,
                    dv.labels = c("AveragedBestModel_SiteEndemics"),
                    file=(file.path(OutputSubdirectory, "BestModel_SiteEndemics.html"))
        )

    #plot with random intercepts 

        textPositions1<-c(max(dat[,"site.endemic"], na.rm=T), min(dat[,"nceas_score"], na.rm=T))
        textPositions2<-textPositions1-c(diff(range(dat[,"site.endemic"], na.rm=T))*0.1, 0)
        textPositions3<-textPositions1-c(diff(range(dat[,"site.endemic"], na.rm=T))*0.2, 0)

        dat %>%
            filter(!is.na(site.endemic))%>%
            cbind("prediction"=predict(mod, type="response"))%>%
            ggplot(aes(y=site.endemic, x=mean_16s))+
                #geom_smooth(method=lm, color="black")+
                geom_point(aes(color=Ecoregion))+
                geom_line(aes(y=prediction, color=Ecoregion))+
                annotate("text", label = paste0("R2m = ", signif( MuMIn::r.squaredGLMM(mod)[1,1], 2)), y = textPositions1[1], x=textPositions1[2], size = 6, color = "black",  hjust = "inward")+
                annotate("text", label = paste0("R2c = ", signif( MuMIn::r.squaredGLMM(mod)[1,2], 2)), y = textPositions2[1], x=textPositions2[2], size = 6, color = "black",  hjust = "inward")+
                annotate("text", label = paste0("p = ", signif(summary(mod)$coefficients[2,4], 2)), y = textPositions3[1], x=textPositions3[2], size = 6, color = "black",  hjust = "inward")

        ggsave("Figures/BestFrequentistModel_SiteEndemics.png", width = 8, height = 7)

  #bayesian version of the best model
        StressBayesianModel <- brms::brm(site.endemic ~mean_16s+(1|Ecoregion),
                    data = df,
                    warmup = 1000, iter = 10000,
                    cores = 2, chains = 2,
                    control = list(adapt_delta = 0.99, max_treedepth = 10))
        summary(StressBayesianModel)
        sjPlot::tab_model(StressBayesianModel,
                    pred.labels = c("Intercept", "nceas_score"),
                    dv.labels = c("BayesianMixedModel_Stress-Metab"),
                    file=(file.path(OutputSubdirectory, "BayesianMixedModel_Stress-Metab.html"))
        )
