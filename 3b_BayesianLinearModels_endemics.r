# 3b_endemics
# Fits Bayesian linear models, first tests relationships with correlations and frequetist linear models


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
    df<-read.csv(file.path("Outputs", "RichnessStats.csv")) %>%
        left_join(read.csv(file.path("Data", "MetaboliteData.csv"))) %>%
        left_join(read.csv(file.path("Data", "nceas_scores.csv"))) %>%
        as_tibble

    metabEndemics<-read.csv(file.path("Data", "Hartmann_metabolite_endemics.csv")) %>%
        select(ARMS, Metabolite_endemic_richness, site.endemic, region.endemic)

    df<-left_join(df, metabEndemics, by=c("Metab_Sample"= "ARMS"))

#2) Correlation matrices
  df %>% 
    select(mean_COI, mean_16s, Metabolite_endemic_richness, nceas_score)%>%
    GGally::ggpairs()

  ggsave("Figures/endemics_correlations.png", width = 8, height = 7)


#3) metab-esv: frequentist linear models (as a sanity check)
  plotList<-c(FitLinearModel(df, "Metabolite_endemic_richness", "mean_16s"),
    FitLinearModel(df, "Metabolite_endemic_richness", "mean_16s", EcoregionRandomIntercept=TRUE),
    FitLinearModel(df, "Metabolite_endemic_richness", "mean_COI"),
    FitLinearModel(df, "Metabolite_endemic_richness", "mean_COI", EcoregionRandomIntercept=TRUE))

  png("Figures/endemics_LinearModels.png", height = 16.6, width = 11.7, units = 'in', res = 300)
    cowplot::plot_grid(plotlist=plotList, nrow=4)
  dev.off()

#4) stress-metab endemics: frequentist linear model (exploratory only)

    lm(Metabolite_endemic_richness~nceas_score, df) %>% summary
    lmerTest::lmer(Metabolite_endemic_richness~nceas_score+(1|Ecoregion), df) %>% summary
    
    #random intercepts model
        OutputSubdirectory<-file.path("Outputs", "endemicmetab~stress_FrequentistMixedModel")
        dir.create(OutputSubdirectory)

        mod<-lmerTest::lmer(Metabolite_endemic_richness~nceas_score+(1|Ecoregion), df)
        png(file.path(OutputSubdirectory, "FrequentistMixedModel_qqplot.png"), height = 8.3, width = 11.7, units = 'in', res = 300)
            car::qqPlot(residuals(mod))
        dev.off()

        sjPlot::tab_model(mod,
                    pred.labels = c("Intercept", "nceas_score"),
                    dv.labels = c("FrquentistMixedModel_Stress-Metab"),
                    file=(file.path(OutputSubdirectory, "FrequentistMixedModel_Stress-Metab.html"))
        )

    #plot with random intercepts 

        textPositions1<-c(max(df[,"Metabolite_endemic_richness"], na.rm=T), min(df[,"nceas_score"], na.rm=T))
        textPositions2<-textPositions1-c(diff(range(df[,"Metabolite_endemic_richness"], na.rm=T))*0.1, 0)
        textPositions3<-textPositions1-c(diff(range(df[,"Metabolite_endemic_richness"], na.rm=T))*0.2, 0)

        df %>%
            filter(!is.na(Metabolite_endemic_richness))%>%
            cbind("prediction"=predict(mod))%>%
            ggplot(aes(y=Metabolite_endemic_richness, x=nceas_score))+
                #geom_smooth(method=lm, color="black")+
                geom_point(aes(color=Ecoregion))+
                geom_line(aes(y=prediction, color=Ecoregion))+
                annotate("text", label = paste0("R2m = ", signif( MuMIn::r.squaredGLMM(mod)[1], 2)), y = textPositions1[1], x=textPositions1[2], size = 6, color = "black",  hjust = "inward")+
                annotate("text", label = paste0("R2c = ", signif( MuMIn::r.squaredGLMM(mod)[2], 2)), y = textPositions2[1], x=textPositions2[2], size = 6, color = "black",  hjust = "inward")+
                annotate("text", label = paste0("p = ", signif(summary(mod)$coefficients[2,5], 2)), y = textPositions3[1], x=textPositions3[2], size = 6, color = "black",  hjust = "inward")

        ggsave("Figures/EndemicMetab~Stress_LinearModel.png", width = 8, height = 7)





# 5) bayesian models
    M5<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_endemic_richness",
        predictor = "mean_16s",
        measurementError= "SD_16s",
        ModelName="M5:16s-MetabEndemics",
        test=FALSE
    )
    saveRDS(M5, file.path("Outputs", "M5_Output.RDS"))

    M6<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_endemic_richness",
        predictor = "mean_COI",
        measurementError= "SD_COI",
        ModelName="M6:COI-MetabEndemics",
        test=FALSE
    )
    saveRDS(M6, file.path("Outputs", "M6_Output.RDS"))

# 6) save bayesian model input data
    saveRDS(df, file.path("Outputs", "endemics_BayesianModelInputData.RDS"))
