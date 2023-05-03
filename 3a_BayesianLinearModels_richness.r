# 3a_BayesianLinearModels
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

#2) Correlation matrices
  df %>% 
    select(mean_COI, mean_16s, Metabolite_diversity, Metabolite_richness, nceas_score)%>%
    GGally::ggpairs()

  ggsave("Figures/correlations.png", width = 8, height = 7)

#3) metab-esv: frequentist linear models (as a sanity check)
  plotList<-c(FitLinearModel(df, "Metabolite_richness", "mean_16s"),
    FitLinearModel(df, "Metabolite_diversity", "mean_16s"),
    FitLinearModel(df, "Metabolite_richness", "mean_COI"),
    FitLinearModel(df, "Metabolite_diversity", "mean_COI"))

  png("Figures/LinearModels.png", height = 8.3, width = 11.7, units = 'in', res = 300)
    cowplot::plot_grid(plotlist=plotList, nrow=2)
  dev.off()


#4) stress-metab: frequentist linear model (exploratory only)

    lm(Metabolite_richness~nceas_score, df) %>% summary
    #tangetial check of microbe macrobe patterns
        lm(mean_16s~nceas_score, df) %>% summary
        lm(mean_COI~nceas_score, df) %>% summary
        lm(mean_16s~mean_COI, df) %>% summary

        lmerTest::lmer(mean_16s~nceas_score+(1|Ecoregion), df) %>% summary
        lmerTest::lmer(mean_COI~nceas_score+(1|Ecoregion), df) %>% summary
        lmerTest::lmer(mean_16s~mean_COI+(1|Ecoregion), df) %>% summary


    #random intercepts model
        OutputSubdirectory<-file.path("Outputs", "FrequentistMixedModel")
        dir.create(OutputSubdirectory)

        mod<-lmerTest::lmer(Metabolite_richness~nceas_score+(1|Ecoregion), df)
        png(file.path(OutputSubdirectory, "FrequentistMixedModel_qqplot.png"), height = 8.3, width = 11.7, units = 'in', res = 300)
            car::qqPlot(residuals(mod))
        dev.off()

        sjPlot::tab_model(mod,
                    pred.labels = c("Intercept", "nceas_score"),
                    dv.labels = c("FrquentistMixedModel_Stress-Metab"),
                    file=(file.path(OutputSubdirectory, "FrequentistMixedModel_Stress-Metab.html"))
        )

    #plot with random intercepts 

        textPositions1<-c(max(df[,"Metabolite_richness"], na.rm=T), min(df[,"nceas_score"], na.rm=T))
        textPositions2<-textPositions1-c(diff(range(df[,"Metabolite_richness"], na.rm=T))*0.1, 0)
        textPositions3<-textPositions1-c(diff(range(df[,"Metabolite_richness"], na.rm=T))*0.2, 0)

        df %>%
            filter(!is.na(Metabolite_richness))%>%
            cbind("prediction"=predict(mod))%>%
            ggplot(aes(y=Metabolite_richness, x=nceas_score))+
                #geom_smooth(method=lm, color="black")+
                geom_point(aes(color=Ecoregion))+
                geom_line(aes(y=prediction, color=Ecoregion))+
                annotate("text", label = paste0("R2m = ", signif( MuMIn::r.squaredGLMM(mod)[1], 2)), y = textPositions1[1], x=textPositions1[2], size = 6, color = "black",  hjust = "inward")+
                annotate("text", label = paste0("R2c = ", signif( MuMIn::r.squaredGLMM(mod)[2], 2)), y = textPositions2[1], x=textPositions2[2], size = 6, color = "black",  hjust = "inward")+
                annotate("text", label = paste0("p = ", signif(summary(mod)$coefficients[2,5], 2)), y = textPositions3[1], x=textPositions3[2], size = 6, color = "black",  hjust = "inward")

        ggsave("Figures/StressLinearModel.png", width = 8, height = 7)


    #bayesian plot exploration
        StressBayesianModel <- brms::brm(Metabolite_richness~nceas_score+(1|Ecoregion),
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

# 5)  metab-esv: bayesian linear models
    M1<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_richness",
        predictor = "mean_16s",
        measurementError= "SD_16s",
        ModelName="M1:16s-MetabRichness",
        test=FALSE
    )
    saveRDS(M1, file.path("Outputs", "M1_Output.RDS"))

    M2<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_diversity",
        predictor = "mean_16s",
        measurementError= "SD_16s",
        ModelName="M2:16s-MetabDiversity",
        test=FALSE
    )
    saveRDS(M2, file.path("Outputs", "M2_Output.RDS"))

    M3<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_richness",
        predictor = "mean_COI",
        measurementError= "SD_COI",
        ModelName="M3:COI-MetabRichness",
        test=FALSE
    )
    saveRDS(M3, file.path("Outputs", "M3_Output.RDS"))

    M4<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_diversity",
        predictor = "mean_COI",
        measurementError= "SD_COI",
        ModelName="M4:COI-MetabDiversity",
        test=FALSE
    )
    saveRDS(M4, file.path("Outputs", "M4_Output.RDS"))


# 6) save bayesian model input data
    saveRDS(df, file.path("Outputs", "BayesianModelInputData.RDS"))
