# 3_BayesianLinearModels
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
FitLinearModel<- function( data, outcome, predictor) {
  model_formula <- formula(paste0(outcome, " ~ ", predictor))
  fit = lm(model_formula, data = data)
  plot(fit, which=2)
  qqplot <- recordPlot() # assign plot

  textPositions1<-c(max(data[,outcome], na.rm=T), min(data[,predictor], na.rm=T))
  textPositions2<-textPositions1-c(diff(range(data[,outcome], na.rm=T))*0.1, 0)

  lmplot <- ggplot (data, aes(y=.data[[outcome]], x=.data[[predictor]])) + 
    geom_point(aes(color = factor(Ecoregion)), show.legend = FALSE) + 
    geom_smooth(method = lm, color = "black") +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab(predictor) +
    scale_color_manual(values = c("Aceh" = "green4", "Bali" = "firebrick2", "LI" = "black", "PaSer" = "gold2", "Samoa" = "darkorange2", "VIP" = "steelblue3"), 
                      labels = c("Aceh", "Bali", "Line Island", "Seribu Island", "Samoa", "Verde Island")) +
    annotate("text", label = paste0("R2 = ", signif(summary(fit)$r.squared,2)), y = textPositions1[1], x=textPositions1[2], size = 6, color = "black",  hjust = "inward")+
    annotate("text", label = paste0("p = ", signif(summary(fit)$coefficients[2,4], 2)), y = textPositions2[1], x=textPositions2[2], size = 6, color = "black",  hjust = "inward")

  return(list("qqplot"=qqplot, "lmplot"=lmplot))
}

FitBayesianLinearModel<-function( data, outcome, predictor, measurementError=FALSE, ModelName="NoName", test=FALSE) {

    #define inputs
        if (test) {
            iter=10000
            max_treedepth=10
        } else {
            iter=30000
            max_treedepth=20
        }

        if (measurementError==FALSE) {
            bform <- bf(paste0(outcome, " ~ ", predictor, " + (1|Ecoregion)"))
        } else {
            bform <- bf(paste0(outcome, " ~ me(", predictor, ", ", measurementError,") + (1|Ecoregion)")) + set_mecor(FALSE)
        }

    #make output directory
        OutputSubdirectory<-paste0("Outputs/ModelDiagnositics_",ModelName)
        dir.create(OutputSubdirectory)
    
    #fit model
        brms::get_prior(bform, data = df) %>%
            write.csv(file=file.path(OutputSubdirectory, paste0("Priors_",ModelName,".csv")))

        model <- brms::brm(bform,
                    data = df,
                    warmup = 1000, iter = iter,
                    cores = 2, chains = 2,
                    control = list(adapt_delta = 0.99, max_treedepth = max_treedepth))

    #save diagnostics and model table
        print(sjPlot::tab_model(model,
                    pred.labels = c("Intercept", predictor),
                    dv.labels = c(ModelName),
                    file=(file.path(OutputSubdirectory, paste0(ModelName,".html")))
        ))

        png(file.path(OutputSubdirectory, paste0(ModelName,"_TraceAndDensityPlots.png")), height = 8.3, width = 11.7, units = 'in', res = 300)
            plot(model)
        dev.off()

            plotList<-list()
            plotList[["PosteriorPredictiveCheckPlot"]]<-brms::pp_check(model)
            plotList[["MeanCapturePlot"]]<-brms::pp_check(model, type='stat', stat='mean')
            plotList[["PredictiveErrorPlot"]]<-brms::pp_check(model, type='error_scatter_avg')
            plotList[["PredictiveErrorByObservation"]]<-brms::pp_check(model, type='intervals')
            plot(brms::loo(model), main="Pareto smoothed importance sampling\n>0.7 a problem")
            plotList[["OutliersPlot"]]<- recordPlot()
            
        png(file.path(OutputSubdirectory, paste0(ModelName,"_OtherDiagnosticPlots.png")), height = 8.3, width = 11.7, units = 'in', res = 300)
            print({
                cowplot::plot_grid(plotlist=plotList, nrow=3)
            })
        dev.off()

    return(model)
}

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
