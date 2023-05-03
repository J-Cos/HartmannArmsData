# 3c_endemics vs endemics
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
    df<-read.csv(file.path("Outputs", "EndemicityStats.csv")) %>%
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

  ggsave("Figures/endemics_v_endemics_correlations.png", width = 8, height = 7)


#3) metab-esv: frequentist linear models (as a sanity check)
  plotList<-c(FitLinearModel(df, "Metabolite_endemic_richness", "mean_16s"),
    FitLinearModel(df, "Metabolite_endemic_richness", "mean_16s", EcoregionRandomIntercept=TRUE),
    FitLinearModel(df, "Metabolite_endemic_richness", "mean_COI"),
    FitLinearModel(df, "Metabolite_endemic_richness", "mean_COI", EcoregionRandomIntercept=TRUE))

  png("Figures/endemics_v_endemics_LinearModels.png", height = 16.6, width = 11.7, units = 'in', res = 300)
    cowplot::plot_grid(plotlist=plotList, nrow=4)
  dev.off()

# 4) bayesian models
    M7<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_endemic_richness",
        predictor = "mean_16s",
        measurementError= "SD_16s",
        ModelName="M7:16sEndemics-MetabEndemics",
        test=FALSE
    )
    saveRDS(M7, file.path("Outputs", "M7_Output.RDS"))

    M8<-FitBayesianLinearModel(
        data = df,
        outcome = "Metabolite_endemic_richness",
        predictor = "mean_COI",
        measurementError= "SD_COI",
        ModelName="M8:COIEndemics-MetabEndemics",
        test=FALSE
    )
    saveRDS(M8, file.path("Outputs", "M8_Output.RDS"))

# 5) save bayesian model input data
    saveRDS(df, file.path("Outputs", "endemics_v_endemics_BayesianModelInputData.RDS"))
