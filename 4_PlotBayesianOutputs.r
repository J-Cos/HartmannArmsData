# 4_PlotBayesianOutputs
# Plots...


#packages
  library(tidyverse)
  library(ggplot2)
  library(tidybayes)
  library(brms)

#functions
MakePredictionPlot<-function(data, model, outcome, predictorGene, endemicLabelling="", gridAndDraws=100){

    sd<-paste0("SD_",predictorGene)
    mean<-paste0("mean_", predictorGene)
    if(predictorGene=="16s"){
        xlabel<-paste0(predictorGene, " rRNA gene ",endemicLabelling, "richness")
    } else if (predictorGene=="COI"){
        xlabel<-paste0(predictorGene, " gene ",endemicLabelling, "richness")
    }

    dataGrid<- data %>%
        modelr::data_grid(mean = modelr::seq_range(!! rlang::sym(mean), n=gridAndDraws), SD = modelr::seq_range(!! rlang::sym(sd), n=gridAndDraws), Ecoregion) 

    names(dataGrid)<-c(mean, sd, "Ecoregion")

    dataGrid %>%
        tidybayes::add_predicted_draws(model, ndraws=gridAndDraws) %>%
        ggplot(aes(y = .data[[outcome]], x = .data[[mean]])) + 
        stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .8, .5), color = "gray85", show.legend = FALSE) +
        scale_fill_grey(start = 0.9, end = 0.3) +
        ggnewscale::new_scale_fill() +
        geom_point(data = data, aes(color = factor(Ecoregion), shape = factor(Ecoregion), fill = factor(Ecoregion)), size = 3) + 
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab(xlabel) +
        ylab(str_replace_all(outcome, "_", " "))+
        customAes[["col"]]+customAes[["fill"]]+customAes[["shape"]]+
        theme(
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
}
MakeEcoregionCredibleIntervalPlot<-function(model, outcomeLabel, predictorLabel){
    model %>%
        spread_draws(b_Intercept, r_Ecoregion[Ecoregion,]) %>%
        mutate(Ecoregion_mean = b_Intercept + r_Ecoregion) %>%
        ggplot(aes(y = Ecoregion, x = Ecoregion_mean, fill = Ecoregion)) +
        scale_y_discrete(breaks = c("Aceh", "Bali", "LI", "PaSer", "Samoa", "VIP"), 
                        labels = (c("Aceh", "Bali", "Line Island", "Seribu Island", "Samoa", "Verde Island"))) +
        scale_x_continuous(name = paste0(outcomeLabel, " posterior distributions")) +
        stat_halfeye() +
        geom_vline(aes(xintercept = mean(b_Intercept)), linetype = "dashed") +
        scale_fill_manual(values = c("green4", "firebrick2", "gray30", "gold2", "darkorange2", "steelblue3"), 
                        labels = c("Aceh", "Bali", "Line Island", "Seribu Island", "Samoa", "Verde Island")) +
        theme(legend.position = "none") +
        ggtitle(paste0(outcomeLabel, " ~ ", predictorLabel)) +
        theme(  plot.title = element_text( size=10),,
                axis.title = element_text(size = 8))
}
#set parameters
    customAes<-list(
        "col"=scale_color_manual(values = c("Aceh" = "green4", "Bali" = "firebrick2", "LI" = "black", "PaSer" = "gold2", "Samoa" = "darkorange2", "VIP" = "steelblue3"), 
                      labels = c("Aceh", "Bali", "Line Island", "Seribu Island", "Samoa", "Verde Island"), name = "Region"),
        "fill"= scale_fill_manual(values = c("Aceh" = "green4", "Bali" = "firebrick2", "LI" = "black", "PaSer" = "gold2", "Samoa" = "darkorange2", "VIP" = "steelblue3"),  
                      labels = c("Aceh", "Bali", "Line Island", "Seribu Island", "Samoa", "Verde Island"), name = "Region"),
        "shape"= scale_shape_manual(values = c("Aceh"=21,"Bali"=24,"LI"=23, "PaSer"=21, "Samoa"=25, "VIP"=22),  
                      labels = c("Aceh", "Bali", "Line Island", "Seribu Island", "Samoa", "Verde Island"), name = "Region")
    )


# 1) load data
    df<-readRDS(file.path("Outputs", "BayesianModelInputData.RDS"))
    df_endemics<-readRDS(file.path("Outputs", "endemics_BayesianModelInputData.RDS"))
    df_endemics_v_endemics<-readRDS(file.path("Outputs", "endemics_v_endemics_BayesianModelInputData.RDS"))
    M1<-readRDS(file.path("Outputs", "M1_Output.RDS"))
    M2<-readRDS(file.path("Outputs", "M2_Output.RDS"))
    M3<-readRDS(file.path("Outputs", "M3_Output.RDS"))
    M4<-readRDS(file.path("Outputs", "M4_Output.RDS"))
    M5<-readRDS(file.path("Outputs", "M5_Output.RDS"))
    M6<-readRDS(file.path("Outputs", "M6_Output.RDS"))
    M7<-readRDS(file.path("Outputs", "M7_Output.RDS"))
    M8<-readRDS(file.path("Outputs", "M8_Output.RDS"))


# 2) plot posterior prediction curves
    MakePredictionPlot(data=df, model=M1, outcome="Metabolite_richness", predictorGene="16s")
    ggsave(filename="Figures/M1_plot.png", width = 6, height = 6)

    MakePredictionPlot(data=df, model=M2, outcome="Metabolite_diversity", predictorGene="16s")
    ggsave(filename="Figures/M2_plot.png", width = 6, height = 6)

    MakePredictionPlot(data=df, model=M3, outcome="Metabolite_richness", predictorGene="COI")
    ggsave(filename="Figures/M3_plot.png", width = 6, height = 6)

    MakePredictionPlot(data=df, model=M4, outcome="Metabolite_diversity", predictorGene="COI")
    ggsave(filename="Figures/M4_plot.png", width = 6, height = 6)

    MakePredictionPlot(data=df_endemics, model=M5, outcome="Metabolite_endemic_richness", predictorGene="16s")
    ggsave(filename="Figures/M5_plot.png", width = 6, height = 6)

    MakePredictionPlot(data=df_endemics, model=M6, outcome="Metabolite_endemic_richness", predictorGene="COI")
    ggsave(filename="Figures/M6_plot.png", width = 6, height = 6)

    MakePredictionPlot(data=df_endemics_v_endemics, model=M7, outcome="Metabolite_endemic_richness", endemicLabelling="endemic ", predictorGene="16s")
    ggsave(filename="Figures/M7_plot.png", width = 6, height = 6)

    MakePredictionPlot(data=df_endemics_v_endemics, model=M8, outcome="Metabolite_endemic_richness", endemicLabelling="endemic ", predictorGene="COI")
    ggsave(filename="Figures/M8_plot.png", width = 6, height = 6)

# 3) plot credible intervals per ecoregion
    MakeEcoregionCredibleIntervalPlot(M1, "Metabolite richness", "16s rRNA gene richness")
    ggsave(filename="Figures/M1_CredibleIntervalPlot.png", width = 6, height = 6)

    MakeEcoregionCredibleIntervalPlot(M2, "Metabolite diversity", "16s rRNA gene richness")
    ggsave(filename="Figures/M2_CredibleIntervalPlot.png", width = 6, height = 6)

    MakeEcoregionCredibleIntervalPlot(M3, "Metabolite richness", "COI gene richness")
    ggsave(filename="Figures/M3_CredibleIntervalPlot.png", width = 6, height = 6)

    MakeEcoregionCredibleIntervalPlot(M4, "Metabolite diversity", "COI gene richness")
    ggsave(filename="Figures/M4_CredibleIntervalPlot.png", width = 6, height = 6)
    
    MakeEcoregionCredibleIntervalPlot(M5, "Metabolite endemic richness", "COI gene richness")
    ggsave(filename="Figures/M5_CredibleIntervalPlot.png", width = 6, height = 6)

    MakeEcoregionCredibleIntervalPlot(M6, "Metabolite endemic richness", "COI gene richness")
    ggsave(filename="Figures/M6_CredibleIntervalPlot.png", width = 6, height = 6)

    MakeEcoregionCredibleIntervalPlot(M7, "Metabolite endemic richness", "16s rRNA gene endemic richness")
    ggsave(filename="Figures/M7_CredibleIntervalPlot.png", width = 6, height = 6)

    MakeEcoregionCredibleIntervalPlot(M8, "Metabolite endemic richness", "COI gene endemic richness")
    ggsave(filename="Figures/M8_CredibleIntervalPlot.png", width = 6, height = 6)