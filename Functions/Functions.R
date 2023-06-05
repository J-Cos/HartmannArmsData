FitLinearModel<- function( data, outcome, predictor, EcoregionRandomIntercept=FALSE) {
    if(!EcoregionRandomIntercept) {
        model_formula <- formula(paste0(outcome, " ~ ", predictor))
        fit = lm(model_formula, data = data)
        GetRsquared<-function(fit){ signif(summary(fit)$r.squared,2)}
        GetPvalue<-function(fit){signif(summary(fit)$coefficients[2,4], 2)}
    } else if (EcoregionRandomIntercept) {
        model_formula <- formula(paste0(outcome, " ~ ", predictor, " + (1|Ecoregion)"))
        fit = lmerTest::lmer(model_formula, data = data)
        GetRsquared<-function(fit){ signif(    MuMIn::r.squaredGLMM(fit)[1], 2) }
        GetPvalue<-function(fit){signif(summary(fit)$coefficients[2,5], 2)}

    }

  car::qqPlot(residuals(fit))
  qqplot <- recordPlot() # assign plot

  textPositions1<-c(max(data[,outcome], na.rm=T), min(data[,predictor], na.rm=T))
  textPositions2<-textPositions1-c(diff(range(data[,outcome], na.rm=T))*0.1, 0)


  lmplot <- ggplot (data, aes(y=.data[[outcome]], x=.data[[predictor]])) + 
    geom_point(aes(color = factor(Ecoregion)), show.legend = FALSE) + 
    #geom_smooth(method = lm, color = "black") +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab(predictor) +
    scale_color_manual(values = c("Aceh" = "green4", "Bali" = "firebrick2", "LI" = "black", "PaSer" = "gold2", "Samoa" = "darkorange2", "VIP" = "steelblue3"), 
                      labels = c("Aceh", "Bali", "Line Island", "Seribu Island", "Samoa", "Verde Island")) +
    annotate("text", label = paste0("R2 = ", GetRsquared(fit)), y = textPositions1[1], x=textPositions1[2], size = 6, color = "black",  hjust = "inward")+
    annotate("text", label = paste0("p = ", GetPvalue(fit)), y = textPositions2[1], x=textPositions2[2], size = 6, color = "black",  hjust = "inward")

    if(!EcoregionRandomIntercept) {
        lmplot<- lmplot + geom_smooth(method = lm, color = "black")
    } else if (EcoregionRandomIntercept) {
        predictionData<-data %>%
            filter(!is.na(get(outcome))& !is.na(get(predictor))) %>%
            mutate(prediction=predict(fit))

        lmplot<- lmplot + 
            geom_line(data=predictionData, aes(y=prediction, x=.data[[predictor]], color=Ecoregion))

    }

  return(list("qqplot"=qqplot, "lmplot"=lmplot))
}

PlotVIF<-function(dat, outcome="Metabolite_richness"){
    ## variance inflation factor
    model_formula <- formula(paste0(outcome, " ~ ."))

    model_all <- lm(model_formula, data=dat)  # with all the independent variables in the dataframe
    summary(model_all)
    vif_values <-car::vif(model_all)
    barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue") #create horizontal bar chart to display each VIF value
    abline(v = 5, lwd = 3, lty = 2)    #add vertical line at 5 as after 5 there is severe correlation
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