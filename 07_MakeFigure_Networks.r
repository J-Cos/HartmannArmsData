# desired clustering at which to generate networks
# either ASV (set to "") or cOTU (set to "_cOTU")
cluster<-c("")
#cluster<-c("_cOTUs")

# packages and functions
    library(tidyverse)
    library(gridExtra)
    library(SpiecEasi)

getMultiNetworkDataframe<-function(se, matrices){
        
    #### MB method ####
    # Symmetrize the MB assymetric graph by taking the maximum value of the edge pairs e_ij / e_ji
    sebeta <- symBeta(getOptBeta(se), mode='maxabs')
    elist <- Matrix::summary(sebeta)

    # Calculate the number of columns - i.e. number of OTUs/metabolites 
    n_bac <- ncol(matrices[[1]])
    n_euk <- ncol(matrices[[2]])
    n_met <- ncol(matrices[[3]])

    # Create a new column that encodes the types of association between two nodes 
    # Nodes are numbered in order from bac to euk to mets 

    elist <- elist %>%
    mutate(interaction.type = case_when(
        i %in% 1:n_bac & 
        j %in% 1:n_bac ~ "Bacteria-Bacteria",

        i %in% (n_bac+1):(n_bac+n_euk) & 
        j %in% (n_bac+1):(n_bac+n_euk) ~ "Eukaryote-Eukaryote",
        
        i %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) & 
        j %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Metabolite-Metabolite",

        i %in% 1:n_bac & 
        j %in% (n_bac+1):(n_bac+n_euk) ~ "Bacteria-Eukaryote",
        j %in% 1:n_bac & 
        i %in% (n_bac+1):(n_bac+n_euk) ~ "Bacteria-Eukaryote",

        i %in% 1:n_bac & 
        j %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Bacteria-Metabolite",
        j %in% 1:n_bac & 
        i %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Bacteria-Metabolite",

        i %in% (n_bac+1):(n_bac+n_euk) & 
        j %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Eukaryote-Metabolite",
        j %in% (n_bac+1):(n_bac+n_euk) & 
        i %in% (n_bac+n_euk+1):(n_bac+n_euk+n_met) ~ "Eukaryote-Metabolite",

        ))
    return(elist)
}

getEdgeDf<-function(se.test, se.rand, matrices){
    df1<-cbind(getMultiNetworkDataframe(se.test, matrices), group="Observed")
    df2<-cbind(getMultiNetworkDataframe(se.rand, matrices), group="Null")
    df<-rbind(df1, df2)
    df<-as.data.frame(df)
    df<-as_tibble(df) %>%
        mutate(value=as.numeric(x))
    return(df)
}

makeNetworkFigure<-function(edge_df){

    means<-edge_df %>%
        group_by(group, interaction.type) %>%
        summarise(mean=mean(value))


    edge_df %>% 
        left_join(means, by= c("group", "interaction.type"))%>%
        ggplot(data=., aes(x=value, fill=group))+
            geom_histogram(alpha=0.5, position="identity",  binwidth=0.01)+
            scale_fill_manual(values=c("black", "blue"))+
            scale_color_manual(values=c("black", "blue"))+
            facet_wrap(~interaction.type, scales = "free_y")+
            geom_vline(aes(xintercept = mean, color=group))+
            #geom_vline(xintercept=0, color="red", linetype = 'dotted')+
            theme_classic()+
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            xlab ("Edge Weights")+ ylab("Frequency")+
            theme(legend.title=element_blank(), strip.background=element_blank(), 
                strip.text = element_text(size = 12))
}


makeDiffFigure<-function(edge_df, asProportion=TRUE){

    edgesInObserved<-edge_df %>%
        filter(group=="Observed") %>%
        group_by(interaction.type) %>%
        summarise(n=n())

    p<-makeNetworkFigure(edge_df)

    g1<-ggplot_build(p)$data[[1]] %>% 
        as_tibble() %>%
        filter(group==1)

    g2<-ggplot_build(p)$data[[1]] %>% 
        as_tibble() %>%
        filter(group==2) %>%
        select(PANEL, xmin, xmax, count) %>%
        rename(count2=count)

    newPlot_df<-left_join(g1, g2, by=c("PANEL", "xmin", "xmax")) %>%
        mutate(diff=count2-count) %>%
        mutate(interaction.type=(factor(PANEL, labels=c("Bacteria-Bacteria", "Bacteria-Eukaryote", "Bacteria-Metabolite", "Eukaryote-Eukaryote", "Eukaryote-Metabolite", "Metabolite-Metabolite"))))    %>%
        mutate(positive=xmin>0) %>%
        left_join(edgesInObserved, by="interaction.type") %>%
        mutate(diff_proportion=diff/n*100)

    if (asProportion) {
        newPlot_df %>%
            ggplot(data=., aes(xmin=xmin,xmax=xmax,ymax=diff_proportion,ymin=0, fill=positive)) + 
                geom_rect(show.legend = FALSE)+
                facet_wrap(~interaction.type)+
                geom_vline(xintercept=0, linetype=2)+
                geom_hline(yintercept=0, linetype=2)+
                theme_classic()+
                xlab ("Edge Weights")+ ylab("Difference in edge frequency\nbetween observed and null networks\n(% of total edges of this type in observed network)")+
                theme(legend.title=element_blank(), 
                    strip.background=element_blank(), 
                    strip.text = element_text(size = 12))
    } else {
        newPlot_df %>%
            ggplot(data=., aes(xmin=xmin,xmax=xmax,ymax=diff,ymin=0, fill=positive)) + 
                geom_rect(show.legend = FALSE)+
                facet_wrap(~interaction.type, scales="free_y")+
                geom_vline(xintercept=0, linetype=2)+
                geom_hline(yintercept=0, linetype=2)+
                theme_classic()+
                xlab ("Edge Weights")+ ylab("Difference in edge frequency\nbetween observed and null networks")+
                theme(legend.title=element_blank(), 
                    strip.background=element_blank(), 
                    strip.text = element_text(size = 12))

    }

}


#load data
matrix_list<-readRDS( paste0("Outputs/NetworkMatrices",cluster,".RDS"))


# make plots

#all
se_all<-readRDS(paste0("Outputs/se_all", cluster, ".RDS"))
se_all_rand<-readRDS(paste0("Outputs/se_all_rand", cluster, ".RDS"))
getStability(se_all)
getStability(se_all_rand)
edge_all_df<-getEdgeDf(se_all, se_all_rand, matrices=matrix_list[["all"]])
p_all<-makeNetworkFigure(edge_all_df)
ggsave(paste0("Figures/Figure_all",cluster,".pdf"), p_all,  width = 12, height = 8)

p_all_diff<-makeDiffFigure(edge_all_df, asProportion=TRUE)
ggsave(paste0("Figures/Figure_all_diff_proportions", cluster, ".pdf"), p_all_diff,  width = 12, height = 8)

p_all_diff<-makeDiffFigure(edge_all_df, asProportion=FALSE)
ggsave(paste0("Figures/Figure_all_diff_absolute", cluster, ".pdf"), p_all_diff,  width = 12, height = 8)



################################################################
# unused coastal and oceanic split (not enough data to fit)
# multiple cluster options therefore are not implemented

#coastal
se_coastal<-readRDS("Outputs/se_coastal.RDS")
se_coastal_rand<-readRDS("Outputs/se_coastal_rand.RDS")
getStability(se_coastal)
getStability(se_coastal_rand)
edge_coastal_df<-getEdgeDf(se_coastal, se_coastal_rand, matrices=matrix_list[["coastal"]])
p_coastal<-makeNetworkFigure(edge_coastal_df)
ggsave("Figures/Figure_coastal.pdf", p_coastal,  width = 12, height = 8)

p_coastal_diff<-makeDiffFigure(edge_coastal_df, asProportion=TRUE)
ggsave("Figures/Figure_coastal_diff_proportions.pdf", p_coastal_diff,  width = 12, height = 8)

p_coastal_diff<-makeDiffFigure(edge_coastal_df, asProportion=FALSE)
ggsave("Figures/Figure_coastal_diff_absolute.pdf", p_coastal_diff,  width = 12, height = 8)


#oceanic
se_oceanic<-readRDS("Outputs/se_oceanic.RDS")
se_oceanic_rand<-readRDS("Outputs/se_oceanic_rand.RDS")
getStability(se_oceanic)
getStability(se_oceanic_rand)
edge_oceanic_df<-getEdgeDf(se_oceanic, se_oceanic_rand, matrices=matrix_list[["oceanic"]])
p_oceanic<-makeNetworkFigure(edge_oceanic_df)
ggsave("Figures/Figure_oceanic.pdf", p_oceanic,  width = 12, height = 8)

p_oceanic_diff<-makeDiffFigure(edge_oceanic_df, asProportion=TRUE)
ggsave("Figures/Figure_oceanic_diff_proportions.pdf", p_oceanic_diff,  width = 12, height = 8)

p_oceanic_diff<-makeDiffFigure(edge_oceanic_df, asProportion=FALSE)
ggsave("Figures/Figure_oceanic_diff_absolute.pdf", p_oceanic_diff,  width = 12, height = 8)
