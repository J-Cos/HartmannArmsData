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

#load data
matrix_list<-readRDS( "Outputs/NetworkMatrices.RDS")


# make plots

#all
se_all<-readRDS("Outputs/se_all.RDS")
se_all_rand<-readRDS("Outputs/se_all_rand.RDS")
getStability(se_all)
getStability(se_all_rand)
edge_all_df<-getEdgeDf(se_all, se_all_rand, matrices=matrix_list[["all"]])
p_all<-makeNetworkFigure(edge_all_df)
ggsave("Figures/Figure_all.pdf", p_all,  width = 12, height = 8)


#coastal
se_coastal<-readRDS("Outputs/se_coastal.RDS")
se_coastal_rand<-readRDS("Outputs/se_coastal_rand.RDS")
getStability(se_coastal)
getStability(se_coastal_rand)
edge_coastal_df<-getEdgeDf(se_coastal, se_coastal_rand, matrices=matrix_list[["coastal"]])
p_coastal<-makeNetworkFigure(edge_coastal_df)
ggsave("Figures/Figure_coastal.pdf", p_coastal,  width = 12, height = 8)

#oceanic
edge_oceanic_df<-getEdgeDf(se_oceanic, se_oceanic_rand, matrices=matrix_list[["oceanic"]])
p_oceanic<-makeNetworkFigure(edge_oceanic_df)
ggsave("Figures/Figure_oceanic.pdf", p_oceanic,  width = 12, height = 8)
