# CIBERSORTx
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(readxl)
library(umap)
library(ggfortify) 
library(ggpubr)
library(rstatix)
library(RecordTest)  ## for fisher.method
# BiocManager::install("metaSeq")
# library(metaSeq)



## input is from input_ms.R
# meatadata information as in condition_ms 
# CIBERSORTx relative abundance of the immune cells in in the cibersort_pt_ms


#faceted box plot 
cibersort_pt_ms %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% 
  filter(!str_detect(Disease, "normal")) %>%
  pivot_longer( !c("Mixture","Grade", "Patient_id","Category","Disease"  ), 
                names_to = "cells",values_to = "pcnt_cells") %>%
  mutate(cells = factor(cells, levels= colnames(cibersort_pt_ms)[2:23]) ) %>%
  arrange(cells) %>%
  group_by(Patient_id , cells) %>%
  ggplot(aes((Disease), pcnt_cells)) +
  geom_boxplot() + 
  facet_wrap(~cells, scales = "free_y", nrow = 6, strip.position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="Patient ID", y="Immune cells in fraction" )



# Boxplot with statistical test
df_plot_ms <- cibersort_pt_ms %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% 
  filter(!str_detect(Disease, "normal")) %>% 
  droplevels() %>%
  mutate(
    Disease = str_wrap(Disease, 10),
    Disease =factor(Disease, levels= c(str_wrap(levels(.$Disease), 10))),
    groups = ifelse(str_detect(Disease, "PDAC"), "PDAC", "non_PDAC"),
    groups = factor(groups, levels = c( "non_PDAC", "PDAC"))
  )
cell_list <- c("B.cells.naive", "B.cells.memory" ,   "Plasma.cells",
               "T.cells.CD8",   "T.cells.CD4.naive", "T.cells.CD4.memory.resting",
               "T.cells.CD4.memory.activated", "T.cells.follicular.helper","T.cells.regulatory..Tregs.",  
               "T.cells.gamma.delta","NK.cells.resting", "NK.cells.activated" ,
               "Monocytes", "Macrophages.M0", "Macrophages.M1", 
               "Macrophages.M2" , "Dendritic.cells.resting" , "Dendritic.cells.activated", 
               "Mast.cells.resting", "Mast.cells.activated", "Eosinophils","Neutrophils" )

## comparison analsysis for a cell
cell_type <- "T.cells.gamma.delta"
df_plot_ms$T.cells.gamma.delta

p_values_ms <- df_plot_ms %>% 
  dplyr::rename("cell_type" = cell_type) %>%
  wilcox_test( cell_type ~ groups ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")
p_values_ms  <- p_values_ms %>% add_xy_position(x = "Disease", group = "groups") %>%
  mutate(y.position = max(attributes(p_values_ms)[[4]]$data$cell_type) * 1.08 )
p_values_ms$xmin <- 1.5
p_values_ms$xmax <- 4
p_values_ms$p.adj.signif[p_values_ms$p.adj.signif =="ns"] <- ''

st_line_ms <- data.frame(xmin = c(1,3),
                      xmax = c(2,5),
                      y.position = c(rep(0.97*p_values_ms$y.position, 2)),
                      p="",
                      group1=1,
                      group2 =2)

color_ms <- c( #"#2F56B5",
           "#2295C7", "#00A86B", "#FCC92A", "#F4691E",  "#AD0000")

ggboxplot(df_plot_ms, # %>% 
                # mutate(Macrophages = as.numeric(Macrophages.M0)+ as.numeric(Macrophages.M1)+ as.numeric(Macrophages.M2)) ,
               x = "Disease", y = cell_type,
               fill = "Disease", palette = color_ms,
               # add = "jitter", 
          legend.title ="") +
  labs(x="MS") +  ylim(c(0, 0.03)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )
        # plot.margin = margin(t=10,b= 10, r =50),
        # legend.margin=margin(r=100)
        ) +
  stat_pvalue_manual( p_values_ms,  label = "{p.adj}{p.adj.signif}"  ) +
  stat_pvalue_manual( st_line_ms, tip.length = 0 )  +
  guides(fill = guide_legend(nrow = 2, byrow=TRUE)) 

## other custom graphs



color_ms <- c( #"#2F56B5",
  #"#2295C7", 
  "#00A86B", #"#FCC92A", 
  "#F4691E" #,  "#AD0000"
  )
df_plot_ms_subset <- df_plot_ms %>%
  filter(Disease %in% c( "isolated\nHGD", "progressor\nHGD w PDAC")) %>%  droplevels()

ggboxplot(df_plot_ms_subset, # %>% 
          # mutate(Macrophages = as.numeric(Macrophages.M0)+ as.numeric(Macrophages.M1)+ as.numeric(Macrophages.M2)) ,
          x = "Disease", y = cell_type,
          fill = "Disease", palette = color_ms,
          add = "jitter", legend.title ="") +
  labs(x="MS") +  #ylim(c(0, 0.125)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )
        # plot.margin = margin(t=10,b= 10, r =50),
        # legend.margin=margin(r=100)
  ) +
  stat_compare_means(label.y = 0.12 )


### combined p values for the MS and BE data set using fisher method
p_test <- RecordTest::fisher.method(c(p_values_ms$p, p_values$p))$p.value
paste0( "Fisher's combined probability test p = ", round(p_test, 4))






###UMAP

ciber_ms <- cibersort_ms%>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>%
  filter(Mixture %in% condition_ms$Sample_ID) %>% arrange(Mixture)
rownames(ciber) <- ciber$Mixture


set.seed(95) #65  134 140 144, 146 148, 151, 153, 411 443
umap_res <- umap(ciber[-1] )
# dim(condition_ms)
umap_df <- data.frame(umap_res$layout) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::left_join(condition_ms, by = "Sample_ID") %>%
  mutate( Disease = factor(Disease, levels = c("normal-LGD-acinar" ,"normal-LGD-duct" ,
                                               "normal-HGD-acinar" ,"normal-HGD-duct" ,
                                               "normal-PDAC-acinar" ,"normal-PDAC-duct" ,
                                               "isolated LGD",
                                               "progressor LGD w HGD" , "progressor LGD w PDAC",
                                               "isolated HGD" ,"progressor HGD w PDAC", "PDAC")),
                   )

plotly::ggplotly(ggplot(umap_df,aes(  x = X1, y = X2,  color= Disease) )+ 
    geom_point(size=3) +
    scale_colour_manual(values = c("normal-LGD-acinar"  = "gray90",
                                   "normal-LGD-duct" = "gray90",
                                   "normal-HGD-acinar" = "gray50",
                                   "normal-HGD-duct" = "gray50",
                                   "normal-PDAC-acinar"  = "gray10",
                                   "normal-PDAC-duct"  = "gray10",
                                   "isolated LGD" ="#DBF3FA", #"gray40",
                                   "progressor LGD w HGD" = "#7AD7F0",
                                   "progressor LGD w PDAC"= "#266CA9", #"#FFA80F",
                                   "isolated HGD"= "#FFA80F",  #"#266CA9",
                                   "progressor HGD w PDAC"="#FE8116",
                                   "PDAC"="#FE5A1D"))+
    labs(title= "Cibersort Immune Cell Proportions: UMAP", x ="UMAP 1", y="UMAP 2") +
    theme_bw()
) 
