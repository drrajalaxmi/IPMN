## relative immune cell abundance from the CIBERSORTx analysis
## input data from input.R
library(rstatix)
cibersort_pt %>% group_by(Disease) %>%
  summarise(count = n())

df_plot <- cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% 
  filter(!Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>% droplevels() %>%
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
cell_subset <- c("Macrophages.M2", "NK.cells.activated", "T.cells.CD4.memory.activated",
                 "T.cells.gamma.delta", "Macrophages.M1", "Monocytes", "Neutrophils")

cell_subset <- c("Plasma.cells")
# box plots with wilcox test for PDAC vs Non-PDAC groups
for( cell_name in cell_subset){
    p_values <- df_plot %>% 
      dplyr::rename("cell_name" = cell_name) %>%
      wilcox_test( cell_name ~ groups ) %>%
      adjust_pvalue( method = "bonferroni" ) %>%
      add_significance("p.adj")
    
    p_values  <- p_values %>% add_xy_position(x = "Disease", group = "groups") %>%
                        mutate(y.position = max(attributes(p_values)[[4]]$data$cell_name) * 1.08 )
    p_values$xmin <- 2
    p_values$xmax <- 5
    p_values$p.adj.signif[p_values$p.adj.signif =="ns"] <- ''
    
    st_line <- data.frame(xmin = c(1,4),
                          xmax = c(3,6),
                          y.position = c(rep(0.97*p_values$y.position, 2)),
                          p="",
                          group1=1,
                          group2 =2)
    
    # color <- c("#48217390","#433E8590", "#38598C90","#2BB07F90","#85D54A90", "#FDE725FF") ##2BB07F50")
    color <- c("#2F56B5", "#2295C7", "#00A86B", "#FCC92A", "#F4691E",  "#AD0000")
    p <- ggboxplot(df_plot , x = "Disease", y = cell_name,
                   fill = "Disease",palette = color,
                   add = "jitter",
                   legend.title ="")
    #  Add p-value
    # Change method
    # p + stat_compare_means(# comparisons = list(c("Isolated LGD", "Isolated HGD") )
    #               #list(c(levels(df_plot$Disease)[c(4,5)]),c(levels(df_plot$Disease)[c(4,6)] ))
    #                  ) 
    
    
    p <- p + #stat_compare_means(label.y = 0.35 ) +
      labs(x="") + # ylim(c(0, 0.14)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) +
      stat_pvalue_manual(  p_values,  label = "{p.adj}{p.adj.signif}" ) +
      stat_pvalue_manual(  st_line, tip.length = 0  ) 
    
    print(p)
}


## 1. Additional custom comparison 
p_values <- df_plot %>% 
  filter( Disease %in% c("Isolated\nLGD", "Isolated\nHGD" )) %>% base::droplevels() %>%
  dplyr::select(Plasma.cells,Disease )

# color <- c("#48217390","#433E8590", "#38598C90","#2BB07F90","#85D54A90", "#FDE725FF") ##2BB07F50")
color <- c("#2F56B5", #"#2295C7", 
           "#00A86B" #, "#FCC92A", "#F4691E",  "#AD0000"
           )
p <- ggboxplot(p_values , x = "Disease", y = "Plasma.cells",
               # color = "Disease",
               fill = "Disease", palette = color,
               add = "jitter",
               legend.title ="")
#  Add p-value
p + stat_compare_means( method = "wilcox",
                        comparisons = list(c("Isolated\nLGD", "Isolated\nHGD") )
              #list(c(levels(df_plot$Disease)[c(4,5)]),c(levels(df_plot$Disease)[c(4,6)] ))
                 )+
  labs(x="") 

## 2. Additional custom comparison 
p_values <- df_plot %>% 
  filter( Disease %in% c("Isolated\nLGD", "Progressor\nLGD w HGD", "Progressor\nLGD w PDAC" )) %>% droplevels() %>%
  dplyr::select(Plasma.cells,Disease )

# color <- c("#48217390","#433E8590", "#38598C90","#2BB07F90","#85D54A90", "#FDE725FF") ##2BB07F50")
color <- c("#2F56B5", "#2295C7", #"#00A86B" #, 
           "#FCC92A" #, "#F4691E",  "#AD0000"
)
p <- ggboxplot(p_values , x = "Disease", y = "Plasma.cells",
               # color = "Disease",
               fill = "Disease", palette = color,
               add = "jitter",
               legend.title ="")
#  Add p-value
p + stat_compare_means( method = "wilcox",
                        # comparisons = list(c("Isolated\nLGD", "Isolated\nHGD") )
                        comparisons = list(c(levels(p_values$Disease)[c(1,2)]),c(levels(p_values$Disease)[c(1,3)] ))
                       ) +  labs(x="") 

