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
library(tibble)
library(grDevices)
library(RColorBrewer)


BE02_case_list <-read_excel("data/BE02_case_list.xlsx", sheet = "FINAL")
column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE02_",BE02_case_list$Sample_ID, "A")
sample_list <- c("BE02_18A",  "BE02_59A",  "BE02_60A",  "BE02_61A",  "BE02_74A" , "BE02_75A",  "BE02_101A", "BE02_102A", "BE02_104A", "BE02_105A", "BE02_113A", "BE02_114A")
condition <- BE02_case_list %>% 
  filter(!Sample_ID %in% sample_list)

# sample file:  
ciber_sample <- condition %>% mutate(Grade = factor(Grade, levels = c("normal", "low grade", "high grade", "invasive") ),
                                     Disease = factor(Disease, levels = c("Normal-LGD" ,"Normal-HGD", "Normal-PDAC","Isolated LGD",
                                                                          "Progressor LGD w HGD" , 
                                                                          "Isolated HGD" ,"Progressor LGD w PDAC",
                                                                          "Progressor HGD w PDAC", "PDAC")),
                                     Patient_id = factor(Path_ID, levels = unique(Path_ID))) %>%
                              arrange(Patient_id , Grade, Disease)

str(ciber_sample)
# crrelation file
cibersort <- read.csv("data/CIBERSORTx_Job4_Results.csv")
# str(cibersort)
cibersort <- cibersort %>% mutate(Mixture = paste0(str_replace(Mixture, "BE", "BE02_"), "A"))

cibersort_pt <- cibersort %>%  #mutate(Mixture = paste0(str_replace(Mixture, "BE", "BE02_"), "A")) %>%
          left_join(ciber_sample %>% dplyr::select(Sample_ID, Grade, Patient_id, Path_diagnosis,Main_duct, Category, Disease), by = c("Mixture" = "Sample_ID")) %>%
                filter(!is.na(Grade))   



cell_list <- c("B.cells.naive","B.cells.memory" ,"Plasma.cells",
           "T.cells.CD8","T.cells.CD4.naive", "T.cells.CD4.memory.resting",
           "T.cells.CD4.memory.activated", "T.cells.follicular.helper","T.cells.regulatory..Tregs.",  
           "T.cells.gamma.delta","NK.cells.resting",
           "NK.cells.activated" ,"Monocytes", "Macrophages.M0", 
           "Macrophages.M1", "Macrophages.M2" , "Dendritic.cells.resting" ,    
           "Dendritic.cells.activated", "Mast.cells.resting", "Mast.cells.activated",
           "Eosinophils","Neutrophils" )
suppressive <- c( #"T.cells.regulatory..Tregs.",  
  "NK.cells.activated" ,"Monocytes", "T.cells.CD4.memory.activated",
                 "T.cells.gamma.delta",
                 "Macrophages.M1" , "Dendritic.cells.activated" ,
  "Neutrophils", "Monocytes")
apcs <-  c("B.cells.memory" , "Macrophages.M0", 
           "Macrophages.M1",   
           "Dendritic.cells.activated" )
inflammatory <- c( "Macrophages.M1","Mast.cells.resting", "Mast.cells.activated",
                  "Eosinophils","Neutrophils" )
toxic <- c("B.cells.memory", "Plasma.cells",
           "T.cells.CD8",
           "T.cells.CD4.memory.activated", "T.cells.follicular.helper",
           "T.cells.gamma.delta","NK.cells.resting",
           "NK.cells.activated" , "Dendritic.cells.resting" ,    
           "Dendritic.cells.activated" )

plot(as.factor(cibersort_pt$Patient_id), cibersort_pt$T.cells.CD8, col = cibersort_pt$Patient_id)

### from Neoantigen
# stage_lists <- neoantigen_data %>%
#   filter(Category == "LGD") %>%
#   group_by(Path_ID, sample_ID) %>%
#   summarise(Neoantigens = list(unique(mt_epitope_seq)))
# 
# pt_list <- stage_lists %>% group_by(Path_ID,sample_ID,Neoantigens) %>%
#   summarize(count = length(unlist(Neoantigens))) %>% arrange(count) 


#####


plotly::ggplotly(cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% 
               # select(colnames(.)[1:(ncol(.)-2)]) %>% 
                pivot_longer( !c("Mixture","Grade", "Patient_id","Path_diagnosis","Main_duct", "Category", "Disease" ), names_to = "cells",values_to = "pcnt_cells") %>%
                # filter(cells %in% c( "Plasma.cells")) %>%   # apcs toxic , suppressive , inflammatory
               filter(Grade != "normal") %>%
               # filter(Patient_id %in% pt_list$Path_ID) %>%
                #filter(Mixture %in% pt_list$sample_ID) %>%
                arrange(Category,as.numeric(pcnt_cells)) %>%
                mutate(#Mixture = factor(Mixture, levels = pt_list$sample_ID),
                       Patient_id = factor(Patient_id , levels = unique(Patient_id)),
                       cells = factor(cells, levels= colnames(cibersort_pt)[2:23]),
                       color = ifelse(as.numeric(Grade)==1, "black", ifelse(as.numeric(Grade)==2, 'blue', ifelse(as.numeric(Grade)==3, "green", "orange")))  ,
                       Category = factor(Category,levels= c("LGD", "HGD", "PDAC"))) %>%
                # arrange(pcnt_cells) %>%
                # group_by(Patient_id ) %>%
               # ggplot(aes(Grade, pcnt_cells)) +
   # ggplot(aes(Patient_id, pcnt_cells)) +
  ggplot(aes(Category, pcnt_cells)) +
                 
                geom_boxplot() +
               geom_jitter(size= 0.5) +
  # geom_col(position= position_dodge(), width = 0.5) + 
  facet_wrap(~cells, scales = "free_y") +
  #  theme_minimal() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +

  # scale_colour_manual(values = c("grey", "blue", "green", "red") ) +
#  coord_flip() #+
  labs(x="", y="Immune cells in fraction" ) )
  #     title=paste0("plot_name") )


#line plot
cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  pivot_longer( !c("Mixture","Grade", "Patient_id",
                   "Path_diagnosis","Main_duct","Category","Disease" 
  ), names_to = "cells",values_to = "pcnt_cells") %>%
  mutate(cells = factor(cells, levels= colnames(cibersort_pt)[2:23]),
         color = ifelse(as.numeric(Grade)==1, "black", ifelse(as.numeric(Grade)==2, 'blue', ifelse(as.numeric(Grade)==3, "green", "orange")))   ) %>%
  arrange(cells) %>%
  # filter( cells %in% unique(.$cells)[1:7]) %>%
  group_by(Patient_id , cells) %>%
  ggplot(aes((Grade), pcnt_cells)) +
  # ggplot(aes(Patient_id, pcnt_cells)) +
    # geom_point(aes(colour = Grade)) + 
   geom_boxplot() + 
  # geom_line(aes(group=cells, colour = cells)) +
  facet_wrap(~cells, scales = "free_y", nrow = 6, strip.position = "top") +
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  
  # scale_colour_manual(values = c("grey", "blue", "green", "red") ) +
  #  coord_flip() #+
  labs(x="Patient ID", y="Immune cells in fraction" )
#     title=paste0("plot_name") )
#####



##### boxplot with p values
colnames(cibersort_pt)
#line plot
cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  filter(Disease %in% c( "Isolated LGD", "Isolated HGD")) %>%
  # pivot_longer( !c("Mixture","Grade", "Patient_id",
  #                  "Path_diagnosis","Main_duct","Category","Disease" 
  # ), names_to = "cells",values_to = "pcnt_cells") %>%
  # mutate(cells = factor(cells, levels= colnames(cibersort_pt)[2:23]),
  #        color = ifelse(as.numeric(Grade)==1, "black", ifelse(as.numeric(Grade)==2, 'blue', ifelse(as.numeric(Grade)==3, "green", "orange")))   ) %>%
  # arrange(cells) %>%
  # filter( cells %in% unique(.$cells)[1:7]) %>%
  #group_by(Patient_id , cells) %>%
  ggplot(aes(Disease, Macrophages.M1)) +
  # ggplot(aes(Patient_id, pcnt_cells)) +
  # geom_point(aes(colour = Grade)) + 
  geom_boxplot() + 
  # geom_line(aes(group=cells, colour = cells)) +
 # facet_wrap(~cells, scales = "free_y", nrow = 6, strip.position = "top") +
  # theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  
  # scale_colour_manual(values = c("grey", "blue", "green", "red") ) +
  #  coord_flip() #+
  labs(x="", y="Immune cells in fraction" )
#     title=paste0("plot_name") )

colnames(cibersort_pt)
levels(cibersort_pt$Disease)

df_plot <- cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
   # filter(Disease %in% c( "Isolated LGD", "Isolated HGD"))
 filter(!Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
  #filter(T.cells.CD8 < 0.3) %>%
# filter(Disease %in% c( "Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC" )) %>%
  droplevels() %>%
   mutate(
         Disease = str_wrap(Disease, 10),
          Disease =factor(Disease, levels= c(str_wrap(levels(.$Disease), 10))),
         groups = ifelse(str_detect(Disease, "PDAC"), "PDAC", "non_PDAC"),
         groups = factor(groups, levels = c( "non_PDAC", "PDAC"))
     # Disease = ifelse(str_detect(Disease, "Progressor"),"Progressor LGD", "Isolated LGD") ,
     #      Disease = factor(Disease , levels = c("Isolated LGD", "Progressor LGD"))
     )


p_values <- df_plot %>% 
  wilcox_test( B.cells.naive ~ groups ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_values  <- p_values %>% add_xy_position(x = "Disease", group = "groups")
p_values$xmin <- 2
p_values$xmax <- 5
#p_values$y.position <- p_values$y.position +0.05
p_values$p.adj.signif[p_values$p.adj.signif =="ns"] <- ''

st_line <- data.frame(xmin = c(1,4),
                      xmax = c(3,6),
                      y.position = c(p_values$y.position-0.01, p_values$y.position-0.01 ),
                      p="",
                      group1=1,
                      group2 =2)
  
color <- c("#48217390","#433E8590", "#38598C90","#2BB07F90","#85D54A90", "#FDE725FF") ##2BB07F50")
p <- ggboxplot(df_plot , x = "Disease", y = "B.cells.naive",
             # color = "Disease",
               fill = "Disease",palette = color,
               add = "jitter",
              legend.title ="")
#  Add p-value

# Change method
# p + stat_compare_means(# comparisons = list(c("Isolated LGD", "Isolated HGD") )
# 
#               #list(c(levels(df_plot$Disease)[c(4,5)]),c(levels(df_plot$Disease)[c(4,6)] ))
#                  ) 


 p + #stat_compare_means(label.y = 0.35 ) +
  labs(x="") + # ylim(c(0, 0.14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) +
   stat_pvalue_manual(
     p_values,  label = "{p.adj}{p.adj.signif}" #, tip.length = 0.02
     #step.increase = 0.05
   ) +
   stat_pvalue_manual(
     st_line, tip.length = 0
     # step.increase = 0.05
   ) 
 # +
 #   scale_y_continuous(
 #     breaks = c(0.025, 0.050, 0.075, 0.1, 0.125 ))
 


####stack bar graph
library(viridis)


str(cibersort_pt)
pt_list <-cibersort_pt %>% arrange(#Disease, 
                                   #T.cells.regulatory..Tregs.
                                    T.cells.gamma.delta
                                  # NK.cells.resting
                                  #NK.cells.activated
                                  #Neutrophils,
                                   # Mast.cells.activated,
                                 # T.cells.CD8
                                  # B.cells.naive,
                                  # desc(Plasma.cells)
                                 #  Macrophages.M1
  
                                   
                                   #,desc(B.cells.naive)
                                   )   %>% dplyr::select(Mixture,Disease,T.cells.gamma.delta )



getPalette <- colorRampPalette(c( #viridis_pal(option = "G")(4)[-1],
                               rev(brewer.pal(3, "Blues")),
                                 brewer.pal(7, "YlGnBu"), # PRGn  YlGnBu
                                #viridis_pal(option = "H")(7),
                                
                               # viridis_pal(option = "D")(8),
                            
                              brewer.pal(5, "BrBG")[2:3],
                               brewer.pal(7, "OrRd")[2:6],
                            # brewer.pal(5, "Blues"),
                              brewer.pal(4, "Purples")))
  #c(brewer.pal(5, "Viridis"),brewer.pal(7, "Reds")[1:7],brewer.pal(5, "GnBu"), brewer.pal(4, "Greys") ) )
cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  filter(Grade !="normal") %>%
  pivot_longer( !c("Mixture","Grade", "Patient_id",
                   "Path_diagnosis","Main_duct","Category","Disease" 
               ), names_to = "cells", values_to = "pcnt_cells" ) %>%
  group_by( Disease, pcnt_cells, cells ) %>% 
  # arrange(Disease, desc(as.numeric(pcnt_cells)) ) %>%
  mutate( 
         Mixture = factor(Mixture, levels = unique(pt_list$Mixture) ),
         x = paste0(Disease, "  ",Mixture),
         x = factor(x , levels = unique(paste0( pt_list$Disease, "  ", pt_list$Mixture))),
         cells = factor(cells, levels= colnames(cibersort_pt)[2:23]),
         color = ifelse(as.numeric(Grade)==1, "black", ifelse(as.numeric(Grade)==2, 'blue', ifelse(as.numeric(Grade)==3, "green", "orange")))  
         ) %>%
  
  # filter( cells %in% unique(.$cells)[1:7]) %>%
  #group_by(Patient_id , cells) %>%
  ggplot(aes(x, pcnt_cells, fill = cells)) +
  # ggplot(aes(Patient_id, pcnt_cells)) +
  # geom_point(aes(colour = Grade)) + 
  geom_col(position = position_stack()) + 
  # geom_line(aes(group=cells, colour = cells)) +
  # facet_wrap(~Disease, scales = "free_x", nrow = 1, strip.position = "top") +
   theme_minimal() +
  # labs(ylim= 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1 )) +
  # scale_fill_viridis(discrete = T) +
  scale_colour_manual(
    values = getPalette(22),
    aesthetics = c("colour", "fill")
  ) +
  # scale_fill_brewer(palette="Spectral") + 
  
  # scale_colour_manual(values = c("grey", "blue", "green", "red") ) +
  #  coord_flip() #+
  labs(x="", y="Immune cells in fraction")
  
#     title=paste0("plot_name") )


### create a medican of this tulip garden
cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
   filter(Grade !="normal") %>%
  pivot_longer( !c("Mixture","Grade", "Patient_id",
                   "Path_diagnosis","Main_duct","Category","Disease" 
                   ), names_to = "cells", values_to = "pcnt_cells" ) %>%
  group_by( Category, Patient_id, cells ) %>% 
  #summarize(pcnt_cells = median(pcnt_cells)) %>% 
  arrange( Category,  desc(as.numeric(pcnt_cells)) ) %>%
  mutate( 
  #  Mixture = factor(Mixture, levels = unique(pt_list$Mixture) ),
    x = paste0(Category, Patient_id, Mixture),
    # x = factor(x ,  levels = c( "Isolated LGD",
    #                             "Isolated HGD" ,
    #                             "Progressor LGD w HGD" , 
    #                             "Progressor LGD w PDAC",
    #                             "Progressor HGD w PDAC", "PDAC") ),
    cells = factor(cells, levels= colnames(cibersort_pt)[2:23])
         )  %>%
  ggplot(aes(x, pcnt_cells, fill = cells )) +
  geom_col(position = position_stack(), color = "black") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1 )) +
  scale_colour_manual(
    values = getPalette(22),
    aesthetics = c("colour", "fill")
         ) +
  labs(x="", y="Immune cells in fraction")


### tulip garden cell enrichment

df_fc <- cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  # filter(Category =="PDAC") %>%
  pivot_longer( !c("Mixture","Grade", "Patient_id",
                   "Path_diagnosis","Main_duct","Category","Disease" 
  ), names_to = "cells", values_to = "pcnt_cells" ) %>%
  mutate( pcnt_cells = pcnt_cells + 0.0001,
           Disease = as.character(Disease) ,
         Disease= ifelse(str_detect(Disease, "Normal"), "Normal", Disease )) %>%
  dplyr::select( -c(Grade, Mixture, Path_diagnosis, Main_duct)) %>% 
  pivot_wider( names_from = Disease, values_from = pcnt_cells ) %>%
  filter(!is.na(Normal))  %>% 
  # mutate(Normal = ifelse(Normal == 0, 0.0001, Normal)) %>%
  group_by( Category, Patient_id, cells ) %>% 
  summarize( #Normal = Normal,
    `Isolated LGD` = (`Isolated LGD` - Normal)/Normal,
    `Isolated HGD` = (`Isolated HGD` - Normal) /Normal,
    `Progressor LGD w HGD` = (`Progressor LGD w HGD`- Normal)/ Normal,
    `Progressor LGD w PDAC` = (`Progressor LGD w PDAC`- Normal)/Normal,
    `Progressor HGD w PDAC` = ( `Progressor HGD w PDAC` - Normal)/Normal,
     PDAC = (PDAC -Normal)/Normal ) %>%
  pivot_longer( !c( Category, Patient_id, cells) , names_to = "Disease", values_to = "fold_change"  )

  


df_fc %>% filter(!is.na(fold_change) ) %>%
  mutate(
  Disease = factor(Disease, 
                   levels = c( "Isolated LGD",
                               "Isolated HGD" ,
                               "Progressor LGD w HGD" ,
                               "Progressor LGD w PDAC",
                                "Progressor HGD w PDAC", "PDAC")
    ),
  cells = factor(cells, levels= colnames(cibersort_pt)[2:23]),
  fold_change = ifelse(fold_change <0, -1*log(abs(fold_change), 2) ,
                ifelse(fold_change ==0 , 0,
                       log(abs(fold_change+0.01), 2) ))
                 ) %>% 
  # group_by(Category, cells, Disease) %>%
  # summarise(fold_change = mean(fold_change, na.rm =TRUE)) %>%
  ggplot(aes(fold_change, Disease )) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ cells )+
  labs(x="Log2 Fold Change from Normal")



 
 
#####
library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()
library(scales)
show_col(viridis_pal()(12))

brewer.pal(22, "Spectral")



cibersort_pt %>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  filter(Grade !="normal") %>%
  pivot_longer( !c("Mixture","Grade", "Patient_id",
                   "Path_diagnosis","Main_duct","Category","Disease" 
  ), names_to = "cells",values_to = "pcnt_cells") %>%
  group_by(Disease,pcnt_cells, cells) %>% 
  # arrange(Disease,desc(as.numeric(pcnt_cells)) ) %>%
  arrange(Category, Patient_id,desc(as.numeric(pcnt_cells)) ) %>%
  mutate(Mixture = factor(Mixture, levels = unique(pt_list$Mixture)),
         # x = paste0(Disease, "  ",Mixture),
         x = paste0(Category, "  ",Patient_id, " ", Mixture),
        # x = factor(x , levels = unique(paste0( pt_list$Disease, "  ", pt_list$Mixture))),
         cells = factor(cells, levels= colnames(cibersort_pt)[2:23]),
         color = ifelse(as.numeric(Grade)==1, "black", ifelse(as.numeric(Grade)==2, 'blue', ifelse(as.numeric(Grade)==3, "green", "orange")))   ) %>%
  
  # filter( cells %in% unique(.$cells)[1:7]) %>%
  #group_by(Patient_id , cells) %>%
  ggplot(aes(x, pcnt_cells, fill = cells)) +
  # ggplot(aes(Patient_id, pcnt_cells)) +
  # geom_point(aes(colour = Grade)) + 
  geom_col() + 
  coord_flip() +
  # geom_line(aes(group=cells, colour = cells)) +
  facet_wrap(~ cells, scales = "free_x", nrow = 1, strip.position = "top",
             labeller = label_wrap_gen(width = 5, multi_line = TRUE)) +
  theme_minimal() +
  # labs(ylim= 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1 )) +
  # scale_fill_viridis(discrete = T) +
  # scale_colour_manual(
  #   values = getPalette(22),
  #   aesthetics = c("colour", "fill")
  # ) +
  # scale_fill_brewer(palette="Spectral") + 
  
  # scale_colour_manual(values = c("grey", "blue", "green", "red") ) +
  #  coord_flip() #+
  labs(x="", y="Immune cells in fraction") +
  guides(fill = guide_legend(ncol = 1, byrow=TRUE)) 

#     title=paste0("plot_name") )


###UMAP

cibersort_pt <- cibersort %>% left_join(ciber_sample %>% dplyr::select(Sample_ID, Grade, Patient_id), by = c("Mixture" = "Sample_ID")) %>%
  filter(!is.na(Grade)) 
ciber <- cibersort%>% dplyr::select(-c("P.value", "Correlation", "RMSE")) %>%
            filter(Mixture %in% condition$Sample_ID) %>% arrange(Mixture)
rownames(ciber) <- ciber$Mixture
ciber[-1]
#for (n in seq(461,500)){
  set.seed(65) #65  134 140 144, 146 148, 151, 153, 411 443
umap_res <- umap(ciber[-1] )
# dim(condition)
umap_df <- data.frame(umap_res$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("Sample_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::left_join(condition, by = "Sample_ID") %>%
  mutate( Disease = factor(Disease, levels = c("Normal-LGD" ,"Normal-HGD", "Normal-PDAC","Isolated LGD",
                                               "Progressor LGD w HGD" , "Progressor LGD w PDAC",
                                               "Isolated HGD" ,"Progressor HGD w PDAC", "PDAC")),
  )

plotly::ggplotly(ggplot(umap_df,aes(
  x = X1,
  y = X2,
  color= Disease) )+ 
    geom_point(size=3) +
   scale_colour_manual(values = c("Normal-LGD"  = "gray70",
                                 "Normal-HGD" = "gray50", #"#DBF3FA",
                                "Normal-PDAC" = "gray30", #"#FFF600",
                               "Isolated LGD" ="#DBF3FA", #"gray40",
                              "Progressor LGD w HGD" = "#7AD7F0",
                             "Progressor LGD w PDAC"= "#266CA9", #"#FFA80F",

                           "Isolated HGD"= "#FFA80F",  #"#266CA9",

                         "Progressor HGD w PDAC"="#FE8116",
                        "PDAC"="#FE5A1D"))+
    
    labs(title= "Cibersort Immune Cell Proportions: UMAP", x ="UMAP 1", y="UMAP 2") +
    theme_bw()
  #scale_colour_brewer(palette = "Greens")
) 

#}
# UMAP dimensionality reduction
# umap_result <- umap(ciber[-1])
umap_embedding <- umap_res$layout

# K-means clustering on UMAP embeddings
k <- 4  # Number of clusters
clusters <- kmeans(umap_embedding, centers = k)$cluster

# Plot
ggplot(data.frame(UMAP1 = umap_embedding[,1], UMAP2 = umap_embedding[,2], Cluster = as.factor(clusters)), 
       aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "UMAP Clustering of Immune Profiles", x = "UMAP1", y = "UMAP2") +
  theme_minimal()




cnt_pca <- prcomp((ciber[-1]))
plotly::ggplotly(autoplot(cnt_pca, data = condition %>% arrange(Sample_ID),color = "Path_diagnosis"))
plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (all Grades)") 


library(Rtsne)

matx <- as.matrix(ciber[-1])
dim(matx)
tsne_out <- Rtsne(matx, perplexity = 10)


tsne_plot <- data.frame(ciber[1],
                        x = tsne_out$Y[,1], 
                        y = tsne_out$Y[,2]) %>% 
  dplyr::left_join(condition, by = c("Mixture" ="Sample_ID"))

str(tsne_out)



str(tsne_plot)

# Plotting the plot using ggplot() function

plotly::ggplotly(ggplot(tsne_plot,aes(x=x, y=y, color = Disease)) +
           geom_point(aes(size = 0.04))
)


heatmap(cor(ciber[-1]),  margins = c(12, 12), main ="Correlation plot (pearson)") 
cor_value <- cor(ciber[-1])
cor_value[cor_value < -0.4]
pairs(ciber[-1])



cor(ciber$NK.cells.resting, ciber$NK.cells.activated)
str(ciber)

 
ciber %>% left_join(condition, by = c("Mixture"= "Sample_ID")) %>%
   filter(Path_diagnosis != "biliary and gastric") %>%
ggplot(aes(x=Grade,y= T.cells.CD4.naive/T.cells.CD4.memory.activated  ) ) + 
          #geom_bar(stat = "identity") +
          geom_boxplot() +
          geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="", y="Immune cells in fraction" , title = "Neutrophils") +
  facet_grid(Path_diagnosis ~ Category ) 






##### group comparission and p-value
levels(cibersort_pt$Disease)
p_table <-  data.frame(cell = colnames(cibersort_pt)[2:23] )







for(cell in p_table$cell) {
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell])
  n2 <-length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" , "Progressor LGD w PDAC"), cell])
  col <- paste0("Isolated_LGD_Vs_Progressor_LGD_all(n = ",n1," & ",n2,")")
  res <-t.test(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" , "Progressor LGD w PDAC"), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
}
      
for(cell in p_table$cell) {
  n1 <-length(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell]) 
  n2 <- length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" ), cell])
  col <- paste0("Isolated_LGD_Vs_Progressor_LGD_wHDG(n = ",n1," & ",n2,")")
  res <-t.test(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
}   

for(cell in p_table$cell) {
  n1<- length(cibersort_pt[cibersort_pt$Disease== "Isolated HGD", cell] )
  n2 <- length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  col <- paste0("Isolated_HGD_Vs_Progressor_HGD_wPDAC(n = ",n1," & ",n2,")")
  res <-t.test(cibersort_pt[cibersort_pt$Disease== "Isolated HGD", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 



for(cell in p_table$cell) {
  
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] )
  n2<-length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  col<- paste0("Progressor_LGD_wPDAC_Vs_Progressor_HGD_wPDAC(n = ",n1," & ",n2,")")
  res <-t.test(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 

for(cell in p_table$cell) {
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] )
  n2<- length(cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell] )
  col <- paste0("Progressor LGD w PDAC_Vs_PDAC(n = ",n1," & ",n2,")")
  res <-t.test(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 

for(cell in p_table$cell) {
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Progressor HGD w PDAC", cell] )
  n2<- length(cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell])
  col<-paste0("Progressor_HGD_wPDAC_Vs_PDAC(n = ",n1," & ",n2,")")
  res <-t.test(cibersort_pt[cibersort_pt$Disease== "Progressor HGD w PDAC", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 


write.xlsx(p_table, "data/ImmuneCell_p_value.xlsx")

p_table <-  data.frame(cell = colnames(cibersort_pt)[2:23] )

for(cell in p_table$cell) {
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell])
  n2 <-length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" , "Progressor LGD w PDAC"), cell])
  col <- paste0("Isolated_LGD_Vs_Progressor_LGD_all(n = ",n1," & ",n2,")")
  res <-wilcox.test(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" , "Progressor LGD w PDAC"), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
}

for(cell in p_table$cell) {
  n1 <-length(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell]) 
  n2 <- length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" ), cell])
  col <- paste0("Isolated_LGD_Vs_Progressor_LGD_wHDG(n = ",n1," & ",n2,")")
  res <-wilcox.test(cibersort_pt[cibersort_pt$Disease== "Isolated LGD", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor LGD w HGD" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
}   

for(cell in p_table$cell) {
  n1<- length(cibersort_pt[cibersort_pt$Disease== "Isolated HGD", cell] )
  n2 <- length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  col <- paste0("Isolated_HGD_Vs_Progressor_HGD_wPDAC(n = ",n1," & ",n2,")")
  res <-wilcox.test(cibersort_pt[cibersort_pt$Disease== "Isolated HGD", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 



for(cell in p_table$cell) {
  
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] )
  n2<-length(cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  col<- paste0("Progressor_LGD_wPDAC_Vs_Progressor_HGD_wPDAC(n = ",n1," & ",n2,")")
  res <-wilcox.test(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("Progressor HGD w PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 

for(cell in p_table$cell) {
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] )
  n2<- length(cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell] )
  col <- paste0("Progressor LGD w PDAC_Vs_PDAC(n = ",n1," & ",n2,")")
  res <-wilcox.test(cibersort_pt[cibersort_pt$Disease== "Progressor LGD w PDAC", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 

for(cell in p_table$cell) {
  n1<-length(cibersort_pt[cibersort_pt$Disease== "Progressor HGD w PDAC", cell] )
  n2<- length(cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell])
  col<-paste0("Progressor_HGD_wPDAC_Vs_PDAC(n = ",n1," & ",n2,")")
  res <-wilcox.test(cibersort_pt[cibersort_pt$Disease== "Progressor HGD w PDAC", cell] ,cibersort_pt[cibersort_pt$Disease %in% c("PDAC" ), cell])
  p_table[p_table$cell==cell, col] <- res$p.value
} 


write.xlsx(p_table, "data/ImmuneCell_p_value_wilcox.test.xlsx")



#### sakey diagram

library(networkD3)
(levels(cibersort_pt$Disease))
print(samkey_df %>% arrange(Disease,desc(values)), n= 22)


# cibersort_pt[cibersort_pt==0 ] <-0.0000001
colnames(cibersort_pt[,c(1:23, 32)])
samkey_df <- cibersort_pt[,c(1:23, 32)] %>% 
  #  mutate(
  #    
  #    Mast_cell = Mast.cells.activated + Mast.cells.resting,
  #    NK_cell = NK.cells.activated + NK.cells.resting,
  #    Den_cell = Dendritic.cells.activated + Dendritic.cells.resting,
  #    M = Macrophages.M0+ Macrophages.M1 +Macrophages.M2,
  #    B_cell = B.cells.naive + B.cells.memory +Plasma.cells,
  #    
  #    
  #    Mast_NK = NK_cell/Mast_cell,
  #    M_neu = Monocytes/Neutrophils,
  #    M_N = Monocytes/M,
  #    M_E = Mast_cell/Eosinophils,
  #    M_K =  NK_cell/M,
  #    
  #    
  # #        M1_M2= Macrophages.M1/Macrophages.M2,
  # #      #  NK_resting_activated =NK.cells.activated/NK.cells.resting,
  # #        
  # #        Treg_memResting = T.cells.regulatory..Tregs./ T.cells.CD4.memory.resting,
  # #      Treg_memActivated = T.cells.regulatory..Tregs./ T.cells.CD4.memory.activated,
  # #      Treg_cd4Naive = T.cells.regulatory..Tregs./T.cells.CD4.naive, 
  # #      Treg_cd8 = T.cells.regulatory..Tregs./T.cells.CD8,
  #       T_cell = T.cells.CD8 +
  #                   T.cells.CD4.naive +
  #                 T.cells.CD4.memory.resting+
  #                   T.cells.CD4.memory.activated +T.cells.follicular.helper + T.cells.regulatory..Tregs.+ T.cells.gamma.delta,
  # #      T.cells.regulatory..Tregs. = T.cells.regulatory..Tregs./T.cells.CD4.naive,
  # #      T.cells.CD8 = T.cells.CD8/ T_cell,
  # #     # T.cells.CD4.naive = T.cells.CD4.naive/ T_cell,
  # #      T.cells.CD4.memory.resting = T.cells.CD4.memory.resting/T.cells.CD4.naive,
  # #      T.cells.CD4.memory.activated = T.cells.CD4.memory.activated/T.cells.CD4.naive,
  # #      T.cells.follicular.helper = T.cells.follicular.helper/ T.cells.CD4.naive,
  # #      T.cells.gamma.delta = T.cells.gamma.delta/ T_cell,
  # #      Tfh_Tgd = T.cells.CD4.naive /T.cells.CD4.memory.activated,
  # #       M0_M2 = Macrophages.M0/Macrophages.M2,
  #         ) %>%
         
  
  pivot_longer(-c("Mixture","Disease"), names_to = 'cell', values_to = 'values' ) %>%
  
  
   # dplyr::select(-Mixture) %>% 
   # filter(cell %in% c("Mast_cell", "B_cell", 
   #                    "Neutrophils", "Eosinophils", "Monocytes",
   #                    "NK_cell", "T_cell", "M","Den_cell"
                     # "Mast_NK", "M_K" #, "M_M" #,"M","M_N"
   #  "M1_M2",
  #                  # "T.cells.regulatory..Tregs.", #"T.cells.CD8",
  #                   # "T.cells.CD4.naive", 
  #                  "T.cells.CD4.memory.resting", 
  #                    "T.cells.CD4.memory.activated","T.cells.follicular.helper"#,"T.cells.gamma.delta"
                    # )) %>%
  #filter(cell %in% c(  "M1_M2", "Treg_memResting" ,"Treg_memActivated", "Treg_cd4Naive", "Treg_cd8")) %>% droplevels() %>%
 # filter(Disease %in% c(  "Progressor LGD w PDAC", "Progressor HGD w PDAC" , "PDAC" )) %>% droplevels() %>%
 # filter(Disease %in% c(  "Progressor LGD w HGD" , "Isolated HGD" )) %>% droplevels() %>%
  filter(Disease %in% c( "Isolated LGD", "Progressor LGD w HGD" , "Progressor LGD w PDAC")) %>% droplevels() %>%
 # filter(Disease %in% c( "Isolated HGD" ,"Progressor HGD w PDAC", "PDAC" )) %>% droplevels() %>%
 # filter(Disease %in% c( "Normal-LGD" ,"Normal-HGD", "Normal-PDAC" )) %>% droplevels() %>%
 # filter(Disease %in% c( "Normal-LGD" , "Isolated LGD", "Progressor LGD w HGD" , "Progressor LGD w PDAC", "PDAC")) %>% droplevels() %>%
  #filter(Disease %in% c( "Normal-HGD" , "Isolated HGD", "Progressor HGD w PDAC", "PDAC")) %>% droplevels() %>% 
  #filter(Disease %in% c( "Normal-HGD" , "Progressor LGD w HGD","Isolated HGD")) %>% droplevels() %>% 
 # filter(Disease %in% c( "Normal-LGD" , "Normal-HGD","Normal-PDAC" , "Progressor LGD w PDAC", "Progressor HGD w PDAC", "PDAC")) %>% droplevels() %>% 
  # group_by(cell, Disease) %>%
  #    summarize(
  #        values= mean(values, na.rm=T) ) %>% ungroup() %>%
  # arrange(desc(values)) %>%
   arrange(Disease, desc(values), cell) %>%
  separate(cell , sep="\\.", c("cell_type"), remove = FALSE, extra = "drop", fill="left") %>%
   mutate(order = 1:nrow(.) ,
         # cell_type = str_split(cell_type, ".", c("cell_type", NA)),
         # order = order-1,
         #values = log(values*100),
         # values = ifelse(values==0, 0.0000001, values),
          cell_order = paste0(Disease,cell)
         ) %>% 
    arrange(Disease, cell)  
  # dplyr::select(-values) %>%
  #   pivot_wider(names_from  = Disease, values_from = "order")
# min(samkey_df$values[samkey_df$values!=0])
 samkey_df %>% arrange(Disease, desc(values)) %>%  print(n=23)
df_prop <- samkey_df %>% filter(cell == "T.cells.CD8")
# 
# print(samkey_df, n= 30)
# 
# samkey_df3 <- data.frame(source = c(rep(0, length(samkey_df$order[samkey_df$Disease == levels(samkey_df$Disease)[1]] )),samkey_df$order[samkey_df$Disease != last(levels(samkey_df$Disease))] ),
# target = c(samkey_df$order),
# value = c(samkey_df$values*100),
# name = c(samkey_df$cell),
# status = c(rep(1, 22), rep(2,22)) )
# 
# 
# samkey_df3 <- samkey_df3 %>% arrange (  name) %>% 
#                mutate(cell_state = paste0(status, name))
# 
# 
# 
# nodes <- data.frame("name"= c("Immune cells",(samkey_df3$cell_state)))
#  
# 
# sankeyNetwork(Links = samkey_df3, Nodes = nodes,
#               Source = "source", Target = "target",
#               Value = "value", NodeID = "name", 
#               fontSize= 10, nodeWidth = 10, iterations = 32)
# 
# URL <- paste0('https://cdn.rawgit.com/christophergandrud/networkD3/',
#               'master/JSONdata/energy.json')
# energy <- jsonlite::fromJSON(URL)
# 
# 
# 
# 
# # Plot
# sankeyNetwork(Links = energy$links, Nodes = energy$nodes, Source = 'source',
#               Target = 'target', Value = 'value', NodeID = 'name',
#               units = 'TWh', fontSize = 12, nodeWidth = 30)
# 
# 
# 
# library(ggalluvial)
# head(as.data.frame(UCBAdmissions), n = 12)
# is_alluvia_form(as.data.frame(UCBAdmissions), axes = 1:3, silent = TRUE)
# ggplot(as.data.frame(UCBAdmissions),
#        aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
#   geom_alluvium(aes(fill = Admit), width = 1/12) +
#   geom_stratum(width = 1/12, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#   scale_fill_brewer(type = "qual", palette = "Set1") +
#   ggtitle("UC Berkeley admissions and rejections, by sex and department")
# 
# data(vaccinations)
# vaccinations <- transform(vaccinations,
#                           response = factor(response, rev(levels(response))))
# ggplot(vaccinations,
#        aes(x = survey, stratum = response, alluvium = subject,
#            y = freq,
#            fill = response, label = response)) +
#   scale_x_discrete(expand = c(.1, .1)) +
#   geom_flow() +
#   geom_stratum(alpha = .5) +
#   geom_text(stat = "stratum", size = 3) +
#   theme(legend.position = "none") +
#   ggtitle("vaccination survey responses at three points in time")
# 
# ggplot(samkey_df %>% arrange(Disease, values),
#        aes(x = Disease, stratum = cell, alluvium = cell,
#            y = values,
#            fill = cell, 
#            label = cell)) +
#   scale_x_discrete(expand = c(.1, .1)) +
#   geom_flow() +
#   geom_stratum(alpha = .5) +
#   geom_text(stat = "stratum", size = 3) +
#   theme(legend.position = "none") +
#   ggtitle("vaccination survey responses at three points in time")

ggplot(data = samkey_df ,#%>% mutate(cell_type = paste0(Mixture, cell)),
       aes(x = Disease, y = (values), alluvium = cell)) +
  geom_alluvium(aes(fill = cell, colour = cell),
                alpha = 0.5, decreasing = FALSE) +
 # scale_x_continuous(breaks = seq(2003, 2013, 2)) +
  theme_bw() +
  #theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = -30, hjust = 0)) +
 scale_fill_manual(values = myPal) +
 # scale_fill_brewer(type = "qual", palette =  "stepped" ) +
  #scale_color_brewer(type = "qual", palette = "stepped") +
 # facet_wrap(~ cell_type, scales = "free") +
  ggtitle("Immune cells over the disease progression (Mean values)")

library(RColorBrewer)
cols <- colorRampPalette(brewer.pal( 8,"Accent"))
myPal <- cols(9)

library(lme4)
library(emmeans)

prop.reg_type <- glm( B.cells.naive ~ Disease , 
                 data = df_prop,
                 family = quasibinomial)
summary(prop.reg_type)

plot(fitted(prop.reg_type), resid(prop.reg_type), main="glm with quasibinomial")
abline(0,0)


betaReg_df <- data.frame(cell = colnames(cibersort_pt)[2:23])
glm_df <- data.frame(cell = colnames(cibersort_pt)[2:23])
glmmTMB_df <- data.frame(cell = colnames(cibersort_pt)[2:23])

colnames(cibersort_pt[,c(1:23, 32)])
levels(cibersort_pt$Disease)


df_prop <- cibersort_pt %>% 
           # filter(Disease %in% c( "Isolated HGD", "Progressor HGD w PDAC" , "PDAC")) %>% 
            filter(Disease %in% c( "Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC" )) %>% 
           # filter(Disease %in% c( "Progressor LGD w PDAC" , "Progressor HGD w PDAC" , "PDAC")) %>% 
            droplevels()  %>%
            mutate(Progressor_LGD_all = ifelse(Disease=="Progressor LGD w PDAC", "Progressor_LGD_all", 
                                               ifelse(Disease == "Progressor LGD w HGD", "Progressor_LGD_all", as.character(Disease))) ,
                   Disease = factor(Progressor_LGD_all, levels = c("Isolated LGD", "Progressor_LGD_all" )))

com1 <- "Isolated_HGD_Vs_ProgressorHGDwPDAC"
com2 <- "Isolated_HGD_Vs_PDAC"

com1 <- "Isolated_LGD_Vs_ProgressorLGDwHGD"
com2 <- "Isolated_LGD_Vs_ProgressorLGDwPDAC"

com1 <- "ProgressorLGDwHGD_Vs_ProgressorLGDwPDAC"
com2 <- ""


com1 <- "Isolated_LGD_Vs_ProgressorLGD_all"
com2 <- ""

com1 <- "ProgressorLGDwPDAC_Vs_ProgressorHGDwPDAC"
com2 <- "ProgressorLGDwPDAC_Vs_PDAC"

com1<- "ProgressorLGDwPDAC_Vs_PDAC"
com2 <- ""





  
cell_type = colnames(df_prop)[23]

df_prop$cell_type <- df_prop[[cell_type]]


plotly::ggplotly(ggplot( df_prop, aes(x= cell_type, fill = paste0( Disease),alpha = 0.01)) +
  geom_density() +
    labs(title =cell_type) )
  # scale_x_continuous( limits = c(0,1))  +
  # labs(title =cell_type) )


 for (cell_type in colnames(cibersort_pt)[c(2:23)] ) {
   df_prop$cell_type <- df_prop[[cell_type]]

# print("Glm with quasibinomial")
prop.reg <- glm( cell_type ~ Disease , 
                      data = df_prop,
                      family = quasibinomial)
p_glm <- coef(summary(prop.reg))[,"Pr(>|t|)"]

(summary(prop.reg) )

glm_df[glm_df$cell==cell_type, com1] <- p_glm[2]
if(com2 !=""){
glm_df[glm_df$cell==cell_type, com2] <- p_glm[3]
}

# qqnorm(residuals(prop.reg), main =  paste0( cell_type,": glm with quasibinomial"))
# qqline(residuals(prop.reg))
# 1-pchisq(prop.reg$deviance, prop.reg$df.residual)


m2 <- MASS::glmmPQL(cell_type/100~ Disease,random=~1|Mixture, family="quasibinomial",data=df_prop)



# qqnorm(residuals(m2), main =  paste0( cell_type,": glmm with glmmPQL \n-residual plot(quasibinomial) "))
# qqline(residuals(m2))
# qqnorm(ranef(m2)[[1]], main =  paste0( cell_type,": glmm with glmmPQL \n- random effect (quasibinomial) "))
# qqline(ranef(m2)[[1]])




print("beta regression")
fit_betareg <- tryCatch({betareg(cell_type ~ Disease , 
                       data = df_prop
                       #control=betareg.control("nlminb")
                       )
},
error=function(error_message) {
  })


if(!is.null(fit_betareg)) {
  if(summary(fit_betareg)$converged){

p_value<-coef(summary(fit_betareg))[[1]][,'Pr(>|z|)']
betaReg_df[betaReg_df$cell==cell_type, com1] <- p_value[2]
if(com2!=""){
betaReg_df[betaReg_df$cell==cell_type, com2] <- p_value[3]
}}}


# qqnorm(residuals(fit_betareg), main =  paste0( cell_type,": Beta distribution"))
# qqline(residuals(fit_betareg))




df_prop$cell_type[df_prop$cell_type ==0]<- 0.0000001
print("glmmTMB with beta distribution and random effect for within subject")
fit_TMB <- glmmTMB::glmmTMB(cell_type ~ Disease + (1| Mixture), 
             data = df_prop,
             family = beta_family(link = "logit")
            # control = glmmTMBControl(optimizer= optim)
             )
# fit_TMB_NR <- glmmTMB::glmmTMB(cell_type ~ Disease ,
#                             data = df_prop,
#                             family = beta_family(link = "logit")
#                              )
p_glmmTMB <-coef(summary(fit_TMB))$cond[,'Pr(>|z|)']

 
glmmTMB_df[glmmTMB_df$cell==cell_type, com1] <- p_glmmTMB[2]
if(com2 != ""){
glmmTMB_df[glmmTMB_df$cell==cell_type, com2] <- p_glmmTMB[3]
}

# 
# qqnorm(residuals(fit_TMB), main = paste0( cell_type,": Beta disribution \nwith fixed effectglmmTMB"))
# qqline(residuals(fit_TMB))
# 
# qqnorm(ranef(fit_TMB)[[1]][[1]][[1]], main = paste0( cell_type,": Beta disribution \nwith random effect glmmTMB"))
# qqline(ranef(fit_TMB)[[1]][[1]][[1]])



}
  
  library(openxlsx)
wb <- createWorkbook("glm")
addWorksheet(wb ,  sheet = "glm")
addWorksheet(wb ,  sheet = "betaReg")
addWorksheet(wb ,  sheet = "glmmTMB")
writeData(wb, betaReg_df,  sheet = "betaReg")
writeData(wb, glmmTMB_df, sheet = "glmmTMB")
writeData(wb, glm_df, sheet= "glm")

saveWorkbook(wb, "data/ImmuneCell_p_value_glm.xlsx", overwrite = TRUE)
  


ratios <- cibersort_pt %>% 
  mutate(Lymphocytes = rowSums(cibersort_pt[2:11] ),
         B_cells = rowSums(cibersort_pt[2:4]),
         T_cells = rowSums(cibersort_pt[5:11]),
         Macrophages = rowSums(cibersort_pt[15:17]),
         APCs = rowSums(cibersort_pt[c(15:19)]),
         Tcells_to_APCs = T_cells/APCs,
         M1_to_M2 = Macrophages.M1/Macrophages.M2,
         plama_to_Bnaive = Plasma.cells/B.cells.naive,
         TMemActiv_To_TMemRest = T.cells.CD4.memory.activated/T.cells.CD4.memory.resting,
         TCD8_to_TMemRest = T.cells.CD8/T.cells.CD4.memory.resting,
         TCD8_to_TMemActiv = T.cells.CD4.memory.activated/T.cells.CD8,
         TCD8_to_Treg = T.cells.CD8/T.cells.regulatory..Tregs.,
         BMem_to_Treg = B.cells.memory/T.cells.regulatory..Tregs.,
         M1_to_DenActiv = Macrophages.M1/Dendritic.cells.activated,
         M1_to_NKActiv= Macrophages.M1/NK.cells.activated
         
  ) %>% 
  filter(Disease %in% c( "Isolated HGD", "Progressor HGD w PDAC" , "PDAC")) %>% 
  droplevels()  %>%
  select(-c(2:24))
colnames(ratios)

ratios[ratios=="Inf"] <- NA
ratio_type = colnames(ratios)[24]


ratios$ratio_type <- ratios[[ratio_type]]


plotly::ggplotly(ggplot( ratios, aes(x= ratio_type, fill = paste0( Disease),alpha = 0.01)) +
                   geom_density() +
                   labs(title =ratio_type) )
# scale_x_continuous( limits = c(0,1))  +
# labs(title =cell_type) )
prop.reg <- lm( ratio_type ~ Disease , 
                 data = ratios )
# p_glm <- coef(summary(prop.reg))[,"Pr(>|t|)"]

(summary(prop.reg) )







library(caret)
library(dplyr)

dfexample <- data.frame(
  subject = c('English', 'English', 'Math', 'Science'),
  enrollment = c(100,200,50,70),
  white = c(0.5,0.5,0.6,0.7),
  black = c(0.25,0.20, 0.10, 0.25),
  hispanic = c(0.25, 0.30, 0.30, 0.05),
  classid = c('1a','3f','4d','5a')
);

# replicating data frame rows to make our example works    
dfexample = dfexample[rep(seq_len(nrow(dfexample)), each = 20), ]
dfexample <- samkey_df %>% select("Mixture", "Disease", "cell", "values") %>% pivot_wider(names_from = cell, values_from = values)
trainIndex <- createDataPartition(dfexample$Disease, p = .6, 
                                  list = FALSE, 
                                  times = 1)
dataTrain <- dfexample[ trainIndex,]
dataTest  <- dfexample[-trainIndex,]   #   colnames(dfexample)
modelFit <- train(Disease ~  B_cell + Den_cell + Eosinophils + M + Mast_cell + Monocytes + 
                    NK_cell + Neutrophils + T_cell  , data = dataTrain, 
                  method = "randomforest",  
                  #method = "gbm",  
                  verbose = FALSE
)

print(modelFit)
predictions <- predict(modelFit, newdata = dataTest)

cm = confusionMatrix(predictions, as.factor(dataTest$Disease))

print(cm)





### MiXCR.. 3X3 analysis

mixcer <- read.table("data/BE_final_postanalysis.diversity.IGK.tsv", sep = "\t", header = TRUE)

mixcer_sample <- mixcer[1:2] %>%
         separate(sample, sep = "A_", c("Sample_ID", NA)) %>%
         mutate(Sample_ID = paste0(Sample_ID, "A")) %>% 
  left_join(ciber_sample %>% dplyr::select(Sample_ID, Grade, Patient_id, Path_diagnosis,Main_duct, Category, Disease), by = "Sample_ID") 

str(mixcer_sample)

df_plot <- mixcer_sample %>%
  filter(!Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
  droplevels() %>%
  mutate(
    Disease = str_wrap(Disease, 10),
    Disease =factor(Disease, levels= c(str_wrap(levels(.$Disease), 10))),
    groups = ifelse(str_detect(Disease, "PDAC"), "PDAC", "non_PDAC"),
    groups = factor(groups, levels = c( "non_PDAC", "PDAC"))
    # Disease = ifelse(str_detect(Disease, "Progressor"),"Progressor LGD", "Isolated LGD") ,
    #      Disease = factor(Disease , levels = c("Isolated LGD", "Progressor LGD"))
  )


p_values <- df_plot %>% 
  wilcox_test( Observed.diversity ~ groups ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_values  <- p_values %>% add_xy_position(x = "Disease", group = "groups")
p_values$xmin <- 2
p_values$xmax <- 5
p_values$y.position <- p_values$y.position +50
p_values$p.adj.signif[p_values$p.adj.signif =="ns"] <- ''

st_line <- data.frame(xmin = c(1,4),
                      xmax = c(3,6),
                      y.position = c(p_values$y.position-25, p_values$y.position-25 ),
                      p="",
                      group1=1,
                      group2 =2)
color <- c( "#2F56B5",  "#2295C7", "#00A86B", "#FCC92A", "#F4691E", "#AD0000")

# color <- c("#48217390","#433E8590", "#38598C90","#2BB07F90","#85D54A90", "#FDE725FF") ##2BB07F50")
p <- ggboxplot(df_plot , x = "Disease", y = "Observed.diversity",
               # color = "Disease",
               fill = "Disease",palette = color,
               add = "jitter",
               legend.title ="")
#  Add p-value

# Change method
# p + stat_compare_means(# comparisons = list(c("Isolated LGD", "Isolated HGD") )
#
#               #list(c(levels(df_plot$Disease)[c(4,5)]),c(levels(df_plot$Disease)[c(4,6)] ))
#                  )


p + #stat_compare_means(label.y = 0.35 ) +
  labs(x="", y = "Observed diversity - IGK") + # ylim(c(0, 0.14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) +
  stat_pvalue_manual(
    p_values,  label = "{p.adj}{p.adj.signif}" #, tip.length = 0.02
    #step.increase = 0.05
  ) +
  stat_pvalue_manual(
    st_line, tip.length = 0
    # step.increase = 0.05
  )


        