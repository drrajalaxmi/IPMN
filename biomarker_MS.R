### Biomarker analysis NCI Cancer Moon Shot

# Box plot comparing the low grade associated IPMN to high grade associated samples
# statistical test wilcox test with p < 0.05 as significant value

## input data condition (sample metadata) and counts (count_matrix) from input.R
## condition has 5 samples excluded: "BE02_65A"  "BE02_77A"  "BE02_78A"  "BE02_80A"  "BE02_81A"  "BE02_120A"

### gene lists
checkpoint_genes <- c( ### biomarker for the HGD and LGD
  # "KRT19", "MUC1", "CEACAM6", "MUC5AC", "CLDN18", "CDH1" ,"CTHRC1"
  # "GNMT", "IGF2BP3",  "CLDN4","HOXB3", "TRIM29",
  # "LAMC2", "CD55","AREG", "PGC", "NKX6-2"
  
  "PARP1", "MUC5AC", "GNMT", "CTHRC1", "IGF2BP3", "NKX6-2"
  
  
  ### immune checkpoints
  # "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1", "VTCN1"
  
  ### antigen presentation
  # "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
  # "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "CALR"
  
  ## Lipid biosynthesis
  # "FA2H", "CERS2", "CERS6",  "UGT8", "GAL3ST1"
)



condition_subset <- condition_ms[condition_ms$Sample_ID %in% colnames(df_count_ms),]

# condition_subset <- condition_subset %>%
#   mutate(Category = factor(Category , levels = c("LGD",  "HGD" , "PDAC")),
#          Disease = factor(Disease, levels = c("Normal-LGD" ,"Normal-HGD", "Normal-PDAC",
#                                               "Isolated LGD",
#                                               "Isolated HGD" ,
#                                               "Progressor LGD w HGD" , 
#                                               "Progressor LGD w PDAC",
#                                               "Progressor HGD w PDAC", "PDAC"))) %>% 
#   arrange(factor(Category),Path_ID,factor(Disease))

counts.DEseq = DESeqDataSetFromMatrix(countData=round(df_count_ms), colData = condition_subset, design = ~ 1)
dds_cp <- estimateSizeFactors(counts.DEseq)
checkpoint_expr <- counts(dds_cp, normalized = TRUE)
# rownames(norm_counts) = unlist(lapply(strsplit(rownames(norm_counts), '\\.'), '[[',1))
# checkpoint_expr <- norm_counts[rownames(norm_counts) %in% ids_IC$ENSEMBL, ]


# 
symbol_df <- checkpoint_expr %>% as.data.frame() %>% 
  filter(rownames(.) %in% checkpoint_genes)

## faceted boxplot for all gbiomarkers
ggplot(t(symbol_df) %>% as.data.frame() %>%
         rownames_to_column("Sample_ID") %>%
         pivot_longer(!Sample_ID,values_to = "counts", names_to = "gene") %>%
         left_join(condition_subset) ,
       aes(Disease, counts ) ) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "MS") +
  theme(axis.text.x = element_text(angle = 45, size= 8,vjust = 1, hjust=1)) +
  facet_wrap(~gene, scales = "free_y", nrow = 2, strip.position = "top")


### for comparison
df_conuts_sample <- t(symbol_df) %>% as.data.frame() %>% rownames_to_column("Sample_ID") %>%
  left_join(condition_subset) %>%            
  filter(!str_detect(Disease ,"normal")) %>%
  droplevels()
df_plot <- df_conuts_sample %>%
  mutate( groups = ifelse(str_detect(Disease, "PDAC"), "PDAC", "non_PDAC"),
          groups = factor(groups, levels = c( "non_PDAC", "PDAC")),
          Disease = str_wrap(Disease, 10),
          Disease =factor(Disease, levels= c(str_wrap(levels(df_conuts_sample$Disease), 10)))
  )




for( gene_name in checkpoint_genes){
  # gene_name = "IGF2BP3"
  p_values <- df_plot %>% 
    dplyr::rename("gene_name" = gene_name) %>%
    wilcox_test( gene_name ~ groups ) %>%
    adjust_pvalue( method = "bonferroni" ) %>%
    add_significance("p.adj")

  
  p_values  <- p_values %>% add_xy_position(x = "Disease", group = "groups") %>%
    mutate(y.position = max(attributes(p_values)[[4]]$data$gene_name) * 1.08 )
  p_values$xmin <- 1.5
  p_values$xmax <- 4
  
  st_line <- data.frame(xmin = c(1,3),
                        xmax = c(2,5),
                        y.position = c(rep(0.97*p_values$y.position, 2)),
                        p="",
                        group1 = 1,
                        group2 = 2)

  
  color <- c(#"#2F56B5", 
             "#2295C7", "#00A86B", "#FCC92A", "#F4691E",  "#AD0000")
  p <- ggboxplot(df_plot %>%
                   dplyr::rename("gene_name" = gene_name) , x = "Disease", y = "gene_name",
                 fill = "Disease", palette = color,
                 add = "jitter", legend.title ="") +
    labs(x="MS", y = gene_name) + # ylim(c(0, 150)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) 
  if(p_values$p.adj.signif !='ns'){
    p <- p +
      stat_pvalue_manual( p_values,  label = "{p.adj}" ) +
      stat_pvalue_manual( st_line, tip.length = 0  ) 
  }
  
  print(p)
  
}



### additional custom plots

gene_name = "MUC5AC"
p_values <- df_plot %>% 
  rename("gene_name" = gene_name) %>%
  wilcox_test( gene_name ~ groups ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_values  <- p_values %>% add_xy_position(x = "Disease", group = "groups") %>%
  mutate(y.position = max(attributes(p_values)[[4]]$data$gene_name) * 1.48 )
p_values$xmin <- 2
p_values$xmax <- 5

st_line <- data.frame(xmin = c(1,4),
                      xmax = c(3,6),
                      y.position = c(rep(0.97*p_values$y.position, 2)),
                      p="",
                      group1 = 1,
                      group2 = 2)

p_additional <- df_plot %>% 
  rename("gene_name" = gene_name) %>%
  filter(Disease != "Isolated\nLGD") %>%
  mutate(additional_group = ifelse(Disease=="Isolated\nHGD", "gr1", "gr2"),
         additional_group = factor(additional_group, levels= c("gr1", "gr2"))) %>%
  wilcox_test( gene_name ~ additional_group ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_additional  <- p_additional %>% add_xy_position(x = "Disease", group = "groups") %>%
  mutate(y.position = max(attributes(p_additional)[[4]]$data$gene_name) * 1.28 )
p_additional$xmin <- 2
p_additional$xmax <- 5
p_additional$p.adj.signif[p_additional$p.adj.signif =="ns"] <- ''

st_line_additional <- data.frame(xmin = c(2,3),
                                 xmax = c(2,6),
                                 y.position = c(rep(0.97*p_additional$y.position, 2)),
                                 p="",
                                 group1 = 1,
                                 group2 = 2)

p_additional_2 <- df_plot %>% 
  rename("gene_name" = gene_name) %>%
  filter(Disease %in% c( "Isolated\nHGD", "Progressor\nLGD w PDAC" )) %>% droplevels() %>%
  # mutate(additional_group = ifelse(Disease=="Isolated\nHGD", "gr1", "gr2"),
  #        additional_group = factor(additional_group, levels= c("gr1", "gr2"))) %>%
  wilcox_test( gene_name ~ Disease ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_additional_2  <- p_additional_2 %>% add_xy_position(x = "Disease", group = "groups") %>%
  mutate(y.position = p_additional$y.position * 0.88 )
p_additional_2$xmin <- 2
p_additional_2$xmax <- 4
p_additional_2$p.adj.signif[p_additional_2$p.adj.signif =="ns"] <- ''

color <- c("#2F56B5", "#2295C7", "#00A86B", "#FCC92A", "#F4691E",  "#AD0000")
p <- ggboxplot(df_plot , x = "Disease", y = gene_name,
               fill = "Disease", palette = color,
               add = "jitter", legend.title ="") +
  labs(x="") + # ylim(c(0, 400)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) 

p +
  stat_pvalue_manual( p_values,  label = "{p.adj}{p.adj.signif}" ) +
  stat_pvalue_manual( st_line, tip.length = 0  )  +
  stat_pvalue_manual( p_additional,  label = "{p.adj}{p.adj.signif}" ) +
  stat_pvalue_manual( st_line_additional, tip.length = 0  )  +
  stat_pvalue_manual( p_additional_2,  label = "{p.adj}{p.adj.signif}" ) +
  stat_pvalue_manual( st_line_additional, tip.length = 0  )  

# stat_compare_means(# comparisons = list(c("Isolated LGD", "Isolated HGD") )
#   position = 10,
#   comparisons =list(c(levels(df_plot$Disease)[c(2,4)]) ))



gene_name <- "CTHRC1"
unique(df_conuts_sample$Disease)
df_plot <- df_conuts_sample %>%            
  # filter(Disease %in% c("Normal-PDAC", "Progressor LGD w PDAC")) %>%
  filter(Disease %in% c("Normal-PDAC", "Progressor HGD w PDAC")) %>%
  droplevels() %>%
  mutate( 
    Disease = str_wrap(Disease, 10),
    Disease =factor(Disease, levels= c(str_wrap(levels(df_conuts_sample$Disease)[c(3,8)], 10)))
  )


p_additional_3 <- df_plot %>% select(Path_ID,Disease, CTHRC1 ) %>%
  pivot_wider(names_from = Disease, values_from =  CTHRC1)%>%
  drop_na() %>%
  pivot_longer(!Path_ID, values_to = "CTHRC1", names_to = "Disease")


color <- c("grey",# "#2F56B5", "#2295C7", "#00A86B",
           # "#FCC92A" ,
           "#F4691E" #  "#AD0000"
)
ggboxplot(p_additional_3 , x = "Disease", y = "CTHRC1", group = "Path_ID",
          fill = "Disease", palette = color,
          add = "jitter", legend.title ="") +
  labs(x="") + # ylim(c(0, 400)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) +
  stat_compare_means( 
    # comparisons = list(c("Normal-PDAC", "Progressor\nLGD w PDAC") ),
    comparisons = list(c("Normal-PDAC", "Progressor\nHGD w PDAC") ),
    paired =TRUE)






## PCA plot with expression color 
ipmn_subset <- condition_ms%>%
  filter(!str_detect(Disease, "normal"))  %>%
  filter(Sample_ID %in% colnames(checkpoint_expr))

tpm_ipmn <- checkpoint_expr %>% as.data.frame() %>%
  dplyr::select(ipmn_subset$Sample_ID)


tpm_ipmn_subset <- checkpoint_expr %>% as.data.frame() %>%
  filter(rownames(.) %in% checkpoint_genes)

# ipmn_subset <- ipmn_subset %>% left_join(tpm_ipmn_subset)

cnt_pca <- prcomp(t(tpm_ipmn))
#summary(cnt_pca)

# plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (all Grades)") 
# pca.plot <- autoplot(cnt_pca,
#                      data = t(counts)) 

colours <- c("Isolated LGD" = "#2F56B5",#115C80", 
             "Isolated HGD" ="#2295C7", #309898" ,#479685", 
             "Progressor LGD w HGD" = "#00A86B", #FF9F00", 
             "Progressor LGD w PDAC" = "#FCC92A", #FF9F00", 
             "Progressor HGD w PDAC" = "#F4691E", 
             "PDAC" = "#AD0000")


colours <- c(#"#2F56B5",
             "#2295C7", "#00A86B", "#FCC92A", "#F4691E", "#AD0000")

autoplot(cnt_pca , data = df_plot , color ="Disease", size = 2.5  ) + 
  geom_point(color = "white") +
  scale_color_manual(values = colours) 



ggplot(as.data.frame(cnt_pca$x[, 1:2]),aes( x= PC1, y = PC2 , size = 2.5)) +
    geom_point(color = "white", size = 2.8) +
   geom_point(data = df_plot , color ="Disease", size = 2.5 )

for(gene_name in checkpoint_genes){
  
  
  p <- autoplot(cnt_pca , data = df_plot %>% mutate(gene = log(df_plot[[gene_name]]+1)), color = "gene", size = 2.5  ) + 
    scale_color_continuous(name = paste0(gene_name, " (log)"),low = "grey", high = "red") 
  
  print(p)
  
}

condition_nonNorm <- condition_subset %>% 
  filter(!str_detect(Disease, "Normal"))
cnt_pca <- prcomp(checkpoint_expr%>% as.data.frame() %>% select(df_plot$Sample_ID) %>% t())


dim(norm_counts)
#label.size = 3)#label.size = 3)fill = 

