# Box plot comparing the low grade associated IPMN to high grade associated samples
# statistical test wilcox test with p < 0.05 as significant value

## input data condition (sample metadata) and counts (count_matrix) from input.R
## condition has 5 samples excluded: "BE02_65A"  "BE02_77A"  "BE02_78A"  "BE02_80A"  "BE02_81A"  "BE02_120A"

### gene lists
checkpoint_genes <- c( ### biomarker for the HGD and LGD
                       # "KRT19", "MUC1", "CEACAM6", "MUC5AC", "CLDN18", "CDH1" ,"CTHRC1"
                           # "GNMT", "IGF2BP3",  "CLDN4","HOXB3", "TRIM29",
                           # "LAMC2", "CD55","AREG", "PGC", "NKX6-2",
                           "PARP1", "MUC5AC", "GNMT", "CTHRC1", "IGF2BP3", "NKX6-2"

                        ### immune checkpoints
                        # "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1", "VTCN1"

                        ### antigen presentation
                        # "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
                          # "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "CALR"

                        ## Lipid biosynthesis
                        # "FA2H", "CERS2", "CERS6",  "UGT8", "GAL3ST1"
                      )

## fetching the gene symbols
ids_IC<-bitr( checkpoint_genes, fromType="SYMBOL",
             toType="ENSEMBL",
             OrgDb="org.Hs.eg.db" ) #NA values drop by default

condition_subset <- condition[condition$Sample_ID %in% colnames(counts),]
condition_subset <- condition_subset %>%
                    mutate(Category = factor(Category , levels = c("LGD",  "HGD" , "PDAC")),
                           Disease = factor(Disease, levels = c("Normal-LGD" ,"Normal-HGD", "Normal-PDAC",
                                                                "Isolated LGD",
                                                                "Isolated HGD" ,
                                                                "Progressor LGD w HGD" , 
                                                                "Progressor LGD w PDAC",
                                                                "Progressor HGD w PDAC", "PDAC"))) %>% 
                    arrange(factor(Category),Path_ID,factor(Disease))
                    

# library(DESeq2)
# set.seed(1000)
counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts), colData = condition_subset, design = ~ 1)
dds_cp <- estimateSizeFactors(counts.DEseq)
norm_counts <- counts(dds_cp, normalized = TRUE)
rownames(norm_counts) = unlist(lapply(strsplit(rownames(norm_counts), '\\.'), '[[',1))
checkpoint_expr <- norm_counts[rownames(norm_counts) %in% ids_IC$ENSEMBL, ]



symbol_df <- checkpoint_expr %>% as.data.frame() %>% rownames_to_column("ENSEMBL")  %>% 
            left_join( ids_IC ) %>% 
            column_to_rownames("SYMBOL") %>%
            dplyr::select(-ENSEMBL) %>%
            dplyr::select(condition_subset$Sample_ID)

## faceted boxplot for all gbiomarkers
ggplot(t(symbol_df) %>% as.data.frame() %>%
         rownames_to_column("Sample_ID") %>%
         pivot_longer(!Sample_ID,values_to = "counts", names_to = "gene") %>%
         left_join(condition_subset) ,
       aes(Disease, counts ) ) +
  geom_boxplot() + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, size= 8,vjust = 1, hjust=1)) +
  facet_wrap(~gene, scales = "free_y", nrow = 2, strip.position = "top")


### for comparison
df_conuts_sample <- t(symbol_df) %>% as.data.frame() %>% rownames_to_column("Sample_ID") %>%
            left_join(condition_subset) 
df_plot <- df_conuts_sample %>%            
            filter(!Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
            droplevels() %>%
            mutate( groups = ifelse(str_detect(Disease, "PDAC"), "PDAC", "non_PDAC"),
                    groups = factor(groups, levels = c( "non_PDAC", "PDAC")),
                    Disease = str_wrap(Disease, 10),
                    Disease =factor(Disease, levels= c(str_wrap(levels(df_conuts_sample$Disease), 10)))
          )




for( gene_name in checkpoint_genes){
  # gene_name = "TRIM29"
  p_values <- df_plot %>% 
              dplyr::rename("gene_name" = gene_name) %>%
              wilcox_test( gene_name ~ groups ) %>%
              adjust_pvalue( method = "bonferroni" ) %>%
              add_significance("p.adj")

  p_values  <- p_values %>% add_xy_position(x = "Disease", group = "groups") %>%
                mutate(y.position = max(attributes(p_values)[[4]]$data$gene_name) * 1.08 )
  p_values$xmin <- 2
  p_values$xmax <- 5

  st_line <- data.frame(xmin = c(1,4),
                        xmax = c(3,6),
                        y.position = c(rep(0.97*p_values$y.position, 2)),
                        p="",
                        group1 = 1,
                        group2 = 2)

  color <- c("#2F56B5", "#2295C7", "#00A86B", "#FCC92A", "#F4691E",  "#AD0000")
  p <- ggboxplot(df_plot %>% dplyr::rename("gene_name" = gene_name) , x = "Disease", y = "gene_name",
               fill = "Disease", palette = color,
               add = "jitter", legend.title ="") +
               labs(x="", y = gene_name) +  
               theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) 
  if(p_values$p.adj.signif !='ns'){
    p <- p +
        ylim(c(0, p_values$y.position*1.02)) +
           stat_pvalue_manual( p_values, # label = "{p.adj}{p.adj.signif}"
                               label = "{p.adj}") +
           stat_pvalue_manual( st_line, tip.length = 0  ) 
     }

  print(p)

}



### additional custom plots

gene_name = "MUC5AC"
p_values <- df_plot %>% 
  dplyr::rename("gene_name" = gene_name) %>%
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
  dplyr::rename("gene_name" = gene_name) %>%
  # filter(Disease != "Isolated\nLGD") %>%
  filter(str_detect(Disease, "LGD")) %>%
  mutate(
         # additional_group = ifelse(Disease=="Isolated\nHGD", "gr1", "gr2"),
         additional_group = ifelse(str_detect(Disease, "PDAC"), "gr2", "gr1"), 
         additional_group = factor(additional_group, levels= c("gr1", "gr2"))) %>%
  wilcox_test( gene_name ~ additional_group ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_additional  <- p_additional %>% add_xy_position(x = "Disease", group = "groups") %>%
  mutate(y.position = max(attributes(p_additional)[[4]]$data$gene_name) * 1.28 )
p_additional$xmin <- 2
p_additional$xmax <- 4
p_additional$p.adj.signif[p_additional$p.adj.signif =="ns"] <- ''

st_line_additional <- data.frame(xmin = c(1,4),
                      xmax = c(3,4),
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
  labs(x="") + ylim(c(0, p_values$y.position*1.02)) + # ylim(c(0, 400)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) 
                 
 p +
    stat_pvalue_manual( p_values,  label = "{p.adj}" ) +
    stat_pvalue_manual( st_line, tip.length = 0  )  +
    stat_pvalue_manual( p_additional,  label = "{p.adj}" ) +
    stat_pvalue_manual( st_line_additional, tip.length = 0.05  )  
   # stat_pvalue_manual( p_additional_2,  label = "{p.adj}{p.adj.signif}" ) +
   # stat_pvalue_manual( st_line_additional, tip.length = 0  )  
   
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
 ipmn_subset <- condition%>%
      filter(!str_detect(Disease, "Normal")) %>%
      mutate(Disease = factor(Disease, levels = c("Isolated LGD" ,
                                                  "Progressor LGD w HGD" ,
                                                  "Isolated HGD" ,
                                                  
                                                  "Progressor LGD w PDAC" ,
                                                  "Progressor HGD w PDAC",
                                                  "PDAC" ) ))
 
tpm_ipmn <- counts_tpm %>%
                dplyr::select(ipmn_subset$Sample_ID)

   
tpm_ipmn_subset <- counts_tpm %>% as.data.frame() %>%
   dplyr::select(ipmn_subset$Sample_ID) %>% 
   rownames_to_column("ENSEMBL")  %>% 
   mutate(ENSEMBL = unlist(lapply(strsplit(ENSEMBL, '\\.'), '[[',1)) ) %>%
   filter(ENSEMBL %in% ids_IC$ENSEMBL) %>%
   left_join( ids_IC ) %>% 
   select(-ENSEMBL) %>% 
   column_to_rownames("SYMBOL") %>% t() %>% as.data.frame() %>%
   rownames_to_column("Sample_ID")  
   
ipmn_subset <- ipmn_subset %>% left_join(tpm_ipmn_subset)

 cnt_pca <- prcomp(t(tpm_ipmn))
 #summary(cnt_pca)
 
 # plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (all Grades)") 
 # pca.plot <- autoplot(cnt_pca,
 #                      data = t(counts)) 

colours <- c("Isolated LGD" = "#2F56B5",#115C80",
             "Progressor LGD w HGD" = "#2295C7",
              "Isolated HGD" ="#00A86B", #309898" ,#479685", 
              "Progressor LGD w PDAC" = "#FCC92A", #FF9F00", 
              "Progressor HGD w PDAC" = "#F4691E", 
              "PDAC" = "#AD0000")
 
 
colours <- c("#2F56B5", "#2295C7", "#00A86B", "#FCC92A", "#F4691E", "#AD0000")

# library(ggfortify) 
for(gene_name in checkpoint_genes[10]){
 

p <- autoplot(cnt_pca , data = df_plot %>% mutate(gene = log(df_plot[[gene_name]]+1)), color = "gene", size = 2.5  ) + 
                    scale_color_continuous(name = paste0(gene_name, " (log)"),low = "grey", high = "red") 
                    
print(p)

}

condition_nonNorm <- condition_subset %>% 
         filter(!str_detect(Disease, "Normal"))
cnt_pca <- prcomp(norm_counts%>% as.data.frame() %>% select(condition_nonNorm$Sample_ID) %>% t())


dim(norm_counts)
#label.size = 3)#label.size = 3)fill = 
autoplot(cnt_pca , data = ipmn_subset , color = "Disease", size = 2.5  ) + 
  scale_color_manual(values = colours) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank() )
  # scale_color_continuous(name = paste0(gene_name, " (log)"),low = "grey", high = "red") 



Ca_cor <- tpm_ipmn%>% as.data.frame() %>% select(condition_nonNorm$Sample_ID) %>% as.matrix() %>% cor()

Ca_cor <- cor(as.matrix(tpm_ipmn))
dim(condition)
y_axis_colors <- rep(c("red"), each = nrow(Ca_cor))  # Example: alternating red and blue




ann_colors <- list(
  Category = c("LGD" = "#2F56B5", "HGD" = "#00A86B", "PDAC" = "#AD0000"),
  
  Disease= c(
             "Normal-LGD" =  "grey90" ,
             "Normal-HGD" =  "grey60" ,
             "Normal-PDAC" =  "grey10" ,
             
             "Isolated LGD" ="#2F56B5",
             "Progressor LGD w HGD" = "#2295C7",
             "Isolated HGD"= "#00A86B" ,
             
             "Progressor LGD w PDAC" ="#FCC92A",
             "Progressor HGD w PDAC" = "#F4691E",
             "PDAC" = "#AD0000" )
)

color_df <- data.frame(Disease= c("Isolated LGD" ,
                                  "Progressor LGD w HGD" ,
                                  "Isolated HGD" ,
                                  
                                  "Progressor LGD w PDAC" ,
                                  "Progressor HGD w PDAC",
                                  "PDAC" ),
                       colors = c("#2F56B5", "#2295C7", "#00A86B", "#FCC92A", "#F4691E", "#AD0000")
) 
color_df <- data.frame(Category= c("LGD", "HGD",
                                  "PDAC" ),
                       colors = c("#2F56B5", "#00A86B",  "#AD0000")
) 
df_color <- data.frame(Sample_ID = rownames(Ca_cor)) %>% left_join(condition %>% dplyr::select(Sample_ID, Disease) %>% 
                                                        left_join(color_df) %>% filter(!is.na(colors)) ) %>%
     mutate(colors = ifelse(is.na(colors), "black", colors))

df_color <- data.frame(Sample_ID = rownames(Ca_cor)) %>% left_join(condition %>% dplyr::select(Sample_ID, Category) %>% 
                                                                     left_join(color_df) %>% filter(!is.na(colors)) ) %>%
  mutate(colors = ifelse(is.na(colors), "black", colors))
# heatmaply::heatmaply(Ca_cor,
                     # heatmap_layers = theme(axis.text.y = element_text(color = c("red",rep(c("red", "blue"), (nrow(Ca_cor)-1)/2)))) #, margins = c(5, 5), scale = "none")

#create data frame for annotations

dfh<-data.frame(Sample_ID=as.character(colnames(Ca_cor)))%>%
  left_join(condition %>% select(Sample_ID,  Disease, Category)) %>%
 column_to_rownames("Sample_ID")

dfh
library(pheatmap)
library(ComplexHeatmap)
stats::heatmap(Ca_cor,RowSideColors = df_color$colors, labCol = NA)
pheatmap::pheatmap(Ca_cor,# cutree_rows = 4, 
                   color = colorRampPalette(
  c( "yellow", "red"))(10),
                   annotation_row = dfh, border_color = NA,drop_levels = TRUE,
                   annotation_colors = ann_colors, angle_col = "45", show_colnames = FALSE,
                   fontsize_row = 7) 

Ca_cor %>% as.data.frame() %>% rownames_to_column("X") %>%
  pivot_longer(!X, names_to = "Y", values_to = "cor") %>% 
  mutate(y_axis_colors = rep(c("red", "blue"), each = n()/2)  ) %>%

ggplot(aes(X, Y, fill = cor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") 
  # theme(axis.text.y = element_text(color = color_vector)) +
  # labs(y = "Y-axis with Colors", x = "X-axis")
 


require(graphics); require(grDevices)
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
stats::heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "specification variables", ylab =  "Car Models",
              main = "heatmap(<Mtcars data>, ..., scale = \"column\")")
