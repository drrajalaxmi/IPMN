## ssgsea

# BiocManager::install("GSVA")
library(GSVA)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggfortify)
library(umap)
# library(biomaRt)
library(clusterProfiler)
library(tidyr)
library(tibble)




### input from input.R

# counts_tpm tmp matrix and condition is the metadata

## generating the ENTREZID and adding it to the tpm matrix
rownames(counts_tpm) = unlist(lapply(strsplit(rownames(counts_tpm), '\\.'), '[[',1))
ids<-bitr(rownames(counts_tpm), fromType="ENSEMBL",
           toType="ENTREZID",
           OrgDb="org.Hs.eg.db") #NA values drop by default


counts_entreID <- counts_tpm %>% 
  rownames_to_column(var ="ENSEMBL") %>% 
  mutate(ENSEMBL = unlist(lapply(strsplit( ENSEMBL, '\\.'), '[[',1)) ) %>%
  left_join(ids) %>%
  filter(!is.na(ENTREZID)) 
counts_entreID <- counts_entreID[!(duplicated(counts_entreID$ENTREZID)),] 
rownames(counts_entreID) <- counts_entreID$ENTREZID
counts_entreID <- counts_entreID %>% dplyr::select(-c("ENSEMBL", "ENTREZID"))

# condition_subset <- condition %>% filter(!str_detect(Disease, "Normal"))
# counts_entreID <- counts_entreID %>% dplyr::select(condition_subset$Sample_ID)

counts_entreID <- counts_entreID %>% dplyr::select(condition$Sample_ID)
#ssgsea pathway enrichment

pathways.hallmark <- fgsea::gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/h.all.v2024.1.Hs.entrez.gmt") # for human

gsvapar <- ssgseaParam(as.matrix(counts_entreID), pathways.hallmark)

ssgsea_hallmark <- gsva(gsvapar )

stats::heatmap(ssgsea_hallmark, main = "ssGSEA: Hallmark Pathways")
        
 
## careting UMAP                                                      
 #27 , 37,  62  , 69 , 80, 86, 61, 95, 96, 99, 108, 109, 114, 151, 168, 179, 181
# 11760  , 16898

set.seed(1876)
umap_res <- umap(t(ssgsea_hallmark)) # %>% as.data.frame() %>% dplyr::select(condition$Sample_ID[condition$Category %in% c("PDAC","HGD")]) ) %>% as.matrix() )
umap_df <- data.frame(umap_res$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("Sample_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::left_join(condition, by = "Sample_ID") %>%
  mutate(Disease = factor(Disease, levels = c("Normal-LGD", "Normal-HGD",
                                               "Normal-PDAC",
                                              "Isolated LGD",
                                              "Progressor LGD w HGD",
                                              "Isolated HGD", 
                                              "Progressor LGD w PDAC",
                                              "Progressor HGD w PDAC",
                                              "PDAC") ))

color <- c( "grey95", "grey80", "grey60", 
           "#2F56B5",#115C80", 
           "#2295C7", #309898" ,#479685", 
           "#00A86B", #FF9F00", 
           "#FCC92A", #FF9F00", 
           "#F4691E", 
           "#AD0000")

plotly::ggplotly(ggplot(umap_df,aes(
  x = X1,
  y = X2,
  label =Sample_ID,
 color= Disease, stroke = 0.2,
   #color= log(Neutrophils)
) )+ 
  geom_point(data= umap_df %>% filter(Grade=="normal"), 
             color = "grey", size = 3.2) +
   geom_point(size = 3) +
  # scale_colour_manual(values = c("Normal-LGD"  = "gray70",
  #                                "Normal-HGD" = "#90c0ce",
  #                                "Normal-PDAC" ="#FFF600",
  #                                "Isolated LGD" ="gray40",
  #                                "Progressor LGD w HGD" = "#7AD7F0",
  #                                "Isolated HGD"= "#266CA9",
  #                                "Progressor LGD w PDAC"= "#FFA80F",
  #                                "Progressor HGD w PDAC"="#FE8116",
  #                                "PDAC"="#FE5A1D"))+
  scale_colour_manual(values = color) +
  theme_bw () +
   labs(title= "ssGSEA: UMAP", x ="UMAP 1", y="UMAP 2") ) 



### boxplot of the ssgsea
hallmark_subcat <-read_excel("data/BE02_case_list.xlsx", sheet = "hallmark")

t(ssgsea_hallmark) %>% as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(!Sample_ID ,names_to = "pathway", values_to = 'score') %>%
  left_join(condition) %>%
  mutate(Disease = factor(Disease, levels = c("Normal-LGD", "Normal-HGD",
                                             "Normal-PDAC", "Isolated LGD",
                                             "Progressor LGD w HGD",
                                             "Isolated HGD", 
                                             "Progressor LGD w PDAC",
                                             "Progressor HGD w PDAC",
                                             "PDAC") ),
          pathway = str_remove(pathway, "HALLMARK_")) %>%
  left_join(hallmark_subcat %>% 
  dplyr::select(Process_category, Description), by= c(pathway = "Process_category")) %>%
  filter(Description %in% unique(Description)[7:8]) %>%
  ggplot(aes(x = score, y = pathway, colour = Disease)) +
        geom_boxplot(position = position_dodge()) +
        facet_wrap(~Description, scales = "free_y",nrow = 1)



### creating the line plot including the normal
df_gsea <- t(ssgsea_hallmark) %>% as.data.frame() %>% 
            
            rownames_to_column('Sample_ID') %>%
          pivot_longer(!Sample_ID, values_to = "score", names_to = 'pathways') %>%
           left_join(condition %>% dplyr::select(Sample_ID, Category, Disease, Grade, Path_ID))
           
summary(df_gsea$score)
# for (path in unique(df_gsea$pathways)[11:15]){
p <- df_gsea %>% #filter(Category !="PDAC") %>%
    # filter(pathways == path) %>%
  mutate(pathways = str_replace(str_remove(pathways, "HALLMARK_"), "_"," "),
     Grade = factor(Grade , levels = c("normal", "low grade", "high grade", "invasive")),
         Category = factor(Category, levels = c("LGD", "HGD", "PDAC"))) %>%
 ggplot(aes(x = Grade, y = score,  group = Path_ID,color =  Category)) +
  geom_line() +
  # labs(title = path) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~pathways, nrow = 8, scales= "free_y", labeller = label_wrap_gen(width = 15, multi_line = TRUE))
print(p)

# }
