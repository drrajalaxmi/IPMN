## differential gene expression analysis and gsea analysis

## input data condition (sample metadata) and counts (raw count data) from input.R
## condition has 5 samples excluded: "BE02_65A"  "BE02_77A"  "BE02_78A"  "BE02_80A"  "BE02_81A"  "BE02_120A"

## outputs
# 1. Volcano plot
# 2. Pathway enrichment graphs
# 3. Table output for DE genes and pathways
# 4. Combined graphs for multiple comparison


library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)



## fetch gene symbols
ensemble_list = unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))
ids<-bitr(ensemble_list, fromType="ENSEMBL",
          toType=c("ENTREZID", "SYMBOL"),
          OrgDb="org.Hs.eg.db") #NA values drop by default



# plot_name <- "Isolated HGD vs Isolated LGD "
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Isolated HGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Isolated LGD", "Isolated HGD")))
# sheet_name <- "IsoHGD_Vs_IsoLGD"

# plot_name <- "Progressor HGD w PDAC vs\n Isolated HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated HGD", "Progressor HGD w PDAC")) %>%
#   mutate( Disease = factor(Disease, levels= c("Isolated HGD", "Progressor HGD w PDAC")))
# sheet_name <- "ProHGDwPDAC_Vs_IsoHGD"
# 
# 
# plot_name <- "Progressor LGD all (combined HGD or PDAC) vs Isolated LGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC")) %>%
#   mutate( Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#           Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD")))
# sheet_name <- "ProLGD(All)_Vs_IsoLGD"
# 
# 
# plot_name <- "Progressor LGD w HGD vs Isolated LGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD")) %>%
#   mutate( #Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#     Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w HGD")))
# sheet_name <- "ProLGDwHGD_Vs_IsoLGD"
# 
# 
# plot_name <- "Progressor LGD w PDAC vs Isolated LGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#     Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w PDAC")))
# sheet_name <- "ProLGDwPDAC_Vs_IsoLGD"

list <- c(unique(condition$Disease)[-9], "Isolated HGD", "Progressor LGD w HGD")
list_dysplasia <- unique(condition$Disease)[unique(condition$Disease)!= "Normal-LGD"]

list_dysplasia <-c( "Progressor LGD w HGD" ,"Isolated HGD")
list_dysplasia <-c("Progressor LGD w PDAC", "Progressor HGD w PDAC", "PDAC")
list_dysplasia <-c( "Isolated HGD")
list_dysplasia <-c( "Progressor HGD w PDAC", "PDAC")
list_dysplasia <-c(  "PDAC")
list_dysplasia <-c( "Progressor LGD w HGD", "Isolated HGD", "Progressor LGD w PDAC","Progressor HGD w PDAC", "PDAC" )

list_dysplasia <-c(  "Progressor LGD w PDAC","Progressor HGD w PDAC", "PDAC" )

counts_mtx <- counts %>% dplyr::select(which(colnames(counts) %in% sample_mtx$Sample_ID)) %>% as.matrix()


for( dysplasia in list_dysplasia ){}
  plot_name <- paste0(dysplasia, " vs Progressor LGD w PDAC")
  sample_mtx <- condition %>% filter(Disease %in% c("Progressor LGD w PDAC", dysplasia)) %>%
    mutate( #Disease = ifelse(Disease =="Progressor LGD w HGD", "Normal-PDAC","Normal-Non-Progressor" ) ,
      Disease = factor(Disease, levels= c("Progressor LGD w PDAC", dysplasia) ))
  sheet_name <- paste0(i, "_Vs_isoLGD")
  
  
  plot_name <- paste0(dysplasia, " vs Isolated HGD")
  sample_mtx <- condition %>% filter(Disease %in% c("Isolated HGD")) %>%
    mutate( #Disease = ifelse(Disease =="Progressor LGD w HGD", "Normal-PDAC","Normal-Non-Progressor" ) ,
      Disease = ifelse(Path_ID %in% c("1-S-22-18334", "SNBU22-14478", "SNBU22-17893", "1-S-23-07849"), "higer grade", "lower grade"),
      Disease = factor(Disease, levels= c( "lower grade", "higer grade") ))
  
  counts_mtx <- counts %>% dplyr::select(which(colnames(counts) %in% sample_mtx$Sample_ID)) %>% as.matrix()
  
levels(sample_mtx$Disease)
dim(sample_mtx)
dim(counts_mtx)


## diffential expression DEseq2

counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts_mtx), colData = sample_mtx, design = ~ Disease)
dds <- DESeq(counts.DEseq)
comparison <- resultsNames(dds)[length(resultsNames(dds))] #lists the coefficients
res <- results(dds)


####
# ## adding gene symbols
res2 <- res %>% as.data.frame() %>%
  rownames_to_column(var ="ENSEMBLID") %>%
  mutate(ENSEMBL = unlist(lapply(strsplit(ENSEMBLID, '\\.'), '[[',1)),
         comparison = plot_name) %>% relocate(comparison) %>%
  left_join(ids)   %>%
  filter(!is.na(ENTREZID) ) %>%
  filter(!is.na(padj) ) %>%     
  filter(!duplicated(SYMBOL)) 
  # filter( !str_detect(SYMBOL, "LOC|RNA" )) %>%
  #  filter(!is.na(ENTREZID) ) %>%
  # arrange(desc(abs(log2FoldChange)))


write.csv(res2 %>% column_to_rownames("SYMBOL") %>%
            mutate(logFC = log2FoldChange,
                   pval = padj,
                   expr = baseMean) , 
           "data/DE_LGDwPDACVsHGDwPDAC.csv")


# 
# # cutoffs
# lfc = 2        # log fold change
# pval = 0.05    # adjusted p value
# 
# 
# # make a data frame with the log2 fold-changes, adjusted p-values and symbol
tab_res = data.frame(logFC = res2$log2FoldChange, negLogPval = -log10(res2$padj),
                     symbol = res2$SYMBOL)
signGenes_res = (abs(tab_res$logFC) > lfc & tab_res$negLogPval > -log10(pval))
tab_res$delabel[signGenes_res] <- tab_res$symbol[signGenes_res]

##volcano plot
ggplot(tab_res, aes(x = logFC, y = negLogPval,  label = delabel )) +
  geom_point(size = 0.8) +
  geom_point(data= tab_res[signGenes_res,], size = 0.8, color = "red") +
  geom_text_repel( size = 3, max.overlaps=11, color = "black") + # Adjust box.padding and max.overlaps as needed
  geom_vline(xintercept=c(-lfc, lfc), linewidth = 0.5,col="blue",linetype= 2) +
  geom_hline(yintercept=-log10(pval), col = "green3", linetype= 2) +
  labs(title = plot_name,
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme_bw() +
  theme(panel.grid = element_blank())


### pathway enrichment
res3 <- res %>% as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL") %>% 
  mutate(ENSEMBL = unlist(lapply(strsplit(ENSEMBL, '\\.'), '[[',1)) ) %>%
  left_join(ids) %>%
  filter(!is.na(ENTREZID), !is.na(stat)) %>%
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat)) 

ranks <- res3$stat
names(ranks) <- res3$SYMBOL

# pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/h.all.v2024.1.Hs.symbols.gmt") # for hallmark pathways
# pathways.reactome <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/c2.cp.reactome.v2024.1.Hs.symbols.gmt") # for reactome pathways


fgseaRes_hallmark <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)
fgseaRes_reactome <- fgseaMultilevel(pathways=pathways.reactome, stats=ranks)

# Tidy the results:
fgseaResTidy_hallmark <- fgseaRes_hallmark %>% 
                            as_tibble() %>% 
                            arrange(desc(NES)) # order by normalized enrichment score (NES)
fgseaResTidy_reactome <- fgseaRes_reactome %>%
                            as_tibble() %>% 
                            arrange(desc(NES)) # order by normalized enrichment score (NES)


### hallmark pathway plot
fgseaRes_hallmark %>%
          filter(padj <= 0.05 )  %>%
          mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
          ggplot(aes(reorder(pathway, NES), NES, fill = "red")) +
          geom_col(show.legend = FALSE) +
          scale_fill_manual(values = "red") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          coord_flip() +
          labs(x="Hallmark Pathways", y="Normalized Enrichment Score",
               title=paste0(plot_name) )


### table output of pathways
fgseaResTidy_table_hallmark <- fgseaResTidy_hallmark %>%
                                  group_by(pathway) %>%
                                  mutate(leadingEdge = paste0(unlist(leadingEdge), collapse = ","),
                                         comparison = plot_name)
fgseaResTidy_table_reactome <- fgseaResTidy_reactome %>%
                                  group_by(pathway) %>%
                                  mutate(leadingEdge = paste0(unlist(leadingEdge), collapse = ","),
                                         comparison = plot_name)


# wb_hallmark <- createWorkbook("BE02_hallmark")
# wb_reactome <- createWorkbook("BE02_reactome")
## Add a worksheets
addWorksheet(wb_hallmark, sheet_name, gridLines = FALSE)
writeData(wb_hallmark, sheet = sheet_name, fgseaResTidy_table_hallmark, rowNames = TRUE)
# saveWorkbook(wb_hallmark, "data/BE02_hallmark_norm.xlsx", overwrite = TRUE)

addWorksheet(wb_reactome, sheet_name, gridLines = FALSE)
writeData(wb_reactome, sheet = sheet_name, fgseaResTidy_table_reactome, rowNames = TRUE)
# saveWorkbook(wb_reactome, "data/BE02_reactome_norm.xlsx", overwrite = TRUE)

### Table output DE genes
# wb <- createWorkbook("BE02_sig_DE_gene")  ## run only for the first time
addWorksheet(wb, sheet_name, gridLines = FALSE)  ## Add a worksheets
writeData(wb, sheet = sheet_name, res2, rowNames = TRUE) ## write data to worksheet
# saveWorkbook(wb, "data/BE02_sig_DE_gene_norm.xlsx", overwrite = TRUE)  ## save wb


### Multiple comparison sets for dot plot

## create list of items and iniiate the list
# hallmark_subcat <-read_excel("data/BE02_case_list.xlsx", sheet = "hallmark") ## file categorizing pathway to functions
# names_list <- c()
# mlist <- list()
# i <- 1

names_list <- c(names_list, plot_name)
mlist[[i]]<-fgseaRes_hallmark
i <- i+1


}

# names_list <- c("Progressor\nLGD w PDAC vs\n Isolated LGD"  , "Progressor\nHGD w PDAC vs\n Isolated HGD")

names(mlist)<- names_list
mresult<-merge_result(mlist)

cluster_list <- names_list

cluster_list <- c(
  "Normal-HGD vs Normal-LGD" ,
  "Isolated LGD vs Normal-LGD",
  
  "Progressor LGD w HGD vs Normal-LGD" ,
  "Isolated HGD vs Normal-LGD",
  "Normal-PDAC vs Normal-LGD" ,
  "Progressor LGD w PDAC vs Normal-LGD",
  "Progressor HGD w PDAC vs Normal-LGD",
  "PDAC vs Normal-LGD"
  
  # 
  # "Progressor\nLGD w PDAC vs\n Isolated LGD",
  # "Progressor\nHGD w PDAC vs\n Isolated HGD"
)

write.xlsx(mresult %>% as.data.frame(), "data/mresult.xlsx")

mresult %>% as.data.frame() %>% 
         mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
         left_join(hallmark_subcat %>% 
         dplyr::select(Process_category, Description), by= c(pathway = "Process_category")) %>%
         filter(Cluster %in% cluster_list) %>% 
         # separate(Cluster, sep = " vs", c("Disease", NA), remove = FALSE) %>%
         filter(padj <= 0.05) %>%
         # filter(abs(NES) >= 1) %>% 
         arrange( NES, Description )%>% 
         mutate(
           Cluster = factor(Cluster, levels = cluster_list),
           Description = str_to_title(Description),
           Pathway = factor(pathway, levels = unique(pathway)),
           Pathways = ifelse(NES < 0,"Down", "Up") ,
           Pathways = factor(Pathways, levels= c("Up","Down" ) ),
           NES = abs(NES))  %>% 
        ggplot(aes(Cluster,Pathway , size = NES, color = Pathways, group = Description)) +
        geom_point() + 
        labs(x ="", y = "") +
        scale_color_manual(values = c("red", "#2F56B5")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, size= 8,vjust = 1, hjust=1)) +
        facet_wrap(~Description, scales = "free_y",nrow = 4)





df_pathway <- mresult %>% as.data.frame() %>% 
  mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
  left_join(hallmark_subcat %>% 
              dplyr::select(Process_category, Description), by= c(pathway = "Process_category")) %>%
  filter(Cluster %in% cluster_list) %>% 
  separate(Cluster, sep = " vs ", c("Disease", "Base"), remove = FALSE) %>%
  filter(padj <= 0.05) %>%
  dplyr::select(Cluster,Disease, pathway, Description, Base, NES) #%>%
  # pivot_wider( names_from = Base, values_from = NES)


color <- c("purple", "grey80", 
           "#2F56B5",#115C80", 
           "#2295C7", #309898" ,#479685", 
           "#00A86B", #FF9F00", 
           "#FCC92A", #FF9F00", 
           "#F4691E", 
           "#AD0000")


df <- df_pathway %>%
  mutate(Disease = factor(Disease, levels = c("Normal-HGD" , "Normal-PDAC" ,"Isolated LGD" ,
                                              "Progressor LGD w HGD",
                                              "Isolated HGD", "Progressor LGD w PDAC", "Progressor HGD w PDAC",
                                              "PDAC" ))) 
  plotly::ggplotly(ggplot(df,aes(x=NES,y = pathway, color = Disease, shape = Base)) +
  geom_point() +
  scale_color_manual(values = color) +
  facet_grid(Description~Disease+Base, scales = "free_y"))
  

  df <- df %>% 
        #filter(Cluster %in% cluster_list[c(1,13)]) %>% 
        arrange(Base)  %>%
        mutate(pathway = factor(pathway, levels = unique(pathway)))
        
  plotly::ggplotly(ggplot(df,aes(x=Disease,y = pathway, color = Disease, shape = Base)) +
                     geom_point() +
  
                     facet_wrap(~Base, scales = "free_y",nrow = 4))




  
### enrichment plot
  fgseaResTidy_reactome$pathway[str_detect(fgseaResTidy_reactome$pathway, "MHC")]
  plotEnrichment(pathways.reactome[["REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION"]],
                 ranks) + labs(title="Class I MHC mediated antigen processing presentation",
                               x = plot_name, y = "Enrichment Score")
  
