library(readxl)
library(ggfortify) 
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(apeglm)
library(biomaRt)
library(clusterProfiler)
library(stringr)
library(tibble)
library(fgsea)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(pathview)
# devtools::install_github("zhangyuqing/sva-devel")


log(1.2,2)





BE02_case_list <-read_excel("data/BE02_case_list.xlsx", sheet = "FINAL")
column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE02_",BE02_case_list$Sample_ID, "A")


# BE_count <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
BE_count <-read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)
colnames(BE_count) <- paste0(stringr::str_replace(colnames(BE_count), "X", "BE02_"), "A")
BE_condition <- BE02_case_list %>% #filter(!Sample_ID %in% c("BE02_65A" , "BE02_77A", "BE02_78A",  "BE02_80A",  "BE02_81A",  "BE02_120A")) %>%
 # filter(Grade !="normal") %>%
  filter(Sample_ID %in% colnames(BE_count))




length(BE_condition$Sample_ID)
BE_counts <- BE_count  %>% dplyr::select(BE_condition$Sample_ID)
length(colnames(BE_counts))

rownames(BE_counts) = unlist(lapply(strsplit(rownames(BE_counts), '\\.'), '[[',1))
ids<-bitr(rownames(BE_counts), fromType="ENSEMBL",
          toType=c("ENTREZID", "SYMBOL"),
          OrgDb="org.Hs.eg.db") #NA values drop by default

BE_mtx <- BE_counts %>% rownames_to_column("ENSEMBL") %>% left_join(ids%>% dplyr::select(-ENTREZID)) 
BE_mtx <- BE_mtx[!duplicated(BE_mtx$SYMBOL),]  
BE_mtx <- BE_mtx[!is.na(BE_mtx$SYMBOL),]  
BE_mtx <- BE_mtx %>% select(-c(ENSEMBL )) %>% relocate(SYMBOL)



##### Cancer Moonshot

MS_case_list <-read_excel("data/cancer_moonshot_case_list_final.xlsx")
str(MS_case_list)
sort(unique(MS_case_list$Sample_ID))

column_name <- stringr::str_replace(colnames(MS_case_list), " ", "_")
# column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
# colnames(MS_case_list) <- column_name
MS_case_list$Sample_ID <- paste0("MS",MS_case_list$Sample_ID)
str(MS_case_list)




MS_count <-read.csv("data/pancreas_counts_min10counts_with_IDs.csv", header=TRUE, row.names = 1)
colnames(MS_count) <- stringr::str_replace(colnames(MS_count), "X", "MS")
colnames(MS_count) <- stringr::str_replace(colnames(MS_count), "\\.", "_")
str(MS_count)

MS_cols <- colnames(MS_count)
MS_case_list <- MS_case_list %>% filter(Sample_ID %in% MS_cols)

MS_tpm <- count2tpm(MS_count, idType = "symbol")

MS_condition <- MS_case_list %>% 
 # filter(!str_detect(Disease, "normal")) %>%
  filter(Sample_ID %in% colnames(MS_count)) %>%
  # filter(Sample_ID %in% colnames(MS_tpm))
  mutate(Disease = ifelse(str_detect(Disease, "^p"), str_replace(Disease, "^p", "P"),
                   ifelse(str_detect(Disease, "^i"), str_replace(Disease,"^i", "I"), 
                          str_replace(Disease,"^n", "N")  )) )



length(MS_condition$Sample_ID)
MS_mtx <- MS_count  %>% dplyr::select(MS_condition$Sample_ID) %>% 
       rownames_to_column("SYMBOL")


length(colnames(MS_mtx))




##### merging the tpms

counts <- MS_mtx %>% inner_join(BE_mtx)  %>% 
        column_to_rownames("SYMBOL")


MS_mtx$SYMBOL[!MS_mtx$SYMBOL %in% BE_mtx$SYMBOL]



dim(MS_mtx)
dim(counts)
dim(BE_mtx)

any(is.na(counts))

##### merging the samples

str(MS_condition)
str(BE_condition)
unique(BE_condition$Category)
condition_all <- MS_condition %>% rename(Grade = "Category") %>%
        mutate(cohort = "MS") %>%
      rbind(BE_condition %>% select(Path_ID, Sample_ID, Disease, Category ) %>%
              mutate(cohort = "BE") )

condition <- condition_all %>%
         filter(!str_detect(Disease, "Normal"))

dim(counts)
dim(condition_all)
counts <- counts %>% select(condition_all$Sample_ID)


###### correcting for batch effect

adjusted_counts <- sva::ComBat_seq(counts, batch=condition_all$cohort, group=NULL)


adjusted_tpm <- count2tpm(adjusted_counts, idType = "symbol")



condition$Disease <- factor(condition$Disease, levels = c("Isolated LGD" ,
                                                          "Isolated HGD" ,
                                                          "Progressor LGD w HGD" ,
                                                          "Progressor LGD w PDAC" ,
                                                          "Progressor HGD w PDAC",
                                                          "PDAC" ) )

cnt_pca <- prcomp(t(log(adjusted_tpm+1)))

counts[is.na(counts)]
df_txt  <- df_txt[!is.na(df_txt$gene_symbol),]

unique(condition$Disease)
colours <- c("Isolated LGD" = "#2F56B5",#115C80", 
             "Isolated HGD" ="#2295C7", #309898" ,#479685", 
             "Progressor LGD w HGD" = "#00A86B", #FF9F00", 
             "Progressor LGD w PDAC" = "#FCC92A", #FF9F00", 
             "Progressor HGD w PDAC" = "#F4691E", 
             "PDAC" = "#AD0000")

autoplot(cnt_pca, data = condition_all,color = "Disease", size = 2) 
  theme_minimal()
plotly::ggplotly(autoplot(cnt_pca, data = condition_all,color = "Disease", size = 2.5, shape = "cohort" ))
                 +
                   scale_color_manual(values = paste0(colours)) )
theme_classic())


