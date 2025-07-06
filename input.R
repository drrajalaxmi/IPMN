## input data count matrix and the sample metadata

library(dplyr)
library(readxl)
### metadata
BE02_case_list <-read_excel("data/BE02_case_list.xlsx", sheet = "FINAL")
column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE02_",BE02_case_list$Sample_ID, "A")


## count matrix
df_count_tpm <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)   ### tpm file
df_count <-read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)   ## raw count
colnames(df_count) <- paste0(stringr::str_replace(colnames(df_count), "X", "BE02_"), "A")
colnames(df_count_tpm) <- paste0(stringr::str_replace(colnames(df_count_tpm), "X", "BE02_"), "A")

condition <- BE02_case_list %>% 
            filter(!Sample_ID %in% c("BE02_65A" , "BE02_77A", "BE02_78A",  "BE02_80A",  "BE02_81A",  "BE02_120A")) %>%  ## samples excluded
            filter(Sample_ID %in% colnames(df_count))
counts <- df_count  %>% dplyr::select(condition$Sample_ID)
counts_tpm <- df_count_tpm  %>% dplyr::select(condition$Sample_ID)

length(unique(condition$Sample_ID))
dim(counts)


### CIBERSORTx output dataset 
cibersort <- read.csv("data/CIBERSORTx_Job4_Results.csv") %>%
                   mutate(Mixture = paste0(str_replace(Mixture, "BE", "BE02_"), "A"))

cibersort_pt <- cibersort %>%  
                    filter(Mixture %in% condition$Sample_ID) %>% 
                    left_join(condition %>% dplyr::select(Sample_ID, Grade, Path_ID, Path_diagnosis,Main_duct, Category, Disease), 
                              by = c("Mixture" = "Sample_ID")) %>%
                    mutate(Grade = factor(Grade, levels = c("normal", "low grade", "high grade", "invasive") ),
                    Disease = factor(Disease, levels = c("Normal-LGD" ,"Normal-HGD", "Normal-PDAC","Isolated LGD",
                                                        "Progressor LGD w HGD" , 
                                                        "Isolated HGD" ,"Progressor LGD w PDAC",
                                                        "Progressor HGD w PDAC", "PDAC")) )


dim(cibersort_pt)
