
### IPMN cancer moonshot project data

MS_case_list <-read_excel("data/cancer_moonshot_case_list_final.xlsx")

condition_ms <- MS_case_list %>%
          mutate( 
            Sample_ID = paste0("MS", Sample_ID),
            Grade = factor(Grade, levels = c("LGD", "HGD", "PDAC") ),
            Disease = factor(Disease, levels = c("normal-LGD-acinar" ,"normal-LGD-duct" ,
                                                 "normal-HGD-acinar" ,"normal-HGD-duct" ,
                                                 "normal-PDAC-acinar" ,"normal-PDAC-duct" ,
                                                 "isolated LGD",
                                                 "progressor LGD w HGD" , 
                                                 "isolated HGD" ,"progressor LGD w PDAC",
                                                 "progressor HGD w PDAC", "PDAC")),
            Patient_id = factor(Path_ID, levels = unique(Path_ID))) %>%
            arrange(Patient_id , Grade, Disease)


## counts matrix
df_count_ms <-read.csv("data/pancreas_counts_min10counts_with_IDs.csv", header=TRUE, row.names = 1)
colnames(df_count_ms) <- stringr::str_replace(colnames(df_count_ms), "X", "MS")
colnames(df_count_ms) <- stringr::str_replace(colnames(df_count_ms), "\\.", "_")
df_tpm_ms <- IOBR::count2tpm(df_count_ms%>% as.matrix(), idType = "symbol")



## CIBERSORTx output data
cibersort_ms <- read.csv("data/CIBERSORTx_Job28_tpm.csv")
cibersort_pt_ms <- cibersort_ms %>%  
  left_join(condition_ms  %>% select(-Path_ID), by = c("Mixture" = "Sample_ID")) %>%
  mutate(Category = Grade) 

