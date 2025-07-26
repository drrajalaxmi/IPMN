### untrimmed count
library(readxl)
library(dplyr)
library(ggplot2)
library(ggfortify)
BE02_case_list <- read_excel("data/BE02 case list.xlsx")

column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE",BE02_case_list$Sample_ID)




df <- read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
# df <- read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)
# df <- read.csv("data/counts_fromBam_untrimmed.csv", header=TRUE, row.names = 1)
# df <- read.csv("data/counts_salmon_fq_untrimmed.csv", header=TRUE, row.names = 1)
str(df)
colnames(df) <- stringr::str_replace(colnames(df), "X", "BE")
BE02_subset <- BE02_case_list %>% filter(Sample_ID %in% colnames(df))
length(BE02_subset$Sample_ID)
counts <- df  %>% select(BE02_subset$Sample_ID)

BE02_subset_test <- BE02_subset  %>% 
                    # filter(!BE02_subset$Sample_ID %in% c("BE120", "BE65")) %>%
                    # filter(!BE02_subset$Notes %in% c("FAIL")) %>%
                    filter(!BE02_subset$Sample_ID %in% outliers$Sample_ID)
colnames(BE02_subset_test)
counts <- df %>% select(BE02_subset_test$Sample_ID) 

dim(counts)
cnt_pca <- prcomp(t(counts))
summary(cnt_pca)


# par(mar=c(4,4,1,1))
# 
# dev.off()
# # Check if it works:
# 
# plot(rnorm(50), rnorm(50))
# plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (all groups)") 
# pca.plot <- autoplot(cnt_pca,
#                      data = t(counts)) 
autoplot(cnt_pca, data = BE02_subset_test,color = "Grade", shape = FALSE, label.size = 3)
autoplot(cnt_pca, data = BE02_subset_test,color = "Grade")
# str(counts)
# unique(BE02_case_list$Grade[BE02_case_list$Color=="green"])


count_mtx <- as.matrix(counts)

Ca_cor <-round(cor(count_mtx), 2)
heatmaply::heatmaply((Ca_cor))
df_val <- df_mtx[!rowSums(df_mtx)==0, ]
dim(df_val)

heatmap(Ca_cor)
# heatmap(Ca, symm = TRUE, margins = c(6,6))

# -------------------------------------------------------------------------

colnames(BE02_case_list)
BE02_case_list_pct <- BE02_case_list %>% left_join(map_df)
write.csv(BE02_case_list_pct, "data/BE02_case_list_pct_softclip")
outliers <- BE02_case_list_pct %>% select(Sample_ID, pcnt_mapped, Grade,Notes, Sample_Quality1) %>% filter(pcnt_mapped <20 & !is.na(pcnt_mapped)) %>% 
  arrange(pcnt_mapped) 
outliers$Sample_ID
