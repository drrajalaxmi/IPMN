### mixcer


tra <- read_tsv("data/BE_final_postanalysis.diversity.TRA.tsv")
igH <- read_tsv("data/BE_final_postanalysis.IsotypeUsage.IGH.tsv")

umap_df <- igH %>% separate(sample, sep = "A_", c("sample_ID", NA), remove = TRUE ) %>%
       mutate(sample_ID = paste0(sample_ID, "A"))# %>% column_to_rownames("sample_ID")

umap_tra <- tra %>% separate(sample, sep = "A_", c("sample_ID", NA), remove = TRUE ) %>%
  mutate(sample_ID = paste0(sample_ID, "A"))# %>% column_to_rownames("sample_ID")

umap_df <- umap_df %>% left_join(umap_tra, by = "sample_ID" ) %>% column_to_rownames("sample_ID")
  
umap_df[umap_df=="null"] <- 0
umap_df <- umap_df %>% type.convert(., as.is = TRUE)

str(umap_df)
umap_res <- umap((umap_df)  )
umap_df <- data.frame(umap_res$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("Sample_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::left_join(condition, by = "Sample_ID") #%>%
# dplyr::left_join(cibersort, by = c("Sample_ID" = "Mixture"))
# str(umap_df)

plotly::ggplotly(ggplot(umap_df,aes(
  x = X1,
  y = X2,
  label =Path_ID,
  color= Category, stroke = 0.2,
  #color= log(Neutrophils)
) ) +
  geom_point(size = 3) )
