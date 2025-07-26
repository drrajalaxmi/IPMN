
library(readxl)

sample <- read_excel("data/GSE278670_Iyer_WTA_GeoMx_2022_GEO_Initial_Dataset.xlsx", sheet = 1)
df <- read_excel("data/GSE278670_Iyer_WTA_GeoMx_2022_GEO_Initial_Dataset.xlsx", sheet = 2)
str(sample)
colnames(df)[c(4,13: ncol(df))]

df_counts <- df %>% select(c(4,13: ncol(df)))
head(sample$ROILabel)
sample_mtx <- sample %>% filter(PanCK == "True" & NMLD =="False") %>%
      mutate(Disease = paste0(tags, "_", class),
             Sample_ID = SegmentDisplayName) %>% select(2,5, Disease, Sample_ID, AOISurfaceArea, AOINucleiCount) 


count_mtx <- df_counts %>% select(1, sample_mtx$SegmentDisplayName) %>% 
  filter(!duplicated(HUGOSymbol), !is.na(HUGOSymbol) ) %>%
  column_to_rownames("HUGOSymbol")
stats::heatmap(as.matrix(count_mtx))




unique(sample_mtx$Disease)
color_df <- data.frame(Disease= c("PanCK_LGD_LGD" ,
                                  "PanCK_LGD_INV" ,
                                   "PanCK_HGD_INV",
                                  "PanCK_INV_INV"),
                       
                       colors = c("#2F56B5",  "#FCC92A", "#F4691E", "#AD0000")
) 
# color_df <- data.frame(Category= c("LGD", "HGD",
#                                    "PDAC" ),
#                        colors = c("#2F56B5", "#00A86B",  "#AD0000")
# ) 
df_color <- data.frame(Sample_ID = colnames(count_mtx)) %>% left_join(sample_mtx %>% dplyr::select(Sample_ID, Disease) %>% 
                                                                     left_join(color_df) %>% filter(!is.na(colors)) ) %>%
  mutate(colors = ifelse(is.na(colors), "black", colors))


# df_color <- data.frame(Sample_ID = rownames(Ca_cor)) %>% left_join(condition %>% dplyr::select(Sample_ID, Category) %>% 
#                                                                      left_join(color_df) %>% filter(!is.na(colors)) ) %>%
#   mutate(colors = ifelse(is.na(colors), "black", colors))
# heatmaply::heatmaply(Ca_cor,
# heatmap_layers = theme(axis.text.y = element_text(color = c("red",rep(c("red", "blue"), (nrow(Ca_cor)-1)/2)))) #, margins = c(5, 5), scale = "none")

#create data frame for annotations

dfh<-data.frame(Sample_ID=as.character(colnames(count_mtx)))%>%
  left_join(sample_mtx %>% select(Sample_ID,  Disease)) %>%
  column_to_rownames("Sample_ID")

dfh
library(pheatmap)
library(ComplexHeatmap)

count_colsum <- colSums(count_mtx)
norm_colsum <- sweep(count_mtx, MARGIN = 2, STATS = count_colsum, FUN = "/")

colnames(sample)
df_area <- data.frame(Sample_ID = colnames(count_mtx)) %>% left_join(sample_mtx %>% dplyr::select(Sample_ID, Disease, AOISurfaceArea, AOINucleiCount) 

norm_area <- sweep(count_mtx, MARGIN = 2, STATS = df_area$AOISurfaceArea, FUN = "/")
norm_nuclei <- sweep(count_mtx, MARGIN = 2, STATS = df_area$AOINucleiCount, FUN = "/")

stats::heatmap(as.matrix(norm_area),ColSideColors = df_color$colors, labCol = NA)


pheatmap::pheatmap(as.matrix(count_mtx),# cutree_rows = 4, 
                   color = colorRampPalette(
                     c( "yellow", "red"))(10),
                   annotation_row = dfh, border_color = NA,drop_levels = TRUE,
                   annotation_colors = ann_colors, angle_col = "45", show_colnames = FALSE,
                   fontsize_row = 7) 

