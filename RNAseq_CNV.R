devtools::install_github(repo = "honzee/RNAseqCNV")
library(RNAseqCNV)
launchApp()

write.table(df_count %>% rownames_to_column("ENSEMBL") %>% 
              mutate(sample = round(BE02_2A),
                     ENSEMBL = unlist(lapply(strsplit(ENSEMBL, '\\.'), '[[',1))) %>%
              select(ENSEMBL,sample) , "data/count_BE02.txt",sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


RNAseqCNV_wrapper(config = "data/config_cnv.txt", metadata = "data/metadata_cnv.txt", snv_format = "vcf") 




cnv <-read.table("data/alteration_matrix.tsv", sep = "\t", header =TRUE)

colnames(cnv) <- str_replace(colnames(cnv), "X", "chr")

cols_list <- c()
for (i  in colnames(cnv)[-1]){
   if (sum(cnv[[i]], na.rm = TRUE)!=0){
     cols_list <- c(cols_list, i)
   }
}


df_cnv <- cnv %>% dplyr::select("sample",cols_list) %>%left_join(condition %>% dplyr::select(Sample_ID, Disease, Path_ID), by = c("sample" = "Sample_ID")) %>% arrange(Path_ID) %>%
  relocate(Path_ID, Disease)

df_cnv %>% filter(str_detect(Disease, "LGD"))


write.xlsx(df_cnv, "data/cnv_matx1.xlsx")





#### using z score to plot the heatmap
library(biomaRt)


  
zScore <- scale(t(as.matrix(counts_tpm)), center = TRUE, scale = TRUE)

boxplot(zScore[1,])
x <- matrix(1:10, ncol = 2)
(scale(x,center = TRUE, scale = FALSE))
cov(centered.scaled.x <- scale(x))

rownames(counts_tpm) = unlist(lapply(strsplit(rownames(counts_tpm), '\\.'), '[[',1))

# Retrieve gene information
gene_info <- getBM(
  attributes = c("chromosome_name", "start_position",#c "end_position", "external_gene_name", 
                 "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = rownames(counts_tpm),
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
)

# Print the retrieved information
counts_chrPos <- t(zScore) %>% as.data.frame() %>%
  rownames_to_column(var ="ENSEMBL") %>% 
  # mutate(ENSEMBL = unlist(lapply(strsplit( ENSEMBL, '\\.'), '[[',1)) ) %>%
  left_join(gene_info %>% filter(
    str_detect(chromosome_name, pattern = "^[1-9]|^[X-Y]") ), by =c("ENSEMBL"= "ensembl_gene_id")) %>%
  filter(!is.na(chromosome_name))  %>%  arrange(chromosome_name,start_position)
counts_chrPos <- counts_chrPos[!(duplicated(counts_chrPos$ENSEMBL)),]  
rownames(counts_chrPos) <- paste0(counts_chrPos$chromosome_name,counts_chrPos$ENSEMBL)
counts_chrPos <- counts_chrPos %>% dplyr::select(-c("ENSEMBL","chromosome_name", "start_position"))

head(counts_chrPos[1:10,])
stats::heatmap(as.matrix(counts_chrPos[1:78639,]),  Rowv = NA)
heatmaply::heatmaply(as.matrix(counts_chrPos)[1:100])
