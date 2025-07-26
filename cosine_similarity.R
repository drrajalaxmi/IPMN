## cosine similarity

lm22 <- read.table("data/LM22.txt", sep = "\t", header = TRUE)



immune_gene_tpm <- df_tpm_ms[rownames(df_tpm_ms) %in% lm22$Gene.symbol,]
subset_lm22 <- lm22[lm22$Gene.symbol %in% rownames(immune_gene_tpm), ]


subset_tpm <- immune_gene_tpm[1:2]

cell = colnames(subset_lm22)[2]


cos_res_mtx <- data.frame(matrix(ncol = 23, nrow = ncol(immune_gene_tpm)))
# Assign column names
colnames(cos_res_mtx) <- c("Sample_ID", colnames(subset_lm22)[-1])
cos_res_mtx$Sample_ID <- colnames(immune_gene_tpm)

for (cell in colnames(subset_lm22)[-1]) {
  cos_res <- c()
for(col in colnames(immune_gene_tpm)){
     cos_vec <- cosine(immune_gene_tpm[[col]], subset_lm22[[cell]])
     cos_res <- c(cos_res, cos_vec)
}
  names(cos_res) <- cell
  cos_res_mtx[[cell]] <- cos_res
  }

rownames(cos_res_mtx) <- cos_res_mtx[1]

heatmaply::heatmaply(t(cos_res_mtx[-1]),cluster_rows= FALSE, cluster_cols= FALSE )









