# To calculate the coorelation between the gene expression and grade of the tumour
# 1/10/2025
# Rajalamxi


gene_list <- read.csv("data/signif_genelist_lcf.csv")
gene_list$corr_v <- c()
condition <- BE02_case_list %>% #filter(!Sample_ID %in% c("BE65" , "BE77", "BE78",  "BE80",  "BE81",  "BE120")) %>%
  filter(Sample_ID %in% colnames(df_count))
sample_grade <- condition %>% mutate(Grade = as.numeric(factor(Grade, levels = c("normal", "low grade", "high grade", "invasive") ))) %>%
                     arrange(Grade) %>%
                       select(Sample_ID, Grade)
# sample_grade %>% filter(Grade==1) %>% tail()
df_count <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
# df_count <-read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)
colnames(df_count) <- stringr::str_replace(colnames(df_count), "X", "BE")

counts_sign_gene <- df_count %>% dplyr::select(sample_grade$Sample_ID) 
str(counts_sign_gene)
dim(sample_grade)
gene_corr <- data.frame(geneID = rownames(counts_sign_gene)) 
gene_corr$corr_v = c()
 for( gene in gene_list$geneID){
  # for (gene in rownames(counts_sign_gene)){
  y = log(as.numeric(counts_sign_gene[rownames(counts_sign_gene)==gene,]) + 1)
  x= sample_grade$Grade
  gene_list$corr_v[gene_list$geneID==gene] <- cor( x,y )
  # gene_corr$corr_v[gene_corr$geneID==gene] <- cor( x,y)
  
  }

head(gene_corr)
gene_corr$geneID[which(abs(gene_corr$corr_v) >0.6)] 

# writeLines(gene_corr$geneID[which(abs(gene_corr$corr_v) >0.6)], "data/all_gene_corr_0.6.txt")
gene_corr_list <- read.table("data/all_gene_corr_0.6.txt")
gene_corr_list <- unlist(lapply(strsplit(gene_corr_list$V1, '\\.'), '[[',1))
(gene_corr_symbol <- getBM(attributes=c('entrezgene_id','uniprot_gn_symbol', 'ensembl_gene_id'), filters = 'ensembl_gene_id',
      values = gene_corr_list, mart = ensembl))

gene_corr_sym <- gene_corr_symbol %>% filter(!is.na(entrezgene_id))
mf_go <- enrichGO(gene_corr_symbol$entrezgene_id, 'org.Hs.eg.db', ont= "MF", pvalueCutoff=0.01) # %>%
# as.data.frame() 
plot_name <- "Correlated genes"
dotplot(mf_go, title = paste0(plot_name, ": GO-MF"))


bp_go <- enrichGO(gene_corr_symbol$entrezgene_id, 'org.Hs.eg.db', ont= "BP", pvalueCutoff=0.01)  #%>%
# as.data.frame() 

dotplot(bp_go, title = paste0(plot_name, ": GO-BP"))

cc_go <- enrichGO(gene_corr_symbol$entrezgene_id, 'org.Hs.eg.db', ont= "CC", pvalueCutoff=0.01) # %>%
# as.data.frame() 
dotplot(cc_go, title = paste0(plot_name, ": GO-CC"))
# emapplot(enrichplot::pairwise_termsim(bp_go))
kegg <- enrichKEGG(gene_corr_symbol$entrezgene_id, 'hsa', pvalueCutoff=0.01)   # for human
# kegg <- enrichKEGG(gene_corr_symbol$entrezgene_id, 'mmu', pvalueCutoff=0.01)   # for mouse
dotplot(kegg,  title = paste0(plot_name, ": KEGG"))



head(gene_list)
gene_list$geneID[which(abs(gene_list$corr_v) >0.6)] 
# writeLines(gene_list$geneID[which(abs(gene_list$corr_v) >0.6)] , "data/signif_gene_corr_0.6.txt")
# write.csv(gene_list , "data/signif_gene_allSets_corr.csv")
# gene_cor_list <- gene_list$geneID[which(abs(gene_list$corr_v) >0.6)] 
# gene_cor_list <- unlist(lapply(strsplit(gene_cor_list, '\\.'), '[[',1))
# converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id','hgnc_symbol','uniprot_gn_symbol'), filters = 'ensembl_gene_id',
#                    values = gene_cor_list, mart = ensembl)
gene_df_count <- counts_sign_gene
rownames(gene_df_count) = unlist(lapply(strsplit(rownames(gene_df_count), '\\.'), '[[',1))

converted_gene_df_count_uni <- getBM(attributes=c('uniprot_gn_symbol', 'ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = rownames(gene_df_count), mart = ensembl)
converted_gene_df_count_uni$uniprot_gn_symbol[which(duplicated(converted_gene_df_count_uni$uniprot_gn_symbol))]
converted_gene_df_count_uni_flt <- converted_gene_df_count_uni[converted_gene_df_count_uni$uniprot_gn_symbol!="",] 
converted_gene_df_count_uni_flt <- converted_gene_df_count_uni_flt[!duplicated(converted_gene_df_count_uni_flt$uniprot_gn_symbol), ]

gene_df_count_symbol_full <- gene_df_count %>% rownames_to_column(var ="ensembl_gene_id") %>% left_join(converted_gene_df_count_uni_flt) %>%
  filter(!is.na(uniprot_gn_symbol))
dim(gene_df_count_symbol_full)
# rownames(gene_df_count_symbol_full) <- gene_df_count_symbol_full$uniprot_gn_symbol
df_count_symbol <- gene_df_count_symbol_full %>% select(-"ensembl_gene_id")%>% select(uniprot_gn_symbol, colnames(.)[1:ncol(.)-1]) 
str(df_count_symbol)
write.table(df_count_symbol, "data/df_count_symbol.txt", sep = '\t')

# write.csv(converted, "data/signif_corGene_symbol_0.6.csv")
logFC_signif <- logFC_signif %>% select("enseble_ID",  "nVsl_logFC","nVsh_logFC", "nVsi_logFC")
logFC_signif$corr_v <- c()
for( gene in logFC_signif$enseble_ID){
  # for (gene in rownames(counts_sign_gene)){
  y = as.numeric(logFC_signif[logFC_signif$enseble_ID==gene,2:4]) 
  x= 1:3
  logFC_signif$corr_v[logFC_signif$enseble_ID==gene] <- cor( x,y )
  # gene_corr$corr_v[gene_corr$geneID==gene] <- cor( x,y)
}

logFC_signif[which(logFC_signif$corr_v >0.8),] %>% dim()

dim(logFC_signif)
plot(x,as.numeric(logFC_signif[1099,2:4]))
plot(sample_grade$Grade,as.numeric(counts_sign_gene[rownames(counts_sign_gene)[5],]))
plot(sample_grade$Grade,log(as.numeric(counts_sign_gene[5,] + 1)))









