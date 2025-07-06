### time series analysis


pdac <- metadata %>% filter(Category %in% c("PDAC", "HGD")) %>% 
      filter(Grade != "invasive")
      # filter(Path_ID %in% c("1-S-20-13193", "S20-9197", "S21-15821")) %>% arrange(Path_ID)
# pdac$ind.n <- factor(rep(rep(1:3,each=4),1))
pdac_counts <- counts %>% select(pdac$Sample_ID)



dds <-  DESeqDataSetFromMatrix(countData=round(pdac_counts), colData = pdac, design = ~   Grade + Category + Grade:Category)


ddsTC <- DESeq(dds, test="LRT", reduced = ~ Grade + Category)
resTC <- results(ddsTC)
resTC
resultsNames(ddsTC)
res30 <- results(ddsTC, name="Grade_invasive_vs_normal", test="Wald")

plot_name <- "Change from normal:\nPDAC"

res30[which(res30$padj < 0.05), ]

### pathway enrichment
res3 <- res30 %>% as.data.frame() %>% 
  filter(!is.na(stat)) %>%
  rownames_to_column(var ="ENSEMBL") %>% 
  mutate(ENSEMBL = unlist(lapply(strsplit(ENSEMBL, '\\.'), '[[',1)) ) %>%
  left_join(ids) %>%
  filter(!is.na(ENTREZID), !is.na(stat)) %>%
  group_by(ENTREZID) %>% 
  summarize(stat=mean(stat)) 

ranks <- res3$stat
names(ranks) <- res3$ENTREZID

# pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/h.all.v2024.1.Hs.symbols.gmt") # for hallmark pathways
# pathways.reactome <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/c2.cp.reactome.v2024.1.Hs.symbols.gmt") # for reactome pathways


fgseaRes_hallmark <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)


# # Tidy the results:
fgseaResTidy_hallmark <- fgseaRes_hallmark %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)



### hallmark pathway plot
fgseaRes_hallmark %>%
   # filter(padj <= 0.05 )  %>%
  mutate(pathway = str_remove(pathway, "HALLMARK_"),
         color = ifelse(padj <= 0.05, "sig", "ns"),
         color = factor(color , levels= c("sig", "ns"))) %>%
  ggplot(aes(reorder(pathway, NES), NES, fill = color)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("red",NA))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Hallmark Pathways", y="Normalized Enrichment Score",
       title=paste0(plot_name) )





fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("Grade", "Category"), returnData = TRUE)
# fiss$Grade <- as.numeric(as.factor(fiss$Grade))
ggplot(fiss,
       aes(x = Grade, y = count, color = Category, group = Category)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()

betas <- coef(ddsTC)
colnames(betas)
topGenes <- head(order(resTC$padj),50)
mat <- betas[topGenes,c(-1)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap::pheatmap(as.matrix(mat) , breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
