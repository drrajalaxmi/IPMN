
### glmmseq
library(glmmSeq)
library(lme4)
library(fgsea)

condition_119 <- condition %>% #filter(!Path_ID %in% c("SNBU23-19368", "SNBU22-14478")) %>%
 # filter(Path_ID %in% c("S21-695", "SNBU22-17893", "S21-3908","1-S-23-07849", "1-S-23-08528"))     %>%
 # filter(!Grade %in% c("invasive", "high grade")) %>%
  #filter(Category%in% c( "PDAC", "HGD")) %>%
  
  mutate(Category= factor(Category, levels= c("LGD", "HGD", "PDAC")),
         Grade= as.numeric(factor(Grade, levels= c("normal", "low grade", "high grade", "invasive"))) )%>%
  arrange(Path_ID, Grade,Category  )  %>% as.data.frame()

rownames(condition_119) <- condition_119$Sample_ID
counts_119 <- counts %>% dplyr::select(condition_119$Sample_ID)


#df_count <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
df_count <-read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)
colnames(df_count) <- stringr::str_replace(colnames(df_count), "X", "BE")


counts <- df_count  %>% dplyr::select(counts_119$Sample_ID)


counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts_119), colData = condition_119, design = ~ 1)
dds <- DESeq(counts.DEseq)
disp <- setNames(dispersions(dds), rownames(counts_119))

# summary(lme4::lmer(as.vector(unlist(counts_119[4,]))~  Grade+ (1 | Path_ID), data = condition_119))

res <- glmmSeq( ~ Category * Grade + (1 | Path_ID),
                countdata = counts_119,
                   metadata = condition_119[c("Category", "Sample_ID", "Grade", "Path_ID")],
               removeSingles = FALSE,
              #returnList= TRUE,
                   dispersion =disp, verbose= TRUE)


glmm_stat <- summary(res)
colnames(glmm_stat)

glmm_stat <- glmm_stat %>% as.data.frame()

head(glmm_stat)

glmm_stat[order(glmm_stat[, 'P_Category:Grade']), ]






library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

str(ddsTC)
model.matrix(fission, ~ strain + minute + strain:minute)
#### timeseries data
str(condition_TC)
condition_TC[is.na(condition$Category),]
condition_TC <- condition %>% mutate(Grades = factor(Grade, levels = c("normal", "low grade", "high grade", "invasive")),
                                     Category = factor(Category, levels=c("LGD", "HGD", "PDAC")))
counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts), colData = condition_TC, design = ~Category + Grades + Category:Grades)

ddsTC <- DESeq(counts.DEseq, test="LRT", reduced = ~ Category + Grade)

resTC <- results(ddsTC)
resTC$symbol <- mcols(dds)$symbol
dim(resTC[order(resTC$padj),])
dim(resTC[which(resTC$padj < 0.05),])
resultsNames(dds)
comparison <- resultsNames(dds)[length(resultsNames(dds))] #lists the coefficients


LFC <- lfcShrink(dds, coef = comparison, res=res, type = "apeglm")

lfc = 2
pval = 0.05

tab = data.frame(logFC = resTC$log2FoldChange, negLogPval = -log10(resTC$padj))#make a data frame with the log2 fold-changes and adjusted p-values


signGenes = abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval)
gene <- rownames(resTC[which(signGenes),])

plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), main = paste0("MLH", "\nSignificant genes from DEseq")) #replace main = with your title
points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")



#Show the cut-off lines
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)

mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.6, line = 0.5)


### continue with Time series
fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("Grade", "Category"), returnData = TRUE)
fiss$Grade <- as.numeric(fiss$Grade)
ggplot(fiss,
       aes(x = Grade, y = count, color = Category, group = Category)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()


betas <- coef(ddsTC)
colnames(betas)
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,3)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
heatmap(mat)
        , breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)


# data(PEAC_minimal_load)
# disp <- apply(tpm, 1, function(x) {
#   (var(x, na.rm = TRUE)-mean(x, na.rm = TRUE))/(mean(x, na.rm = TRUE)**2)
# })
# MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#                      countdata = tpm[1:2, ],
#                      metadata = metadata,
#                      dispersion = disp,
#                      verbose = TRUE)
# names(attributes(MS4A1glmm))



#hallmarkpathway enrichment

resTC$row <- rownames(resTC)
dim(resTC)


rownames(resTC) = unlist(lapply(strsplit(rownames(resTC), '\\.'), '[[',1))

converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = rownames(resTC), mart = ensembl)




glmm_stat$row <- rownames(glmm_stat)
dim(stats)


rownames(glmm_stat) = unlist(lapply(strsplit(rownames(glmm_stat), '\\.'), '[[',1))

glmm_stat_sig <- glmm_stat %>% filter(`P_Category:Grade` <= 0.05)

res2 <- glmm_stat_sig %>% 
  rownames_to_column(var ="ensembl_gene_id") %>% left_join(converted) %>%
  filter(!is.na(entrezgene_id)) %>%
  group_by(entrezgene_id) %>% 
  summarize(stat=mean(as.numeric(Chisq_Category:Grade), na.rm =TRUE)) 

as.numeric(glmm_stat$`Chisq_Category:Grade`)

ranks <- res2$stat
names(ranks) <- res2$entrezgene_id



### downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
# pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/mm_data_240814/mm_data_240814/data/mh.all.v2024.1.Mm.entrez.gmt")   # for mouse
# head(pathways.hallmark)
pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/h.all.v2024.1.Hs.entrez.gmt") # for human

#Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

# To see what genes are in each of these pathways:
gene.in.pathway <- pathways.hallmark %>% 
  enframe("pathway", "entrezgene_id") %>% 
  unnest(cols = c(entrezgene_id)) %>% 
  mutate(entrezgene_id = as.numeric(entrezgene_id)) %>%
  inner_join( res2, by="entrezgene_id")

fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy %>% mutate(pathway = str_remove(pathway, "HALLMARK_")), aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Hallmark Pathways", y="Normalized Enrichment Score",
       title=paste0("Time Series differential gene hallmark pathway enrichment") )








