
##______________________________
## Rajalaxmi 
# 240825 
# mm_data_240814
##-------------------

#library
library(readxl)
library(ggfortify) 
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(apeglm)
library(biomaRt)
library(clusterProfiler)
library(stringr)
library(tibble)
library(fgsea)
# library(EnhancedVolcano)
library(org.Hs.eg.db)
library(pathview)
#BiocManager::install("biomaRt")
# BiocManager::install('EnhancedVolcano')

mb_df <- read_excel("/Users/rajalaxmi/Downloads/gene_counts (1).xlsx")
write.table(df_txt,"/Users/rajalaxmi/Downloads/genecounts.txt", row.names = FALSE, sep ="\t", quote=FALSE)
rownames(df_txt)
mouse_df <- read.csv("/Users/rajalaxmi/Desktop/41598_2017_BFsrep40508_MOESM313_ESM.csv")
write.table(mouse_df,"/Users/rajalaxmi/Downloads/signatureMatrixMouse.txt", row.names = FALSE, sep ="\t", quote=FALSE)

df_count <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
##...... a below till counts
 
counts$ENSEMBL = unlist(lapply(strsplit(rownames(df_count), '\\.'), '[[',1))


dim(counts)
dim(counts_txt)
counts_txt <- counts %>% rownames_to_column("ENSEMBL") %>% left_join(ids%>% dplyr::select(-ENTREZID)) 
counts_txt <- counts_txt[!duplicated(counts_txt$SYMBOL),]  
counts_txt <- counts_txt[!is.na(counts_txt$SYMBOL),]  


counts_txt <- counts_txt %>% select(-c(ENSEMBL )) %>% relocate(SYMBOL)

write.table(counts_txt,"data/tpm_untrimmed.txt", row.names = FALSE, sep ="\t", quote=FALSE)

read.csv("data/tpm_untrimmed.txt", header=TRUE, sep ="\t")


df_txt <- df_txt[!duplicated(df_txt$gene_symbol),]


  
  
df_txt  <- df_txt[!is.na(df_txt$gene_symbol),]

df_count$ensembl_gene_id <- rownames(df_count)
df_count <- df_count %>% separate(ensembl_gene_id, c("ensembl_gene_id", NA) , remove = TRUE)
df_gene <- left_join(df_count,mb_df[1:2] )
dim(df_gene[146])

df_txt <-df_gene %>% select(gene_symbol, 1: 144) 

BE02_case_list <-read_excel("data/BE02_case_list.xlsx", sheet = "FINAL")
column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE02_",BE02_case_list$Sample_ID, "A")
str(BE02_case_list)

# write.table(BE02_case_list %>% select(Category, Sample_ID, Path_ID), "data/case_lits.tsv",col.names = FALSE, sep = "\t")
# sample_list <- c("BE02_18A",  "BE02_59A",  "BE02_60A",  "BE02_61A",  "BE02_74A" , "BE02_75A",  "BE02_101A", "BE02_102A", "BE02_104A", "BE02_105A", "BE02_113A", "BE02_114A")

 #df_count <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
df_count <-read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)
colnames(df_count) <- paste0(stringr::str_replace(colnames(df_count), "X", "BE02_"), "A")
condition <- BE02_case_list %>% filter(!Sample_ID %in% c("BE02_65A" , "BE02_77A", "BE02_78A",  "BE02_80A",  "BE02_81A",  "BE02_120A")) %>%
                           # filter(Grade !="normal") %>%
                       filter(Sample_ID %in% colnames(df_count))

length(unique(condition$Path_ID))

for(col in colnames(df_count)){ 
write.table(df_count %>% rownames_to_column("ENSEMBL") %>% 
              mutate(sample = round(df_count[[col]]),
                     ENSEMBL = unlist(lapply(strsplit(ENSEMBL, '\\.'), '[[',1))) %>%
              select(ENSEMBL,sample) , paste0("RNAseqCNV/count_", col,".txt"),sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
BE02_case_list$Disease[BE02_case_list$Sample_ID %in% c("BE02_65A" , "BE02_77A", "BE02_78A",  "BE02_80A",  "BE02_81A",  "BE02_120A")]

length(unique(condition$Path_ID))
condition %>% group_by(Disease) %>%
              summarize(n = length(Sample_ID))            
  
counts <- df_count  %>% dplyr::select(condition$Sample_ID)
# colnames(counts) <- paste0(condition$Path_ID,"(",condition$Sample_ID,")",condition$Disease
#                            )
length(colnames(counts))
# counts <-read.csv("data/counts.csv", header = TRUE,stringsAsFactors=FALSE)
# rownames(counts) <- counts[[1]]
# counts <- counts[-1]
# head(counts[,1:5])

gene_ensemble = unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))

symbols <- mapIds(org.Hs.eg.db, keys = gene_ensemble,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(gene_ensemble, names(symbols))]

set.seed(1200)
condition$Disease <- factor(condition$Disease, levels = c("Isolated LGD" ,
                                                          "Isolated HGD" ,
                                                          "Progressor LGD w HGD" ,
                                                          "Progressor LGD w PDAC" ,
                                                          "Progressor HGD w PDAC",
                                                          "PDAC" ) )

cnt_pca <- prcomp(t(counts_tpm))
#summary(cnt_pca)

# plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (all Grades)") 
# pca.plot <- autoplot(cnt_pca,
#                      data = t(counts)) 
unique(condition$Disease)
colours <- c("Isolated LGD" = "#2F56B5",#115C80", 
             "Isolated HGD" ="#2295C7", #309898" ,#479685", 
             "Progressor LGD w HGD" = "#00A86B", #FF9F00", 
             "Progressor LGD w PDAC" = "#FCC92A", #FF9F00", 
             "Progressor HGD w PDAC" = "#F4691E", 
             "PDAC" = "#AD0000")

autoplot(cnt_pca, data = condition,color = "Disease", size = 2) +
  scale_color_manual(values = paste0(colours)) +
  theme_minimal()
plotly::ggplotly(autoplot(cnt_pca, data = condition,color = "Disease", size = 2.5 ) ) 
                   scale_color_manual(values = paste0(colours)) )
                   theme_classic())
         #label.size = 3)#label.size = 3)fill = 

count_mtx <- as.matrix(counts)

Ca_cor <-round(cor(count_mtx), 2)
# heatmaply::heatmaply((Ca_cor))
# df_val <- df_mtx[!rowSums(df_mtx)==0, ]
# dim(df_val)

heatmap(Ca_cor, margins = c(5, 5), scale = "none")

# colnames(condition)
#autoplot(kmeans(t(counts), 5), data = t(counts))
# Inspired by https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/01c_RNAseq_count_distribution.html
mean_counts <- apply(counts, 1, mean)
variance_counts <- apply(counts, 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()

# write.csv(condition,"data/sample_mtx.csv")
unique(condition$Disease)



# LGD
# -3 way comparison: isolated LGD vs. progressor LGD w HGD vs. progressor LGD w PDAC
# -2 way comparisons: isolated LGD vs. progressor LGD w HGD
# -2 way comparisons: isolated LGD vs. progressor LGD w PDAC
# -2 way comparisons: isolated LGD vs. progressor LGD all (combined HGD or PDAC)
# HGD
# -2 way comparisons: isolated HGD vs. progressor HGD w PDAC
# Normal
# -3 way comparison: normal-LGD vs. normal-HGD vs. normal-PDAC
# -2 way comparisons: normal-LGD vs. (normal-HGD and normal-PDAC)
# -2 way comparisons: (normal-LGD and normal-HGD) vs. normal-PDAC

# LGD
# -3 way comparison: isolated LGD vs. progressor LGD w HGD vs. progressor LGD w PDAC
# plot_name <- "Isolated LGD vs. Progressor LGD w HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#           Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w HGD")))
# plot_name <- "Isolated LGD vs. Progressor LGD w HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#           Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w HGD")))


# -2 way comparisons: isolated LGD vs. progressor LGD w HGD


# -2 way comparisons: isolated LGD vs. progressor LGD w PDAC

# plot_name <- "Isolated LGD vs. Progressor LGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#           Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w PDAC")))


# -2 way comparisons: isolated LGD vs. progressor LGD all (combined HGD or PDAC)
# plot_name <- "isolated LGD vs. progressor LGD all (combined HGD or PDAC)"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC")) %>%
#   mutate( Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#     Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD")))

# plot_name <- "Non PDACs vs. PDACs"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD","Isolated HGD", "Progressor LGD w HGD","Progressor HGD w PDAC", "Progressor LGD w PDAC", "PDAC")) %>%
#                mutate( Disease = ifelse(Category =="PDAC", "PDAC","non_PDAC" ) ,
#                        Disease = factor(Disease, levels= c("non_PDAC", "PDAC")))
# sheet_name <- "NonPDAC_Vs_PDAC"

# plot_name <- "Progressor LGD w HGD Vs Progressor LGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Progressor LGD w HGD", "Progressor LGD w PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#     Disease = factor(Disease, levels= c("Progressor LGD w HGD", "Progressor LGD w PDAC")))

# plot_name <- "isolated LGD vs. progressor LGD w HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD")) %>%
#   mutate( #Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#     Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w HGD")))

# plot_name <- "progressor LGD w HGD Vs progressor HGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Progressor LGD w HGD", "Progressor HGD w PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
    # Disease = factor(Disease, levels= c("Progressor LGD w HGD", "Progressor HGD w PDAC")))

# HGD
# -2 way comparisons: isolated HGD vs. progressor HGD w PDAC
# plot_name <- "Isolated LGD vs Isolated HGD "
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Isolated HGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Isolated LGD", "Isolated HGD")))
# sheet_name <- "IsoLGD_Vs_IsoHGD"

# plot_name <- "isolated HGD vs progressor HGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated HGD", "Progressor HGD w PDAC")) %>%
#   mutate( Disease = factor(Disease, levels= c("Isolated HGD", "Progressor HGD w PDAC")))

# Normal
# -3 way comparison: normal-LGD vs. normal-HGD vs. normal-PDAC

# plot_name <- "Normal-LGD vs Normal-HGD vs Normal-PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Normal-HGD","Normal-PDAC")) %>%
#   mutate( Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#     Disease = factor(Disease, levels= c("Normal-Non-Progressor", "Normal-PDAC"))
#     )

# plot_name <- "Normal-LGD vs. Normal-HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Normal-HGD")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#     Disease = factor(Disease, levels= c("Normal-LGD", "Normal-HGD")))


# plot_name <- "Normal-HGD Vs Normal-PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-HGD", "Normal-PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#     Disease = factor(Disease, levels= c("Normal-HGD", "Normal-PDAC")))


# plot_name <- "Normal-PDAC Vs progressor LGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-PDAC", "progressor LGD w PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#     Disease = factor(Disease, levels= c("Normal-PDAC", "progressor LGD w PDAC")))


# plot_name <- "Normal-PDAC Vs progressor HGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-PDAC", "progressor HGD w PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#     Disease = factor(Disease, levels= c("Normal-PDAC", "progressor HGD w PDAC")))

# -2 way comparisons: normal-LGD vs. (normal-HGD and normal-PDAC)
# plot_name <- "normal-LGD vs. (normal-HGD and normal-PDAC)"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
#   mutate( Disease = ifelse(Disease =="Normal-LGD", "Normal-LGD","Normal-Progressor" ) ,
#           Disease = factor(Disease, levels= c("Normal-LGD","Normal-Progressor")))


# -2 way comparisons: (normal-LGD and normal-HGD) vs. normal-PDAC
# plot_name <- "(normal-LGD and normal-HGD) vs. normal-PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
#   mutate( Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#           Disease = factor(Disease, levels= c("Normal-Non-Progressor","Normal-PDAC")))


### normal to higer grades
# plot_name <- "Normal LGD Vs Isolated LGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Isolated LGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Normal-LGD","Isolated LGD")))


# plot_name <- "Normal HGD Vs Isolated HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-HGD", "Isolated HGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Normal-HGD","Isolated HGD")))

# plot_name <- "Normal HGD Vs Progressor LGD w HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-HGD", "Progressor LGD w HGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Normal-HGD","Progressor LGD w HGD")))

# 
# plot_name <- "Normal LGD Vs Progressor LGD w HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Progressor LGD w HGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Normal-LGD","Progressor LGD w HGD")))



# plot_name <- "Normal PDAC Vs PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Normal-PDAC", "PDAC")) %>%
#   mutate( Disease = factor(Disease, levels= c("Normal-PDAC","PDAC")))


# plot_name <- "Progressor LGD w HGD Vs Isolated HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Progressor LGD w HGD", "Isolated HGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Progressor LGD w HGD", "Isolated HGD")))

# plot_name <- "Progressor LGD w PDAC Vs Progressor HGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Progressor LGD w PDAC", "Progressor HGD w PDAC")) %>%
#   mutate( Disease = factor(Disease, levels= c("Progressor LGD w PDAC", "Progressor HGD w PDAC")))

# plot_name <- "Progressor HGD w PDAC Vs PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Progressor HGD w PDAC", "PDAC")) %>%
#   mutate( Disease = factor(Disease, levels= c("Progressor HGD w PDAC", "PDAC")))

# plot_name <- "Progressor LGD w PDAC Vs PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Progressor LGD w PDAC", "PDAC")) %>%
#   mutate( Disease = factor(Disease, levels= c("Progressor LGD w PDAC", "PDAC")))

# levels(sample_mtx$Disease)
# 
# plot_name <- "Normal Vs Low grade"
# sample_mtx <- condition %>% filter(Grade %in% c("normal", "low grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("normal", "low grade")))

# plot_name <- "Normal Vs Invasive"
# sample_mtx <- condition %>% filter(Grade %in% c("normal", "invasive")) %>%
#   mutate( Grade = factor(Grade, levels= c("normal", "invasive")))
#
# plot_name <- "Normal Vs High grade"
# sample_mtx <- condition %>% filter(Grade %in% c("normal", "high grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("normal", "high grade")))
# # 
# plot_name <- "Low grade Vs Invasive"
# sample_mtx <- condition %>% filter(Grade %in% c("low grade", "invasive")) %>%
#   mutate( Grade = factor(Grade, levels= c("low grade", "invasive")))
# 
# plot_name <- "Low grade Vs High grade"
# sample_mtx <- condition %>% filter(Grade %in% c("low grade", "high grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("low grade", "high grade")))
# 
# plot_name <- "Invasive Vs High grade"
# sample_mtx <- condition %>% filter(Grade %in% c("invasive", "high grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("invasive", "high grade")))

# c("normal", "invasive", "low grade", "high grade")


### 


plot_name <- "Isolated LGD vs Isolated HGD "
sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Isolated HGD")) %>%
  mutate( Disease = factor(Disease, levels= c("Isolated LGD", "Isolated HGD")))
sheet_name <- "IsoLGD_Vs_IsoHGD"

plot_name <- "isolated HGD vs progressor HGD w PDAC"
sample_mtx <- condition %>% filter(Disease %in% c("Isolated HGD", "Progressor HGD w PDAC")) %>%
  mutate( Disease = factor(Disease, levels= c("Isolated HGD", "Progressor HGD w PDAC")))
sheet_name <- "IsoHGD_Vs_ProHGDwPDAC"


plot_name <- "isolated LGD vs. progressor LGD all (combined HGD or PDAC)"
sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC")) %>%
  mutate( Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
    Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD")))
sheet_name <- "IsoLGD_Vs_ProLGD(All)"


plot_name <- "isolated LGD vs. progressor LGD w HGD"
sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD")) %>%
  mutate( #Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
    Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w HGD")))
sheet_name <- "IsoLGD_Vs_ProLGDwHGD"


plot_name <- "Isolated LGD vs. Progressor LGD w PDAC"
sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w PDAC")) %>%
  mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
          Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w PDAC")))
sheet_name <- "IsoLGD_Vs_ProLGDwPDAC"

levels(sample_mtx$Disease)




# size <- sample_mtx %>% #mutate(Diseases = as.character(Disease)) %>% 
#                       group_by(Disease) %>% summarise(Diseases = unique(Diseases) ,
#                                                     n = length(Sample_ID))
counts_mtx <- counts %>% dplyr::select(which(colnames(counts) %in% sample_mtx$Sample_ID)) %>% as.matrix()

dim(sample_mtx)
dim(counts_mtx)
#Import to DEseq2
# sample_mtx$Category = factor(sample_mtx$Category , levels = c("LGD", "HGD","PDAC"))
# m1 <- model.matrix(~  Category , sample_mtx)
# all.zero <- apply(m1, 2, function(x) all(x==0))
# all.zero
# 
# idx <- which(all.zero)
# m1 <- m1[,-idx]
# unname(m1)

counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts_mtx), colData = sample_mtx, design = ~ Disease)
dds <- DESeq(counts.DEseq)
# dds<- DESeq(counts.DEseq, full = m1)
resultsNames(dds)
comparison <- resultsNames(dds)[length(resultsNames(dds))] #lists the coefficients

# results(dds, contrast=list("CategoryHGD", "CategoryPDAC"))

# cnt_pca <- prcomp(t(counts_mtx))
# summary(cnt_pca)
# 
# plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (Normal Vs Low grade)") 
# pca.plot <- autoplot(cnt_pca,
#                      data = t(counts_mtx)) 
# autoplot(cnt_pca, data = sample_mtx,color = "Grade", shape = FALSE, label.size = 3)


# par(mar=c(4,4,1,1))
# plotDispEsts(dds)
# 
# rld <- rlog(dds, blind = TRUE)
# 
# #plotPCA from DEseq2 plots uses the top 500 genes:
# data = plotPCA(rld, intgroup = c("Grade"), returnData = TRUE)
# p <- ggplot(data, aes(x = PC1, y = PC2, color =  Grade))
# p <- p + geom_text(aes(label = name)) + 
#   labs(title = paste0("PCA ", plot_name))+ theme ()
# print(p)



res <- results(dds) #, contrast = list(resultsNames(dds1)[2], resultsNames(dds1)[3]))
# rownames(res) = unlist(lapply(strsplit(rownames(res), '\\.'), '[[',1))
res2 <- res %>% as.data.frame() %>% 
  rownames_to_column(var ="ENSEMBLID") %>% 
 
 # filter(!is.na(ENTREZID), !is.na(stat)) %>%
  # filter( ENSEMBLID  != "ENSG00000273730.1") %>%
  # filter(abs(log2FoldChange) > 2 & padj < 0.05 ) %>%
  mutate(ENSEMBL = unlist(lapply(strsplit(ENSEMBLID, '\\.'), '[[',1)),
         comparison = plot_name) %>% relocate(comparison) %>%
  left_join(ids)   %>%
  filter(!is.na(ENTREZID), !is.na(padj)) %>%
  filter( !str_detect(SYMBOL, "LOC|RNA" )) %>%
#  filter(!is.na(ENTREZID) ) %>% 
   arrange(desc(abs(log2FoldChange)))

print(plot_name)

# wb <- createWorkbook("BE02_sig_DE_gene")
## Add a worksheets
# addWorksheet(wb, "IsoLGD_Vs_IsoHGD", gridLines = FALSE)
# addWorksheet(wb, "IsoHGD_Vs_ProHGDwPDAC", gridLines = FALSE)
# addWorksheet(wb, "IsoLGD_Vs_ProLGD(All)", gridLines = FALSE)
# addWorksheet(wb, "IsoLGD_Vs_ProLGDwHGD", gridLines = FALSE)
# addWorksheet(wb, "IsoLGD_Vs_ProLGDwPDAC", gridLines = FALSE)
addWorksheet(wb, sheet_name, gridLines = FALSE)
## write data to worksheet 1
writeData(wb, sheet = sheet_name, res2, rowNames = TRUE)
saveWorkbook(wb, "data/BE02_sig_DE_gene.xlsx", overwrite = TRUE)


tab_res = data.frame(logFC = res2$log2FoldChange, negLogPval = -log10(res2$padj))  #make a data frame with the log2 fold-changes and adjusted p-values
  
signGenes_res = (abs(tab_res$logFC) > lfc & tab_res$negLogPval > -log10(pval))

#Substitute the '????' with a comparison, selected from the resultsNames(dds) shown above
# LFC <- lfcShrink(dds, contrast = c('Disease',"Isolated LGD", "Progressor LGD w HGD"), res=res, type = 'normal')

# LFC <- lfcShrink(dds, coef = comparison, res=res, type = "normal")

# #Add a title to reflect your comparison 
# plotMA(LFC, main = plot_name, cex = 0.5, ylim=c(-2,2))
# #The distribution of p-values
# hist(LFC$pvalue, breaks = 50, col = 'grey', main = plot_name, xlab = 'p-value')
# 
# #The false-discovery rate distribution
# hist(LFC$padj, breaks = 50, col = 'grey', main = plot_name, xlab = 'Adjusted p-value')
# 
# #Allow for more space around the borders of the plot
# par(mar = c(5, 4, 4, 4))

#Set your log-fold-change and p-value thresholds
lfc = 2
pval = 0.05

tab_res = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj)) #make a data frame with the log2 fold-changes and adjusted p-values

           
 
tab = data.frame(logFC = LFC$log2FoldChange, negLogPval = -log10(LFC$padj))#make a data frame with the log2 fold-changes and adjusted p-values


# ensembl <- useEnsembl(biomart = 'genes',
#            dataset = 'hsapiens_gene_ensembl', host ="useast.ensembl.org",
#            version = 113)



# (which(str_detect(listAttributes(ensembl)$name, "symbol")))
# listAttributes(ensembl)$name[63]

# converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id','uniprot_gn_symbol', 'hgnc_symbol'), filters = 'ensembl_gene_id',
#                    values = gene_ensemble, mart = ensembl)



# 
# rownames(res)[ which((!(is.na(converted$hgnc_symbol) | converted$hgnc_symbol=='' ) )& abs(LFC$log2FoldChange) >1)] <- converted$hgnc_symbol[which((!(is.na(converted$hgnc_symbol) | converted$hgnc_symbol=='' ) )& abs(LFC$log2FoldChange)>1) ]
# converted[which(!(is.na(converted$hgnc_symbol) | converted$hgnc_symbol=='' ) ) ,]
# names_1 <- converted$hgnc_symbol[which(!(is.na(converted$hgnc_symbol) | converted$hgnc_symbol=='' ) ) ]
# 
# EnhancedVolcano(LFC,
#                 lab = rownames(res),
#                 x = 'log2FoldChange',
#                 y = 'pvalue')

#Genes with a fold-change greater than 2 and p-value<0.05:
# signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
# 
signGenes_res = (abs(tab_res$logFC) > lfc & tab_res$negLogPval > -log10(pval))
# 
# gene_res <- rownames(res[which(signGenes_res),])
# gene <- rownames(LFC[which(signGenes),])


# gene_list <- unique(c(gene,  gene_list))



# logFC <- LFC %>% as.data.frame() %>% select(log2FoldChange) %>% rownames_to_column("enseble_ID") %>% 
#           rename(nVsh_logFC = log2FoldChange) %>% #head()
#           left_join(logFC)
# 
# 
# 
# 
# logFC_signif <- logFC %>% filter(enseble_ID %in% gene_list)
# head(logFC_signif)
# 
# write.csv(logFC_signif,  "data/signif_gene_log2FC.csv")
# write.csv(data.frame("geneID" = c(gene_list[,2])),  "data/signif_genelist_lcf.csv")


# 
# length(gene_list)
# length(gene)
# sum(!gene %in% gene_res)

# plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
#      ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nSignificant genes from lcfshrink")) #replace main = with your title

# gene_ensemble <- rownames(LFC)[which(LFC$log2FoldChange < -5)]
# gene_ensemble = unlist(lapply(strsplit(gene_ensemble, '\\.'), '[[',1))
# 
# converted <- getBM(attributes=c('entrezgene_id','hgnc_symbol'), filters = 'hgnc_symbol',
#                    values = gene_ensemble, mart = ensembl)
#Colour these red
# points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")
# plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
#      ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nSignificant genes from DEseq")) #replace main = with your title
# points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")

plot(tab_res, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nSignificant genes from DEseq")) #replace main = with your title
points(tab_res[signGenes_res, ], pch = 16, cex = 0.5, col = "red")



#Show the cut-off lines
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)

mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.6, line = 0.5)

# plot(tab_res, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
#      ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nColored for significant gene from lcfshrink")) #replace main = with your title
# points(tab_res[signGenes, ], pch = 16, cex = 0.5, col = "red")


# plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
#      ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nColored for significant gene from DEseq")) #replace main = with your title
# points(tab[signGenes_res, ], pch = 16, cex = 0.5, col = "red")




# #increased expression
# attach(as.data.frame(res))
# 
# #The total number of DEGs with an adjusted p-value<0.05
# summary(res, alpha=0.05)
# 
# 
# #The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
# a <-sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
# 
# #Decreased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange < 0) #any fold-change
# b <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange > (2)) #fold-change greater than 2
# 
# 
# #Increased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
# c <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange < -2) #fold-change greater than 2
# 
# cat("Sample Size: \n",paste0(size[1,]), "\n",paste0(size[2,]), "\n \nTotal DE genes: ", a, "\nUp regulated: ", b, "\nDown regulated: ", c)

# #increased expression
# attach(as.data.frame(LFC))
# 
# #The total number of DEGs with an adjusted p-value<0.05
# summary(LFC, alpha=0.05)
# 
# 
# #The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
# a <-sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
# 
# #Decreased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange < 0) #any fold-change
# b <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange > (2)) #fold-change greater than 2
# 
# 
# #Increased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
# c <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange < -2) #fold-change greater than 2
# 
# cat("Sample Size: \n",paste0(size[1,]), "\n",paste0(size[2,]), "\n \nTotal DE genes: ", a, "\nUp regulated: ", b, "\nDown regulated: ", c)

# attach(as.data.frame(LFC))
attach(as.data.frame(res))
#At this stage it may be useful to create a copy of the results with the gene version removed from the gene name, to make it easier for you to search for the gene name etc. 
#The rownames currently appear as 'ENSG00000175197.12, ENSG00000128272.15' etc.
#To change them to 'ENSG00000175197, ENSG00000128272'
# LFC.gene = as.data.frame(LFC)
LFC.gene = as.data.frame(res)

#Some gene names are repeated if they are in the PAR region of the Y chromosome. Since dataframes cannot have duplicate row names, we will leave these gene names as they are and rename the rest.
whichgenes = which(!grepl('PAR', rownames(LFC.gene)))
rownames(LFC.gene)[whichgenes] = unlist(lapply(strsplit(rownames(LFC.gene)[whichgenes], '\\.'), '[[',1))

#subset the significant genes
LFC.sig = LFC.gene[LFC.gene$padj < 0.05 & !is.na(LFC.gene$padj) & abs(LFC.gene$log2FoldChange) > 2,]#subset the significant genes
# LFC.sig %>% arrange(log2FoldChange) %>% t
# rownames(LFC.sig)
#We can add a column with the HGNC gene names
# mart = useMart('ensembl')
# listDatasets(mart)
# listEnsembl()
# ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
# mmusculus_gene_ensembl  - mouse
# hsapiens_gene_ensembl - human
# 
# ensembl <- useEnsembl(biomart = 'genes',
#            dataset = 'hsapiens_gene_ensembl', host ="useast.ensembl.org",
#            version = 113)
# symbol <-read.csv("data/df_count_symbol.csv")
# str(symbol)

# (which(str_detect(listAttributes(ensembl)$name, "symbol")))
# listAttributes(ensembl)$name[98]

# converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), filters = 'ensembl_gene_id',
#                    values = rownames(LFC.sig), mart = ensembl)

# converted <- converted %>% filter(ensembl_gene_id %in% gene_corr_symbol$ensembl_gene_id)

#Add gene names to the LFC.sig data-frame
# LFC.sig_gene <- LFC.sig %>% rownames_to_column(var ="ensembl_gene_id") %>% left_join(converted) %>%
#            filter(!is.na(entrezgene_id))

LFC.sig_gene <- LFC.sig %>% as.data.frame() %>% 
  rownames_to_column(var ="ENSEMBL") %>% left_join(ids) %>%
  filter(!is.na(ENTREZID) )


points(tab_res[signGenes_res, ], pch = 16, cex = 0.5, col = "red")


#LFC.sig$hgnc = converted[converted[,1] %in% rownames(LFC.sig),2]

#View the top 10 genes with the most significant (adjusted) p-values
# head(LFC.sig, n = 10)

#The largest fold-changes with a significant p-value
# LFC.sig[order(abs(LFC.sig$log2FoldChange), decreasing = TRUE),][1:10,] #add the [1:10,] to see the top 10 rows



#org.Mm.eg.db - mouse
#org.Hs.eg.db - Human
mf_go <- enrichGO(LFC.sig_gene$ENTREZID, 'org.Hs.eg.db', ont= "MF", pvalueCutoff=0.05) # %>%
  # as.data.frame() 
dotplot(mf_go, title = paste0(plot_name, ": GO-MF"))


bp_go <- enrichGO(LFC.sig_gene$ENTREZID, 'org.Hs.eg.db', ont= "BP", pvalueCutoff=0.05)  #%>%
  # as.data.frame() 

dotplot(bp_go, title = paste0(plot_name, ": GO-BP"))

cc_go <- enrichGO(LFC.sig_gene$ENTREZID, 'org.Hs.eg.db', ont= "CC", pvalueCutoff=0.05) # %>%
  # as.data.frame() 
dotplot(cc_go, title = paste0(plot_name, ": GO-CC"))
# emapplot(enrichplot::pairwise_termsim(bp_go))
kegg <- enrichKEGG(LFC.sig_gene$ENTREZID, 'hsa', pvalueCutoff=0.05)   # for human
# kegg <- enrichKEGG(LFC.sig_gene$entrezgene_id, 'mmu', pvalueCutoff=0.01)   # for mouse
dotplot(kegg,  title = paste0(plot_name, ": KEGG"))




for (id in kegg$ID){
  pathview <- pathview(gene.data  = c(str_split(kegg$geneID[kegg$ID==id],"/", simplify = TRUE)) ,
                     pathway.id = id,
                     species    = "hsa")
                    # limit = list(gene=max(c(str_split(kegg$geneID[kegg$ID==id],"/", simplify = TRUE))), cpd=1))

print(pathview)
}

#### to be removed

# 
# plot_name <- "Isolated LGD vs Isolated HGD "
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Isolated HGD")) %>%
#   mutate( Disease = factor(Disease, levels= c("Isolated LGD", "Isolated HGD")))
# sheet_name <- "IsoLGD_Vs_IsoHGD"

# plot_name <- "isolated HGD vs progressor HGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated HGD", "Progressor HGD w PDAC")) %>%
#   mutate( Disease = factor(Disease, levels= c("Isolated HGD", "Progressor HGD w PDAC")))
# sheet_name <- "IsoHGD_Vs_ProHGDwPDAC"


# plot_name <- "isolated LGD vs. progressor LGD all (combined HGD or PDAC)"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC")) %>%
#   mutate( Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#           Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD")))
# sheet_name <- "IsoLGD_Vs_ProLGD(All)"


# plot_name <- "isolated LGD vs. progressor LGD w HGD"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD")) %>%
#   mutate( #Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#     Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w HGD")))
# sheet_name <- "IsoLGD_Vs_ProLGDwHGD"


# plot_name <- "Isolated LGD vs. Progressor LGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w PDAC")) %>%
#   mutate( #Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
#     Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD w PDAC")))
# sheet_name <- "IsoLGD_Vs_ProLGDwPDAC"

levels(sample_mtx$Disease)

counts_mtx <- counts %>% dplyr::select(which(colnames(counts) %in% sample_mtx$Sample_ID)) %>% as.matrix()

counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts_mtx), colData = sample_mtx, design = ~ Disease)
dds <- DESeq(counts.DEseq)
resultsNames(dds)
comparison <- resultsNames(dds)[length(resultsNames(dds))] #lists the coefficients

res <- results(dds)
res$row <- rownames(res)


#### to be removed

#######

res$row <- rownames(res)
# dim(res)
# res.sig.gene <- res[which(abs(res$log2FoldChange) >2 & res$padj >0.05),]
# dim(res.sig.gene)
# ens2symbol <- AnnotationDbi::select(org.Mm.eg.db,
#                                     key=res$row, 
#                                     columns="SYMBOL",
#                                     keytype="ENSEMBL")
# names(ens2symbol)[1] <- "row"

# dim(LFC.sig_gene)

 
# min(res$stat[max(res$log2FoldChange, na.rm =T)], na.rm = T)

# res[order(res$stat,decreasing = TRUE ), ]
# res %>% as_tibble() %>% arrange((stat)) %>% 
#      mutate(test = log2FoldChange/lfcSE)

# rownames(res) = unlist(lapply(strsplit(rownames(res), '\\.'), '[[',1))
# 
# converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), filters = 'ensembl_gene_id',
#                    values = rownames(res), mart = ensembl)
# 
# res2 <- res %>% as.data.frame() %>% 
#   rownames_to_column(var ="ensembl_gene_id") %>% left_join(converted) %>%
#   filter(!is.na(entrezgene_id), !is.na(stat)) %>%
#   group_by(entrezgene_id) %>% 
#   summarize(stat=mean(stat)) 
# 
 # rownames(counts) = unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))
ids<-bitr(rownames(counts), fromType="ENSEMBL",
          toType=c("ENTREZID", "SYMBOL"),
          OrgDb="org.Hs.eg.db") #NA values drop by default




res$row <- rownames(res)

rownames(res) = unlist(lapply(strsplit(rownames(res), '\\.'), '[[',1))
res2 <- res %>% as.data.frame() %>% 
  rownames_to_column(var ="ENSEMBL") %>% left_join(ids) %>%
  filter(!is.na(ENTREZID), !is.na(stat)) %>%
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat)) 


# res %>% as.data.frame() %>% 
#   rownames_to_column(var ="ENSEMBL") %>% left_join(ids) %>%
#   filter(!is.na(ENTREZID), !is.na(stat)) %>%
#   filter(SYMBOL %in% c("FA2H", "CERS2", "CERS6",  "UGT8", "GAL3ST1"))


# gene_name <- bitr(rownames(counts), fromType="ENSEMBL",
#                   toType="GENENAME",
#                   OrgDb="org.Hs.eg.db") #NA values drop by default
# res %>%  as.data.frame() %>% arrange(padj) %>% head(n=5) %>% 
#   rownames_to_column(var ="ENSEMBL") %>% left_join(gene_name)


nrow(res2)
ranks <- res2$stat
names(ranks) <- res2$SYMBOL



### downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
# pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/mm_data_240814/mm_data_240814/data/mh.all.v2024.1.Mm.entrez.gmt")   # for mouse
# head(pathways.hallmark)
# pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/h.all.v2024.1.Hs.entrez.gmt") # for human

# pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/h.all.v2024.1.Hs.symbols.gmt") # for hallmark pathways
# pathways.reactome <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/c2.cp.reactome.v2024.1.Hs.symbols.gmt") # for reactome pathways

# antigen_assembly <- pathways.reactome[ str_detect(names(pathways.reactome), "ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY")]
#Running fgsea algorithm:

fgseaRes_hallmark <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)
fgseaRes_reactome <- fgseaMultilevel(pathways=pathways.reactome, stats=ranks)
# fgseaRes <- fgseaMultilevel(pathways=antigen_assembly, stats=ranks)
# Tidy the results:
fgseaResTidy_hallmark <- fgseaRes_hallmark %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
fgseaResTidy_reactome <- fgseaRes_reactome %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)


fgseaResTidy_table_hallmark <- fgseaResTidy_hallmark %>%
  group_by(pathway) %>%
  mutate(leadingEdge = paste0(unlist(leadingEdge), collapse = ","))
fgseaResTidy_table_reactome <- fgseaResTidy_reactome %>%
  group_by(pathway) %>%
  mutate(leadingEdge = paste0(unlist(leadingEdge), collapse = ","))


fgseaResTidy_table_reactome[str_detect(fgseaResTidy_table_reactome$pathway, "MHC"),]


# wb_hallmark <- createWorkbook("BE02_hallmark")
# wb_reactome <- createWorkbook("BE02_reactome")
## Add a worksheets
addWorksheet(wb_hallmark, sheet_name, gridLines = FALSE)

## write data to worksheet 1
writeData(wb_hallmark, sheet = sheet_name, fgseaResTidy_table_hallmark, rowNames = TRUE)

addWorksheet(wb_reactome, sheet_name, gridLines = FALSE)

## write data to worksheet 1
writeData(wb_reactome, sheet = sheet_name, fgseaResTidy_table_reactome, rowNames = TRUE)



saveWorkbook(wb_hallmark, "data/BE02_hallmark.xlsx", overwrite = TRUE)
saveWorkbook(wb_reactome, "data/BE02_reactome.xlsx", overwrite = TRUE)


# fgseaResTidy$leadingEdge
# tail(fgseaResTidy)
# To see what genes are in each of these pathways:
# gene.in.pathway <- pathways.hallmark %>% 
#   enframe("pathway", "entrezgene_id") %>% 
#   unnest(cols = c(entrezgene_id)) %>% 
#   mutate(entrezgene_id = as.numeric(entrezgene_id)) %>%
#   inner_join( res2, by="entrezgene_id")

fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy %>% 
  filter(padj <= 0.05 & abs(NES) >= 2)  %>%       
  mutate(pathway = str_remove(pathway, "HALLMARK_")), aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Hallmark Pathways", y="Normalized Enrichment Score",
       title=paste0(plot_name) )
# ggsave(paste0("docs/",str_remove(comparison, "Grade_"),"_hallmark.png"), width = 4.5, height = 7)
# dir()
# 
# plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes, 
#               gseaParam=0.5)
# 
# # "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
# 
# 
# plotEnrichment(pathways.hallmark[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]],
#                ranks) + labs(title="wnt signaling")
# 
# 
# 
# ### for Dr Shen
# rownames(counts) = unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))
# 
# converted <- getBM(attributes=c('ensembl_gene_id','mgi_symbol' ), filters = 'ensembl_gene_id',
#                    values = rownames(counts), mart = ensembl)
# 
# converted <- converted %>% rename( "gene_symbol" = "mgi_symbol")
# 
# counts_name <- converted %>% left_join( counts %>% rownames_to_column("ensembl_gene_id") ) 
# write.xlsx(counts_name, "data/gene_counts.xlsx")                          
# 
# 


# Git clone <url from github>
#   Copy over your files to the new folder created from the command
# Then do git add . To add all the files to the GitHub repo
# Then git commit -m “<your message here>”
# Git push -u origi
# 
# Git push -u origin
# Git reset —hard HEAD
# Git reset —hard <commit_id



#mlist <- list()
# i <- 1
# names_list <- c()
# names_list <- names_list[-c(18,19)]

names_list <- c(names_list, plot_name)
mlist[[i]]<-fgseaRes
i <- i+1


names_list_new <- str_replace(str_replace(names_list, "Vs|vs", "Vs\n"), "\\(", "\n\\(")

names(mlist)<- names_list_new
mresult<-merge_result(mlist)

cluster_list <- c(
  "Normal-LGD Vs\n. Normal-HGD",
  "Normal-HGD Vs\n Normal-PDAC",
  
  "Normal LGD Vs\n Isolated LGD",
  "Normal LGD Vs\n Progressor LGD w HGD"  ,  
  "Normal HGD Vs\n Progressor LGD w HGD",
  "Normal-PDAC Vs\n progressor LGD w PDAC",
  
  "isolated LGD Vs\n. progressor LGD w HGD",
  "isolated LGD Vs\n isolated HGD ",
  "Progressor LGD w HGD Vs\n Isolated HGD",
  "Progressor LGD w HGD Vs\n Progressor LGD w PDAC",
  "progressor LGD w HGD Vs\n progressor HGD w PDAC" ,
  "Progressor LGD w PDAC Vs\n Progressor HGD w PDAC",
  "Normal-PDAC Vs\n progressor HGD w PDAC" ,
  
  "isolated HGD Vs\n progressor HGD w PDAC",
  "Normal HGD Vs\n Isolated HGD",
  "Normal PDAC Vs\n PDAC",
  "Progressor LGD w PDAC Vs\n PDAC" ,
  "Progressor HGD w PDAC Vs\n PDAC"
  )

cluster_list <- c(

  # 
  # "Normal-LGD Vs\n. Normal-HGD",
  # "Normal-HGD Vs\n Normal-PDAC",
  # "Normal LGD Vs\n Isolated LGD",
  # "Normal LGD Vs\n Progressor LGD w HGD"  ,
  # "Normal HGD Vs\n Progressor LGD w HGD",
  # "Normal-PDAC Vs\n progressor LGD w PDAC",
  # 
   "Isolated LGD Vs\n. Progressor LGD w HGD",
   "isolated LGD Vs\n. progressor LGD all \n(combined HGD or PDAC)",
   
   "Isolated LGD Vs\n Isolated HGD ",
   "Progressor LGD w HGD Vs\n Isolated HGD",
  # "Progressor LGD w HGD Vs\n Progressor LGD w PDAC",
  # "progressor LGD w HGD Vs\n progressor HGD w PDAC" ,
  # 
  # "Normal PDAC Vs\n PDAC",
  # "Progressor LGD w PDAC Vs\n PDAC" ,
  # 
  # "Normal-PDAC Vs\n progressor LGD w PDAC",
  # "Normal-PDAC Vs\n progressor HGD w PDAC" ,
  # "Progressor LGD w PDAC Vs\n Progressor HGD w PDAC",
  # "Progressor HGD w PDAC Vs\n PDAC"
  

  "isolated HGD Vs\n progressor HGD w PDAC"
  # "Normal HGD Vs\n Isolated HGD",
  # "Normal PDAC Vs\n PDAC",
  # 

)

# hallmark_subcat <-read_excel("data/BE02_case_list.xlsx", sheet = "hallmark")

ggplot(mresult %>% as.data.frame() %>% 
         mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%
         left_join(hallmark_subcat %>% dplyr::select(Process_category, Description), by= c(pathway = "Process_category")) %>%
         filter(Cluster %in% cluster_list) %>% 
         filter(padj <= 0.05) %>% 
         mutate(Cluster = factor(Cluster, levels = cluster_list)) %>%
         arrange( NES, Description )%>% 
         mutate(
                
                pathway = factor(pathway, levels = unique(pathway)),
                color = ifelse(NES < 0,"negative", "positive") ,
                color = factor(color, levels= c("positive","negative" ) )) , 
        aes(Cluster, pathway , size = abs(NES), color = color)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, size= 8,vjust = 1, hjust=1)) +
  facet_wrap(vars(Description), scales = "free_y",nrow = 4)



heatmap_df <- mresult %>% as.data.frame() %>% filter(padj <= 0.05) %>% 
  mutate(pathway = str_remove(pathway, "HALLMARK_")) %>%  filter(Cluster %in% cluster_list) %>% 
   dplyr::select(pathway,NES, Cluster) %>% 
  pivot_wider(id_cols = pathway,names_from = Cluster, values_from = NES, values_fill =0)

heatmap_df <- column_to_rownames(heatmap_df, 'pathway') 

# heatmap(as.matrix(heatmap_df[-1]), symmetric = TRUE ,scale ="none")

pheatmap(as.matrix(heatmap_df[-1]))





#### immune check point
checkpoint_genes <- c("KRT19", "MUC1", "CEACAM6", "MUC5AC", "CLDN18", "CDH1")  ### biomarker for the HGD and LGD
checkpoint_genes <- c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1", "VTCN1")

immune_genes <- c(
  # Checkpoints
  # "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1", "VTCN1"
  
  # # MHC & antigen presentation
  # "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
  # "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "CALR"
  
  # Lipid biosynthesis
  "FA2H", "CERS2", "CERS6",  "UGT8", "GAL3ST1"
  
  
)
ids_IC<-bitr(checkpoint_genes, fromType="SYMBOL",
             toType="ENSEMBL",
             OrgDb="org.Hs.eg.db") #NA values drop by default
# ids_IC <- ids_IC[-2,]
condition_subset <- condition[condition$Sample_ID %in% colnames(counts),]
counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts), colData = condition_subset, design = ~ 1)
dds_cp <- estimateSizeFactors(counts.DEseq)
norm_counts <- counts(dds_cp, normalized = TRUE)
rownames(norm_counts) = unlist(lapply(strsplit(rownames(norm_counts), '\\.'), '[[',1))
checkpoint_expr <- norm_counts[rownames(norm_counts) %in% ids_IC$ENSEMBL, ]
unique(condition_subset$Category)
condition_subset <- condition_subset %>%
      mutate(Category = factor(Category , levels = c("LGD",  "HGD" , "PDAC")),
        Disease = factor(Disease, levels = c("Normal-LGD" ,"Normal-HGD", "Normal-PDAC","Isolated LGD",
                                     "Progressor LGD w HGD" , 
                                     "Isolated HGD" ,
                                     "Progressor LGD w PDAC",
                                    "Progressor HGD w PDAC", "PDAC"))) %>% 
  
                        arrange(factor(Category),Path_ID,factor(Disease))

heatmap_df <-checkpoint_expr %>% as.data.frame() %>% rownames_to_column("ENSEMBL")  %>% left_join( ids_IC ) %>% 
           mutate(#ENSEMBL = paste0(ENSEMBL, "_", SYMBOL)
                  ENSEMBL = SYMBOL)%>% 
            dplyr::select(-SYMBOL) %>%
            column_to_rownames("ENSEMBL") %>%
          dplyr::select(condition_subset$Sample_ID)

ggplot(t(heatmap_df) %>% as.data.frame() %>%
         rownames_to_column("Sample_ID") %>%
         pivot_longer(!Sample_ID,values_to = "counts", names_to = "gene") %>%
         left_join(condition_subset) ,
       aes(Disease, counts ) ) +
  geom_boxplot() + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, size= 8,vjust = 1, hjust=1)) +
  facet_wrap(~gene, scales = "free_y", nrow = 2, strip.position = "top")
  


colnames(heatmap_df) <- paste0(condition_subset$Sample_ID, "_", condition_subset$Path_ID, "_",condition_subset$Disease)

pheatmap::pheatmap(heatmap_df  , scale = "column")




dds_IC <- estimateSizeFactors(counts.DEseq)
norm_counts <- counts(dds_IC, normalized = TRUE)


checkpoint_expr <- norm_counts[rownames(norm_counts) %in% ids_IC$ENSEMBL, ]

rownames(checkpoint_expr) <- ids_IC$SYMBOL[match(rownames(checkpoint_expr), ids_IC$ENSEMBL)]



pheatmap(checkpoint_expr, scale = "row", show_rownames = TRUE)


msigdbr_species <- "Homo sapiens"
m_df <- msigdbr::msigdbr(species = msigdbr_species, category = "C2", subcategory = "CP:REACTOME")
m_df <- m_df[str_detect(m_df$gs_name, "ANTIGEN"), ]
gene_sets <- split(m_df$ensembl_gene, m_df$gs_name)

fgsea_res <- fgsea(
  pathways = gene_sets,
  stats    = ranks,
  eps      = 0.0,  # Avoid p=0
  minSize  = 10,   # Min genes in pathway
  maxSize  = 500   # Max genes in pathway
)

# Plot enrichment plot
plotEnrichment(gene_sets[[1]], ranks) +
  labs(title = "GSEA: Antigen Presentation Pathway")





#### boxplot with p-value

#str(heatmap_df)
# str(condition_subset) 

df_plot <- t(heatmap_df) %>% as.data.frame() %>% rownames_to_column("Sample_ID") %>%
      left_join(condition_subset) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  # filter(Disease %in% c( "Isolated LGD", "Isolated HGD"))
  filter(!Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
  #filter(T.cells.CD8 < 0.3) %>%
  # filter(Disease %in% c( "Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC" )) %>%
  droplevels() %>%
  mutate(
    Disease = factor(Disease, levels = c("Isolated LGD" ,
                                         "Progressor LGD w HGD" ,
                                         "Isolated HGD",
                                         "Progressor LGD w PDAC",
                                         "Progressor HGD w PDAC", 
                                         "PDAC" )),
    # Disease = str_wrap(Disease, 10),
    # Disease =factor(Disease, levels= c(str_wrap(levels(.$Disease), 10))),
    groups = ifelse(str_detect(Disease, "PDAC"), "PDAC", "non_PDAC"),
    groups = factor(groups, levels = c( "non_PDAC", "PDAC"))
    
    # Disease = ifelse(str_detect(Disease, "Progressor"),"Progressor LGD", "Isolated LGD") ,
    #      Disease = factor(Disease , levels = c("Isolated LGD", "Progressor LGD"))
  )

df_plot <- df_plot %>% mutate(
  Disease = str_wrap(Disease, 10),
  Disease =factor(Disease, levels= c(str_wrap(levels(df_plot$Disease), 10)))
)

colnames(df_plot)
p_values <- df_plot %>% 
  wilcox_test( MUC1 ~ groups ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_values  <- p_values %>% add_xy_position(x = "Disease", group = "groups")
p_values$xmin <- 2
p_values$xmax <- 5
# p_values$y.position <- p_values$y.position +15
p_values$p.adj.signif[p_values$p.adj.signif =="ns"] <- ''

st_line <- data.frame(xmin = c(1,4),
                      xmax = c(3,6),
                      y.position = c(p_values$y.position-3, p_values$y.position-3),
                      p="",
                      group1=1,
                      group2 =2)

color <- c("#2F56B5",#115C80", 
             "#2295C7", #309898" ,#479685", 
             "#00A86B", #FF9F00", 
          "#FCC92A", #FF9F00", 
             "#F4691E", 
            "#AD0000")

# color <- c("#48217390","#433E8590", "#38598C90","#2BB07F90","#85D54A90", "#FDE725FF") ##2BB07F50")
p <- ggboxplot(df_plot , x = "Disease", y = "MUC1",
               # color = "Disease",
               fill = "Disease",palette = color,
               add = "jitter",
               legend.title ="")

p + #stat_compare_means(label.y = 0.35 ) +
  labs(x="") + # ylim(c(0, 0.14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) +
  stat_pvalue_manual(
    p_values,  label = "{p.adj}{p.adj.signif}" #, tip.length = 0.02
    #step.increase = 0.05
  ) +
  stat_pvalue_manual(
    st_line, tip.length = 0
    # step.increase = 0.05
  ) 





#### ensemble id to posiition
library(biomaRt)

# Select the Ensembl human gene dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ensemble <- unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))
# Get gene information
results <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = gene_ensemble , # Replace with your Ensembl gene ID
                 mart = ensembl)

print(results)

# counts
# condition

sorted_counts <-counts %>% rownames_to_column("ENSEMBLID") %>%
  mutate(ENSEMBL = unlist(lapply(strsplit(ENSEMBLID, '\\.'), '[[',1)) ) %>% 
  left_join(results %>%
              filter(!str_detect(chromosome_name, "GL|KI|MT" )) %>%
              mutate(chromosome_name = paste0("chr", chromosome_name),
                     chromosome_name = factor(chromosome_name , levels =paste0("chr", c(seq(1:22), "X", "Y"))) ) %>%
              arrange(chromosome_name, start_position) %>%
              group_by(chromosome_name) %>%
              mutate(xaxis = seq(1, length(start_position)),
                     start_position = start_position,
                     ensembl_gene_id = ensembl_gene_id),
               by = c("ENSEMBL" = "ensembl_gene_id")) %>%
  filter(!is.na(chromosome_name)) %>%
  # filter(!str_detect(chromosome_name, "GL|KI|MT" )) %>%
  # mutate(chromosome_name = paste0("chr", chromosome_name),
  #        chromosome_name = factor(chromosome_name , levels =paste0("chr", c(seq(1:22), "X", "Y"))) ) %>%
  arrange(chromosome_name, start_position) %>% 
   mutate(x = seq(1, nrow(.))) 
  
 # dplyr::select(-c("chromosome_name","start_position", "end_position" )) %>%
  # column_to_rownames("ENSEMBL")
 

   

v_line <- data.frame(v_line = sorted_counts$x[sorted_counts$xaxis==1])


for(col in colnames(sorted_counts)[70:75]){
  p <- sorted_counts %>% 
    mutate(col_name = sorted_counts[[col]] ) %>% 
     filter(col_name > 0) %>%
ggplot( aes(x = x, y = log(col_name,2), color = as.character(chromosome_name), group = chromosome_name) )+
  geom_point()  +
 # geom_vline(aes(xintercept = v_line-1), v_line) +
  labs(title = col) +
  theme(legend.position = "none") 
  
  print(p)
  
  
}



cosine_sim_matrix_vector <- function(mat, vec) {
  # Normalize matrix rows
  mat_norm <- mat / sqrt(rowSums(mat^2))
  
  # Normalize vector
  vec_norm <- vec / sqrt(sum(vec^2))
  
  # Compute cosine similarity
  sim <- mat_norm %*% vec_norm
  
  return(as.vector(sim))
}

# Example usage
if (interactive()) {
  mat <- matrix(c(1, 2, 3,
                  4, 5, 6,
                  7, 8, 9), nrow = 3, byrow = TRUE)
  
  vec <- c(1, 0, 1)
  
  result <- cosine_sim_matrix_vector(mat, vec)
  print(result)
}



library(ggplot2)

# Sample data and colors
x <- c("A", "B", "C", "D", "E")
y <- c("Var1", "Var2", "Var3", "Var4", "Var5")
data <- expand.grid(X = x, Y = y)
data$Z <- runif(25, 0, 5)
color_vector <- c("red", "blue", "green", "purple", "orange")

# Create the heatmap
ggplot(data, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.y = element_text(color = color_vector)) +
  labs(y = "Y-axis with Colors", x = "X-axis")
