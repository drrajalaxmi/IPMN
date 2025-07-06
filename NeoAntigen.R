# Neoantigen Comparison Across IPMN Stages
# Requires: neoantigen data with columns: Sample_ID, Neoantigen_sequence, Stage, Binding_affinity, etc.
library(readxl)
library(tidyverse)
library(ggpubr)
library(UpSetR)
library(ComplexHeatmap)
library(dplyr)

# AntiGen <- read.table("data/BE02_67A_mhc-I_neoepitopes_all.txt", header = TRUE) %>% 
#            mutate(TPM = as.numeric(TPM),
#                   agretopicity = as.numeric(agretopicity),
#                   mt_immunogenicity = as.numeric(mt_immunogenicity)
#                   )



BE02_case_list <-read_excel("data/BE02_case_list.xlsx", sheet = "FINAL")
column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE02_",BE02_case_list$Sample_ID, "A")
sample_list_not <- c("BE02_18A",  "BE02_59A",  "BE02_60A",  "BE02_61A",  "BE02_74A" , "BE02_75A",  "BE02_101A", "BE02_102A", "BE02_104A", "BE02_105A", "BE02_113A", "BE02_114A")
condition <- BE02_case_list %>%
  filter(!Sample_ID %in% sample_list_not) %>%
  # filter(Sample_ID %in% colnames(df_count)) %>%
  mutate(Category= factor(Category, levels=c("LGD", "HGD", "PDAC")),
         Grade = factor( Grade, levels = c("normal", "low grade", "high grade", "invasive")),
         Disease = factor(Disease , levels = c("Normal-LGD", "Normal-HGD",
                                               "Normal-PDAC", "Isolated LGD",
                                               "Isolated HGD", "Progressor LGD w HGD",
                                               "Progressor LGD w PDAC",
                                               "Progressor HGD w PDAC",
                                               "PDAC") )) %>%
          arrange( Disease, Category)
dim(condition)
thresholds <- list( 
  # mt_immunogenicity = 0.5,     # Keep top 50%  (Higher = more immunogenic)
  self.similarity = 0.5,       # Exclude peptides similar to self (lower = safer) 0.3
  # pathogen_similarity = 0.6,   # Prefer peptides resembling pathogens (Higher = microbial mimicry) 
  proteome_similarity = 0.3,   # Exclude peptides matching human proteome (Lower = more unique to tumor)
  vaf = 0.2,                   # Minimum clonality (Higher = clonal mutation)
  TPM = 10,                    # Minimum expression (better expression)
  agretopicity = 0.3 ,          # ≥ 0.7 = High (prioritize) , 0.4–0.7 = Moderate ,< 0.4 = Low (exclude) 
  dp = 10 ,                     ## AF , Depth (DP = 5 )  at least 2 reads 
  ic50 = 100
 ) 

AntiGen_flt <- data.frame()


file_list <- list.files("data/")[grep(list.files("data/"), pattern="BE02_4[0-9]A_mhc-I_neoepitopes_all.txt")]

for (file in file_list){
  print(file)
  
  AntiGen <- read.table(paste0("data/",file), header = TRUE, sep ="\t") %>% 
    mutate(TPM = as.numeric(TPM),
           agretopicity = as.numeric(agretopicity),
           mt_immunogenicity = as.numeric(mt_immunogenicity),
           sample_ID = unlist(str_split(file, "_m"))[1]
    )
  
  AntiGen_flt <- rbind( AntiGen_flt, AntiGen)
  
  
}
AntiGen_flt_old <- AntiGen_flt


#######
## neoAnti_norm_ag_ic50.flt.tsv   ag < 0.3 and ic50 < 100
AntiGen_flt <- read.table('data/neoAnti_norm_ag_ic50_dp.flt.tsv', sep ="\t", header = TRUE)
# dim(AntiGen_flt)
length(unique(AntiGen_flt$sample_ID))
AntiGen_flt$TPM_count <- AntiGen_flt$TPM
AntiGen_flt$TPM[str_detect(AntiGen_flt$TPM, "\\|")] <- -1

AntiGen_flt <- AntiGen_flt %>% 
                 mutate(TPM = as.numeric(TPM))

# str(AntiGen_flt)
# unique(AntiGen_flt$chrom[AntiGen_flt$PTC_exon_number!="chr1"])
######


neo_flt <- AntiGen_flt %>%  
     # filter(mt_epitope_ic50 < 100) #%>%
     filter(
            # mt_immunogenicity >= thresholds$mt_immunogenicity,
             self.similarity <= thresholds$self.similarity,   ## 0.5
           # ( pathogen_similarity >= thresholds$pathogen_similarity | pathogen =="." ),
             proteome_similarity <= thresholds$proteome_similarity ,   ### 0.3
            ( vaf >= thresholds$vaf | vaf == -1 ),   ## 0.2
            ( TPM >= thresholds$TPM | TPM == -1 ),    ### 10
            # (agretopicity >= thresholds$agretopicity | agretopicity == -1)    ## < 0.3
            ( dp >= thresholds$dp | is.na(dp)) ## 10
           ) %>%
     mutate(
      Composite_score = (
        0.3 * mt_immunogenicity +
        0.2 * agretopicity +
        0.2 * pathogen_similarity +
        0.15 * (1 - self.similarity) +
        0.1 * (1 - proteome_similarity) +
        0.05 * vaf +               # Clonality
        0.05 * log10(TPM + 1)      # Log-transformed expression
             )
        ) %>%
     arrange(desc(Composite_score))


   

neoantigen_data <- neo_flt %>% #left_join(condition, by = c(sample_ID = "Sample_ID") )
       mutate(Grade = factor(Grade, levels = c("low grade", "high grade", "invasive")),
       Disease = factor( Disease , levels = c("Isolated LGD" , "Progressor LGD w HGD",
                                              "Isolated HGD" , "Progressor LGD w PDAC",
                                              "Progressor HGD w PDAC", "PDAC")))
# unique(AntiGen_flt$sample_ID)[!unique(AntiGen_flt$sample_ID) %in% unique(neoantigen_data$sample_ID)]

ggplot(neoantigen_data, aes(x = agretopicity, y = mt_immunogenicity)) +
  geom_point(aes(color = Composite_score), size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Aggretopicity vs. Immunogenicity")

ggplot(neoantigen_data %>% filter(agretopicity > -1), aes(x = agretopicity)) +
  geom_histogram(bins = 20, fill = "purple", alpha = 0.7) +
  labs(title = "Aggretopicity Distribution in Filtered Neoantigens")




### 2. Neoantigen Burden Analysis ------------------------
# Calculate burden (neoantigens per sample)
unique(neoantigen_data$Disease)
burden <- neoantigen_data %>%
  mutate(Grade = factor(Grade, levels = c("low grade", "high grade", "invasive")),
         Disease = factor( Disease , levels = c("Isolated LGD" , "Progressor LGD w HGD",
                                                "Isolated HGD" ,"Progressor LGD w PDAC",
                                                "Progressor HGD w PDAC", "PDAC")),
         Disease = str_wrap(Disease, 10),
         Disease =factor(Disease, levels= c(str_wrap(levels(cibersort_pt$Disease), 10)))) %>%
  group_by(sample_ID, Disease) %>%
  summarise(Neoantigen_count = n_distinct(mt_epitope_seq),
            .groups = 'drop')

# Statistical comparison
burden_stats <- compare_means(
  Neoantigen_count ~ Disease, 
  data = burden,
  method = "kruskal.test"
)

# Visualization
ggplot(burden, aes(x = Disease, y = log(Neoantigen_count), fill = Disease)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  stat_compare_means(label = "p.format", method = "kruskal.test") +
  labs(title = "Neoantigen Burden Across IPMN Stages",
       y = "Neoantigens per Sample (log)",
       x = "Dysplasia Stage") +
 # scale_fill_manual(values = c("normal" = "gray","low grade" = "#4DBBD5", "high grade" = "#E64B35", "invasive" = "#00A087")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

### adding with ggpubr
p_values <- burden %>% 
  mutate(groups = ifelse(str_detect(Disease, "PDAC"), "PDAC", "non_PDAC"),
         groups = factor(groups, levels = c( "non_PDAC", "PDAC"))) %>%
  t_test( Neoantigen_count ~ groups ) %>%
  adjust_pvalue( method = "bonferroni" ) %>%
  add_significance("p.adj")

p_values  <- p_values %>% add_xy_position(x = "groups", group = "groups")
p_values$xmin <- 2
p_values$xmax <- 5
p_values$y.position <- 11

st_line <- data.frame(xmin = c(1,4),
                      xmax = c(3,6),
                      y.position = c(p_values$y.position-0.5,p_values$y.position-0.5),
                      p="",
                      group1=1,
                      group2 =2)

color <- c("#48217390","#433E8590", "#38598C90","#2BB07F90","#85D54A90", "#FDE725FF") ##2BB07F50")
p <- ggboxplot(burden , x = "Disease", y  = "log(Neoantigen_count)",
               # color = "Disease",
               fill = "Disease",palette = color,
               add = "jitter",
               legend.title ="") +
  labs(#title = "Neoantigen Burden Across IPMN Stages",
       y = "Neoantigens per Sample (log)")
#  Add p-value

# Change method
# p + stat_compare_means(# comparisons = list(c("Isolated LGD", "Isolated HGD") )
#                       
#               #list(c(levels(df_plot$Disease)[c(4,5)]),c(levels(df_plot$Disease)[c(4,6)] ))
#                  ) +
p+ #stat_compare_means(label.y = 0.35 ) +
  labs(x="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1 )) +
  
  stat_pvalue_manual(
    p_values,  label = "{p.adj}" #{p.adj.signif}" #, tip.length = 0.02
    # step.increase = 0.05
  ) +
  stat_pvalue_manual(
    st_line, tip.length = 0
    # step.increase = 0.05
  )



### 3. Shared Neoantigen Analysis -----------------------
# Create binary presence/absence matrix
shared_matrix <- neoantigen_data %>%
  distinct(sample_ID, mt_epitope_seq, .keep_all = TRUE) %>%
  mutate(Present = 1) %>% dplyr::select(sample_ID, mt_epitope_seq, Present) %>%
  pivot_wider(names_from = mt_epitope_seq, values_from = Present, values_fill = 0)

shared_matrix_chrom <- neoantigen_data %>%
  distinct(chrom, start, mt_epitope_seq, .keep_all = TRUE) %>%
  mutate(Present = 1) %>% dplyr::select(sample_ID, mt_epitope_seq, Present) %>%
  pivot_wider(names_from = mt_epitope_seq, values_from = Present, values_fill = 0)

# UpSet plot for shared neoantigens
stage_lists <- neoantigen_data %>%
  filter(Category == "PDAC") %>%
  group_by(sample_ID) %>%
  summarise(Neoantigens = list(unique(mt_epitope_seq)))




upset(fromList(stage_lists$Neoantigens %>% setNames(stage_lists$sample_ID)),
      nsets = 20,
      order.by = "freq",
      mainbar.y.label = "Shared Neoantigens (PDAC)",
      sets.x.label = "Neoantigens per Stage")

### 4. Neoantigen Quality Comparison --------------------
# Binding affinity comparison
ggplot(neoantigen_data, aes(x = factor(Disease), y = mt_epitope_ic50, fill = Grade)) +
  geom_violin(alpha = 0.7) +
  #geom_jitter() +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Neoantigen Binding Affinity by Grade",
       y = "Predicted Binding Affinity (nM)",
       x = "") +
  #scale_fill_manual(values = c("normal" = "gray","low grade" = "#4DBBD5", "high grade" = "#E64B35", "invasive" = "#00A087")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  stat_compare_means(comparisons = list(c("low grade", "high grade"), 
                                        c("high grade", "PDAC"),
                                        c("low grade", "PDAC")))


# Immunogenicity comparison
neoantigen_data %>%
 # group_by(Disease, mt_epitope_seq) %>%
 # summarise(Mean_immunogenicity = mean(mt_immunogenicity), .groups = 'drop') %>%
  ggplot(aes(x = Disease, y = mt_immunogenicity, fill = Grade)) +
  geom_boxplot() +
  labs(title = "Neoantigen Immunogenicity by Grade",
       y = "Predicted Immunogenicity Score",
       x = "") +
 # stat_compare_means(method = "kruskal.test", label.y = 1.1) +
 # scale_fill_manual(values = c("normal" = "gray","low grade" = "#4DBBD5", "high grade" = "#E64B35", "invasive" = "#00A087")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

### 5. Stage-Specific Neoantigen Patterns ---------------
# Heatmap of top neoantigens (example for 50 neoantigens)
str(neoantigen_data)
top_neos <- neoantigen_data %>%
  #count( mt_epitope_seq, Grade) %>%
  group_by(mt_epitope_seq, Grade) %>%
  summarise(n = n()) %>%
  filter(n() > 1) %>% # Present in >1 stage
  arrange(desc(n)) %>%
  slice_head(n = 20) %>%
  pull(mt_epitope_seq)

heatmap_data <- neoantigen_data %>%
  filter(mt_epitope_seq %in% top_neos) %>%
  group_by(mt_epitope_seq, Grade) %>%
  summarise(Mean_immunogenicity = mean(mt_immunogenicity),
            .groups = 'drop') %>%
  pivot_wider(names_from = Grade, values_from = Mean_immunogenicity) %>%
  column_to_rownames("mt_epitope_seq")

heatmap(as.matrix(heatmap_data),
        name = "Immunogenicity",
        col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 8),
        column_names_rot = 45,
        show_row_names = FALSE,
        column_title = "Stage-Specific Neoantigen Immunogenicity")

### 6. Statistical Reporting ----------------------------
# Generate summary tables
neoantigen_summary <- neoantigen_data %>%
  group_by(Disease) %>%
  summarise(
    Median_count = median(n_distinct(mt_epitope_seq)),
    Mean_affinity = mean(mt_epitope_ic50),
    Mean_immunogenicity = mean(mt_immunogenicity),
    Shared_neos = n_distinct(mt_epitope_seq[duplicated(mt_epitope_seq)])
  )

# Save results
write_csv(neoantigen_summary, "neoantigen_stage_comparison.csv")







library(umap)
library(proxy)  # For distance matrices

# Create a binary matrix (samples x neoantigens)
binary_matrix <- table(neoantigen_data$sample_ID, neoantigen_data$mt_epitope_seq)
str(binary_matrix)
# Jaccard similarity
jaccard_dist <- proxy::dist(binary_matrix, method = "Jaccard")
hclust_result <- hclust(jaccard_dist)
plot(hclust_result, main = "Hierarchical Clustering by Neoantigen Similarity")

# UMAP (if high-dimensional)
umap_result <- umap(binary_matrix)
plot(umap_result$layout, col = as.factor(rownames(binary_matrix)), pch = 19)





