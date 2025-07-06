
## neoantigen output from the scannneo2


#setting a thresholds for the quality neoantigen
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


####combining predicted neoantigen from individual sample are merged into a matrix

# AntiGen_flt <- data.frame()
# file_list <- list.files("data/")[grep(list.files("data/"), pattern="BE02_4[0-9]A_mhc-I_neoepitopes_all.txt")]
# for (file in file_list){
#   print(file)
#   
#   AntiGen <- read.table(paste0("data/",file), header = TRUE, sep ="\t") %>% 
#     mutate(TPM = as.numeric(TPM),
#            agretopicity = as.numeric(agretopicity),
#            mt_immunogenicity = as.numeric(mt_immunogenicity),
#            sample_ID = unlist(str_split(file, "_m"))[1]
#     )
#   
#   AntiGen_flt <- rbind( AntiGen_flt, AntiGen)
#   
#   
# }


## saved neoantigen matrix file neoAnti_norm_ag_ic50.flt.tsv with filters ag < 0.3 and ic50 < 100
AntiGen_flt <- read.table('data/neoAnti_norm_ag_ic50_dp.flt.tsv', sep ="\t", header = TRUE)

length(unique(AntiGen_flt$sample_ID))
AntiGen_flt$TPM_count <- AntiGen_flt$TPM
AntiGen_flt$TPM[str_detect(AntiGen_flt$TPM, "\\|")] <- -1

AntiGen_flt <- AntiGen_flt %>% 
  mutate(TPM = as.numeric(TPM))


#filtering for the quality neoantigen
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

ggplot(neoantigen_data, aes(x = agretopicity, y = mt_immunogenicity)) +
  geom_point(aes(color = Composite_score), size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Aggretopicity vs. Immunogenicity")

ggplot(neoantigen_data %>% filter(agretopicity > -1), aes(x = agretopicity)) +
  geom_histogram(bins = 20, fill = "purple", alpha = 0.7) +
  labs(title = "Aggretopicity Distribution in Filtered Neoantigens")




### 1. Neoantigen Burden Analysis
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

color <- c( "#2F56B5", "#2295C7", "#00A86B", "#FCC92A", "#F4691E",  "#AD0000")
p <- ggboxplot(burden , x = "Disease", y  = "log(Neoantigen_count)",
               # color = "Disease",
               fill = "Disease",palette = color,
               add = "jitter",
               legend.title ="") +
  labs(#title = "Neoantigen Burden Across IPMN Stages",
    y = "Neoantigens per Sample (log)")
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


### 2. Binding affinity comparison
 ggplot(neoantigen_data %>%
          mutate(
            Disease = str_wrap(Disease, 10),
            Disease =factor(Disease, levels= c(str_wrap(levels(.$Disease), 10)))), aes(x = factor(Disease), y = mt_epitope_ic50, fill = Disease)) +
   geom_violin(alpha = 0.7) +
   #geom_jitter() +
   geom_boxplot(width = 0.1, fill = "white") +
   labs(#title = "Neoantigen Binding Affinity by Grade",
        y = "Predicted Binding Affinity (nM)",
        x = "") +
   scale_fill_manual(values = c(color)) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

### 3. Immunogenicity comparison
 neoantigen_data %>%
   mutate(
     Disease = str_wrap(Disease, 10),
     Disease =factor(Disease, levels= c(str_wrap(levels(.$Disease), 10)))) %>%
   ggplot(aes(x = Disease, y = mt_immunogenicity, fill = Disease)) +
   geom_boxplot() +
   labs(#title = "Neoantigen Immunogenicity by Grade",
        y = "Predicted Immunogenicity Score",
        x = "") +
   scale_fill_manual(values = c(color)) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
 