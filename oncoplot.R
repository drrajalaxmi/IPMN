# BiocManager::install("maftools", force = TRUE)
library(maftools)
library(readxl)

###example
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml_ex = read.maf(maf = laml.maf, clinicalData = laml.clin)



BE02_case_list <-read_excel("data/BE02_case_list.xlsx", sheet = "FINAL")
column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE02_",BE02_case_list$Sample_ID, "A")
sample_list_not <- c("BE02_18A",  "BE02_59A",  "BE02_60A",  "BE02_61A",  "BE02_74A" , "BE02_75A",  "BE02_101A", "BE02_102A", "BE02_104A", "BE02_105A", "BE02_113A", "BE02_114A")
condition <- BE02_case_list %>% 
  filter(!Sample_ID %in% sample_list_not)
condition$Tumor_Sample_Barcode <- condition$Sample_ID
condition$Category[condition$Category=="HGD"] <- paste0("2 HGD")
condition$Category[condition$Category=="LGD"] <- paste0("1 LGD")
condition$Category[condition$Category=="PDAC"] <- paste0("3 PDAC")
condition$Grade[condition$Grade=="normal"] <- paste0("1 Normal")
condition$Grade[condition$Grade=="low grade"] <- paste0("2 low grade")
condition$Grade[condition$Grade=="high grade"] <- paste0("3 high grade")
condition$Grade[condition$Grade=="invasive"] <- paste0("4 invasive")


write.table(condition, "data/maf_sample.tsv", sep ="\t")
concat_barcode_pass.maf
cat_PDAC_concat_barcode_pass_p.maf
laml = read.maf(maf = 'data/concat_barcode_pass.maf.gz', clinicalData = 'data/maf_sample.tsv')
laml

# str(laml)

oncoplot(maf = laml, top = 100, clinicalFeatures = 'Category', sortByAnnotation = T, sortByMutation = T,showTumorSampleBarcodes = TRUE, titleText = "All Samples")

df <-lollipopPlot(
  maf = laml,
  gene = 'GNAS',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  printCount = TRUE
)

df$gene <- "GNAS"
df$category <- "All"

table <- rbind(df, table)

chip <- table

library(openxlsx)
write.xlsx(chip , "data/chip_list.xlsx")



rainfallPlot(maf = laml, detectChangePoints = TRUE, pointSize = 0.4)
pws = pathways(maf = laml, plotType = 'treemap')
plotPathways(maf = laml, pathlist = pws, showTumorSampleBarcodes =TRUE)

plotVaf(maf = laml)


getSampleSummary(laml)
str(laml)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml = read.maf(maf = laml.maf)
laml

laml$Tumor_Sample_Barcode
varient <- read.table("data/cosmic_all_varient.txt", sep = "\t")
head(varient)
colnames(varient) = c("chrom", "Position", "ID", "REF", "ALT")
str(condition)
varient_df <- varient %>% separate(chrom , c("path", "chrom"), sep =":") %>% separate(path, sep= "\\/", c(NA, NA, NA,"file")) %>%
 separate(file, sep= "A_cosmic_", c("Sample_ID", NA)) %>% 
     mutate(Sample_ID = str_replace(Sample_ID,"BE02_", "BE")) %>% 
     left_join(BE02_case_list %>%
                  select(Sample_ID, Path_ID,Grade, Path_ID )) %>% 
      group_by( Path_ID, Sample_ID,chrom, Position) %>%
       summarise(Sample_count = length(Sample_ID)) %>% arrange(desc(Sample_count)) %>% ungroup()

as.data.frame(varient_df)

varient_count <-varient_df %>% 
  mutate(Sample_ID = paste0(Path_ID,"_",Sample_ID)) %>% select(-c("Path_ID")) %>%
  pivot_wider(names_from = "Sample_ID", values_from = "Sample_count", names_sort=TRUE) %>% as.data.frame() %>%
  arrange(chrom, Position)
 varient_count[is.na(varient_count) ] <- 0
rownames(varient_count) <- paste0(varient_count$chrom,"_", varient_count$Position) 
varient_count_mtx <- varient_count %>% select(-c("chrom","Position")) %>% as.matrix()

heatmap(varient_count_mtx, Rowv=NA, Colv = "Rowv", scale = NULL)

