library(dplyr)

map_df <- read.table("data/PRI_2pass_pcnt_mapped.txt", col.names = c("sample", "pcnt_mapped"))
map_df$Sample_ID <- stringr::str_replace(stringr::str_replace(map_df$sample , "BE02_","BE" ), "A",'')

map_df <- map_df %>% arrange(pcnt_mapped) %>% mutate(Sample_ID = factor(Sample_ID, levels = .$Sample_ID))
(map_df <- map_df %>% select(Sample_ID, pcnt_mapped))
hist(map_df$pcnt_mapped, main = "Uniquely mapped reads %", xlab = "Uniquely mapped reads %")


boxplot(pcnt_mapped ~ Grade , data = BE02_case_list_pct )
