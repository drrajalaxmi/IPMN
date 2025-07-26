

heatmap_df<- cibersort_pt %>% select(-c("P.value", "Correlation", "RMSE")) 
rownames(heatmap_df) <- paste0(heatmap_df$Patient_id, "_",heatmap_df$Grade)
heatmaply::heatmaply(as.matrix(heatmap_df[2:23]), Colv=FALSE)

#pairwise test

#line plot
grade_pivot <-cibersort_pt %>% select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  pivot_longer( !c("Mixture","Grade", "Patient_id"), names_to = "cells",values_to = "pcnt_cells") %>%
  pivot_wider(id_cols =!Mixture, names_from = Grade, values_from = pcnt_cells) #%>% arrange(cells) 

df_pValue <- data.frame(cells= unique(grade_pivot$cells) )

for ( cell in unique(grade_pivot$cells)){
  df_test <- grade_pivot %>% filter(cells==cell)
  
  df_pValue$LCD_ttest[df_pValue$cells==cell] <-t.test(df_test$normal, df_test$`low grade`, paired = FALSE)$p.value
  df_pValue$LCD_wilcox[df_pValue$cells==cell] <-wilcox.test(df_test$normal, df_test$`low grade`, paired = FALSE)$p.value

  df_pValue$HCD_ttest[df_pValue$cells==cell] <-t.test(df_test$normal, df_test$`high grade`, paired = FALSE)$p.value
  df_pValue$HCD_wilcox[df_pValue$cells==cell] <-wilcox.test(df_test$normal, df_test$`high grade`, paired = FALSE)$p.value
  
  df_pValue$INV_ttest[df_pValue$cells==cell] <-t.test(df_test$normal, df_test$`invasive`, paired = FALSE)$p.value
  df_pValue$INV_wilcox[df_pValue$cells==cell] <-wilcox.test(df_test$normal, df_test$`invasive`, paired = FALSE)$p.value
  
  
}

rownames(df_pValue) <- df_pValue$cells
# df_pValue[df_pValue > 0.05] <- 1
heatmap(as.matrix(df_pValue[2:7]))
heatmaply::heatmaply(as.matrix(df_pValue[2:7]), Colv=FALSE)


long_pivot <- cibersort_pt %>% select(-c("P.value", "Correlation", "RMSE")) %>% #filter(Patient_id== cibersort_pt$Patient_id[1] ) %>%
  # select(colnames(.)[1:(ncol(.)-2)]) %>% 
  pivot_longer( !c("Mixture","Grade", "Patient_id"), names_to = "cells",values_to = "pcnt_cells")
for ( cell in unique(grade_pivot$cells)){
  val <- long_pivot$pcnt_cells[long_pivot$cells==cell]
  grp <- long_pivot$Grade[long_pivot$cells==cell]
print(pairwise.t.test(val, c(grp, , p.adjust.method = "BH"))
  
}

print(pairwise.t.test(long_pivot$pcnt_cells, paste0(long_pivot$Grade, long_pivot$cells) , p.adjust.method = "BH"))

