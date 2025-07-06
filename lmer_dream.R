### linear mixed effect model
# BiocManager::install("variancePartition")
library("variancePartition")
library("edgeR")
library("BiocParallel")
data(varPartDEdata)
dim(dge)
# filter genes by number of counts

isexpr <- rowSums(cpm(counts) > 0.1) >= 5

# Standard usage of limma/voom
dge <- DGEList(counts[isexpr, ])
dge <- calcNormFactors(dge)

# make this vignette faster by analyzing a subset of genes
dge <- dge[1:100, ]


metadata <- metadata %>% 
  mutate(Category= factor(Category , levels= c("LGD", "HGD", "PDAC")),
         Grade = factor(Grade, levels = c("normal", "low grade", "high grade", "invasive")),
         Disease = factor(Disease , levels = c( "Normal-LGD", "Normal-HGD", "Normal-PDAC",
                                                "Isolated LGD", "Progressor LGD w HGD" , "Isolated HGD",
                                                "Progressor LGD w PDAC", "Progressor HGD w PDAC","PDAC" )))


# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)


form <- ~ 0 + Category   * Grade + (1 | Path_ID)
design <- model.matrix(~ 0 + Category , metadata)

L <- makeContrastsDream(form, metadata,
                        contrasts = c(
                         compareHGD_LGD =  "CategoryHGD - CategoryLGD",
                         comparePADC_LGD ="CategoryPDAC - CategoryLGD",
                         comparePDAC_HGD = "CategoryPDAC - CategoryHGD"
                          
                          # comparelow_norm = "`Gradelow grade` - Gradenormal",
                          # comparehigh_norm = "`Gradehigh grade` - Gradenormal",
                          # compareinvasive_norm = "Gradeinvasive - Gradenormal",
                          # comparePDAC_HGD = "CategoryPDAC - CategoryHGD"
                          # compare1 = "`Gradelow grade:CategoryHGD` - `Gradehigh grade:CategoryHGD`",
                          # compare2 = "`Gradelow grade:CategoryPDAC` - `Gradehigh grade:CategoryPDAC`"
                        )
)

# Visualize contrast matrix
plotContrasts(L)
# fit dream model with contrasts
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)
fit <- dream(vobjDream, form, metadata, L)
fit <- eBayes(fit)

# get names of available coefficients and contrasts for testing
colnames(fit)
# extract results from first contrast
topTable(fit, coef = "compareHGD_LGD", number = 3)
str(fit)





# Load necessary packages
library(lme4)
library(ggplot2)
library(modelr)
library(dplyr)

# Load the sleepstudy data
data(sleepstudy)

# Create a fake "condition" variable for demonstration
# Recode Subject IDs into two groups
sleepstudy$condition[sleepstudy$Subject %in% c('308','309','310','330','331','332','333','334')] <- 'control'
sleepstudy$condition[!sleepstudy$Subject %in% c('308','309','310','330','331','332','333','334')] <- 'treatment'


# Fit an lmer model with Time, Condition, and their interaction,
# including random intercepts and slopes for Subject
model <- lmer(Reaction ~ 0  +condition +Days + (1 | ( Subject)), data = sleepstudy)
summary(model)
# Generate a grid of data to predict on
# Predictions are wanted for each Day and each Condition, for a single Subject (to show individual variability)
# You can adjust this to show average trends across subjects if preferred
predicted_values <- data_grid(sleepstudy, Days, condition, Subject = "308") %>%
  add_predictions(model, var = "pred")

newdata <-sleepstudy  #, expand.grid(condition = unique(condition), Days = unique(Days), Subject = unique(Subject)))
newdata$pred <- predict(model, newdata ) #, re.form = ~ condition + Days +(1 | Subject))

ggplot(newdata, aes(x = Days, y = pred, color = condition, group = Subject)) +
  geom_line() +
  ggtitle("Predicted Reaction Time over Days by Subject")

# Plotting the data and model predictions
ggplot(sleepstudy, aes(x = Days, y = Reaction, color = condition)) +
  geom_point(alpha = 0.5) + # Plot the raw data points
  geom_line(data = predicted_values, aes(y = pred, group = condition), size = 1) + # Add predicted lines
  labs(x = "Days of sleep deprivation",
       y = "Average reaction time (ms)",
       color = "Condition") +
  theme_minimal()



