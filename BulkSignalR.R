# install.packages("devtools")
# devtools::install_github("jcolinge/BulkSignalR")
# BiocManager::install("SingleCellSignalR")
library(BulkSignalR)
library(igraph)
library(dplyr)

ensemble_list = unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))
ids<-bitr(ensemble_list, fromType="ENSEMBL",
          toType=c("SYMBOL"),
          OrgDb="org.Hs.eg.db") #NA values drop by default
counts_sym <- counts %>% as.data.frame() %>%
  rownames_to_column(var ="ENSEMBL") %>%
  mutate(ENSEMBL = unlist(lapply(strsplit(ENSEMBL, '\\.'), '[[',1))) %>% 
  left_join(ids)   %>%
  filter(!is.na(SYMBOL)) %>%
  filter(!duplicated(SYMBOL)) %>%
  column_to_rownames("SYMBOL") %>% dplyr::select(-ENSEMBL)

dfh<-data.frame(Sample_ID=as.character(colnames(counts_sym)[1:2]))%>%
  left_join(condition %>% select(Sample_ID,  Disease, Category))


dfh

# write.csv(counts_sym , "data/count_symbol.csv")

# browseVignettes("BulkSignalR")

bsrdm <- prepareDataset(counts = counts_sym[,1:2])
# print object
str(bsrdm)

bsrdm <- learnParameters(bsrdm)

bsrdm <- learnParameters(bsrdm, 
                         plot.folder = "data/", filename = "LigandRS")
bsrinf <- initialInference(bsrdm)
bsrinf.redBP    <- reduceToBestPathway(bsrinf)

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=0.001)

scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP,
                                  name.by.pathway=FALSE)

ha = ComplexHeatmap::HeatmapAnnotation(bar = dfh$Disease,
                        col = list(bar = c("PDAC" = "red", "Progressor HGD w PDAC" = "green", "Progressor LGD w PDAC" = "blue")))

simpleHeatmap(scoresLR[1:20,], 
              path = "data/",
              filename = "LigandRS_scores",
              column.names = TRUE, 
              height = 5, width = 9,
              pointsize = 10,
              hcl.palette = "Cividis",
              bottom.annotation = ha
)
bsrinf.redP  <- reduceToPathway(bsrinf)  
bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 
bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP, qval.thres=0.001)
scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig.redPBP,
                                       name.by.pathway=TRUE)

getPathwayStats(bsrinf)
bsrinf.redP[str_detect(bsrinf.redP$pw.name, "PD-1"),]

simpleHeatmap(scoresPathway, 
              path = "data/",
              filename = "scoresPathway",
              column.names = TRUE, 
              width = 9, 
              height = 12, 
              pointsize = 12,
              hcl.palette = "Blue-Red 2"
)

chordDiagramLR (bsrinf,
                path = "data/",
                filename = "sdc_chord",
                pw.name = "R-HSA-202733",
                limit = 20,
                width = 5, 
                height = 4.5
)

subset <- c("REACTOME_BASIGIN_INTERACTIONS",
            "REACTOME_SYNDECAN_INTERACTIONS",
            "REACTOME_ECM_PROTEOGLYCANS",
            "REACTOME_CELL_JUNCTION_ORGANIZATION")

reactSubset <- BulkSignalR:::.SignalR$BulkSignalR_Reactome[
  BulkSignalR:::.SignalR$BulkSignalR_Reactome$`Reactome name` %in% subset,]
gLR <- getLRNetwork(bsrinf.redBP, qval.thres=1e-8)

# save to file
# write.graph(gLR,file="SDC-LR-network.graphml",format="graphml")
gLR <- read_graph("data/SDC-LR-network.graphml", format="graphml")
# As an alternative to Cytoscape, you can play with igraph package functions.
plot(gLR,
     edge.arrow.size=0.1,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.1)
