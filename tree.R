# BiocManager::install("ggtree")
# BiocManager::install("treeio")

library(ggtree)
library(treeio)
tree <- read.iqtree("/Users/rajalaxmi/Desktop/msa_acn_mycovirus.afa.contree")
ggtree(tree,layout="daylight") 



ggtree(tree, layout="circular") + geom_tiplab(aes(angle=angle), color='blue')
p <- ggtree(tree)  
 
p+ geom_tiplab(size=3, color="purple")

ggtree(tree) + geom_treescale() 



data <- read.csv("/Users/rajalaxmi/Desktop/pfam_hits.csv", ",", header =T )
ggplot(data, aes(x=query.name, y=description.of.target, color=-log10(as.numeric(E.value)))) + geom_point()
str(data)
