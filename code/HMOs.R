if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#Paquetes
list.of.packages = c("tximeta", "tximport", "limma", "edgeR", "tidyverse", "org.Mm.eg.db", "statmod", "pheatmap", "ggplotify", "ggrepel", "SummarizedExperiment", "patchwork", "xlsx", "ragg", "OmnipathR", "clusterProfiler")
#Instalación por CRAN o Bioconductor
new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
if(length(new.packages)> 0) {
  for(i in new.packages) {
    if(i %in% available.packages()[,1]){ #Chequea si el paquete está en un repositorio CRAN y lo instala
      install.packages(i,dependencies=TRUE)
    }else {BiocManager::install(i, update = TRUE, ask = FALSE, version = BiocManager::version()) #Instala por BiocManager si no está en un repositorio CRAN
    }}
}
invisible(lapply(list.of.packages, FUN=library, character.only=TRUE))
rm(list.of.packages, new.packages)
files = data.frame(read.table(file="metadata/metadata.txt", header = TRUE, stringsAsFactors = F))
suppressPackageStartupMessages(library(SummarizedExperiment))
se = tximeta(files)                   
se <- addExons(se)
gse <- summarizeToGene(se, assignRanges="abundant") #Asigna rangos según la isoforma más abundante de los transcritos en vez de desde la isoforma más en 5' a la isoforma más en 3'. Los creadores lo recomiendan.
y <- makeDGEList(gse)
rm(se)
sampleinfo <- read.delim(file ="metadata/sampleinfo.txt", sep = ",", header = T)
sampleinfo
group = sampleinfo[,2]
group = factor(group)
group = relevel(group, ref = "CONTROL")

y$samples$group = group
y$samples
design = model.matrix(~0 + group)
colnames(design) = levels(group)
print(design)
keep <- filterByExpr(y, design)
print(summary(keep))
y <- y[keep, ]
points <- c(1:10) # Formas
colors <- c(1:5) #Colores
mds = plotMDS(y, col=colors[group], pch=points[group])
mds_data <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  Sample = colnames(mds$distance.matrix.squared),
  Group = as.factor(y$samples$group)
)
axis_label = round(mds$var.explained[1:2]*100)
library(ggConvexHull)
mds_data %>% ggplot(aes(x = Dim1, y = Dim2, col = Group, shape = Group)) +
  geom_point(size = 2) +
  labs(title = "MDS Analysis", x = paste0("Dimension 1: ", axis_label[1], "%"), y = paste0("Dimension 2: ", axis_label[2], "%"), ) +
  theme_minimal() +
  scale_shape_manual(values = rep(15:17, len = 10)) + 
  geom_convexhull(aes(fill=Group, 
                      colour= Group),
                  alpha = 0.2) + 
  theme(plot.title = element_text(hjust = 0.5))
y <- estimateDisp(y, design, method = "ROBUST")
plotBCV(y)
fit <- glmQLFit(y, design, method = "ROBUST")
plotQLDisp(fit)
pairwisecomb = combn(levels(group), 2, function(x) paste(x[2], "-", x[1], sep = "")); pairwisecomb
contrast = makeContrasts(contrasts = pairwisecomb, levels=design)
RESs = vector(length = 9L)
names(RESs) = colnames(contrast[,c(1:9)])
for (i in c(1:9)){
  RESs[i] = topTags(glmQLFTest(fit, contrast = contrast[,i]), n = Inf)
}
DEGs = lapply(RESs, function(i) i %>% dplyr::filter(FDR <=0.05 ) %>% dplyr::select(symbol, gene_name, entrezid, logFC, PValue, FDR))
DEGs = DEGs[lapply(DEGs,nrow)>0]
data = lapply(RESs, function(i) i %>%
         dplyr::mutate(DE = case_when(logFC > 0 & FDR <0.05 ~ "UP",
                                      logFC < 0 & FDR <0.05 ~ "DOWN",
                                      .default = "NO")))
volcanoplot = function(data, name){
  ggplot(data = data, aes(x=logFC, y =-log10(PValue), col=DE, label = ifelse(abs(logFC) >= 1 & FDR <0.05, as.character(symbol),  ''))) + geom_point() +
  geom_text_repel(hjust = 0, nudge_x = 0.1, color = "black") +
  scale_color_manual(values = c("DOWN" = "firebrick", "UP" = "dodgerblue", "NO" = "grey")) +
  labs(title = name) +
  theme_minimal()
}
head(data)
plots <- lapply(names(data), function(name) {
  volcanoplot(data[[name]], name)
})
combined_plot <- wrap_plots(plots)
print(combined_plot)
