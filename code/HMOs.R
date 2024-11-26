if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#Paquetes
list.of.packages = c("tximeta", "tximport", "limma", "edgeR", "tidyverse", "org.Mm.eg.db", "statmod", "pheatmap", "ggplotify", "ggrepel", "SummarizedExperiment", "patchwork", "xlsx", "ragg", "OmnipathR", "clusterProfiler", "eulerr", "scales", "enrichplot", "GOSemSim")
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
RESs = vector(length = 13L)
names(RESs) = colnames(contrast[,c(1:9, 10, 34, 41, 32)])
for (i in c(1:length(RESs))){
  RESs[i] = topTags(glmQLFTest(fit, contrast = contrast[,i]), n = Inf)
}
DEGs = lapply(RESs, function(i) i %>% dplyr::filter(FDR <=0.05 ) %>% dplyr::select(description, gene_name, gene_id, entrezid, logFC, PValue, FDR))
DEGs = DEGs[lapply(DEGs,nrow)>0]
data = lapply(RESs, function(i) i %>%
         dplyr::mutate(DE = case_when(logFC > 0 & FDR <0.05 ~ "UP",
                                      logFC < 0 & FDR <0.05 ~ "DOWN",
                                      .default = "NO")))
volcanoplot = function(data, name){
  ggplot(data = data, aes(x=logFC, y =-log10(PValue), col=DE)) + geom_point() +
  scale_color_manual(values = c("DOWN" = "firebrick", "UP" = "dodgerblue", "NO" = "grey")) +
  labs(title = name) +
  theme_minimal()
}
plots <- lapply(names(data), function(name) {
  volcanoplot(data[[name]], name)
})
combined_plot <- wrap_plots(plots)
print(combined_plot)
ggsave("results/volcano.png", scale = 2)

logcpm = cpm(y, log=TRUE)
rownames(logcpm) = y$genes$gene_name #Nombre de genes
colnames(logcpm) =  paste(y$samples$group, y$samples$names, sep = "-")
head(logcpm)

DEG_selection = lapply(DEGs, function(i)
  logcpm[na.omit(i$gene_name),])


heatmap = function(data, samples, rows, title){
  as.ggplot(pheatmap(data[c(1:min(rows, 100)),samples], scale = "row", 
                     clustering_method = "complete",
                     display_numbers = F,
                     border_color = NA, cluster_cols = T, cutree_cols = 2, cutree_rows = 2, show_rownames = F,
                     #annotation_col = annotation, show_rownames = F, annotation_names_col = F,
                     #annotation_row = setNames(data.frame(Cluster = as.factor(cutree(q$tree_row, k=2))), "Cluster"), 
                     annotation_names_row = F, main = title,
                     legend_labels = F,))
}

heatplots = lapply(names(DEG_selection[unlist(lapply(DEG_selection, is.matrix))]), function(name) {
  string = str_split_1(name, pattern = "-")
  cols = grep(paste0("\\b", string, "\\b", collapse = "|"), y$samples$group)
  heatmap(DEG_selection[[name]], cols, nrow(DEG_selection[[name]]), name)
  
})
combined_heatplot <- wrap_plots(heatplots)
print(combined_heatplot)
ggsave("results/heatmaps.png", scale = 2)

write.xlsx2(x = res_corrected[!is.na(res_corrected$table$gene_name) & res_corrected$table$PValue <= 0.05 ,c(1,7,13,16,17)], file = paste0(experimento, "/",experimento, ".xlsx"), col.names = T, row.names = F, append = TRUE)

lapply(names(RESs), function(name) {
  write.xlsx2(x = RESs[[name]][!is.na(RESs[[name]]$gene_name) & RESs[[name]]$PValue <= 0.05, c(1,5,7,13,16,17)],
              file = "results/genes.xlsx", sheetName = name ,col.names = T, row.names = F, append = TRUE)
  gc()
})

lapply(names(DEGs), function(name) {
  write.xlsx2(x = DEGs[[name]],
              file = "results/degs.xlsx", sheetName = name ,col.names = T, row.names = F, append = TRUE)
  gc()
})

comparisons = tail(names(DEGs),4)
genes = lapply(DEGs[comparisons], function(i) i[,"gene_name"])
alldegs = unique(unlist(genes))
matrix = matrix(nrow = length(alldegs), ncol = 5)
colnames(matrix) = c("genes", comparisons)
matrix[,"genes"] = alldegs
for (i in c(2:5)){
  matrix[,i] = matrix[,1] %in% genes[[i-1]]
}
matrix = as.data.frame(matrix)
matrix <- matrix %>% mutate(across(c(2:5), as.logical))
as.ggplot(plot(eulerr::euler(matrix[,2:5]),  quantities = TRUE))
ggsave("results/eulerr.png")
common = matrix[which(rowSums(matrix[,c(2:5)]) >1 ),]
write.xlsx2(x = common, file = "results/comunes.xlsx", col.names = T, row.names = F)
single =  matrix[which(rowSums(matrix[,c(2:5)]) == 1 ),]
write.xlsx2(x = single, file = "results/singulares.xlsx", col.names = T, row.names = F)

#Enrichment Analysis
sig_down = lapply(DEGs, function(i) i %>% dplyr::filter(logFC < 0 ))
sig_up = lapply(DEGs, function(i) i %>% dplyr::filter(logFC > 0 ))
GSEA_data <- lapply(RESs, function(i) i %>% filter(!is.na(FDR)) %>%
  arrange(desc(logFC)))
log2FC = lapply(GSEA_data, function(i) i %>%
  pull(logFC, name = gene_name))
rm(GSEA_data)
names = make.names(names(RESs))
lapply(paste0("results/",names), FUN = dir.create)
#Remember which one you're doing
names[10]
sig_up[[10]]
sig_down[[10]]
#GO enrichment on up-regulated genes. Genes that underwent differential expression testing. Up 1st Down 2nd
ego <- enrichGO(gene= sig_down[[10]]$gene_name, #OVER-REPRESENTATION TEST
                OrgDb= org.Mm.eg.db,
                keyType = "SYMBOL",
                ont= "BP",
                universe=RESs[[10]]$gene_name)
s_ego<-clusterProfiler::simplify(ego)
s_ego@result %>% head()
s_ego@result %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  write.xlsx2(x = . , file = paste0("results/", names[10], "/DOWNGOBP.xlsx"), col.names = T, row.names = F)
s_ego %>% filter(p.adjust < 0.05) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  scale_y_discrete(labels = wrap_format(60)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment of down-regulated genes")
ggsave(paste0("results/", names[10], "/GOUPdotplot.pdf"), bg = "white", scale = 1)
pairwise_termsim(s_ego, method = "Wang", semData = godata('org.Mm.eg.db', ont="BP")) %>%
  treeplot(showCategory = 20, label_format_tiplab = 40)
ggsave(paste0("results/", names[10], "/GOUPtreeplot.pdf"), bg = "white", scale = 1.5)
#KEGG enrichment
kgo <- enrichKEGG(gene = sig_down[[10]]$entrezid,
                  organism = 'mmu',
                  pvalueCutoff  = 0.05,
                  universe = as.character(RESs[[10]]$entrezid)) %>%
  mutate(., Description = sub(pattern = " [-] Mus musculus [(]house mouse[)]", "",  x = Description))
kgo %>% filter(p.adjust < 0.05) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  scale_y_discrete(labels = wrap_format(60)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("KEGG Enrichment of down-regulated genes")
ggsave(paste0("results/", names[10], "/KEGGUPdotplot.pdf"), bg = "white", scale = 1)
kgo@result %>%
  dplyr::select(ID, subcategory, Description, GeneRatio, pvalue, p.adjust) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  write.xlsx2(x = . , file = paste0("results/", names[10], "/UPkegg.xlsx"), col.names = T, row.names = F)
##Merge Up & Down
comparelist<-list(sig_down[[9]]$gene_id,sig_up[[9]]$gene_id)
names(comparelist)<-c("down-regulated","up-regulated") 
cclust<-compareCluster(geneCluster = comparelist, 
                       fun = enrichGO,
                       OrgDb= org.Mm.eg.db,
                       keyType = "ENSEMBL",
                       ont= "BP",
                       universe=RESs[[9]]$gene_id)
#options(enrichplot.colours = c("firebrick","dodgerblue"))
dotplot(cclust,showCategory=10, label_format = 50)
ggsave(paste0("results/", names[9], "/comparedplot.pdf"), bg = "white", scale = 1)
#GSEA with clusterProfiler
##all of the data is used regardless of arbitrary cut-offs like p-values
set.seed(123) #RANDOM ORDERING => USE THE SAME SEED  
gsea_go <- gseGO(geneList = log2FC[[10]],
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 keyType = "SYMBOL",
                 seed=TRUE)
gsea_go@result %>% head()
gsea_go@result %>%
  dplyr::select(ID, Description, enrichmentScore, NES, pvalue, p.adjust) %>%
  write.xlsx2(x = . , file = paste0("results/", names[10], "/GSEA.xlsx"), col.names = T, row.names = F)
ggplot(gsea_go, showCategory=20, aes(NES, fct_reorder(Description, NES),
                                      fill=qvalue)) +
  geom_col() +
  geom_vline(xintercept=0, linetype="dashed", color="black",size=1)+
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE))+
  scale_x_continuous(expand=c(0,0))+
  theme_bw() + 
  #theme(text=element_text(size=8))+
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA with GO")
ggsave(paste0("results/", names[10], "/GSEAGO.pdf"), bg = "white", scale = 1)
pairwise_termsim(gsea_go, method = "Wang", semData = godata('org.Mm.eg.db', ont="BP")) %>%
  treeplot(showCategory = 20, label_format_tiplab = 40)
ggsave(paste0("results/", names[10], "/GSEAtree.pdf"), bg = "white", scale = 1.5)
