library(decoupleR)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(OmnipathR)
library(progeny)
#Inputs
#Logtransformed data
Logcpm = cpm(y, log=TRUE)
rownames(Logcpm) = y$genes$gene_name
rownames(Logcpm) = toupper(rownames(Logcpm))
colnames(Logcpm)
#Design metadata
designmetadata = y$samples[,c(4,1)]
designmetadata
#DEG
deg = DEGs[[1]][!duplicated(DEGs[[1]]$gene_name),] %>%
  dplyr::select(gene_name, logFC, PValue) %>%
  filter(!is.na(gene_name)) %>%
  mutate(., t = -log10(PValue) * logFC)
rownames(deg) = toupper(deg$gene_name)
deg$symbol = toupper(deg$gene_name)
deg
#CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources. 
#TFs and their transcriptional targets
net <- get_collectri(organism='mouse', split_complexes=FALSE)
#To infer TF enrichment scores we will run the Univariate Linear Model (ulm) method.
#For each sample in our dataset (mat) and each TF in our network (net),
#it fits a linear model that predicts the observed gene expression based solely on the TF’s TF-Gene interaction weights. 
#Once fitted, the obtained t-value of the slope is the score.
#If it is positive, we interpret that the TF is active and if it is negative we interpret that it is inactive.
#To run decoupleR methods, we need an input matrix (mat),
#an input prior knowledge network/resource (net),
#and the name of the columns of net that we want to use.
sample_acts <- run_ulm(mat=Logcpm, net=net, .source='source', .target='target',
                       .mor='mor', minsize = 5)

#From the obtained results we will observe the most variable activities across samples in a heat-map:
n_tfs <- 25

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()
sample_acts
# Get top tfs with more variable means across clusters
tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 

#Also infer TF activities from the t-values of the DEGs between KO and WT:
contrast_acts <- run_ulm(mat=deg[, 't', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)

# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)
# Plot
ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "dodgerblue", high = "firebrick", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("TFs")
ggsave(paste0("results/", names[1], "/TFscore.pdf"), bg = "white", scale = 1)
contrast_acts %>%
  arrange(p_value) %>%
  dplyr::select(source, score, p_value) %>%
  write.xlsx2(x = . , file = paste0("results/", names[1], "/TFs.xlsx"), col.names = T, row.names = F)
# visualize the most differential target genes in each TF along their p-values to interpret the results.
#For example, let’s see the genes that are belong to SP1:
tf <- 'CDX2'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('gene', 'logfc', 'p_value', 't_value')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & t_value > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & t_value < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & t_value > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & t_value < 0, '1', color))

ggplot(df, aes(x = logfc, y = -log10(p_value), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(tf)

#Here blue means that the sign of multiplying the mor and t-value is negative, meaning that these genes are “deactivating” the TF,
#and red means that the sign is positive, meaning that these genes are “activating” the TF.

#PATHWAYS ESTO SON PATHWAYS RECUERDA ESTO SON PATHWAYS ASÍ QUE ATIENDE A QUE ESTO SON PATHWAYS
designmetadata
rownames(Logcpm) = y$genes$gene_name
degpath = deg %>%
  select(t)
rownames(degpath) = str_to_title(rownames(degpath))
degpath
#PROGENy is a comprehensive resource containing a curated collection of pathways and their target genes, with weights for each interaction.
netpath <- get_progeny(organism = 'mouse', top = 500)
#predicts the observed gene expression based on all pathways’ Pathway-Gene interactions weights.
#Once fitted, the obtained t-values of the slopes are the scores. If it is positive,
#we interpret that the pathway is active and if it is negative we interpret that it is inactive.
sample_acts1 <- run_mlm(mat=Logcpm, net=netpath, .source='source', .target='target',
                        .mor='weight', minsize = 5)
# Transform to wide matrix
sample_acts_mat1 <- sample_acts1 %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Scale per feature
sample_acts_mat1 <- scale(sample_acts_mat1)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(sample_acts_mat1, border_color = NA, color=my_color, breaks = my_breaks) 

#We can also infer pathway activities from the t-values of the DEGs between KO and WT:
contrast_acts1 <- run_mlm(mat=degpath, net=netpath, .source='source', .target='target',
                          .mor='weight', minsize = 5)
contrast_acts1

# Plot
ggplot(contrast_acts1, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "dodgerblue", high = "firebrick", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")
ggsave(paste0("results/", names[1], "/Pathwayscore.pdf"), bg = "white", scale = 1)
contrast_acts1 %>%
  arrange(p_value) %>%
  dplyr::select(source, score, p_value) %>%
  write.xlsx2(x = . , file = paste0("results/", names[1], "/pathways.xlsx"), col.names = T, row.names = F)

#Especificas
pathway <- 'NFkB'

df1 <- netpath %>%
  dplyr::filter(source == pathway) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter1 <- sort(intersect(rownames(degpath),rownames(df1)))
df1 <- df1[inter1, ]
df1['t_value'] <- degpath[inter1, ]
df1
df1 <- df1 %>%
  mutate(color = if_else(weight > 0 & t_value > 0, '1', color)) %>%
  mutate(color = if_else(weight > 0 & t_value < 0, '2', color)) %>%
  mutate(color = if_else(weight < 0 & t_value > 0, '2', color)) %>%
  mutate(color = if_else(weight < 0 & t_value < 0, '1', color))

ggplot(df1, aes(x = weight, y = t_value, color = color)) + geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(pathway)
