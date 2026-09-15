######################################
# Figure 3A 집컴터.
######################################
dir="C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_collecTRI = "C:/Dropbox/PNU/시스템생물학연구실/DB/collecTRI"

load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl
load(file = sprintf("%s/collectri_set_mouse.Rdata", dir_collecTRI)) #set_m
load(file = sprintf("%s/bulkRNA/Rdata/Figure3A.Rdata", dir)) #gluconeogenesis, gluconeogenesis.genes

library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(cluster)
library(ComplexHeatmap)
library(colorRamp2)
library(ggfortify)
library(openxlsx)
library(gridExtra)
library(circlize)
library(msigdbr)
library(fgsea)
library(igraph)
library(RCy3)
#BiocManager::install('decoupleR')
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(OmnipathR)
library(ggrepel)
library(GOfuncR)
library(biomaRt)
library(ggvenn)
library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')


### target 만 gluconeogenesis gene
deg = ddfl$`Starvation  vs  Control`
deg = deg[!is.na(deg$pvalue),]

deg.tf = deg %>%
  dplyr::select(log2FoldChange, stat, pvalue)
rownames(deg.tf) = rownames(deg)

deg.mtx = as.matrix(deg.tf)
deg.mtx = na.omit(deg.mtx)

gluconeogenesis = get_anno_genes('GO:0006094', database = "org.Mm.eg.db") 
gluconeogenesis.genes = gluconeogenesis$gene #106

set_m_glu = set_m[set_m$target %in% gluconeogenesis.genes, ]

sample.acts = run_ulm(mat = deg.mtx[, 'stat', drop=FALSE], net=set_m_glu, .source='source', .target='target', .mor='mor', minsize = 5)
sample.acts = sample.acts %>% arrange(desc(score))

f.acts = sample.acts %>% mutate(rank = NA)
msk = f.acts$score > 0
f.acts[msk, 'rank'] = rank(-f.acts[msk, 'score'])
f.acts[!msk, 'rank'] = rank(-abs(f.acts[!msk, 'score']))

top10_pos = f.acts %>%
  filter(score > 0) %>%
  arrange(rank) %>%
  slice(1:10)
top10_neg = f.acts %>%
  filter(score < 0) %>%
  arrange(rank) %>%
  slice(1:10)

top10 = rbind(top10_pos, top10_neg)
top3.neg = head(top10_neg$source, 3)


tiff(filename = sprintf("%s/figure/Figure3A.tiff", dir), width = 12, height = 7, units = 'cm', res = 300)
ggplot(top10, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  geom_text(aes(label = ifelse(score > 0 | source %in% top3.neg, source, ""), angle = 90), 
            position = position_stack(vjust = 0),  
            hjust = 0, size = 5)+
  geom_hline(yintercept = 0, color = "black", size = 0.2)+
  geom_vline(xintercept = 10.5, color = "black", size = 0.2) + 
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "white", midpoint = 0) + 
  scale_y_continuous(breaks = seq(-9, 9, by = 3), limits = c(-9, 9)) + 
  scale_x_discrete(limits = rev(levels(reorder(top10$source, top10$score)))) + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks.y = element_blank(), 
        panel.grid.major = element_line(color = "gainsboro", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gainsboro", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))
dev.off()


#save(gluconeogenesis, gluconeogenesis.genes, file = sprintf("%s/bulkRNA/Rdata/Figure3A.Rdata", dir))

bulk.top10 = bulk.top10[1:10,]

########################
# Figure 3B 집컴터.
########################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_collecTRI = "C:/Dropbox/PNU/시스템생물학연구실/DB/collecTRI"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/scRNA/Rdata/5-1_PTcell.Rdata", dir)) #pt.deg
load(file = sprintf("%s/scRNA/Rdata/5-1_PTvsOtherCell.Rdata", dir)) # stv.deg, con.deg
load(file = sprintf("%s/collectri_set_mouse.Rdata", dir_collecTRI)) #set_m
load(file = sprintf("%s/scRNA/Rdata/supple_Figure5B_TFactivity.Rdata", dir)) #acts 


library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(cluster)
library(ComplexHeatmap)
library(colorRamp2)
library(ggfortify)
library(openxlsx)
library(gridExtra)
library(circlize)
library(msigdbr)
library(fgsea)
library(igraph)
library(RCy3)
#BiocManager::install('decoupleR')
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(OmnipathR)
library(ggrepel)
library(GOfuncR)
library(biomaRt)
library(ggvenn)
library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')


### target 만 gluconeogenesis gene
pt.deg = pt.deg[!is.na(pt.deg$p_val),]
pt.deg$p_val[pt.deg$p_val == 0] = min(pt.deg$p_val[pt.deg$p_val != 0])
pt.deg$stat = -log10(pt.deg$p_val)*sign(pt.deg$avg_log2FC)

deg.tf = pt.deg %>%
  dplyr::select(avg_log2FC, stat, p_val)
rownames(deg.tf) = rownames(pt.deg)

deg.mtx = as.matrix(deg.tf)
deg.mtx = na.omit(deg.mtx)

gluconeogenesis = get_anno_genes('GO:0006094', database = "org.Mm.eg.db") 
gluconeogenesis.genes = gluconeogenesis$gene #106

set_m_glu = set_m[set_m$target %in% gluconeogenesis.genes, ]

sample.acts = run_ulm(mat = deg.mtx[, 'stat', drop=FALSE], net=set_m_glu, .source='source', .target='target', .mor='mor', minsize = 5)
sample.acts = sample.acts %>% arrange(desc(score))

f.acts = sample.acts %>% mutate(rank = NA)
msk = f.acts$score > 0
f.acts[msk, 'rank'] = rank(-f.acts[msk, 'score'])
f.acts[!msk, 'rank'] = rank(-abs(f.acts[!msk, 'score']))

top10_pos = f.acts %>%
  filter(score > 0) %>%
  arrange(rank) %>%
  slice(1:10)
top10_neg = f.acts %>%
  filter(score < 0) %>%
  arrange(rank) %>%
  slice(1:10)

top10 = rbind(top10_pos, top10_neg)
top3.neg = head(top10_neg$source, 3)

tiff(filename = sprintf("%s/figure/Figure5B.tiff", dir), width = 12, height = 7, units = 'cm', res = 300)
ggplot(top10, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  geom_text(aes(label = ifelse(score > 0 | source %in% top3.neg, source, ""), angle = 90), 
            position = position_stack(vjust = 0),  
            hjust = 0, size = 5)+
  geom_hline(yintercept = 0, color = "black", size = 0.2)+
  geom_vline(xintercept = 10.5, color = "black", size = 0.2) + 
  scale_fill_gradient2(low = "blue", high = "red", 
                       mid = "white", midpoint = 0) + 
  scale_y_continuous(breaks = seq(-9, 9, by = 3), limits = c(-9, 9)) + 
  scale_x_discrete(limits = rev(levels(reorder(top10$source, top10$score)))) + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks.y = element_blank(), 
        panel.grid.major = element_line(color = "gainsboro", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gainsboro", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))
dev.off()


sc.top10 = sc.top10[1:10,]

#save(bulk.top10, sc.top10, file = sprintf("%s/scRNA/Rdata/Figure3_common_rank.Rdata", dir))




###########################
# Figure 3C
###########################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_collecTRI = "C:/Dropbox/PNU/시스템생물학연구실/DB/collecTRI"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/scRNA/Rdata/5-1_PTcell.Rdata", dir)) #pt.deg
load(file = sprintf("%s/scRNA/Rdata/5-1_PTvsOtherCell.Rdata", dir)) # stv.deg, con.deg
load(file = sprintf("%s/collectri_set_mouse.Rdata", dir_collecTRI)) #set_m
load(file = sprintf("%s/scRNA/Rdata/supple_Figure5B_TFactivity.Rdata", dir)) #acts 

library('Seurat')
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(cluster)
library(ComplexHeatmap)
library(colorRamp2)
library(ggfortify)
library(openxlsx)
library(gridExtra)
library(circlize)
library(msigdbr)
library(fgsea)
library(igraph)
library(RCy3)
#BiocManager::install('decoupleR')
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(OmnipathR)
library(ggrepel)
library(GOfuncR)
library(biomaRt)
library(ggvenn)
library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(ggpubr)


# scRNA collecTri
stv.mat = as.matrix(s.integrated@assays$RNA$data.2)
stv.mat[1:10, 1:10]
row.sum = rowSums(stv.mat)
sum(row.sum == 0)
stv.mat = stv.mat[row.sum != 0, ]

acts = run_ulm(mat=stv.mat, net=set_m, .source='source', .target='target', .mor='mor', minsize = 5)

#save(acts, file = sprintf("%s/scRNA/Rdata/supple_Figure5B_TFactivity.Rdata", dir))

s.integrated[['tfsulm']] = acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(s.integrated) = "tfsulm"
s.integrated = ScaleData(s.integrated)
s.integrated@assays$tfsulm@data = s.integrated@assays$tfsulm@scale.data

DefaultAssay(s.integrated) = "tfsulm"
p1 = (FeaturePlot(s.integrated, features = 'Foxo1', min.cutoff = 0, max.cutoff = 0.8) & 
        scale_colour_gradient2(low = 'lightgrey', mid = 'lightgrey', high = 'violetred')) &
  theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 15), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

DefaultAssay(s.integrated) = "RNA"
p2 = FeaturePlot(s.integrated, features = 'Foxo1', max.cutoff = 0.8) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 15), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())


tiff(filename = sprintf("%s/figure/Figure3C.tiff", dir), width = 18, height = 8, units = 'cm', res = 300)
p1 | p2
dev.off()



######## PT cell 
DefaultAssay(s.integrated) = "RNA"
pt = subset(s.integrated, subset = cellType == "PT")
grcol = c('steelblue','darkorange')

tiff(filename = sprintf("%s/figure/Figure3C_PT_vln.tiff", dir), width = 7, height = 5, units = 'cm', res = 300)
VlnPlot(pt, features = 'Foxo1', split.by = 'sample', group.by = "cellType", cols = grcol, pt.size = 0, assay = 'RNA') + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none") +
  stat_compare_means(method = "wilcox.test", na.rm = TRUE)+
  geom_boxplot(width=0.1, position = position_dodge(0.9), outlier.shape = NA)
dev.off()





##################
# Figure 3D
##################
dir="E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/scRNA/Rdata/5-1_PTcell.Rdata", dir)) #pt.deg
load(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/omnipath_TF/mouse_tft.Rdata") #m.tft
load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl
load(file = sprintf("%s/bulkRNA/Rdata/gluconeogenesis_genes_of_학교버전.Rdata", dir)) #gluconeogenesis.genes

library(GOfuncR)
library(org.Mm.eg.db)
library(igraph)
library(RCy3)
library(OmnipathR)


gluconeogenesis = get_anno_genes('GO:0006094', database = "org.Mm.eg.db") 
gluconeogenesis.genes = gluconeogenesis$gene #94
#save(gluconeogenesis.genes, file = sprintf("%s/bulkRNA/Rdata/gluconeogenesis_genes_of_학교버전.Rdata", dir))

stv.gluconeo.upDEG = pt.deg[pt.deg$avg_log2FC > 0 & pt.deg$p_val < 0.05 & rownames(pt.deg) %in% gluconeogenesis.genes & !is.na(pt.deg$p_val_adj), ] 
stv.glu.up.genes = rownames(stv.gluconeo.upDEG) 

bulk.deg = ddfl$`Starvation  vs  Control`


################## TF-target
glu.tft = m.tft[m.tft$target_genesymbol %in% gluconeogenesis.genes,]
glu.tft = glu.tft[glu.tft$source_genesymbol %in% stv.glu.up.genes, ]
glu.tft = glu.tft[order(glu.tft$source_genesymbol),]

graph = graph_from_edgelist(as.matrix(glu.tft[,1:2]), directed = T)
graph = igraph::simplify(graph)

V(graph)[names(V(graph)) %in% glu.tft$source_genesymbol]$tft = "source"
V(graph)[names(V(graph)) %in% glu.tft$target_genesymbol]$tft = "target"  
V(graph)[names(V(graph)) %in% stv.glu.up.genes]$upDEG.gluconeo = "upDEG.gluconeo"  

V(graph)$sc_log2FC = pt.deg$avg_log2FC[match(names(V(graph)), rownames(pt.deg))]
V(graph)$sc_log2FC[is.na(V(graph)$sc_log2FC)] = 0
V(graph)$bulk_log2FC = bulk.deg$log2FoldChange[match(names(V(graph)), bulk.deg$Genes)]
V(graph)$bulk_log2FC[is.na(V(graph)$bulk_log2FC)] = 0

E(graph)$stimulation = glu.tft$is_stimulation
E(graph)$inhibition = glu.tft$is_inhibition
edge_attr(graph)

createNetworkFromIgraph(graph, "Figure3D")



