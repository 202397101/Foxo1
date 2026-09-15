##################
# Figure 6A
##################

dir = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_geo = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO"

load(file = sprintf("%s/GSE199711_snRNA/Rdata/1_integration.Rdata", dir_geo)) #s.integrated, markers
load(file = sprintf("%s/GSE199711_snRNA/Rdata/2_DEGs.Rdata", dir_geo)) #markers.li, cluster.degl

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(randomcoloR)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)



s.integrated$group = factor(s.integrated$group, levels = c("healthy", "CKD"))

mycol = c("purple","mediumpurple1","violetred","pink","deeppink","thistle","darkgreen","chartreuse1","cadetblue","royalblue","steelblue1","slategray2")
grcol = c("khaki","darkseagreen")
scol = c("moccasin","lightgoldenrod","palegoldenrod","palegreen","lightgreen","mediumspringgreen")

tiff(filename = sprintf("%s/GSE199711_snRNA/figure/1_dimplot_celltype.tiff", dir_geo), width = 15, height = 15, units = 'cm', res = 300)
DimPlot(s.integrated, reduction='umap', group.by='cellType', cols = mycol) + ggtitle('') + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
dev.off()



##################
# Figure 6B
##################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_geo = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO"

load(file = sprintf("%s/GSE199711_snRNA/Rdata/1_integration.Rdata", dir_geo)) #s.integrated, markers
load(file = sprintf("%s/GSE199711_snRNA/Rdata/2_DEGs.Rdata", dir_geo)) #markers.li, cluster.degl

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(randomcoloR)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(igraph)
library(RCy3)



pt = cluster.degl$PT
updeg = pt[pt$p_val_adj < 0.01 & pt$avg_log2FC > 0 & !is.na(pt$p_val_adj),]

ppi = read.csv(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/string/homo_sapiens/9606.protein.links.v12.0.txt", sep = " ", header = T, stringsAsFactors = F, quote = "")
pinfo = read.csv(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/string/homo_sapiens/9606.protein.info.v12.0.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
ppi$protein1 = pinfo$preferred_name[match(ppi$protein1, pinfo$X.string_protein_id)]
ppi$protein2 = pinfo$preferred_name[match(ppi$protein2, pinfo$X.string_protein_id)]
ppi1 = ppi[ppi$protein1 %in% rownames(updeg) & ppi$protein2 %in% rownames(updeg), ]
ppi1 = ppi1[ppi1$combined_score > 400, ]

graph = graph_from_edgelist(as.matrix(ppi1[,1:2]), directed = F)
graph = igraph::simplify(graph)
cluster = cluster_walktrap(graph)
cluster.size = sizes(cluster)
cluster.size

clusters.to.keep = which(cluster.size > 10)
nodes.to.remove = V(graph)$name[!(membership(cluster) %in% clusters.to.keep)]
graph.filtered = delete.vertices(graph, nodes.to.remove)
cluster.filtered = cluster_walktrap(graph.filtered)
sizes(cluster.filtered) 
max(sizes(cluster.filtered)) #1번 모듈.


module.kegg = c()
i=2
for (i in 1:length(cluster.filtered)){
  g = cluster.filtered[[i]]
  entrez = bitr(g, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
  kegg = as.data.frame(enrichKEGG(gene = entrez$ENTREZID, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500))
  kegg = kegg[order(kegg$pvalue),]
  module.kegg[[i]] = kegg
}
lapply(module.kegg, function(i) i[grep('senescence', i$Description),])

mo = module.kegg[[1]]
mo.filtered = mo[mo$category %in% c("Cellular Processes"),]
mo.filtered$p.adjust = -log10(mo.filtered$p.adjust)

ggplot(data = mo.filtered, aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
  coord_flip() +
  #geom_hline(yintercept = 2, color = "red", linewidth = 0.2, linetype = "dashed") +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(low = 'mediumpurple1', mid = 'mediumpurple1', high = 'mediumpurple1',
                       guide = FALSE, midpoint = 2) +
  ylab("-log10 FDR") + xlab("") +
  ggtitle("Module 1") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    plot.title  = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.1, color = "grey"),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 0.2, r = 0.5, b = 0.3, l = 0, unit = "cm")
  ) +
  geom_vline(xintercept = 0, color = "dimgrey", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "dimgrey", linewidth = 0.3) +
  geom_text(aes(label = Description, y = 0.01),
            hjust = 0, size = 4, color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 1))


##################
# Figure 6C
##################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_geo = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO"

load(file = sprintf("%s/GSE199711_snRNA/Rdata/1_integration.Rdata", dir_geo)) #s.integrated, markers
load(file = sprintf("%s/GSE199711_snRNA/Rdata/2_DEGs.Rdata", dir_geo)) #markers.li, cluster.degl

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(randomcoloR)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)



c2 = read.gmt("E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
c2 = split(x = c2$gene, f = c2$term )
names(c2)[grep('SENESCENCE', names(c2))]
inpGS = c2[grep('KAMMINGA_SENESCENCE', names(c2))]

DefaultAssay(s.integrated) = 'RNA'

genes = rownames(s.integrated)
inpGS.filtered = lapply(inpGS, function(i) {intersect(i, genes)})
new = JoinLayers(s.integrated, overwrite = T)
new = AddModuleScore(new, features = inpGS.filtered, name = "signaling.path", assay = 'RNA')
colnames(new@meta.data)[grep("signaling.path", colnames(new@meta.data))] = names(inpGS)

tiff(filename = sprintf("%s/GSE199711_snRNA/figure/5_senescence_module-score.tif", dir_geo), width = 20, height = 8, units = 'cm', res = 300)
FeaturePlot(new, reduction = "umap", pt.size = 0.1, features = names(inpGS), order = TRUE, split.by = 'group', max.cutoff = 0.35, min.cutoff = -0.2) & scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "RdGy"))) & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10)) 
dev.off()


##################
# Figure 6D
##################
sub = subset(new, subset = cellType %in% c('PT'))
grcol = c("khaki","darkseagreen")

tiff(filename = sprintf("%s/GSE199711_snRNA/figure/5_senescence_module-score_PT.tif", dir_geo), width = 7, height = 10, units = 'cm', res = 300)
VlnPlot(sub, pt.size = 0, features ="KAMMINGA_SENESCENCE", group.by = "cellType", split.by = 'group', assay = 'RNA', cols = grcol, y.max=0.48)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none") +
  stat_compare_means(method = "wilcox.test",na.rm = TRUE) +
  geom_boxplot(width=0.4, position = position_dodge(0.9), outlier.shape = NA)
dev.off()





##################
# Figure 6E
##################
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(randomcoloR)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggpubr)



dir_data = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO/GSE197266_FA_mouse_scNRA"

load(file=file.path(dir_data, "Rdata/1_sobj.Rdata")) #s.merged, s.integrated



mycol = c("purple","mediumpurple1","deeppink","thistle","chartreuse1","cyan","cadetblue","yellowgreen","royalblue","slategray2")
names(mycol) = c("Podocyte","PT","ATL/TAL/DCT/CNT/CD-PC","CD-IC","T cell","B cell","Neutrophil","Macrophage","Endothelia","Pericyte")
condcol = c("Ctrl"="slategrey","FA"="orangered")


s.sub = subset(
  s.integrated,
  subset = condition %in% c("Ctrl", "FA")
)

table(s.sub$condition)
Idents(s.sub)


tiff(filename = file.path(dir_data, "/figure/1_dimplot_celltype.tiff"), width = 15, height = 15, units = 'cm', res = 300)
DimPlot(s.sub, reduction='umap', group.by='cellType', cols = mycol) + ggtitle('') + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
dev.off()





##################
# Figure 6F
##################
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(randomcoloR)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggpubr)


dir_data = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO/GSE197266_FA_mouse_scNRA"

load(file=file.path(dir_data, "Rdata/1_sobj.Rdata")) #s.merged, s.integrated


c2 = read.gmt("E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
c2 = split(x = c2$gene, f = c2$term)
names(c2)[grep('GLUCONEO', names(c2))]
gluco_genes  = c2[["MOOTHA_GLUCONEOGENESIS"]]
gluco_mouse  = paste0(toupper(substr(gluco_genes, 1, 1)),
                      tolower(substr(gluco_genes, 2, nchar(gluco_genes))))

DefaultAssay(s.sub) = "RNA"
available    = rownames(s.sub)
gluco_keep   = gluco_mouse[gluco_mouse %in% available]
s.sub = JoinLayers(s.sub, assay="RNA")
s.sub = AddModuleScore(s.sub, features = list(gluconeogenesis = gluco_keep), name = "gluconeogenesis", assay = "RNA")

sene_genes  = c2[["KAMMINGA_SENESCENCE"]]
sene_mouse  = paste0(toupper(substr(sene_genes, 1, 1)),
                     tolower(substr(sene_genes, 2, nchar(sene_genes))))

DefaultAssay(s.sub) = "RNA"
available    = rownames(s.sub)
sene_keep   = sene_mouse[sene_mouse %in% available]
s.sub = JoinLayers(s.sub, assay="RNA")
s.sub = AddModuleScore(s.sub, features = list(senescence = sene_keep), name = "senescence", assay = "RNA")

s.pt = subset(s.sub, subset = cellType == "PT")
s.pt = JoinLayers(s.pt, assay="RNA")
s.pt = AddModuleScore(s.pt, features = list(gluconeogenesis = gluco_keep), name = "gluconeogenesis", assay = "RNA")
colnames(s.pt@meta.data)[grep("^gluconeogenesis", colnames(s.pt@meta.data))] = "gluconeogenesis_score"

condcol = c("Ctrl"="slategrey","FA"="orangered")
s.pt$condition = factor(s.pt$condition, levels=c("Ctrl","FA"))

tiff(filename = file.path(dir_data, "figure/2_PT_gluconeogenesis_score.tif"), width = 7, height = 7, units = 'cm', res = 300)
VlnPlot(s.pt, pt.size = 0, features = "gluconeogenesis_score", group.by = "condition", split.by = 'condition', assay = 'RNA', cols = condcol) + 
  #scale_y_continuous(limits = c(-0.38, 0.8)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, colour = "black", angle = 0, hjust = 0.5),
        plot.title = element_blank(),
        legend.position = "right") +
  geom_boxplot(width=0.3, position = position_dodge(0.9), outlier.shape = NA, linewidth = 0.3)
#stat_compare_means(method = "wilcox.test", na.rm = TRUE)
dev.off()


s.pt = AddModuleScore(s.pt, features = list(senescence = sene_keep), name = "senescence", assay = "RNA")
colnames(s.pt@meta.data)[grep("^senescence", colnames(s.pt@meta.data))] = "senescence_score"

tiff(filename = file.path(dir_data, "figure/2_PT_senescence_score.tif"), width = 7, height = 7, units = 'cm', res = 300)
VlnPlot(s.pt, pt.size = 0, features = "senescence_score", group.by = "condition", split.by = 'condition', assay = 'RNA', cols = condcol) + 
  #scale_y_continuous(limits = c(-0.38, 0.8)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, colour = "black", angle = 0, hjust = 0.5),
        plot.title = element_blank(),
        legend.position = "right") +
  geom_boxplot(width=0.3, position = position_dodge(0.9), outlier.shape = NA, linewidth = 0.3)
#stat_compare_means(method = "wilcox.test", na.rm = TRUE)
dev.off()




##################
# Figure 6H
##################

dir = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/cisplatin"

load(file = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/bulkRNA/Rdata/1_annotation.Rdata") # ginfo, bm
load(file = "E:/Dropbox/PNU/시스템생물학연구실/data/aging/results/Rdata/1_annotation_aging.Rdata") #gtf

load(file = sprintf("%s/Rdata/1_bulkRNA_rawCount_TPM.Rdata", dir)) #m.rcm, m.tpm, r.rcm, r.tpm
load(file = sprintf("%s/Rdata/2_mouse_deg.Rdata", dir)) #m.deg
load(file = sprintf("%s/Rdata/2_rat_deg.Rdata", dir)) #r.degl
load(file = "E:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/MsigDB_humanTOmouse_symbol.Rdata") #wpl, msig
load(file = sprintf("%s/Rdata/Foxo-associated-genes.Rdata", dir)) #foxo


library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(factoextra)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(colorRamp2)


max = apply(r.tpm, 1, max)
summary(max) # 3사분위수 기준. 75%는 최대 발현값이 1.5 이하.
filteredTpm = r.tpm[max >= 1.5,] #7670
group.map = c("CON1" = "CON", "CON2" = "CON", "CON3" = "CON",
              "Cis_48h_1" = "Cis_48h", "Cis_48h_2" = "Cis_48h", "Cis_48h_3" = "Cis_48h",
              "Cis_72h_1" = "Cis_72h", "Cis_72h_2" = "Cis_72h", "Cis_72h_3" = "Cis_72h")
sample.gr = group.map[colnames(filteredTpm)]
avg.expr = sapply(unique(sample.gr), function(v) {
  rowMeans(filteredTpm[, sample.gr == v, drop = FALSE])
})
avg.expr = log2(avg.expr+1)
expr.mat = as.matrix(avg.expr)
var = apply(expr.mat, 1, var, na.rm = TRUE)
expr.mat.filtered = expr.mat[var > 0, ]
scaled.expr = t(scale(t(expr.mat)))

gene.dist = dist(scaled.expr)
hc = hclust(gene.dist, method = "ward.D2")

ep = fviz_nbclust(scaled.expr, FUN = hcut, method = "wss", k.max = 10, linecolor = "black") + 
  ggtitle("Elbow method for hierarchical clustering") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        title = element_blank()) 
clu6 = ep$data[6,]

clusters = cutree(hc, k = 6)
cluster.genes = split(names(clusters), clusters)
cluster.genes = lapply(cluster.genes, function(i) {
  symbol = gtf$gene_name[match(i, gtf$gene_id)]
  symbol = symbol[!is.na(symbol) & symbol != ""]
  unique(symbol)   
})




#################### pathway
oral = list()
i=3
for (i in 1:length(cluster.genes)) {  
  genes = cluster.genes[[i]]
  entrez = bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Rn.eg.db")
  ora = as.data.frame(enrichKEGG(gene = entrez$ENTREZID, organism = 'rno', pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500))
  ora = ora[ora$category %in% c("Metabolism", "Environmental Information Processing"),]
  oral[[names(cluster.genes)[i]]] = ora
}



top10 = lapply(oral, function(i) head(i, 10))

clu = top10$'4'
clu$p.adjust = -log10(clu$p.adjust)

ggplot(data = clu, aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
  coord_flip() +
  #geom_hline(yintercept = 2, color = "red", linewidth = 0.2, linetype = "dashed") +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(low = 'cyan', mid = 'cyan', high = 'cyan',
                       guide = FALSE, midpoint = 2) +
  ylab("-log10 FDR") + xlab("") +
  ggtitle("Cluster 4") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    plot.title  = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.1, color = "grey"),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 0.2, r = 0.5, b = 0.3, l = 0, unit = "cm")
  ) +
  geom_vline(xintercept = 0, color = "dimgrey", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "dimgrey", linewidth = 0.3) +
  geom_text(aes(label = Description, y = 0.01),
            hjust = 0, size = 4, color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 1))



clu = top10$'5'
clu$p.adjust = -log10(clu$p.adjust)

ggplot(data = clu, aes(x = reorder(Description, p.adjust), y = p.adjust, fill = p.adjust)) +
  coord_flip() +
  #geom_hline(yintercept = 2, color = "red", linewidth = 0.2, linetype = "dashed") +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient2(low = 'cyan', mid = 'cyan', high = 'cyan',
                       guide = FALSE, midpoint = 2) +
  ylab("-log10 FDR") + xlab("") +
  ggtitle("Cluster 5") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    plot.title  = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.1, color = "grey"),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 0.2, r = 0.5, b = 0.3, l = 0, unit = "cm")
  ) +
  geom_vline(xintercept = 0, color = "dimgrey", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "dimgrey", linewidth = 0.3) +
  geom_text(aes(label = Description, y = 0.01),
            hjust = 0, size = 4, color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 1))



##################
# Figure 6J,K
##################

dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/cisplatin"

load(file = sprintf("%s/Rdata/1_bulkRNA_rawCount_TPM.Rdata", dir)) #m.rcm, m.tpm, r.rcm, r.tpm
load(file = sprintf("%s/Rdata/2_mouse_deg.Rdata", dir)) #m.deg
load(file = sprintf("%s/Rdata/2_rat_deg.Rdata", dir)) #r.degl
load(file = "C:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/MsigDB_humanTOmouse_symbol.Rdata") #wpl, msig

library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(colorRamp2)


m.deg <- m.deg %>%
  group_by(Genes) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()

rownames(m.deg) <- m.deg$Genes

idx = !is.na(m.deg$stat)
use = m.deg[idx,]
rk = use$stat
names(rk) = use$Genes
rk = sort(rk, decreasing = T)
sum(duplicated(names(rk)))

gsea=as.data.frame(fgsea(pathways=msig,stats=rk,minSize=10,maxSize=500,nproc=1))
mootha = gsea[grep('MOOTHA_', gsea$pathway),]
mootha$score = -log10(mootha$pval) * sign(mootha$NES)
mootha$pathway = gsub("MOOTHA_", "", mootha$pathway)
mootha$pathway = gsub("_", " ", mootha$pathway)
mootha$pathway = sapply(mootha$pathway, function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})
mootha = mootha[order(mootha$score),]
mootha = mootha[mootha$padj < 0.01, ]


ggplot(data=mootha,aes(x=reorder(pathway,-NES),y=NES,fill=NES))+
  geom_bar(stat="identity",width=0.8)+
  coord_flip()+
  scale_y_continuous(limits=c(-3,0))+
  scale_fill_gradient(low="slategray2",high="slategray2",guide=FALSE)+
  geom_hline(yintercept=0,color="black",linewidth=0.5)+
  ggtitle("Cis vs CON")+
  ylab("NES")+xlab("")+
  theme_minimal()+
  theme(axis.text.x=element_text(size=9),
        plot.title=element_text(size=12,hjust=0.5),
        axis.text.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_line(size=0.1,color="grey"),
        panel.grid.minor.x=element_blank(),
        plot.margin=margin(t=0.2,r=0,b=0.3,l=-0.5,unit="cm"))+
  geom_text(aes(label=pathway,y=-0.2),hjust=1,size=3.5,color="black")+
  geom_segment(aes(x=0.5,xend=0.5,y=-3,yend=0),inherit.aes=FALSE,color="dimgrey",size=0.3)+
  geom_segment(aes(x=0.5,xend=length(unique(mootha$pathway))+0.5,y=0,yend=0),color="dimgrey",size=0.3)



#############################################################

gset = msig[grep('MOOTHA_GLUCONEOGENESIS', names(msig))]
labs = list(mt="Gluconeogenesis", redgroup.lab="Cis", bluegroup.lab="CON", mlab="")
length(rk)
xmax = 32300
hcol=c("blue", "white", "red")

mt=labs$mt
redgroup.lab=labs$redgroup.lab
bluegroup.lab=labs$bluegroup.lab
class.name=labs$mlab
a = plotEnrichment(gset[[1]],rk)
fres=a$data
colnames(fres)=c("x","y")
es.range = c(min(fres$y), max(fres$y))
gsea.rnk = data.frame(names(rk),rk)
colnames(gsea.rnk) = c("id", "metric")
gsea.rnk=gsea.rnk[order(- gsea.rnk$metric),]
metric.range = c(min(gsea.rnk$metric), max(gsea.rnk$metric))
ss=fgsea(pathways = gset, stats = rk,nperm=10000,nproc=0)
nes=round(ss$NES,digits=2)
pv=round(ss$pval,digits=2)

gsea.layout= layout(matrix(c(1, 2, 3)), heights = c(1.5,0.2,0.6))
#layout.show(gsea.layout)
par(mar = c(0,3.5,0.1,0.8))
ln1=round(seq(es.range[1], es.range[2], 0.1), digits = 1)
plot(fres$x,fres$y, type = "l", col = "green", lwd = 2.5, xaxt = "n", xlab = "", ylab = ""
     , xaxs = "i", yaxt = "n", yaxs = "i", main = "", xpd=T
     , ylim = c(es.range[1]-0.05,es.range[2]+0.05), xlim=c(0,xmax)
     , panel.first = { abline(h = ln1 , col = "gray95", lty = 2)
       abline(h = 0, col = "black", lty = 2)
       abline(h = es.range[ifelse(nes<0,1,2)], col = "red", lty = 2)} )
ln2=round( seq(es.range[1], es.range[2]+0.05, 0.2), digits = 1)
axis(2,at = ln2,labels=rep("",length(ln2)), cex.axis=1, pos=0,tck=-0.03)
axis(2,at = ln2,labels=ln2, las=2, cex.axis=1, pos=-length(gsea.rnk$metric)*0.0001, tick=FALSE)
axis(2,at = es.range[1]+(es.range[2]-es.range[1])/2, labels ="Enrichment score" ,tick=FALSE, cex.axis=1.2 ,pos=-length(gsea.rnk$metric)*0.068)
mtext(mt,line=0.2, cex=0.8)
plot.coordinates = par("usr")
text(length(gsea.rnk$metric)*0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.15)
     , paste("p-value:", pv,"\nNES:", nes), adj = c(1,1), cex=1)
par(mar = c(0,3.5, 0.1,0.8))
plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n", ylab = "", xlim = c(1, xmax))
abline(v =fres$x, lwd = 0.8)
mx=min(abs(quantile(gsea.rnk$metric,c(0.01,0.99))))
metric.range = c(-mx, mx)
rank.colors = gsea.rnk$metric
rank.colors[rank.colors>=mx] = mx
rank.colors[rank.colors<= -mx]= -mx
rank.colors = (rank.colors - metric.range[1])/ (metric.range[2] - metric.range[1])
rank.colors = ceiling(rank.colors * 255 + 1)
rank.colors = colorRampPalette(hcol)(256)[rank.colors]
rank.colors = rle(rank.colors)

par(mar = c(3,3.5, 0.5,0.8)) # 아래, 왼, 위, 오른:: 아래 숫자 조절
bp=barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n"
           , xlim = c(0, xmax), width = 0.2)
box()
text(length(gsea.rnk$metric)/2, 0.1, labels = ifelse(!missing(class.name), class.name, ""), cex=1)
text(length(gsea.rnk$metric)*0.03, 0.1, redgroup.lab, adj = c(0, NA), cex=1.5) # 2번째가 y축 위치 조절.
text(length(gsea.rnk$metric)*0.97, 0.1, bluegroup.lab, adj = c(1, NA),cex=1.5)
#ln3=c(seq(0,xmax, ifelse(xmax>15000,5000,5000))) 
ln3=seq(0,length(rk),5000)
axis(1,at = ln3 ,labels= rep("",length(ln3)), cex.axis=0.7, pos=-0.05,tck=-0.08,xpd=T) # 컬러바 아래 눈금 위치 pos 조절
axis(1,at = ln3, labels=ln3, cex.axis=1, pos=0.05, tick=FALSE, xpd=T ) # 컬러바 아래 숫자 위치 pos 조절
axis(1,at = xmax/2 ,labels=  "Gene rank", pos=-0.05,tick=FALSE, cex.axis=1) # x축 위치 조절 pos 





#############################################################

load(file = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/bulkRNA/Rdata/1_annotation.Rdata") # ginfo, bm

genes = unique(unlist(mootha[7,]$leadingEdge))
m.deg[m.deg$Genes %in% genes,]

symbol = ginfo$gene_name[match(rownames(m.tpm), ginfo$gene_id)]
rownames(m.tpm) = ifelse(is.na(symbol), rownames(m.tpm), symbol)
m.tpm = m.tpm[,-1]

tpm = m.tpm[rownames(m.tpm) %in% genes, ]
tpm = log2(tpm+1)
ztpm = data.frame(t(apply(tpm, 1, function(x) scale(x, center=T, scale=T))))
colnames(ztpm) = colnames(tpm)
ztpm = ztpm[!is.na(ztpm$Cis2),]


sp = rep(c("1", "2"), times=c(3,4))  
color.ht=colorRamp2(c(-2,-1,0,1,2), c('midnightblue', 'mediumblue','white', 'firebrick', 'darkred'))

tiff(filename = sprintf("%s/figure/3_mouse_Ztpm.tiff", dir), width = 15, height = 12, units = 'cm', res = 300)
Heatmap(t(ztpm),  
        show_heatmap_legend = T,
        col = color.ht, border = T,
        cluster_columns = T, cluster_rows = F, 
        show_column_names = T, column_names_rot = 90, 
        row_names_gp = gpar(fontsize = 15, col = "black"), column_names_gp = gpar(fontsize = 18, fontface = "italic"),
        width = nrow(ztpm) * unit(0.6, "cm"), height = ncol(ztpm) * unit(0.6, "cm"),
        row_split = sp, row_gap = unit(0.15, "cm"))
dev.off()



##################
# Figure 6L
##################


library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(randomcoloR)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggpubr)



dir_data = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO/GSE197266_FA_mouse_scNRA"

load(file=file.path(dir_data, "Rdata/1_sobj.Rdata")) #s.merged, s.integrated

levels(s.integrated$condition)[
  levels(s.integrated$condition) == "CP"
] = "CIS"

mycol = c("purple","mediumpurple1","deeppink","thistle","chartreuse1","cyan","cadetblue","yellowgreen","royalblue","slategray2")
names(mycol) = c("Podocyte","PT","ATL/TAL/DCT/CNT/CD-PC","CD-IC","T cell","B cell","Neutrophil","Macrophage","Endothelia","Pericyte")
condcol = c("Ctrl"="slategrey","CIS"="gold")

s.sub = subset(
  s.integrated,
  subset = condition %in% c("Ctrl", "CIS")
)

table(s.sub$condition)
Idents(s.sub)


tiff(filename = file.path(dir_data, "/figure/1_CIS_dimplot_celltype.tiff"), width = 15, height = 15, units = 'cm', res = 300)
DimPlot(s.sub, reduction='umap', group.by='cellType', cols = mycol) + ggtitle('') + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
dev.off()




##################
# Figure 6M
##################

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)
library(dplyr)
library(SingleCellExperiment)
library(randomcoloR)
library(GOfuncR)
library(stringr)
library(homologene)
library(fgsea)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggpubr)



dir_data = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO/GSE197266_FA_mouse_scNRA"

load(file=file.path(dir_data, "Rdata/1_sobj.Rdata")) #s.merged, s.integrated

levels(s.integrated$condition)[
  levels(s.integrated$condition) == "CP"
] = "CIS"

mycol = c("purple","mediumpurple1","deeppink","thistle","chartreuse1","cyan","cadetblue","yellowgreen","royalblue","slategray2")
names(mycol) = c("Podocyte","PT","ATL/TAL/DCT/CNT/CD-PC","CD-IC","T cell","B cell","Neutrophil","Macrophage","Endothelia","Pericyte")
condcol = c("Ctrl"="slategrey","CIS"="gold")

s.sub = subset(
  s.integrated,
  subset = condition %in% c("Ctrl", "CIS")
)




c2 = read.gmt("C:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/c2.all.v2024.1.Hs.symbols.gmt")
c2 = split(x = c2$gene, f = c2$term)
names(c2)[grep('GLUCONEO', names(c2))]
gluco_genes  = c2[["MOOTHA_GLUCONEOGENESIS"]]
# human SYMBOL → mouse Symbol (첫 글자만 대문자)
gluco_mouse  = paste0(toupper(substr(gluco_genes, 1, 1)),
                      tolower(substr(gluco_genes, 2, nchar(gluco_genes))))

DefaultAssay(s.sub) = "RNA"
available    = rownames(s.sub)
gluco_keep   = gluco_mouse[gluco_mouse %in% available]
s.sub = JoinLayers(s.sub, assay="RNA")
s.sub = AddModuleScore(s.sub, features = list(gluconeogenesis = gluco_keep), name = "gluconeogenesis", assay = "RNA")

names(c2)[grep('SENESCENCE', names(c2))]
sene_genes  = c2[["KAMMINGA_SENESCENCE"]]
sene_mouse  = paste0(toupper(substr(sene_genes, 1, 1)),
                     tolower(substr(sene_genes, 2, nchar(sene_genes))))

DefaultAssay(s.sub) = "RNA"
available    = rownames(s.sub)
sene_keep   = sene_mouse[sene_mouse %in% available]
s.sub = JoinLayers(s.sub, assay="RNA")
s.sub = AddModuleScore(s.sub, features = list(senescence = sene_keep), name = "senescence", assay = "RNA")

names(c2)[grep('SENESCENCE', names(c2))]
sene_genes  = c2[["KAMMINGA_SENESCENCE"]]
sene_mouse  = paste0(toupper(substr(sene_genes, 1, 1)),
                     tolower(substr(sene_genes, 2, nchar(sene_genes))))

s.pt = subset(s.sub, subset = cellType == "PT")
s.pt = JoinLayers(s.pt, assay="RNA")
s.pt = AddModuleScore(s.pt, features = list(gluconeogenesis = gluco_keep), name = "gluconeogenesis", assay = "RNA")
colnames(s.pt@meta.data)[grep("^gluconeogenesis", colnames(s.pt@meta.data))] = "gluconeogenesis_score"

condcol = c("Ctrl"="slategrey","CIS"="gold")
s.pt$condition = factor(s.pt$condition, levels=c("Ctrl","CIS"))

tiff(filename = file.path(dir_data, "figure/2_CIS_PT_gluconeogenesis_score.tif"), width = 7, height = 7, units = 'cm', res = 300)
VlnPlot(s.pt, pt.size = 0, features = "gluconeogenesis_score", group.by = "condition", split.by = 'condition', assay = 'RNA', cols = condcol) + 
  #scale_y_continuous(limits = c(-0.38, 0.8)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, colour = "black", angle = 0, hjust = 0.5),
        plot.title = element_blank(),
        legend.position = "right") +
  geom_boxplot(width=0.3, position = position_dodge(0.9), outlier.shape = NA, linewidth = 0.3)
#stat_compare_means(method = "wilcox.test", na.rm = TRUE)
dev.off()


s.pt = AddModuleScore(s.pt, features = list(senescence = sene_keep), name = "senescence", assay = "RNA")
colnames(s.pt@meta.data)[grep("^senescence", colnames(s.pt@meta.data))] = "senescence_score"

tiff(filename = file.path(dir_data, "figure/2_CIS_PT_senescence_score.tif"), width = 7, height = 7, units = 'cm', res = 300)
VlnPlot(s.pt, pt.size = 0, features = "senescence_score", group.by = "condition", split.by = 'condition', assay = 'RNA', cols = condcol) + 
  #scale_y_continuous(limits = c(-0.38, 0.8)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, colour = "black", angle = 0, hjust = 0.5),
        plot.title = element_blank(),
        legend.position = "right") +
  geom_boxplot(width=0.3, position = position_dodge(0.9), outlier.shape = NA, linewidth = 0.3)
#stat_compare_means(method = "wilcox.test", na.rm = TRUE)
dev.off()



