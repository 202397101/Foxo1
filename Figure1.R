##################
# Figure 1B
##################
dir="E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl
load(file = sprintf("%s/bulkRNA/Rdata/9_gse_GO.Rdata", dir)) #gsego
load(file = sprintf("%s/bulkRNA/Rdata/Figure1B.Rdata", dir)) #stv.gsego, total.lipid, total.glc, total.nucleotide

library(stringr)

# total.lipid 34
# total.glc 11
# total.nucleotide 7

#reg = total.lipid[grep("regulation", total.lipid$Description), "Description"] #14
#total.lipid1 = total.lipid[!total.lipid$Description %in% reg, ] #20 x 11

#total1 = list(total.lipid1, total.glc, total.nucleotide)
#names(total1) = c("Lipid & Fatty acid", "Carbohydrate", "Nucleotide")

#total2 = c()
#i=1
#for (i in 1:length(total1)){
#  df = total1[[i]]
#  df$Description = gsub("metabolic process", "metabolism", df$Description)
#  df$Description = gsub("catabolic process", "catabolism", df$Description)
#  df$Description = gsub("biosynthetic process", "biosynthesis", df$Description)
#  df$score = -log10(df$pvalue) * sign(df$NES)
#  df1 = df[df$score > 0,]
#  df2 = df1[order(df1$score, decreasing = T),]
#  total2[[names(total1[i])]] = df2
#}

#save(total2, file = sprintf("%s/Figure1B.Rdata", dir))

total2 = lapply(total2, function(i) {
  i$pvalue = -log10(i$pvalue)
  colnames(i)[colnames(i) == "pvalue"] = "logPval"
  ora = i[order(i$logPval),]
  return(ora)
})

par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(total2$Carbohydrate$logPval, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', xlab = "", names.arg = NA, width = 0.7, border = NA, col = "lightskyblue", main = " ", cex.main = 1.5)
axis(1, at=seq(0,5,1), labels = seq(0,5,1), cex.axis=1, las=1)
abline(v=seq(1,5,1), lty=3, col= "dimgrey")
bp = barplot(tail(total2$Carbohydrate$logPval, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', xlab = "", names.arg = NA, width = 0.7, border = NA, col = "lightskyblue", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels = str_to_sentence(tail(total2$Carbohydrate$Description, 5)), col = "black", xpd=T, cex=1.3, adj=0)


par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(total2$`Lipid & Fatty acid`$logPval, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', xlab = "", names.arg = NA, width = 0.7, border = NA, col = "khaki1", main = " ", cex.main = 1.5)
axis(1, at=seq(0,5,1), labels = seq(0,5,1), cex.axis=1, las=1)
abline(v=seq(1,5,1), lty=3, col= "dimgrey")
bp = barplot(tail(total2$`Lipid & Fatty acid`$logPval, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', xlab = "", names.arg = NA, width = 0.7, border = NA, col = "khaki1", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels = str_to_sentence(tail(total2$`Lipid & Fatty acid`$Description, 5)), col = "black", xpd=T, cex=1.3, adj=0)


par(mai = c(0.5,0.1,0.1,1))
bp = barplot(tail(total2$Nucleotide$logPval, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', xlab = "", names.arg = NA, width = 0.7, border = NA, col = "pink", main = " ", cex.main = 1.5)
axis(1, at=seq(0,5,1), labels = seq(0,5,1), cex.axis=1, las=1)
abline(v=seq(1,5,1), lty=3, col= "dimgrey")
bp = barplot(tail(total2$Nucleotide$logPval, 5), xlim = c(0, 5), horiz = T, xaxt = 'n', yaxt = 'n', xlab = "", names.arg = NA, width = 0.7, border = NA, col = "pink", main = " ", cex.main = 1.5, add = T)
abline(v=0, lty=1)
text(x=0.2, y=bp ,labels = str_to_sentence(tail(total2$Nucleotide$Description, 5)), col = "black", xpd=T, cex=1.3, adj=0)


##################
# Figure 1C
##################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
library(dplyr)
library(stringr)
library(ggfortify)
library(DESeq2)
library(edgeR)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(GOfuncR)
library(biomaRt)
#library(devEMF)
library(igraph)
library(RCy3)


mycol = c("purple","mediumpurple1","violetred","pink","deeppink","palevioletred","thistle","violet","chartreuse1","cyan","cadetblue","yellowgreen","royalblue")
grcol = c('steelblue','darkorange')

tiff(filename = sprintf("%s/figure/Figure1C_celltype.tiff", dir), width = 15, height = 15, units = 'cm', res = 300)
DimPlot(s.integrated, group.by = 'cellType', label=F, cols = mycol) + ggtitle('')+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
dev.off()




##################
# Figure 1D
##################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/scRNA/Rdata/2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
library(dplyr)
library(stringr)
library(ggfortify)
library(DESeq2)
library(edgeR)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(GOfuncR)
library(biomaRt)
#library(devEMF)
library(igraph)
library(RCy3)

pop= s.integrated@meta.data %>% 
  group_by(cellType) %>% 
  summarise(counts = n()) %>% 
  mutate(total.counts = sum(counts)) %>% 
  mutate(percent = counts/total.counts*100) 

celltypes = unique(s.integrated$cellType)
mycol = c(
  "PT"          = "mediumpurple1",
  "DCT/CNT"     = "pink",
  "Macrophage"  = "yellowgreen",
  "CD-ICB"      = "violet",
  "CNT/CD-PC"   = "deeppink",
  "CD-PC"       = "palevioletred",
  "ATL/TAL"     = "violetred",
  "Endothelia"  = "royalblue",
  "T cell"      = "chartreuse1",
  "CD-ICA"      = "thistle",
  "B cell"      = "mediumorchid1",
  "Podocyte"    = "cyan",
  "Neutrophil"  = "cadetblue"
)

mycol = mycol[celltypes]
grcol = c('steelblue','darkorange')

tiff(filename = sprintf("%s/figure/Figure1D.tiff", dir), width = 8, height = 10, units = 'cm', res = 300)
ggplot(pop, aes(x = reorder(cellType, percent), y = percent, fill = cellType, width = 0.8)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(pop$percent, 1), "%")), position = position_stack(vjust = 1), size = 3) +
  scale_fill_manual(values = mycol) + 
  labs(y = "cell counts (%)d") + 
  geom_vline(xintercept = 0.4, color = "black", size=0.5) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(color = "black", size=10), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x =  element_blank(), 
        axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        legend.position = "none")+
  coord_flip()
dev.off()



##################
# Figure 1E
##################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/scRNA/Rdata/5-2_PC_vs_others_DEG.Rdata", dir)) #markers.li
load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl

library('Seurat')
library('tidyverse')
library('ggplot2')
library(stringr)
library(dplyr)
library(ggplot2)
library(dbplyr)
library(tidyr)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(ggrepel)
library(gridExtra)
library(RColorBrewer)
#install.packages("ComplexUpset")
library(ComplexUpset)
library(ComplexHeatmap)

bulk = ddfl$`Starvation  vs  Control`
bulk.updeg = bulk[bulk$log2FoldChange > 2 & bulk$padj < 0.01 & !is.na(bulk$padj), ]

stv.li = markers.li$stv
scRNA.updeg = list()
i=1
for (i in 1:length(stv.li)){
  deg = stv.li[[i]]  
  updeg = deg[deg$avg_log2FC > 2 & deg$p_val_adj < 0.01 & !is.na(deg$p_val_adj), ]
  scRNA.updeg[[names(stv.li)[i]]] = updeg
}

gl = list()
i=1
for (i in 1:length(scRNA.updeg)){
  celltype.name = names(scRNA.updeg)[i]
  celltype.gene = rownames(scRNA.updeg[[i]])
  gl[[names(scRNA.updeg)[i]]] = celltype.gene
}
gl[["Bulk"]] = bulk.updeg$Genes

Bulk = bulk.updeg$Genes
composite_gl = list()
for (cell in names(gl)) {  
  if (cell != "Bulk") {
    composite_gl[[cell]] = unique(intersect(Bulk, gl[[cell]]))  
  }
}
composite_gl[["Bulk"]] = Bulk
binary_matrix_composite = make_comb_mat(composite_gl)
comb_size(binary_matrix_composite)
comb_degrees = comb_degree(binary_matrix_composite)
filtered_binary_matrix = binary_matrix_composite[comb_degrees == 2]

UpSet(filtered_binary_matrix)


###### fisher test - upDEGs

fisher.pval = list()
i=1
background = unique(c(bulk.updeg$Genes, unlist(lapply(scRNA.updeg, rownames))))

for (i in 1:length(scRNA.updeg)){
  celltype = scRNA.updeg[[i]]
  
  a = length(intersect(bulk.updeg$Genes, rownames(celltype)))
  b = length(setdiff(rownames(celltype), bulk.updeg$Genes)) #scrna만 up.
  c = length(setdiff(bulk.updeg$Genes, rownames(celltype))) #bluk만 up.
  d = length(setdiff(background, c(rownames(celltype), bulk.updeg$Genes)))
  
  contingency.table = matrix(c(a, b, c, d), nrow = 2, dimnames = list(
    "scRNA_DEGs" = c("Yes", "No"),
    "Bulk_DEGs" = c("Yes", "No")
  ))
  fisher.res = fisher.test(contingency.table, alternative = "greater")
  fisher.pval[[names(scRNA.updeg)[i]]] = fisher.res$p.value
}
unlist(fisher.pval)


######## plot
pval = unlist(fisher.pval)
log_pval = -log10(pval)

binary_to_rownames = apply(filtered_binary_matrix, 2, function(x) {
  first_match = rownames(filtered_binary_matrix)[which(x == 1)][1] 
})
binary_mapping = data.frame(
  Binary =  names(comb_size(filtered_binary_matrix)),  # Binary 코드
  Rowname = binary_to_rownames  # 실제 rownames 
)

ordered_logp = log_pval[match(binary_mapping$Rowname, names(log_pval))]
mycol = colorRamp2(c(0, 1, max(log_pval, na.rm = TRUE)), c("grey", "yellow","red"))

UpSet(
  filtered_binary_matrix,
  bottom_annotation = HeatmapAnnotation(
    '-logPvalue' = ordered_logp, 
    col = list('-logPvalue' = mycol),  
    border = TRUE, 
    height = unit(1, "cm"),
    gp = gpar(col = "black", lwd=1),
    annotation_name_gp = gpar(fontsize = 0)
  )
)





##################
# Figure 1F
##################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/Figure1B.Rdata", dir)) #total2


library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
library(dplyr)
library(stringr)
library(ggfortify)
library(DESeq2)
library(edgeR)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(GOfuncR)
library(biomaRt)
#library(devEMF)
library(igraph)
library(RCy3)
library(data.table)
library(ggpubr)



gluco = unique(unlist(total2$Carbohydrate$core_enrichment[4]))
gluco = unlist(strsplit(gluco, "/"))

boxi = unique(unlist(total2$`Lipid & Fatty acid`$core_enrichment[5]))
boxi = unlist(strsplit(boxi, "/"))

rntp = unique(unlist(total2$Nucleotide$core_enrichment[1]))
rntp = unlist(strsplit(rntp, "/"))

gluconeo = unique(unlist(total2$Carbohydrate$core_enrichment[5]))
gluconeo = unlist(strsplit(gluconeo, "/"))

inpGS = list(glucose=gluco, betaOxidation=boxi, ribonucleotide=rntp, gluconeogenesis=gluconeo)
names(inpGS)

DefaultAssay(s.integrated) = 'RNA'
genes = rownames(s.integrated)
inpGS.filtered = lapply(inpGS, function(i) {
  intersect(i, genes)
})
new = JoinLayers(s.integrated, overwrite = T)
new = AddModuleScore(new, features = inpGS.filtered, name = "signaling.path", assay = 'RNA')
colnames(new@meta.data)[grep("signaling.path", colnames(new@meta.data))] = names(inpGS)
new.pc = subset(new, subset = cellType == "PT")

grcol = c('steelblue','darkorange')

DefaultAssay(new.pc)

tiff(filename = sprintf("%s/figure/figure1D_glucose.tiff", dir), width = 6, height = 8, units = 'cm', res = 300)
VlnPlot(new.pc, pt.size = 0, features = names(inpGS)[1], group.by = "cellType", split.by = 'sample', assay = 'RNA', cols = grcol) + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none") +
  #stat_compare_means(method = "wilcox.test", na.rm = TRUE)+
  geom_boxplot(width=0.2, position = position_dodge(0.9), outlier.shape = NA)
dev.off()


tiff(filename = sprintf("%s/figure/figure1D_betaOxidation.tiff", dir), width = 6, height = 8, units = 'cm', res = 300)
VlnPlot(new.pc, pt.size = 0, features = names(inpGS)[2], group.by = "cellType",split.by = 'sample', assay = 'RNA', cols = grcol) + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none") +
  #stat_compare_means(method = "wilcox.test", na.rm = TRUE)+
  geom_boxplot(width=0.2, position = position_dodge(0.9), outlier.shape = NA)
  
dev.off()


tiff(filename = sprintf("%s/figure/figure1D_rNTPmetabolism.tiff", dir), width = 6, height = 8, units = 'cm', res = 300)
VlnPlot(new.pc, pt.size = 0, features = names(inpGS)[3], group.by = "cellType",split.by = 'sample', assay = 'RNA', cols = grcol) + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none") +
  #stat_compare_means(method = "wilcox.test", na.rm = TRUE)+
  geom_boxplot(width=0.2, position = position_dodge(0.9), outlier.shape = NA)
dev.off()

tiff(filename = sprintf("%s/figure/figure1D_gluconeogenesis.tiff", dir), width = 6, height = 8, units = 'cm', res = 300)
VlnPlot(new.pc, pt.size = 0, features = names(inpGS)[4], group.by = "cellType", split.by = 'sample', assay = 'RNA', cols = grcol) + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none") +
  stat_compare_means(method = "wilcox.test", na.rm = TRUE)+
  geom_boxplot(width=0.2, position = position_dodge(0.9), outlier.shape = NA)
dev.off()


##################
# Figure 1G
##################
dir = "E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_string = "E:/Dropbox/PNU/시스템생물학연구실/DB/string"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/scRNA/Rdata/5-1_PTcell.Rdata", dir)) #pt.deg
load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
library(dplyr)
library(stringr)
library(ggfortify)
library(DESeq2)
library(edgeR)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(GOfuncR)
library(biomaRt)
#library(devEMF)
library(igraph)
library(RCy3)
library(ggvenn)
library(RCy3)


updeg.pt = pt.deg[pt.deg$avg_log2FC > 2 & pt.deg$p_val_adj < 0.01 & !is.na(pt.deg$p_val_adj),]

updeg.bulk = ddfl$`Starvation  vs  Control`
updeg.bulk = updeg.bulk[updeg.bulk$log2FoldChange > 2 & updeg.bulk$padj < 0.01 & !is.na(updeg.bulk$padj), ]

intersect.updeg = unique(intersect(rownames(updeg.pt), updeg.bulk$Genes)) #24
inter.li = list(scRNA = rownames(updeg.pt), bulkRNA = updeg.bulk$Genes)

tiff(filename = sprintf("%s/figure/figure1F_venn.tiff", dir), width = 10, height = 10, units = 'cm', res = 300)
ggvenn(inter.li, show_percentage = F, 
       fill_color = c("lightgoldenrod1","turquoise"), fill_alpha = 0.6, stroke_color = "black", 
       stroke_size = 0.3, set_name_size = 3, text_color = "black", text_size = 9) 
dev.off()



##################
# Figure 1I
##################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers

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


DefaultAssay(s.integrated) = "RNA"
gene1 = c('Pck1')
s.integrated@meta.data$sample = factor(s.integrated@meta.data$sample, levels = c("con", "stv"))

tiff(filename = sprintf("%s/figure/Figure1G.tiff", dir), width = 20, height = 8.5, units = 'cm', res = 300)
FeaturePlot(s.integrated, features = gene1, order = TRUE, split.by = 'sample', min.cutoff = 0, max.cutoff = 'q90') & theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()





##################
# Figure 1H
##################
dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_string = "C:/Dropbox/PNU/시스템생물학연구실/DB/string"

load(file = sprintf("%s/scRNA/Rdata/3-2_integrated_seurat_object.Rdata", dir)) #s.integrated, markers
load(file = sprintf("%s/scRNA/Rdata/5-1_PTcell.Rdata", dir)) #pt.deg
load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl

library('Seurat')
library('tidyverse')
library('ggplot2')
library('patchwork')
library(stringr)
library(dplyr)
library(ggplot2)
#remotes::install_version("Matrix", "1.6.1")
library(Matrix)
library(dbplyr)
library(dittoSeq)
library(colorRamp2)
library(viridis)
library(randomcoloR)
library(GPTCelltype)
library(openai)
library(cowplot)
library(dplyr)
library(stringr)
library(ggfortify)
library(DESeq2)
library(edgeR)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(gridExtra)
library(stringr)
library(ggvenn)
library(msigdbr)
library(fgsea)
library(GOfuncR)
library(biomaRt)
#library(devEMF)
library(igraph)
library(RCy3)
library(ggvenn)
library(RCy3)


updeg.pt = pt.deg[pt.deg$avg_log2FC > 2 & pt.deg$p_val_adj < 0.01 & !is.na(pt.deg$p_val_adj),]

updeg.bulk = ddfl$`Starvation  vs  Control`
updeg.bulk = updeg.bulk[updeg.bulk$log2FoldChange > 2 & updeg.bulk$padj < 0.01 & !is.na(updeg.bulk$padj), ]

intersect.updeg = unique(intersect(rownames(updeg.pt), updeg.bulk$Genes)) #24
# PPI network centrality
ppi = read.csv(file = sprintf("%s/mus_musculus/10090.protein.links.v12.0.txt", dir_string), sep = " ", header = T, stringsAsFactors = F, quote = "")
pinfo = read.csv(file = sprintf("%s/mus_musculus/10090.protein.info.v12.0.txt", dir_string), sep = "\t", header = T, stringsAsFactors = F, quote = "")
ppi$protein1 = pinfo$preferred_name[match(ppi$protein1, pinfo$X.string_protein_id)]
ppi$protein2 = pinfo$preferred_name[match(ppi$protein2, pinfo$X.string_protein_id)]
ppi1 = ppi[ppi$protein1 %in% intersect.updeg & ppi$protein2 %in% intersect.updeg & ppi$combined_score>500, ]
ppi1[grep('Pck1', ppi1$protein1),]

graph = graph_from_edgelist(as.matrix(ppi1[,1:2]), directed = F)
graph = igraph::simplify(graph)

V(graph)$vcount=vcount(graph)
V(graph)$degree = igraph::degree(graph)                 
V(graph)$eig = evcent(graph)$vector            
V(graph)$hubs = hub.score(graph)$vector            
V(graph)$authorities = authority.score(graph)$vector   
V(graph)$closeness = closeness(graph)               
V(graph)$betweenness = betweenness(graph)  
centrality = data.frame(row.names   = V(graph)$name,
                        betweenness = V(graph)$betweenness,
                        degree      = V(graph)$degree,
                        closeness   = V(graph)$closeness,
                        eigenvector = V(graph)$eig)
centrality = centrality[order(centrality$betweenness, decreasing = T),]

V(graph)$sc_log2FC = pt.deg$avg_log2FC[match(names(V(graph)), rownames(pt.deg))]
V(graph)$bulk_log2FC = updeg.bulk$log2FoldChange[match(names(V(graph)), updeg.bulk$Genes)]
range(c(V(graph)$bulk_log2FC, V(graph)$sc_log2FC))

createNetworkFromIgraph(graph, "Figure1F")


top8 = head(centrality, 8)
ggplot(top8, aes(x = reorder(rownames(top8), -betweenness), y = betweenness, fill = betweenness)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient2(low = "rosybrown1", mid = "violetred",high = "darkmagenta", midpoint = 20) + 
  labs(x = "Gene", y = "Betweenness") +  
  theme_minimal() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = "grey", size = 0.3),
        panel.grid.minor.y = element_line(color = "grey", size = 0.3),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none")




##################
# Figure 1J
##################
dir="C:/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl

library(ComplexHeatmap)
library(colorRamp2)


stv.deg = ddfl$`Starvation  vs  Control`

pcx = stv.deg[stv.deg$Genes == "Pcx",]
pck1 = stv.deg[stv.deg$Genes == "Pck1",]
eno1 = stv.deg[stv.deg$Genes == "Eno1",]
pgam1 = stv.deg[stv.deg$Genes == "Pgam1",]
pgk1 = stv.deg[stv.deg$Genes == "Pgk1",]
gapdh = stv.deg[stv.deg$Genes == "Gapdh",]
tpi1 = stv.deg[stv.deg$Genes == "Tpi1",]
aldoa = stv.deg[stv.deg$Genes == "Aldoa",]
fbp1 = stv.deg[stv.deg$Genes == "Fbp1",]
gpi1 = stv.deg[stv.deg$Genes == "Gpi1",]
g6p = stv.deg[stv.deg$Genes == "G6pc",]


stv.glu.deg = rbind(pcx, pck1, eno1, pgam1, pgk1, gapdh, tpi1, aldoa, fbp1, gpi1, g6p)

stv.glu.deg.m = stv.glu.deg[,3]

stv.glu.deg.m = as.matrix(stv.glu.deg.m)
rownames(stv.glu.deg.m) = rownames(stv.glu.deg)

color.ht=colorRamp2(c(-3,-2,-1,0,1,2,3), c('blue2', 'blue1','dodgerblue1','white','salmon1', 'red1', 'red2'))
lg.ht=Legend(title="", at=c(-3,-2,-1,0,1,2,3), col_fun=color.ht, border='black', title_position="topcenter", title_gp = gpar(fontsize = 8, fontface = "bold"), legend_height = unit(5, "cm"), grid_width = unit(1, "cm"), labels_gp = gpar(fontsize = 15))

hm = Heatmap(stv.glu.deg.m, show_heatmap_legend=F, col=color.ht, border=T, 
             cluster_columns=T, cluster_rows=F, 
             show_column_names=T, cluster_column_slices=T, column_names_rot=90, 
             width=ncol(stv.glu.deg.m)*unit(1,"cm"), height=nrow(stv.glu.deg.m)*unit(1,"cm"),
             left_annotation = rowAnnotation(foo = anno_text((stv.glu.deg$pvalue), location = 1, just = "right", gp = gpar(fontsize = 10))))


draw(hm, annotation_legend_list=packLegend(list=list(lg.ht)), merge_legend=T, heatmap_legend_side="left", annotation_legend_side='right', padding=unit(c(0,0,0,0), "cm"))
