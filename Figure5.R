##############################
# Figure 5G
##############################
dir="E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"
dir_msigdb_c5="E:/PNU/시스템생물학연구실/data/gluconeogenesis/results/msigdb_c5"

load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl
load(file = sprintf("%s/bulkRNA/Rdata/9_gse_GO.Rdata", dir)) #gsego
source("E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/bulkRNA/Rscripts/11-2_GSEAplot.R")

library(fgsea)
library(msigdbr)
library(qusage)
library(org.Mm.eg.db)
library(GOfuncR)


fastv.ddfl = ddfl$`FA_STV  vs  Starvation`
fastv.gsego = gsego$`FA_STV  vs  Starvation`

fastv.gsego.gluconeo = fastv.gsego[fastv.gsego@result$ID == "GO:0006094"]
gluconeogenesis = get_anno_genes('GO:0006094', database = "org.Mm.eg.db") #94
gluconeogenesis.genes = gluconeogenesis$gene

gene.sets = as.data.frame(msigdbr(species = "Mus musculus", category = "C5")) #1250307
gluco.gsets = gene.sets[grep("GO:0006111", gene.sets$gs_exact_source), c("gs_name", "gene_symbol", "gs_exact_source")] #82 x 3
gluco.gsets1 = split(x=gluco.gsets$gene_symbol, f=gluco.gsets$gs_name)
gluco.gene = unique(unlist(gluco.gsets1))

idx = !is.na(fastv.ddfl$stat)
use = fastv.ddfl[idx, ]
rk = use$stat 
names(rk) = use$Genes
rk = sort(rk, decreasing = T) #29378

gset = gluco.gsets1
fname = sprintf("%s/figure/Figure6G_fastv_gluconeogeneis.tif", dir)
labs = list(mt="Gluconeogenesis", redgroup.lab="Starvation", bluegroup.lab="FA+Starvation", mlab="")
gseaPlot(fname, rk, gset, labs, xmax = 30000)


##############################
# Figure 5H
##############################
dir="E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/bulkRNA/Rdata/5_DEGs.Rdata", dir)) #ddf, ddfl
load(file = sprintf("%s/bulkRNA/Rdata/2_expression_matrix.Rdata", dir)) #ginfo, sinfo, tpm, rcm, ercm, etpm
load(file = sprintf("%s/bulkRNA/Rdata/9_gse_GO.Rdata", dir)) #gsego


library(stringr)
library(ComplexHeatmap)
library(colorRamp2)
library(gridExtra)
library(ggplot2)
library(GOfuncR)
library(org.Mm.eg.db)

gluconeogenesis = get_anno_genes('GO:0006094', database = "org.Mm.eg.db") 
gluconeogenesis.genes = gluconeogenesis$gene #94

ddfl1 = c()
i=1
for (i in 1:length(ddfl)){
  df = ddfl[[i]]
  df = df[match(gluconeogenesis.genes, df$Genes), ]
  df1 = df[!is.na(df$padj) & df$padj<0.01 & df$log2FoldChange < -0.5, ]
  ddfl1[[names(ddfl[i])]] = df1
}

glu.tpm = c()
i=1
for (i in 1:length(ddfl1)) {
  df = ddfl1[[i]]
  glu.g = df$Genes
  tpm1 = tpm[which(rownames(tpm) %in% glu.g), ]
  tpm1 = log2(tpm1+1)
  glu.tpm[[names(ddfl1[i])]] = tpm1
}

fastv.glu.tpm = glu.tpm[[2]][,c(7,8,9,10,11,12)]
fastv.glu.z = t(apply(fastv.glu.tpm, 1 ,function(x) scale(x, center = T, scale = T)))
colnames(fastv.glu.z) = colnames(fastv.glu.tpm)

color.ht=colorRamp2(c(-2,-1,0,1,2), c('midnightblue', 'mediumblue','white', 'firebrick', 'darkred'))
lg.ht=Legend(title="Tansformed\nZ-scale", at=c(-2,0,2), col_fun=color.ht, border='black', title_position="topcenter", title_gp = gpar(fontsize = 8, fontface = "bold"), legend_height = unit(1.5, "cm"), grid_width = unit(0.3, "cm"))
sp = rep(c(" ", ""),c(3,3))

hm = Heatmap(
  t(fastv.glu.z), 
  show_heatmap_legend = F,
  col = color.ht,
  border = T, 
  cluster_columns = T,
  cluster_rows = F,
  show_column_names = T,
  cluster_column_slices = T, 
  column_names_rot = 90,
  column_title = "Gluconeogenesis", 
  width = ncol(fastv.glu.z) * unit(2.5, "cm"),
  height = nrow(fastv.glu.z) * unit(0.1, "cm"), 
  row_names_gp = grid::gpar(fontsize = 9),
  row_split = sp,
  row_gap = unit(0.15, "cm"),
  column_names_gp = grid::gpar(
    col = "black",
    fontsize = 12,
    fontface = "italic"
  )
)

hm.plot = draw(
  hm,
  annotation_legend_list = packLegend(list = list(lg.ht)),
  merge_legend = T,
  heatmap_legend_side = "left",
  annotation_legend_side = "right",
  padding = unit(c(0,0,0,0), "cm")
)
tiff(filename = sprintf("%s/figure/supple_figure9C-2.tiff", dir), width = 22, height = 8, units = 'cm', res = 300)
hm.plot
dev.off()


##############################
# Figure 5I
##############################

dir="E:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis"

load(file = sprintf("%s/bulkRNA/Rdata/2_expression_matrix.Rdata", dir)) #ginfo, sinfo, tpm, rcm, ercm, etpm

#BiocManager::install('GSVA')
library(GSVA)
library(GOfuncR)
library(msigdbr)
library(ggplot2)
#packageVersion("GSVA") 1.50.5

tpm1 = tpm[,c(1,2,3,10,11,12,7,8,9)]
tpm2 = as.matrix(log2(tpm1+1))


gluconeogenesis = get_anno_genes('GO:0006094', database = "org.Mm.eg.db") #94
gluconeogenesis.genes = gluconeogenesis$gene

gset.li = list(gluconeo = gluconeogenesis.genes)
gp = gsvaParam(tpm2, gset.li)
ssgsea = as.data.frame(t(gsva(gp, verbose=TRUE)))
ssgsea[,"sample"] = "1"
ssgsea$sample[grep("Starvation", rownames(ssgsea))] = "2"
ssgsea$sample[grep("FA_STV", rownames(ssgsea))] = "3"

tiff(filename = sprintf("%s/figure/Figure3I.tiff", dir), width = 8, height = 10, units = 'cm', res = 300)
par(mar=c(5, 2, 1.5, 0.8), cex.main=0.8, cex.axis=0.5, cex.lab=0.6, tck = -0.02)
boxplot(gluconeo ~ sample, data = ssgsea, col = "white", xlab="", ylab="ssGSSEA score", lwd = 1.2, main="Gluconeogenesis", boxwex=0.6, xaxt = "n", yaxt = "n", horizontal=F, outline=F, ylim = c(-0.5, 0.5))
stripchart(gluconeo ~ sample, vertical = T, data = ssgsea, method = "jitter", jitter=0.23, add = T, pch = 19, lwd=1.8, cex = 1.2, col = c("lightblue", "blue", "red"))
grid(nx=NA, ny=NULL, lty=3, lwd=1.5, col="lightgrey")
axis(side = 2, at=seq(-0.5, 0.5, 0.5), labels=seq(-0.5, 0.5, 0.5), cex.axis=0.8, las=1, mgp=c(1, 0.5, 0), tck = -0.02) 
axis(side = 1, at=1:4, labels=c("CON", "STARVATION", "FA+STV"), cex.axis=0.5, mgp=c(1, 0.5, 0))
dev.off()

tiff(filename = sprintf("%s/figure/Figure3I_legend.tiff", dir), width = 15, height = 4, units = 'cm', res = 300)
par(mar=c(5,0,0,0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("bottom", legend=c("CON", "Starvation", "FA+Starvation"), col= c("lightblue", "blue", "red"), pch=16,pt.cex=1.5, horiz=T, x.intersp=0.5) 
dev.off()

ssgsea[,"sample"] = "1CON"
ssgsea$sample[grep("Starvation", rownames(ssgsea))] = "2Starvation"
ssgsea$sample[grep("FA_STV", rownames(ssgsea))] = "3FA+Starvation"
g.anova = aov(ssgsea$gluconeo~ssgsea$sample, data = ssgsea)
summary(g.anova)
TukeyHSD(g.anova)



