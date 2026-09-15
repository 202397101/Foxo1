
##################
# Figure 9P,Q
##################
library(clusterProfiler)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ppcor)
library(pheatmap)
library(ggfortify)



# ── Paths ──────────────────────────────────────────────────────────────────────
rdata_dir  = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO/GSE137570_CKD_eGFR/Rdata"
figure_dir = "C:/Dropbox/PNU/시스템생물학연구실/data/gluconeogenesis/GEO/GSE137570_CKD_eGFR/figure"
gmt_file   = "C:/Dropbox/PNU/시스템생물학연구실/DB/msigdb/m2.all.v2025.1.Mm.symbols.gmt"

# ── Load expression matrices ───────────────────────────────────────────────────
load(file.path(rdata_dir, "cohort1_normalized_expr.RData"))  # expr_mat_c1, meta_c1
load(file.path(rdata_dir, "cohort2_normalized_expr.RData"))  # expr_mat_c2, meta_c2

# ── Load MSigDB gene sets ──────────────────────────────────────────────────────
message("Loading MSigDB gene sets...")
c2_raw  = read.gmt(gmt_file)
c2_list = split(x = c2_raw$gene, f = c2_raw$term)

# Use specific pathways only
selected_pathways = c("REACTOME_GLUCONEOGENESIS","REACTOME_CELLULAR_SENESCENCE")
stopifnot(all(selected_pathways %in% names(c2_list)))

# Build gene set list; convert Mm symbols → uppercase to match human symbols
pathway_list = lapply(c2_list[selected_pathways], toupper)
cat("\nSelected pathways:\n")
for (nm in names(pathway_list))
  cat(sprintf("  %s  (%d genes)\n", nm, length(pathway_list[[nm]])))

# Gene overlap check
cat("\n--- Gene overlap (pathway genes vs Cohort 1 row names) ---\n")
for (nm in names(pathway_list)) {
  n_overlap = sum(pathway_list[[nm]] %in% rownames(expr_mat_c1))
  n_total   = length(pathway_list[[nm]])
  cat(sprintf("  %-60s %d / %d matched\n", nm, n_overlap, n_total))
}

# ── Exclude G5_failure from Cohort 1 ──────────────────────────────────────────
excl_c1      = meta_c1$CKD_stage != "G5_failure"
meta_c1      = meta_c1[excl_c1, ]
expr_mat_c1  = expr_mat_c1[, excl_c1]
cat(sprintf("\nCohort 1 after G5_failure exclusion: %d samples\n", ncol(expr_mat_c1)))

# ── 수동 아웃라이어 제거 (S29, S7) ────────────────────────────────────────────
manual_excl  = c("S29", "S7")
keep_manual  = !meta_c1$sample_id %in% manual_excl
meta_c1      = meta_c1[keep_manual, ]
expr_mat_c1  = expr_mat_c1[, keep_manual]
cat(sprintf("Cohort 1 after manual exclusion (%s): %d samples\n",
            paste(manual_excl, collapse = ", "), ncol(expr_mat_c1)))

ckd_colors = c(G1_normal           = "#4575b4",
               G2_mild             = "#91bfdb",
               G3a_mild_moderate   = "#fee090",
               G3b_moderate_severe = "#fc8d59",
               G4_severe           = "#d73027")




# ── GSVA scoring ───────────────────────────────────────────────────────────────
run_gsva = function(expr_mat, gset_list) {
  # GSVA >= 2.0: new param-based API only (legacy API defunct)
  param = gsvaParam(expr_mat, gset_list, kcdf = "Gaussian")
  gsva(param, verbose = FALSE)
}

message("\nCalculating GSVA scores - Cohort 1 (", ncol(expr_mat_c1), " samples)...")
scores_c1 = run_gsva(expr_mat_c1, pathway_list)


# ── Tidy helper ────────────────────────────────────────────────────────────────
scores_to_long = function(scores_mat, meta_df, id_col = "sample_id") {
  df      = as.data.frame(t(scores_mat))
  df[[id_col]] = rownames(df)
  df_long = pivot_longer(df, cols = -all_of(id_col),
                         names_to = "pathway", values_to = "score")
  left_join(df_long, meta_df, by = id_col)
}

long_c1 = scores_to_long(scores_c1, meta_c1)

# pathway label 정리 (strip에서 줄바꿈)
pathway_labels = c(
  REACTOME_GLUCONEOGENESIS     = "Gluconeogenesis",
  REACTOME_CELLULAR_SENESCENCE = "Cellular Senescence"
)

pathway_order = c(
  "Gluconeogenesis",
  "Cellular Senescence"
)

long_c1$pathway = factor(
  pathway_labels[as.character(long_c1$pathway)],
  levels = pathway_order
)


# ══════════════════════════════════════════════════════════════════════════════
# Cohort 1 Visualization (eGFR / CKD stage)
# ══════════════════════════════════════════════════════════════════════════════
ckd_colors = c(G1_normal           = "#4575b4",
               G2_mild             = "#91bfdb",
               G3a_mild_moderate   = "#fee090",
               G3b_moderate_severe = "#fc8d59",
               G4_severe           = "#d73027")

# 1a. Boxplot by CKD stage
p1a = ggplot(long_c1, aes(x = CKD_stage, y = score, fill = CKD_stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = CKD_stage), width = 0.15, size = 2.5, alpha = 0.7) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 2) +
  scale_fill_manual(values  = ckd_colors) +
  scale_color_manual(values = ckd_colors) +
  #stat_compare_means(method = "kruskal.test", label = "p.format",
  #                   label.y.npc = 0.97, label.x.npc = 0.80, size = 3.5) +
  labs(x = "CKD Stage", y = "GSVA Score") +
  theme_bw(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, colour = "black", size = 12),
        axis.text.y  = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.minor.y = element_blank())
p1a
ggsave(file.path(figure_dir, "2_C1_pathway_by_CKD_stage.png"), p1a, width = 6, height = 4, dpi = 150)



################# supple 
p1b = ggplot(long_c1, aes(x = GFR, y = score, color = CKD_stage)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", color = "black", se = TRUE,
              linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 2) +
  scale_color_manual(values = ckd_colors) +
  stat_cor(method = "spearman", color = "black", size = 3.5,
           label.x.npc = 0.05, label.y.npc = 0.97) +
  labs(x = "eGFR (mL/min/1.73m²)", y = "GSVA Score",
       color = "CKD Stage") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        strip.text = element_text(size = 11, face = "bold"),
        axis.text.x  = element_text(angle = 0, hjust = 0.5, colour = "black", size = 12),
        axis.text.y  = element_text(colour = "black", size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank())
p1b
ggsave(file.path(figure_dir, "2_C1_pathway_vs_GFR_scatter.png"), p1b, width = 9, height = 4, dpi = 150)



###############################################

wide_c1 = long_c1 %>%
  pivot_wider(names_from = pathway, values_from = score)

p_cor_c1 = ggplot(wide_c1,
                  aes(x = `Cellular Senescence`, y = Gluconeogenesis,
                      color = CKD_stage)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", color = "black", se = TRUE,
              linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = ckd_colors) +
  #stat_cor(method = "spearman", color = "black", size = 4,
  #         label.x.npc = 0.4, label.y.npc = 0.97) +
  labs(title = "",
       x = "Cellular Senescence Score",
       y = "Gluconeogenesis Score",
       color = "CKD State") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", 
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 13, colour = "black"),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 15, colour = "black", face = "bold"),
        panel.grid.minor = element_blank())

p_cor_c1

ggsave(file.path(figure_dir, "2_C1_senescence_vs_gluconeogenesis_cor.png"),
       p_cor_c1, width = 6, height = 3.5, dpi = 150)

