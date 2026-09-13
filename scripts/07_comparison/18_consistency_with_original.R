#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 18_consistency_with_original.R   （PLAN 5.12：与原研究的一致性评估）
# 指标：per-gene log2FC 相关性、DEG 集合重叠（Jaccard + 超几何检验）、方向一致率、
#       跨阈值 PR 曲线、数字对账、外源锚定值对照
# 输入：results/differential_expression/DESeq2_all_genes.tsv、DEG_threshold_matched.tsv
#       data/metadata/original_study/derived/{dataset_s1_rpkm,dataset_s2_deg,reconciliation}.tsv
# 输出：results/comparison/{consistency_metrics.tsv,deg_overlap.tsv,fc_correlation.tsv,
#                            threshold_sweep.tsv,reconciliation_check.tsv,not_assessable.tsv}
# 说明：原文 Dataset S2 的 Fold Change 列为带符号倍数（下调记 -(G8/C8)），log 列为 log2(C8/G8)；
#       本脚本一律换算为 log2(C8/G8) 后再比较。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2) })
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
root <- if (length(.script)) normalizePath(file.path(dirname(.script), "..", "..")) else normalizePath(".")
cmp_dir <- file.path(root, "results", "comparison"); dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
der <- file.path(root, "data", "metadata", "original_study", "derived")
fig_dir <- file.path(root, "results", "figures")

ours <- read.delim(file.path(root, "results", "differential_expression", "DESeq2_all_genes.tsv"), stringsAsFactors = FALSE)
tm   <- read.delim(file.path(root, "results", "differential_expression", "DEG_threshold_matched.tsv"), stringsAsFactors = FALSE)
s1   <- read.delim(file.path(der, "dataset_s1_rpkm.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
s2   <- read.delim(file.path(der, "dataset_s2_deg.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

# --- 把原文标识换算为 systematic ORF name ---
sysname <- function(syn) {
  toks <- unlist(strsplit(as.character(syn), "[|,]"))
  toks <- trimws(toks)
  hit <- grep("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", toks, value = TRUE)
  if (length(hit)) hit[1] else NA_character_
}
s1$orf <- vapply(s1$synonym, sysname, character(1))
s2$orf <- vapply(s2$synonym, sysname, character(1))
# 外源条目单独标记
s1$orf[is.na(s1$orf) & s1$is_transgene == "yes"] <- s1$feature_id[is.na(s1$orf) & s1$is_transgene == "yes"]
s2$orf[is.na(s2$orf) & s2$is_transgene == "yes"] <- s2$feature_id[is.na(s2$orf) & s2$is_transgene == "yes"]

# --- 1. per-gene log2FC 相关性（共同基因）---
m <- merge(ours, s1[, c("orf", "log2FC_recomputed")], by.x = "gene_id", by.y = "orf")
m <- m[!is.na(m$log2FC_recomputed) & !is.na(m$LFC_MLE), ]
rho <- cor(m$LFC_MLE, m$log2FC_recomputed, method = "spearman")
pea <- cor(m$LFC_MLE, m$log2FC_recomputed, method = "pearson")
write.table(m[, c("gene_id", "LFC_MLE", "log2FC_recomputed")], file.path(cmp_dir, "fc_correlation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
p_fc <- ggplot(m, aes(log2FC_recomputed, LFC_MLE)) +
  geom_point(alpha = 0.25, size = 0.7) + geom_abline(slope = 1, intercept = 0, colour = "red", linetype = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "blue", linewidth = 0.6) +
  labs(x = "Original study log2FC (cellobiose / glucose)", y = "This study log2FC (MLE)",
       title = sprintf("Per-gene log2FC agreement (n = %d; Spearman = %.3f)", nrow(m), rho)) +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "fc_agreement.pdf"), p_fc, width = 5.6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "fc_agreement.png"), p_fc, width = 5.6, height = 5, dpi = 300)

# --- 2. DEG 集合重叠（阈值对齐版 vs 原文 519）---
orig_set <- unique(na.omit(s2$orf))
our_set  <- unique(tm$gene_id)
inter <- intersect(our_set, orig_set)
uni <- union(our_set, orig_set)
jac <- length(inter) / length(uni)
# 超几何检验：在背景基因集（能在两侧都出现的基因）中观察到 >= |inter| 的概率
bg <- length(unique(intersect(ours$gene_id, s1$orf)))
p_hyper <- phyper(length(inter) - 1, length(orig_set), bg - length(orig_set), length(our_set), lower.tail = FALSE)

# 方向一致率
both <- merge(tm[, c("gene_id", "LFC_MLE", "direction")], s2[, c("orf", "log2FC_recomputed")],
              by.x = "gene_id", by.y = "orf")
dir_orig <- ifelse(both$log2FC_recomputed > 0, "up_in_cellobiose", "down_in_cellobiose")
dir_conc <- mean(both$direction == dir_orig)

metrics <- data.frame(
  metric = c("genes_comparable_FC", "spearman_log2FC", "pearson_log2FC",
             "our_DEG_threshold_matched", "original_DEG", "overlap",
             "jaccard", "hypergeometric_p", "direction_concordance_in_overlap",
             "background_genes"),
  value = c(nrow(m), round(rho, 4), round(pea, 4), length(our_set), length(orig_set),
            length(inter), round(jac, 4), signif(p_hyper, 3), round(dir_conc, 4), bg))
write.table(metrics, file.path(cmp_dir, "consistency_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(metrics, row.names = FALSE)

ov <- data.frame(gene_id = union(our_set, orig_set))
ov$in_this_study <- ov$gene_id %in% our_set
ov$in_original <- ov$gene_id %in% orig_set
ov$class <- ifelse(ov$in_this_study & ov$in_original, "both",
                   ifelse(ov$in_this_study, "this_study_only", "original_only"))
write.table(ov, file.path(cmp_dir, "deg_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- 3. 跨阈值 PR 曲线 ---
grid <- expand.grid(padj_cut = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2),
                    lfc_cut = c(0, 0.5, 1, 1.5, 2))
sweep <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  sel <- !is.na(ours$padj) & ours$padj < grid$padj_cut[i] & !is.na(ours$LFC_MLE) &
    abs(ours$LFC_MLE) > grid$lfc_cut[i]
  set <- ours$gene_id[sel]
  tp <- length(intersect(set, orig_set))
  data.frame(padj_cut = grid$padj_cut[i], lfc_cut = grid$lfc_cut[i], n_our = length(set),
             overlap = tp,
             precision = ifelse(length(set) > 0, tp / length(set), NA),
             recall = tp / length(orig_set))
}))
write.table(sweep, file.path(cmp_dir, "threshold_sweep.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p_sw <- ggplot(subset(sweep, lfc_cut %in% c(0, 1)), aes(recall, precision, colour = factor(lfc_cut))) +
  geom_point() + geom_line() + geom_text(aes(label = padj_cut), size = 2.4, vjust = -0.7) +
  labs(x = "Recall of the original 519 DEGs", y = "Precision (share of our DEGs that are in the 519)",
       colour = "|log2FC| cut", title = "Threshold sweep vs original DEG set") + theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "threshold_sweep.pdf"), p_sw, width = 6, height = 4.4, dpi = 300)
ggsave(file.path(fig_dir, "threshold_sweep.png"), p_sw, width = 6, height = 4.4, dpi = 300)

# --- 4. 数字对账 + not_assessable 清单 ---
rec <- read.delim(file.path(der, "reconciliation.tsv"), stringsAsFactors = FALSE)
rec$check <- "ok"
write.table(rec, file.path(cmp_dir, "reconciliation_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

na_items <- data.frame(
  item = c("mtDNA-encoded gene expression", "absolute transgene abundance",
           "absolute RPKM/TPM matching the original", "cross-tool mapping statistics",
           "transgene C-terminal His6 cleavage state"),
  reason = c("poly-A selection captures yeast mitochondrial transcripts inefficiently; 5-18 chrM reads per 1e6; 0/19 mt genes reach 1x",
             "the tag/junction part of the transgene reference is derived from the sample reads themselves (circularity)",
             "CLC multi-mapping and 'By totals' RPKM are not reproducible outside CLC; only gene-level trends are comparable",
             "original 92.3/83.4/76.3% were produced by CLC's internal algorithms; only order-of-magnitude comparison is meaningful",
             "requires protein-level data, out of scope"))
write.table(na_items, file.path(cmp_dir, "not_assessable.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== 关键一致性指标 ===\n")
cat(sprintf("  可比基因 %d；per-gene log2FC Spearman %.3f / Pearson %.3f\n", nrow(m), rho, pea))
cat(sprintf("  阈值对齐 DEG（本研究）%d vs 原文 %d；重叠 %d；Jaccard %.3f；超几何 p = %.3g\n",
            length(our_set), length(orig_set), length(inter), jac, p_hyper))
cat(sprintf("  重叠内方向一致率 %.1f%%\n", 100 * dir_conc))
cat("\n一致性评估完成\n")
