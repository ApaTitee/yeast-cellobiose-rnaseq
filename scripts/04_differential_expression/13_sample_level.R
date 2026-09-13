#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_sample_level.R   （PLAN 5.7：样本水平分析）
# 目的：评估整体结构与重复一致性（PCA / 相关性 / 层次聚类）
# 输入：results/quantification/gene_counts.tsv + config/samplesheet.tsv
# 输出：results/differential_expression/vst_matrix.tsv
#       results/differential_expression/sample_correlation.tsv
#       results/figures/{pca,sample_correlation,clustering}.pdf|png
#       results/differential_expression/sample_level_summary.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap); library(RColorBrewer)
})
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
root <- if (length(.script)) normalizePath(file.path(dirname(.script), "..", "..")) else normalizePath(".")
fig_dir <- file.path(root, "results", "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
de_dir <- file.path(root, "results", "differential_expression"); dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)

ss <- read.delim(file.path(root, "config", "samplesheet.tsv"), stringsAsFactors = FALSE)
cnt <- read.delim(file.path(root, "results", "quantification", "gene_counts.tsv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
mat <- as.matrix(cnt[, -1]); rownames(mat) <- cnt$gene_id
mat <- mat[, ss$sample_id, drop = FALSE]
col_data <- data.frame(row.names = ss$sample_id,
                       condition = factor(ss$condition, levels = c("glucose", "cellobiose")),
                       replicate = factor(ss$replicate))
stopifnot(identical(colnames(mat), rownames(col_data)))

dds <- DESeqDataSetFromMatrix(round(mat), col_data, design = ~ condition)
# 低表达过滤（仅用于样本水平结构展示，与 5.10 的 universe 定义一致）
keep <- rowSums(counts(dds) >= 10) >= 3
cat(sprintf("基因总数 %d；count>=10 且 >=3 样本的基因 %d\n", nrow(dds), sum(keep)))
vsd <- vst(dds[keep, ], blind = TRUE)

# --- 1. VST 矩阵 ---
v <- assay(vsd)
write.table(data.frame(gene_id = rownames(v), v, check.names = FALSE),
            file.path(de_dir, "vst_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- 2. 样本相关性 ---
sp <- cor(v, method = "spearman"); pe <- cor(v, method = "pearson")
write.table(data.frame(metric = "spearman", sample = rownames(sp), sp, check.names = FALSE),
            file.path(de_dir, "sample_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pdf(file.path(fig_dir, "sample_correlation.pdf"), width = 6.5, height = 5.8)
pheatmap(sp, display_numbers = TRUE, number_format = "%.4f", cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         main = "Sample correlation (Spearman, VST)", fontsize = 9)
dev.off()

# --- 3. PCA ---
pca <- prcomp(t(v), scale. = FALSE)
pv <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
pd <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 condition = col_data$condition, replicate = col_data$replicate,
                 sample = rownames(pca$x))
p <- ggplot(pd, aes(PC1, PC2, colour = condition, shape = replicate, label = sample)) +
  geom_point(size = 3.5) + geom_text(vjust = -1.1, size = 2.6, show.legend = FALSE) +
  labs(x = sprintf("PC1 (%.1f%% variance)", pv[1]), y = sprintf("PC2 (%.1f%% variance)", pv[2]),
       title = "PCA of VST-normalised expression", colour = "Carbon source", shape = "Replicate") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "pca.pdf"), p, width = 6.2, height = 4.6)
ggsave(file.path(fig_dir, "pca.png"), p, width = 6.2, height = 4.6, dpi = 300)

# --- 4. 层次聚类 ---
samp_dist <- as.dist(1 - sp)
pdf(file.path(fig_dir, "clustering.pdf"), width = 6.2, height = 5.2)
plot(hclust(samp_dist, method = "average"), main = "Hierarchical clustering of samples (1 - Spearman)",
     xlab = "", sub = "", ylab = "1 - rho")
dev.off()

# --- 5. 汇总结论 ---
within <- function(cond) {
  ids <- ss$sample_id[ss$condition == cond]
  sp[ids, ids][upper.tri(sp[ids, ids])]
}
gm <- within("glucose"); cm <- within("cellobiose")
between <- sp[ss$sample_id[ss$condition == "glucose"], ss$sample_id[ss$condition == "cellobiose"]]
summ <- data.frame(
  item = c("n_genes_used", "PC1_var_pct", "PC2_var_pct",
           "within_glucose_spearman_min", "within_glucose_spearman_max",
           "within_cellobiose_spearman_min", "within_cellobiose_spearman_max",
           "between_condition_spearman_min", "between_condition_spearman_max",
           "separation_verdict"),
  value = c(sum(keep), pv[1], pv[2], min(gm), max(gm), min(cm), max(cm),
            min(between), max(between),
            ifelse(min(gm, cm) > max(between), "conditions separated on PC1", "no clean separation")))
write.table(summ, file.path(de_dir, "sample_level_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)
cat("\n样本水平分析完成\n")
