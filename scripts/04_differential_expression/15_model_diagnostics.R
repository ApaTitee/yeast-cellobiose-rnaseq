#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 15_model_diagnostics.R   （PLAN 5.9：模型诊断）
# 内容：p 值直方图、dispersion 拟合、MA plot（shrunken LFC）、Cook's distance 与离群检查
# 输出：results/figures/model_diagnostics/*.pdf|png
#       results/differential_expression/model_diagnostics.tsv（摘要数字）
# 说明：为保持脚本独立可重跑，此处按与 14 相同的设定重新拟合（universe 同 5.10）
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(DESeq2); library(tximport); library(ggplot2) })
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
root <- if (length(.script)) normalizePath(file.path(dirname(.script), "..", "..")) else normalizePath(".")
fig_dir <- file.path(root, "results", "figures", "model_diagnostics")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
de_dir <- file.path(root, "results", "differential_expression")

ss <- read.delim(file.path(root, "config", "samplesheet.tsv"), stringsAsFactors = FALSE)
tx2gene <- read.delim(file.path(root, "refs", "custom", "tx2gene.tsv"), stringsAsFactors = FALSE)[, 1:2]
files <- file.path(root, "results", "quantification", "salmon",
                   paste0(ss$run, "_", ss$sample_id, "_", ss$condition), "quant.sf")
names(files) <- ss$sample_id
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
col_data <- data.frame(row.names = ss$sample_id,
                       condition = factor(ss$condition, levels = c("glucose", "cellobiose")))
dds <- DESeqDataSetFromTximport(txi, col_data, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- DESeq(dds[keep, ], quiet = TRUE)
res <- results(dds, alpha = 0.05, independentFiltering = TRUE)
resL <- lfcShrink(dds, coef = "condition_cellobiose_vs_glucose", type = "apeglm", quiet = TRUE)
cat(sprintf("检验基因 %d；size factors: %s\n", nrow(dds),
            paste(sprintf("%.3f", sizeFactors(dds)), collapse = ", ")))

# --- 1. p 值直方图 ---
pv <- data.frame(p = res$pvalue[!is.na(res$pvalue)])
g1 <- ggplot(pv, aes(p)) + geom_histogram(bins = 50, boundary = 0, fill = "grey35", colour = "white") +
  labs(x = "Wald test p-value", y = "Number of genes",
       title = "p-value distribution (primary model)") + theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "pvalue_histogram.pdf"), g1, width = 5.4, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "pvalue_histogram.png"), g1, width = 5.4, height = 4, dpi = 300)
frac_small <- mean(pv$p < 0.05)

# --- 2. dispersion 拟合 ---
pdf(file.path(fig_dir, "dispersion.pdf"), width = 5.6, height = 4.6)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# --- 3. MA plot（shrunken LFC）---
pdf(file.path(fig_dir, "ma_plot.pdf"), width = 5.6, height = 4.6)
plotMA(resL, ylim = c(-8, 8), main = "MA (apeglm-shrunken LFC; display only)")
abline(h = c(-1, 1), col = "blue", lty = 2)
dev.off()

# --- 4. Cook's distance ---
cd <- assays(dds)[["cooks"]]
cd_max <- apply(cd, 1, max, na.rm = TRUE)
cook_cut <- qf(0.99, df1 = 3, df2 = ncol(dds) - 3)
n_out <- sum(cd_max > cook_cut, na.rm = TRUE)
df_cd <- data.frame(gene_id = names(cd_max), cook_max = cd_max,
                    is_outlier = cd_max > cook_cut)
write.table(df_cd[order(-df_cd$cook_max), ], file.path(de_dir, "cooks_distance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
g4 <- ggplot(df_cd, aes(x = seq_along(cook_max), y = cook_max)) +
  geom_point(size = 0.4, alpha = 0.4) + scale_y_log10() +
  geom_hline(yintercept = cook_cut, colour = "red", linetype = 2) +
  labs(x = "Gene rank", y = "max Cook's distance (log scale)",
       title = sprintf("Cook's distance (cutoff %.2f; %d genes above)", cook_cut, n_out)) +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "cooks_distance.pdf"), g4, width = 5.6, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "cooks_distance.png"), g4, width = 5.6, height = 4, dpi = 300)

# --- 5. 摘要 ---
# 反保守性检查：小 p 值是否富集（若分布平坦或左端凹陷则有问题）
bins <- hist(pv$p, breaks = seq(0, 1, 0.05), plot = FALSE)$counts
anti <- bins[1] / mean(bins)
summ <- data.frame(
  item = c("genes_tested", "frac_p_below_0.05", "p_hist_first_bin_over_mean",
           "cook_cutoff", "genes_above_cook_cutoff",
           "size_factor_min", "size_factor_max", "size_factor_ratio_max_min",
           "verdict"),
  value = c(nrow(dds), round(frac_small, 4), round(anti, 3), round(cook_cut, 3), n_out,
            round(min(sizeFactors(dds)), 3), round(max(sizeFactors(dds)), 3),
            round(max(sizeFactors(dds)) / min(sizeFactors(dds)), 3),
            ifelse(anti > 1.5 && n_out < 0.05 * nrow(dds), "no systematic anti-conservatism", "check manually")))
write.table(summ, file.path(de_dir, "model_diagnostics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)
cat("\n模型诊断完成\n")
