#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 14_deseq2_de.R   （PLAN 5.8：差异表达）
# 主分析   ：~ condition（glucose 为 reference），independentFiltering 开启
# 敏感性   ：~ replicate + condition（配对）
# 阈值     ：主 padj<0.05 & |log2FC|>1；放宽 padj<0.1；收紧 padj<0.01 & |log2FC|>1.5
#            阈值对齐原文：|log2FC|>=1 & FDR<=0.001（results(alpha=0.001)）
# LFC 口径 ：**阈值判定一律使用未收缩（MLE）LFC**；apeglm 收缩值仅用于排序与作图
# 输入：results/quantification/salmon/*/quant.sf、refs/custom/tx2gene.tsv、config/samplesheet.tsv
# 输出：results/differential_expression/
#         DESeq2_all_genes.tsv（全基因，含 MLE 与 apeglm 两套 LFC）
#         DEG_primary.tsv / DEG_relaxed.tsv / DEG_stringent.tsv /
#         DEG_paired.tsv / DEG_threshold_matched.tsv
#         de_summary.tsv（各口径的上下调计数 + 外源条目计数）
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(DESeq2); library(tximport) })
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
root <- if (length(.script)) normalizePath(file.path(dirname(.script), "..", "..")) else normalizePath(".")
out_dir <- file.path(root, "results", "differential_expression"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ss <- read.delim(file.path(root, "config", "samplesheet.tsv"), stringsAsFactors = FALSE)
tx2gene <- read.delim(file.path(root, "refs", "custom", "tx2gene.tsv"), stringsAsFactors = FALSE)[, 1:2]
files <- file.path(root, "results", "quantification", "salmon",
                   paste0(ss$run, "_", ss$sample_id, "_", ss$condition), "quant.sf")
names(files) <- ss$sample_id
stopifnot(all(file.exists(files)))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

col_data <- data.frame(row.names = ss$sample_id,
                       condition = factor(ss$condition, levels = c("glucose", "cellobiose")),
                       replicate = factor(ss$replicate))
TRANSGENES <- c("cdt-1EGFP", "gh1-1")

fit_model <- function(design, label) {
  dds <- DESeqDataSetFromTximport(txi, col_data, design = design)
  # universe（PLAN 5.10）：在 >=3 个样本中 count >= 10
  keep <- rowSums(counts(dds) >= 10) >= 3
  cat(sprintf("[%s] 过滤后基因数 %d / %d\n", label, sum(keep), nrow(dds)))
  dds <- DESeq(dds[keep, ], quiet = TRUE)
  dds
}

de_table <- function(dds, coef_name, alpha = 0.05) {
  r <- results(dds, name = coef_name, alpha = alpha, independentFiltering = TRUE)
  rs <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  df <- data.frame(gene_id = rownames(r),
                   baseMean = r$baseMean,
                   LFC_MLE = r$log2FoldChange,
                   lfcSE = r$lfcSE,
                   stat = r$stat,
                   pvalue = r$pvalue,
                   padj = r$padj,
                   LFC_apeglm = rs$log2FoldChange[rownames(r)],
                   is_transgene = rownames(r) %in% TRANSGENES,
                   row.names = NULL)
  # 方向以 MLE LFC 为准
  df$direction <- ifelse(is.na(df$LFC_MLE), NA,
                         ifelse(df$LFC_MLE > 0, "up_in_cellobiose", "down_in_cellobiose"))
  df[order(df$padj, -abs(df$LFC_MLE)), ]
}

sel <- function(df, padj_cut, lfc_cut) {
  ok <- !is.na(df$padj) & df$padj < padj_cut & !is.na(df$LFC_MLE) & abs(df$LFC_MLE) > lfc_cut
  df[ok, ]
}
write_out <- function(df, name) {
  write.table(df, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  %-28s %4d 行（上调 %d / 下调 %d；外源 %d）\n", name, nrow(df),
              sum(df$direction == "up_in_cellobiose", na.rm = TRUE),
              sum(df$direction == "down_in_cellobiose", na.rm = TRUE),
              sum(df$is_transgene, na.rm = TRUE)))
  df
}

cat("=== 主分析 ~ condition ===\n")
dds_main <- fit_model(~condition, "main")
coef_main <- "condition_cellobiose_vs_glucose"
all_main <- de_table(dds_main, coef_main, alpha = 0.05)
write.table(all_main, file.path(out_dir, "DESeq2_all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p  <- write_out(sel(all_main, 0.05, 1.0),  "DEG_primary.tsv")
rl <- write_out(sel(all_main, 0.10, 0.0),  "DEG_relaxed.tsv")
st <- write_out(sel(all_main, 0.01, 1.5),  "DEG_stringent.tsv")

# 阈值对齐原文：FDR<=0.001 且 |log2FC|>=1（alpha=0.001 重算 padj，MLE LFC）
all_0001 <- de_table(dds_main, coef_main, alpha = 0.001)
tm <- write_out(sel(all_0001, 0.001, 1.0), "DEG_threshold_matched.tsv")

cat("=== 敏感性分析 ~ replicate + condition ===\n")
dds_paired <- fit_model(~replicate + condition, "paired")
coef_paired <- "condition_cellobiose_vs_glucose"
all_paired <- de_table(dds_paired, coef_paired, alpha = 0.05)
write.table(all_paired, file.path(out_dir, "DESeq2_all_genes_paired.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pp <- write_out(sel(all_paired, 0.05, 1.0), "DEG_paired.tsv")

# --- 汇总结论 ---
summ <- data.frame(
  model = c("~condition (primary)", "~condition relaxed", "~condition stringent",
            "~condition FDR<=0.001 & |LFC|>=1", "~replicate+condition (paired)"),
  n_tested = c(nrow(all_main), nrow(all_main), nrow(all_main), nrow(all_0001), nrow(all_paired)),
  n_deg = c(nrow(p), nrow(rl), nrow(st), nrow(tm), nrow(pp)),
  n_up = c(sum(p$direction == "up_in_cellobiose"), sum(rl$direction == "up_in_cellobiose"),
           sum(st$direction == "up_in_cellobiose"), sum(tm$direction == "up_in_cellobiose"),
           sum(pp$direction == "up_in_cellobiose")),
  n_down = c(sum(p$direction == "down_in_cellobiose"), sum(rl$direction == "down_in_cellobiose"),
             sum(st$direction == "down_in_cellobiose"), sum(tm$direction == "down_in_cellobiose"),
             sum(pp$direction == "down_in_cellobiose")),
  n_transgene = c(sum(p$is_transgene), sum(rl$is_transgene), sum(st$is_transgene),
                  sum(tm$is_transgene), sum(pp$is_transgene)))
write.table(summ, file.path(out_dir, "de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)

# --- 主 vs 配对的一致性 ---
common <- intersect(p$gene_id, pp$gene_id)
cat(sprintf("\n主分析与配对模型的 DEG 交集 %d / %d（主）%d（配对），Jaccard = %.3f\n",
            length(common), nrow(p), nrow(pp),
            length(common) / length(union(p$gene_id, pp$gene_id))))
both <- merge(p[, c("gene_id", "LFC_MLE", "direction")], pp[, c("gene_id", "LFC_MLE", "direction")],
              by = "gene_id", suffixes = c("_main", "_paired"))
if (nrow(both)) {
  cat(sprintf("  两模型共同 DEG 中方向一致率 = %.1f%%\n",
              100 * mean(both$direction_main == both$direction_paired)))
  cat(sprintf("  共同 DEG 的 LFC 相关（Pearson）= %.4f\n", cor(both$LFC_MLE_main, both$LFC_MLE_paired)))
}
# --- 外源条目在主分析中的表现 ---
cat("\n外源条目（主分析）：\n")
print(all_main[all_main$is_transgene, c("gene_id", "baseMean", "LFC_MLE", "LFC_apeglm", "padj", "direction")], row.names = FALSE)
cat("\n差异表达分析完成\n")
