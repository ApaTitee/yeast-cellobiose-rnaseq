#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 11_tximport_gene_matrix.R   （PLAN 5.4）
# 目的：把 salmon 的转录本水平定量汇总为基因水平矩阵（主用 k=31；另出 k=25 供 5.5 敏感性分析）
# 输入：results/quantification/salmon*/<run>_<sample>_<condition>/quant.sf
#       refs/custom/tx2gene.tsv
# 输出：results/quantification/gene_counts.tsv        （lengthScaledTPM 计数，供 DESeq2）
#       results/quantification/gene_tpm.tsv           （TPM）
#       results/quantification/gene_avg_length.tsv
#       results/quantification/gene_counts_k25.tsv    （k=25 敏感性对照）
# 说明：countsFromAbundance="lengthScaledTPM" 的语义写入 README 与 docs/decisions；
#       本文件输出的是 length-scaled TPM 计数，不是原始 read 计数。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(tximport))
# 使用 base R 写表，避免额外依赖（readr 未纳入 environment.yml）
write_tsv <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

# 通过 --file= 参数定位脚本位置，从而稳健地解析项目根目录（不依赖 sys.frame）
.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
root <- if (length(.script)) normalizePath(file.path(dirname(.script), "..", "..")) else normalizePath(".")
stopifnot(file.exists(file.path(root, "config", "samplesheet.tsv")))
sample_sheet <- read.delim(file.path(root, "config", "samplesheet.tsv"), stringsAsFactors = FALSE)
tx2gene <- read.delim(file.path(root, "refs", "custom", "tx2gene.tsv"), stringsAsFactors = FALSE)[, 1:2]

run_one <- function(quant_dir, tag) {
  files <- file.path(root, quant_dir, paste0(sample_sheet$run, "_", sample_sheet$sample_id, "_", sample_sheet$condition), "quant.sf")
  names(files) <- sample_sheet$sample_id
  stopifnot(all(file.exists(files)))
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
  out_dir <- file.path(root, "results", "quantification")
  counts <- as.data.frame(round(txi$counts, 3))
  counts <- cbind(gene_id = rownames(counts), counts)
  write_tsv(counts, file.path(out_dir, sprintf("gene_counts%s.tsv", tag)))
  tpm <- as.data.frame(txi$abundance)
  tpm <- cbind(gene_id = rownames(tpm), tpm)
  write_tsv(tpm, file.path(out_dir, sprintf("gene_tpm%s.tsv", tag)))
  len <- as.data.frame(txi$length)
  len <- cbind(gene_id = rownames(len), len)
  write_tsv(len, file.path(out_dir, sprintf("gene_avg_length%s.tsv", tag)))
  cat(sprintf("[%s] 基因数 %d；样本 %s\n", tag, nrow(counts), paste(names(files), collapse = ",")))
  invisible(list(counts = counts, tpm = tpm))
}

res31 <- run_one(file.path("results", "quantification", "salmon"), "")
res25 <- run_one(file.path("results", "quantification", "salmon_k25"), "_k25")

# --- 自检 ---
tg <- c("cdt-1EGFP", "gh1-1", "YEL021W", "YCR012W")
cat("\n自检（k31，TPM）：\n")
for (g in tg) if (g %in% res31$tpm$gene_id) {
  cat(sprintf("  %-10s %s\n", g, paste(sprintf("%.1f", as.numeric(res31$tpm[res31$tpm$gene_id == g, -1])), collapse = "  ")))
}
cat("\n列顺序：", paste(colnames(res31$counts)[-1], collapse = ", "), "\n")
cat("计数矩阵是否为整数以外的值（lengthScaledTPM 的正常现象）：",
    any(abs(res31$counts[, -1] - round(res31$counts[, -1])) > 1e-6), "\n")
