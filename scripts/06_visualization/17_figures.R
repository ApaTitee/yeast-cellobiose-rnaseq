#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 17_figures.R   （PLAN 5.7/5.8/5.11 的图件汇总）
# 目的：生成报告中使用的核心图件（图注与标签统一英文）
# 输出：results/figures/{volcano,deg_heatmap,transgene_anchor}.pdf|png
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(pheatmap); library(RColorBrewer) })
.cmdargs <- commandArgs(trailingOnly = FALSE)
.script_path <- base::sub("^--file=", "", .cmdargs[grep("^--file=", .cmdargs)])
root <- if (length(.script_path)) normalizePath(file.path(dirname(.script_path), "..", "..")) else normalizePath(".")
fig_dir <- file.path(root, "results", "figures")
de_dir <- file.path(root, "results", "differential_expression")
interp <- file.path(root, "results", "interpretation")

all_g <- read.delim(file.path(de_dir, "DESeq2_all_genes.tsv"), stringsAsFactors = FALSE)
ss <- read.delim(file.path(root, "config", "samplesheet.tsv"), stringsAsFactors = FALSE)

# ---------- 1. Volcano ----------
d <- all_g[!is.na(all_g$padj) & !is.na(all_g$LFC_MLE), ]
d$class <- ifelse(d$padj < 0.05 & abs(d$LFC_MLE) > 1,
                  ifelse(d$LFC_MLE > 0, "up in cellobiose", "down in cellobiose"), "not significant")
d$class <- factor(d$class, levels = c("up in cellobiose", "down in cellobiose", "not significant"))
lab <- rbind(head(d[order(d$padj), ][d$class[order(d$padj)] == "up in cellobiose", ], 4),
             head(d[order(d$padj), ][d$class[order(d$padj)] == "down in cellobiose", ], 4))
p1 <- ggplot(d, aes(LFC_MLE, -log10(padj), colour = class)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_colour_manual(values = c("up in cellobiose" = "#C0392B", "down in cellobiose" = "#2471A3",
                                 "not significant" = "grey75")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40") +
  labs(x = expression(log[2]~"fold change (cellobiose / glucose, MLE)"),
       y = expression(-log[10]~"(BH-adjusted p)"), colour = NULL,
       title = sprintf("Volcano plot (n = %d genes; %d DEGs at padj < 0.05 and |log2FC| > 1)",
                       nrow(d), sum(d$class != "not significant"))) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(fig_dir, "volcano.pdf"), p1, width = 6.4, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "volcano.png"), p1, width = 6.4, height = 5, dpi = 300)

# ---------- 2. DEG 热图（按 padj 取前 50）----------
vst <- read.delim(file.path(de_dir, "vst_matrix.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
m <- as.matrix(vst[, -1]); rownames(m) <- vst$gene_id
deg <- read.delim(file.path(de_dir, "DEG_primary.tsv"), stringsAsFactors = FALSE)
top <- head(deg$gene_id[order(deg$padj)], 50)
top <- intersect(top, rownames(m))
ann <- data.frame(`Carbon source` = factor(ss$condition, levels = c("glucose", "cellobiose")),
                  row.names = ss$sample_id, check.names = FALSE)
pdf(file.path(fig_dir, "deg_heatmap.pdf"), width = 6.4, height = 7.6)
pheatmap(m[top, ss$sample_id], scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
         annotation_col = ann, show_rownames = TRUE, fontsize_row = 5.2, fontsize_col = 8,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Top 50 DEGs (row-scaled VST)")
dev.off()

# ---------- 3. 外源基因锚定图（多估计对照）----------
x <- read.delim(file.path(interp, "h7_transgene_crosscheck.tsv"), stringsAsFactors = FALSE)
x <- x[x$estimator %in% c("salmon_CPM", "salmon_TPM", "reads_1M_subsample", "original_RPKM",
                          "URA3_normalized", "URA3_proxy"), ]
x$gene <- ifelse(grepl("gh1-1", x$gene), "gh1-1", ifelse(grepl("cdt-1", x$gene), "cdt-1EGFP", "URA3 (proxy)"))
x$estimator <- factor(x$estimator,
                      levels = c("reads_1M_subsample", "salmon_CPM", "salmon_TPM",
                                 "original_RPKM", "URA3_proxy", "URA3_normalized"))
p3 <- ggplot(x, aes(estimator, as.numeric(FC), fill = gene)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
  geom_text(aes(label = sprintf("%.2f", as.numeric(FC))), position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.6) +
  labs(x = NULL, y = "Fold change (cellobiose / glucose)", fill = NULL,
       title = "Transgene and plasmid-copy-number proxy: independent estimators") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  coord_cartesian(ylim = c(0, max(as.numeric(x$FC), na.rm = TRUE) * 1.15))
ggsave(file.path(fig_dir, "transgene_anchor.pdf"), p3, width = 7, height = 4.4, dpi = 300)
ggsave(file.path(fig_dir, "transgene_anchor.png"), p3, width = 7, height = 4.4, dpi = 300)

# ---------- 4. 为只输出 PDF 的图件补齐 PNG（HTML 报告需要；使用 poppler 的 pdftoppm）----------
pdfs <- list.files(fig_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
for (f in pdfs) {
  png <- base::sub("\\.pdf$", ".png", f)
  if (!file.exists(png)) {
    system2("pdftoppm", c("-png", "-r", "200", "-singlefile", shQuote(f), shQuote(base::sub("\\.pdf$", "", f))),
            stdout = FALSE, stderr = FALSE)
    if (file.exists(png)) cat(sprintf("  生成 %s\n", basename(png)))
  }
}
cat("图件已生成：volcano、deg_heatmap、transgene_anchor（+ PCA/相关性/聚类/诊断/富集/一致性图），并补齐 PDF->PNG\n")
