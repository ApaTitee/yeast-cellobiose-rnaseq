#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 16_enrichment.R   （PLAN 5.10：功能分析）
# ORA ：clusterProfiler::enricher，TERM2GENE 自建
#         · GO slim  —— 来自 SGD go_slim_mapping.tab（按 SGD 规范）
#         · 全 GO     —— 来自 SGD GAF（gene_association 的官方发布版本）
#       universe = 检出表达的基因（count >= 10 且 >= 3 个样本）——与 5.10 定义一致
# GSEA：以 apeglm 收缩 LFC 排序（clusterProfiler::GSEA）
# KEGG：enrichKEGG(organism = "sce")，原始结果缓存入库；若网络不可用则标记 not_assessable
# 输出：results/enrichment/*.tsv、results/enrichment/kegg_raw.rds、results/figures/enrichment/
# 与原文方法学差异：原文 FunSpec + Bonferroni(0.01)；本研究 clusterProfiler + BH（须在报告中注明）
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler); library(ggplot2); library(GO.db)
})
# 注意：clusterProfiler 的依赖会把 sub() 覆盖为 S4 泛型，必须显式用 base::sub；
# 且变量名不可用 args（base R 的函数）
.cmdargs <- commandArgs(trailingOnly = FALSE)
.script_path <- base::sub("^--file=", "", .cmdargs[grep("^--file=", .cmdargs)])
root <- if (length(.script_path)) normalizePath(file.path(dirname(.script_path), "..", "..")) else normalizePath(".")
enr_dir <- file.path(root, "results", "enrichment"); dir.create(enr_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(root, "results", "figures", "enrichment"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ann <- file.path(root, "refs", "host", "annotation")
de_dir <- file.path(root, "results", "differential_expression")

ORS <- "^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$"

# ---------- TERM2GENE（由 16a_build_term2gene.R 预先生成，保证口径一致）----------
# GO slim：SGD go_slim_mapping.tab（SGD 规范）
# 全 GO  ：SGD GAF 经 GO 层级传播（true path rule）到祖先术语
#          —— SGD GAF 只含最具体术语，不传播会使宽泛术语的基因集严重偏小
t2g_slim_df <- read.delim(file.path(ann, "term2gene_goslim.tsv"), stringsAsFactors = FALSE)
t2g_go_df   <- read.delim(file.path(ann, "term2gene_go_propagated.tsv"), stringsAsFactors = FALSE)
t2g_slim <- unique(data.frame(term = paste0(t2g_slim_df$term_name, " [", t2g_slim_df$term, "]"),
                              gene = t2g_slim_df$gene))
t2g_go   <- unique(data.frame(term = t2g_go_df$term_name, gene = t2g_go_df$gene))
cat(sprintf("GO slim TERM2GENE：%d 个 slim 术语，覆盖 %d 个基因\n",
            length(unique(t2g_slim$term)), length(unique(t2g_slim$gene))))
cat(sprintf("全 GO TERM2GENE（含层级传播）：%d 个术语，覆盖 %d 个基因\n",
            length(unique(t2g_go$term)), length(unique(t2g_go$gene))))

# ---------- 输入：universe 与 DEG ----------
all_genes <- read.delim(file.path(de_dir, "DESeq2_all_genes.tsv"), stringsAsFactors = FALSE)
universe <- all_genes$gene_id
cat(sprintf("universe（检出表达的基因）= %d\n", length(universe)))
read_deg <- function(f) {
  d <- read.delim(file.path(de_dir, f), stringsAsFactors = FALSE)
  list(up = d$gene_id[d$direction == "up_in_cellobiose"],
       down = d$gene_id[d$direction == "down_in_cellobiose"], all = d)
}
prim <- read_deg("DEG_primary.tsv"); tmat <- read_deg("DEG_threshold_matched.tsv")

ora <- function(genes, t2g, label, pcut = 0.05, qcut = 0.05) {
  if (length(genes) < 5) return(NULL)
  e <- enricher(genes, TERM2GENE = t2g, universe = universe, pvalueCutoff = 1, qvalueCutoff = 1,
                minGSSize = 5, maxGSSize = 2000)
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  df <- as.data.frame(e)
  df$comparison <- label
  df <- df[df$p.adjust < qcut & df$pvalue < pcut, ]
  df[order(df$p.adjust), ]
}

res <- list(
  ora(prim$up,   t2g_slim, "primary_UP_GOslim"),
  ora(prim$down, t2g_slim, "primary_DOWN_GOslim"),
  ora(prim$up,   t2g_go,   "primary_UP_GO"),
  ora(prim$down, t2g_go,   "primary_DOWN_GO"),
  ora(tmat$up,   t2g_slim, "threshold_matched_UP_GOslim"),
  ora(tmat$down, t2g_slim, "threshold_matched_DOWN_GOslim")
)
res <- res[!vapply(res, is.null, logical(1))]
if (length(res)) {
  ora_all <- do.call(rbind, res)
  write.table(ora_all, file.path(enr_dir, "ora_go.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("ORA 结果：%d 行（已写入 ora_go.tsv）\n", nrow(ora_all)))
  slim_res <- ora_all[grepl("GOslim", ora_all$comparison), ]
  if (nrow(slim_res)) {
    top <- do.call(rbind, lapply(split(slim_res, slim_res$comparison), function(d) head(d, 10)))
    g <- ggplot(top, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) +
      geom_point(size = 2) + facet_wrap(~comparison, scales = "free_y") +
      labs(x = "-log10(BH-adjusted p)", y = NULL, title = "GO slim enrichment (top 10 per set)") +
      theme_bw(base_size = 9)
    ggsave(file.path(fig_dir, "go_slim_ora.pdf"), g, width = 10, height = 6, dpi = 300)
    ggsave(file.path(fig_dir, "go_slim_ora.png"), g, width = 10, height = 6, dpi = 300)
  }
} else {
  cat("ORA 无显著结果\n")
  write.table(data.frame(), file.path(enr_dir, "ora_go.tsv"))
}

# ---------- GSEA（以 apeglm 收缩 LFC 排序） ----------
rnk <- all_genes$LFC_apeglm
names(rnk) <- all_genes$gene_id
rnk <- sort(rnk[!is.na(rnk)], decreasing = TRUE)
set.seed(1)
gsea <- try(GSEA(rnk, TERM2GENE = as.data.frame(t2g_go), pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, eps = 0, verbose = FALSE), silent = TRUE)
if (!inherits(gsea, "try-error") && !is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
  gd <- as.data.frame(gsea)
  write.table(gd, file.path(enr_dir, "gsea_go.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("GSEA：%d 条显著通路（已写入 gsea_go.tsv）\n", nrow(gd)))
} else {
  cat("GSEA 未能完成或无显著结果\n")
  write.table(data.frame(), file.path(enr_dir, "gsea_go.tsv"))
}

# ---------- KEGG（在线，结果缓存） ----------
kegg_cache <- file.path(enr_dir, "kegg_raw.rds")
kk <- try({
  if (file.exists(kegg_cache)) readRDS(kegg_cache) else {
    k <- enrichKEGG(gene = prim$all$gene_id, organism = "sce", universe = universe,
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    saveRDS(k, kegg_cache); k
  }
}, silent = TRUE)
if (!inherits(kk, "try-error") && !is.null(kk) && nrow(as.data.frame(kk)) > 0) {
  kd <- as.data.frame(kk); kd$comparison <- "primary_all"
  write.table(kd, file.path(enr_dir, "kegg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("KEGG：%d 条通路（原始结果已缓存 %s）\n", nrow(kd), basename(kegg_cache)))
} else {
  cat("KEGG 富集不可用（网络/接口问题）-> 标记为 not_assessable，写入 results/comparison/not_assessable.tsv 由 18 脚本维护\n")
  write.table(data.frame(note = "KEGG enrichment unavailable in this run"), file.path(enr_dir, "kegg.tsv"))
}

# ---------- 与原文 Dataset S3 的对照数据准备 ----------
s3 <- read.delim(file.path(root, "data", "metadata", "original_study", "derived", "dataset_s3_goslim.tsv"),
                 stringsAsFactors = FALSE)
if (nrow(s3)) {
  slim_terms <- unique(data.frame(term = paste0(s3$go_term, " [", s3$go_id, "]"), direction = s3$direction,
                                  p_orig = s3$p_bonferroni, genes_orig = s3$genes))
  write.table(slim_terms, file.path(enr_dir, "original_goslim_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("原文 Dataset S3 的 slim 术语：UP %d / DOWN %d（已整理供对照）\n",
              sum(slim_terms$direction == "UP"), sum(slim_terms$direction == "DOWN")))
}
cat("\n功能分析完成\n")
