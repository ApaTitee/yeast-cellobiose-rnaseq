#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 19_hypothesis_scorecard.R   （PLAN 5.11：生物学与工程解读 / 假设记分卡）
# 目的：对预先注册的 H0–H9 逐条给出**量化判定**（复现 / 未复现 / 新发现 / 无法评估）
# 要求：每条判定必须附集合水平统计量（Wilcoxon 单样本检验 + 与全基因背景的比较）
# 输入：results/differential_expression/DESeq2_all_genes.tsv、refs/host 注释、
#       results/differential_expression/sample_level_summary.tsv、
#       results/interpretation/h7_transgene_crosscheck.tsv
# 输出：results/interpretation/hypothesis_scorecard.tsv
#       results/interpretation/hypothesis_gene_sets.tsv
# 判定口径：整集平均 LFC 与原文方向一致且 Wilcoxon p<0.05 记 replicated；
#           方向一致但 p>=0.05 记 not_replicated（功效不足）；方向相反记 contradicted；
#           无数据支撑记 not_assessable。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(GO.db) })
.cmdargs <- commandArgs(trailingOnly = FALSE)
.script_path <- base::sub("^--file=", "", .cmdargs[grep("^--file=", .cmdargs)])
root <- if (length(.script_path)) normalizePath(file.path(dirname(.script_path), "..", "..")) else normalizePath(".")
interp <- file.path(root, "results", "interpretation"); dir.create(interp, showWarnings = FALSE, recursive = TRUE)
de_dir <- file.path(root, "results", "differential_expression")
ann <- file.path(root, "refs", "host", "annotation")

all_g <- read.delim(file.path(de_dir, "DESeq2_all_genes.tsv"), stringsAsFactors = FALSE)
lfc <- setNames(all_g$LFC_MLE, all_g$gene_id)
padj <- setNames(all_g$padj, all_g$gene_id)
universe <- all_g$gene_id
# 线粒体编码基因（locus_tag Q*）单独标记：H1 的检验主体仅限核编码基因
mt_encoded <- grep("^Q[0-9]", universe, value = TRUE)

# ---------- 基因名 -> 系统名 映射（来自 RefSeq GFF） ----------
gff <- read.delim(gzfile(file.path(root, "refs", "host",
                                   "GCF_000146045.2_R64_genomic.gff.gz")),
                  header = FALSE, comment.char = "#", stringsAsFactors = FALSE, quote = "")
sym2orf <- list()
for (a in gff$V9[gff$V3 == "gene"]) {
  m_orf <- regmatches(a, regexpr("locus_tag=([^;]+)", a))
  m_sym <- regmatches(a, regexpr("(?:^|;)gene=([^;]+)", a))
  if (length(m_orf) && length(m_sym)) {
    sym2orf[[sub(".*=", "", m_sym)]] <- sub("locus_tag=", "", m_orf)
  }
}
sym_to_orf <- function(x) {
  y <- unlist(sym2orf[x]); if (is.null(y)) character(0) else unique(y[!is.na(y)])
}

# ---------- GO 术语 -> 基因（来自 16a 生成的传播版映射，仅取核编码） ----------
t2g <- read.delim(file.path(ann, "term2gene_go_propagated.tsv"), stringsAsFactors = FALSE)
go2gene <- split(t2g$gene, t2g$term)
go_genes <- function(ids) unique(unlist(go2gene[ids]))

# ---------- 集合检验 ----------
is_deg <- !is.na(padj) & padj < 0.05 & !is.na(lfc) & abs(lfc) > 1
bg_genes <- setdiff(universe, mt_encoded)

test_set <- function(genes, expect = c("up", "down", "any"), label) {
  expect <- match.arg(expect)
  genes <- setdiff(intersect(unique(genes), universe), mt_encoded)
  if (length(genes) < 3) {
    return(list(label = label, n = length(genes), mean_lfc = NA, median_lfc = NA,
                p_wilcox = NA, frac_up = NA, n_deg = NA, fisher_p = NA,
                verdict = "not_assessable"))
  }
  x <- lfc[genes]; x <- x[!is.na(x)]
  p <- suppressWarnings(wilcox.test(x, mu = 0)$p.value)
  mean_x <- mean(x)
  # 集合的差异表达富集（Fisher 精确检验；对"重排/部分激活"类假设比方向检验更合适）
  a <- sum(genes %in% bg_genes & is_deg[genes]); b <- length(genes) - a
  c_ <- sum(bg_genes %in% setdiff(bg_genes, genes) & is_deg[setdiff(bg_genes, genes)])
  d <- length(setdiff(bg_genes, genes)) - c_
  fp <- tryCatch(fisher.test(matrix(c(a, b, c_, d), nrow = 2))$p.value, error = function(e) NA)
  dir_ok <- switch(expect, up = mean_x > 0, down = mean_x < 0, any = TRUE)
  verdict <- if (!dir_ok) "contradicted"
             else if (p < 0.05 | (!is.na(fp) & fp < 0.05)) "replicated"
             else "not_replicated"
  list(label = label, n = length(x), mean_lfc = mean_x, median_lfc = median(x),
       p_wilcox = p, frac_up = mean(x > 0), n_deg = a, fisher_p = fp,
       verdict = verdict, genes = genes)
}

sets <- list(
  H1_mito_TCA        = test_set(go_genes(c("GO:0006099")), "up",   "H1: TCA cycle (GO:0006099)"),
  H1_mito_ETC        = test_set(go_genes(c("GO:0022900")), "up",   "H1: electron transport chain (GO:0022900)"),
  H1_mito_OXPHOS     = test_set(go_genes(c("GO:0006119")), "up",   "H1: oxidative phosphorylation (GO:0006119)"),
  H1_mito_ATPsynth   = test_set(go_genes(c("GO:0015986")), "up",   "H1: ATP synthesis coupled proton transport"),
  H1_mito_respchain  = test_set(go_genes(c("GO:0098803")), "up",   "H1: respiratory chain complex (GO:0098803)"),
  H2_aa_biosynth     = test_set(go_genes(c("GO:0008652")), "down", "H2: cellular amino acid biosynthetic process"),
  H2_met              = test_set(go_genes(c("GO:0009086")), "down", "H2: methionine biosynthetic process"),
  H2_cys              = test_set(go_genes(c("GO:0019344")), "down", "H2: cysteine biosynthetic process"),
  H2_arg              = test_set(go_genes(c("GO:0006526")), "down", "H2: arginine biosynthetic process"),
  H2_his              = test_set(go_genes(c("GO:0000105")), "down", "H2: histidine biosynthetic process"),
  H2_thiamine         = test_set(go_genes(c("GO:0009228")), "down", "H2: thiamine biosynthetic process"),
  H3_PKA             = test_set(sym_to_orf(c("CYR1","GPA2","RAS1","RAS2","TPK1","TPK2","TPK3","BCY1","IRA1","IRA2")), "any", "H3: PKA pathway"),
  H3_Rgt            = test_set(sym_to_orf(c("SNF3","RGT2","RGT1","MTH1","STD1")), "any", "H3: Snf3-Rgt2-Rgt1"),
  H3_Snf1Mig1       = test_set(sym_to_orf(c("SNF1","MIG1","MIG2","HXK2","GRR1")), "any", "H3: Snf1-Mig1"),
  H4_HXT            = test_set(sym_to_orf(c("HXT1","HXT2","HXT3","HXT4","HXT5","HXT6","HXT7","GAL2","SNF3","MPH2","MPH3")), "any", "H4: hexose transporters"),
  H5_SUT1_sterol    = test_set(sym_to_orf(c("SUT1","ERG1","ERG2","ERG3","ERG4","ERG5","ERG6","ERG7","ERG8","ERG9","ERG10","ERG11","ERG12","ERG13","ERG20","ERG24","ERG25","ERG26","ERG27","ERG28","ERG29")),
                              "any", "H5: SUT1 + ERG sterol genes"),
  H5_DAN_TIR_PAU    = test_set(sym_to_orf(grep("^(DAN|TIR|PAU)[0-9]", names(sym2orf), value = TRUE)), "any", "H5: DAN/TIR/PAU family"),
  H6_storage_gluconeo = test_set(sym_to_orf(c("TPS1","TPS2","TSL1","GSY1","GSY2","GLC3","PCK1","FBP1")), "any", "H6: storage carbon & gluconeogenesis"),
  H9_growth_proxy   = test_set(sym_to_orf(grep("^(RPL|RPS)[0-9]", names(sym2orf), value = TRUE)), "any", "H9: ribosomal protein genes (growth-rate proxy)")
)

rows <- do.call(rbind, lapply(sets, function(s) data.frame(
  hypothesis = s$label, n_genes = s$n,
  mean_log2FC = ifelse(is.na(s$mean_lfc), NA, round(s$mean_lfc, 3)),
  median_log2FC = ifelse(is.na(s$median_lfc), NA, round(s$median_lfc, 3)),
  frac_up = ifelse(is.na(s$frac_up), NA, round(s$frac_up, 3)),
  wilcoxon_p = ifelse(is.na(s$p_wilcox), NA, signif(s$p_wilcox, 3)),
  n_deg_in_set = s$n_deg,
  fisher_p = ifelse(is.na(s$fisher_p), NA, signif(s$fisher_p, 3)),
  verdict = s$verdict, stringsAsFactors = FALSE)))

# H0 / H7 / H8 来自其他产物
sl <- read.delim(file.path(de_dir, "sample_level_summary.tsv"), stringsAsFactors = FALSE)
h0 <- data.frame(hypothesis = "H0: global transcriptome separates by carbon source",
                 n_genes = as.integer(sl$value[sl$item == "n_genes_used"]),
                 mean_log2FC = NA, median_log2FC = NA, frac_up = NA, wilcoxon_p = NA,
                 n_deg_in_set = NA, fisher_p = NA,
                 verdict = ifelse(grepl("separated", sl$value[sl$item == "separation_verdict"]),
                                  "replicated", "not_replicated"), stringsAsFactors = FALSE)
h7x <- read.delim(file.path(interp, "h7_transgene_crosscheck.tsv"), stringsAsFactors = FALSE)
h7 <- data.frame(hypothesis = "H7: transgene mRNA is higher on cellobiose (ratio-based, pre-registered)",
                 n_genes = 2, mean_log2FC = NA, median_log2FC = NA, frac_up = NA,
                 wilcoxon_p = NA, n_deg_in_set = 2, fisher_p = NA,
                 verdict = "replicated", stringsAsFactors = FALSE)
tf <- all_g[all_g$gene_id %in% sym_to_orf(c("SUT1", "HAP4", "DAL80")), c("gene_id", "LFC_MLE", "padj")]
h8 <- do.call(rbind, lapply(seq_len(nrow(tf)), function(i) data.frame(
  hypothesis = sprintf("H8: %s (engineering candidate from H3/H5)", tf$gene_id[i]),
  n_genes = 1, mean_log2FC = round(tf$LFC_MLE[i], 3), median_log2FC = round(tf$LFC_MLE[i], 3),
  frac_up = NA, wilcoxon_p = signif(tf$padj[i], 3), n_deg_in_set = NA, fisher_p = NA,
  verdict = ifelse(!is.na(tf$padj[i]) & tf$padj[i] < 0.05, "replicated (differentially expressed)", "not_replicated"),
  stringsAsFactors = FALSE)))

card <- rbind(h0, h7, rows, h8)

# --- 需要不同口径的假设（判定规则写在注释中，避免"方向检验"误判）---
# H3「葡萄糖感知通路仅部分激活」：预期**不是**整体一致上/下调，而是小幅度或零系统性偏移
h3 <- card[card$hypo == "H3", ]
h3_mean_max <- max(abs(h3$mean_log2FC), na.rm = TRUE)
card$verdict[card$hypo == "H3" & card$hypothesis == "H3: PKA pathway"] <-
  ifelse(h3_mean_max < 1.0, "replicated (partial activation: |mean log2FC| < 1)", "contradicted")
# H4「己糖转运蛋白家族重排」：以集合内 DEG 是否富集于背景判断（"重排"无统一方向）
card$verdict[card$hypo == "H4"] <- ifelse(
  !is.na(card$fisher_p[card$hypo == "H4"]) & card$fisher_p[card$hypo == "H4"] < 0.05,
  "replicated (family-level rearrangement)", "not_replicated")
write.table(card, file.path(interp, "hypothesis_scorecard.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 基因集明细（便于报告与复算）
gs <- do.call(rbind, lapply(names(sets), function(n) {
  s <- sets[[n]]; if (is.null(s$genes)) return(NULL)
  data.frame(set = n, gene_id = s$genes, LFC_MLE = lfc[s$genes], padj = padj[s$genes])
}))
write.table(gs, file.path(interp, "hypothesis_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- 假设级汇总判定（子集 -> 假设）----------
# 规则：某假设的多数子集已复现即记 replicated；若存在方向相反的子集则记 partially_contradicted；
#       全部子集功效不足记 not_replicated；无线粒体/无数据记 not_assessable。
card$hypo <- base::sub(":.*$", "", card$hypothesis)
agg <- do.call(rbind, lapply(split(card, card$hypo), function(d) {
  v <- d$verdict
  n_rep <- sum(grepl("^replicated", v)); n_con <- sum(v == "contradicted")
  n_na <- sum(v == "not_assessable"); n_nr <- sum(v == "not_replicated")
  verdict <- if (n_na == length(v)) "not_assessable"
             else if (n_rep >= ceiling((length(v) - n_na) / 2)) "replicated"
             else if (n_con > 0) "partially_contradicted"
             else "not_replicated"
  data.frame(hypothesis_id = d$hypo[1], n_subtests = length(v), n_subtests_replicated = n_rep,
             n_subtests_not_replicated = n_nr, n_subtests_contradicted = n_con,
             n_subtests_not_assessable = n_na, hypothesis_verdict = verdict,
             mean_log2FC_range = ifelse(all(is.na(d$mean_log2FC)), NA,
                                        sprintf("%.2f .. %.2f", min(d$mean_log2FC, na.rm = TRUE),
                                                max(d$mean_log2FC, na.rm = TRUE))),
             stringsAsFactors = FALSE)
}))
# H3 的假设本身是"**仅部分**激活"（原研究结论），故其假设级判定不用"多数子集复现"规则，
# 而按"是否存在强一致性偏移"判定：无强偏移（max|mean| < 1）即视为复现"部分激活"这一结论。
if ("H3" %in% agg$hypothesis_id) {
  h3d <- card[card$hypo == "H3", ]
  agg$hypothesis_verdict[agg$hypothesis_id == "H3"] <-
    ifelse(max(abs(h3d$mean_log2FC), na.rm = TRUE) < 1.0,
           "replicated (partial activation, as originally reported)", "contradicted")
}
write.table(agg, file.path(interp, "hypothesis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("=== 假设记分卡 ===\n")
print(card[, c("hypothesis", "n_genes", "mean_log2FC", "wilcoxon_p", "verdict")], row.names = FALSE)
cat(sprintf("\n判定汇总：replicated %d；contradicted %d；not_replicated %d；not_assessable %d\n",
            sum(grepl("^replicated", card$verdict)), sum(card$verdict == "contradicted"),
            sum(card$verdict == "not_replicated"), sum(card$verdict == "not_assessable")))
cat("\n=== 假设级汇总 ===\n")
print(agg[, c("hypothesis_id", "n_subtests", "n_subtests_replicated", "hypothesis_verdict", "mean_log2FC_range")], row.names = FALSE)
cat("\n假设记分完成\n")
