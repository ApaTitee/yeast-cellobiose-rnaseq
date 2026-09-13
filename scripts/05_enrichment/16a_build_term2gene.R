#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 16a_build_term2gene.R   （PLAN 5.10 的前置：注释映射表构建）
# 目的：构建两套 TERM2GENE 映射表，供富集分析（16）与假设记分卡（19）共用，保证口径一致
#   (1) GO slim    —— SGD go_slim_mapping.tab（SGD 规范；本身即层级化的 slim 映射）
#   (2) 全 GO      —— SGD GAF 经 **GO 层级传播**（true path rule）到全部祖先术语
# 为何需要传播：SGD GAF 只记录最具体术语的注释（例如 GO:0008652「细胞氨基酸生物合成」
#   直接注释仅 3 行），不传播会使宽泛术语的基因集严重偏小，导致集合检验与 ORA 失真。
# 输出：refs/host/annotation/term2gene_goslim.tsv
#       refs/host/annotation/term2gene_go_propagated.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(GO.db); library(AnnotationDbi) })
.cmdargs <- commandArgs(trailingOnly = FALSE)
.script_path <- base::sub("^--file=", "", .cmdargs[grep("^--file=", .cmdargs)])
root <- if (length(.script_path)) normalizePath(file.path(dirname(.script_path), "..", "..")) else normalizePath(".")
ann <- file.path(root, "refs", "host", "annotation")
ORS <- "^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$"

# ---------- 1. GO slim ----------
slim <- read.delim(file.path(ann, "go_slim_mapping.tab"), header = FALSE, stringsAsFactors = FALSE)
slim <- slim[grepl(ORS, slim$V1, perl = TRUE), c("V1", "V5", "V6")]
t2g_slim <- unique(data.frame(term = slim$V6, term_name = slim$V5, gene = slim$V1,
                              aspect = "slim", stringsAsFactors = FALSE))
write.table(t2g_slim, file.path(ann, "term2gene_goslim.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("GO slim：%d 术语 / %d 基因 -> term2gene_goslim.tsv\n",
            length(unique(t2g_slim$term)), length(unique(t2g_slim$gene))))

# ---------- 2. 全 GO + 层级传播 ----------
gaf <- read.delim(gzfile(file.path(ann, "sgd.gaf.gz")), header = FALSE, comment.char = "!",
                  stringsAsFactors = FALSE, quote = "")
orf <- vapply(strsplit(gaf$V11, "|", fixed = TRUE), function(x) {
  h <- grep(ORS, x, value = TRUE, perl = TRUE)
  if (length(h)) h[1] else NA_character_
}, character(1))
direct <- unique(data.frame(go = gaf$V5, aspect = gaf$V9, gene = orf, stringsAsFactors = FALSE))
direct <- direct[!is.na(direct$gene), ]
cat(sprintf("GAF 直接注释：%d 个 GO 术语 / %d 基因\n",
            length(unique(direct$go)), length(unique(direct$gene))))

anc <- list(BP = as.list(GOBPANCESTOR), MF = as.list(GOMFANCESTOR), CC = as.list(GOCCANCESTOR))
expanded <- do.call(rbind, lapply(split(direct, direct$aspect), function(d) {
  key <- switch(d$aspect[1], P = "BP", F = "MF", C = "CC", NULL)
  if (is.null(key)) return(NULL)
  a <- anc[[key]]
  do.call(rbind, lapply(unique(d$gene), function(g) {
    gos <- unique(d$go[d$gene == g])
    data.frame(go = unique(c(gos, unlist(a[gos], use.names = FALSE))), gene = g,
               stringsAsFactors = FALSE)
  }))
}))
expanded <- unique(expanded[!is.na(expanded$go), ])
names_tbl <- suppressMessages(AnnotationDbi::select(GO.db, keys = unique(expanded$go),
                                                   columns = "TERM", keytype = "GOID"))
expanded$term_name <- names_tbl$TERM[match(expanded$go, names_tbl$GOID)]
t2g_go <- unique(data.frame(term = expanded$go,
                            term_name = ifelse(is.na(expanded$term_name), expanded$go, expanded$term_name),
                            gene = expanded$gene, aspect = "propagated", stringsAsFactors = FALSE))
write.table(t2g_go, file.path(ann, "term2gene_go_propagated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("全 GO（传播后）：%d 术语 / %d 基因 -> term2gene_go_propagated.tsv\n",
            length(unique(t2g_go$term)), length(unique(t2g_go$gene))))

# 关键术语的规模抽查（用于验证传播是否生效）
for (g in c("GO:0008652", "GO:0006099", "GO:0009228", "GO:0022900")) {
  nm <- names_tbl$TERM[match(g, names_tbl$GOID)]
  cat(sprintf("  %s %-46s 基因数 %d\n", g, ifelse(is.na(nm), "", substr(nm, 1, 46)),
              length(unique(t2g_go$gene[t2g_go$term == g]))))
}
