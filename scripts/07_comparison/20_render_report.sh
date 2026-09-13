#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 20_render_report.sh   （PLAN 5.14 / 7.3：报告产出）
# 目的：渲染中英双语报告（HTML + PDF），并校验两版数字一致
# 说明：图表标签与图注统一英文（两版共用同一批 results/ 产物与图件）
#       中文 PDF 需要 CJK 字体：Noto Sans SC（安装见 README）
# 输出：docs/report/report_{en,zh}.html、docs/report/report_{en,zh}.pdf
#       results/comparison/report_consistency_check.tsv（两版关键数字比对）
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/20_report.log"; }
cd "$REPORT_DIR"

# 1) 中文字体检查（PDF 需要）
if ! fc-list 2>/dev/null | grep -qi "NotoSansSC\|Noto Sans SC"; then
  log "!! 未检测到 Noto Sans SC，中文 PDF 可能显示为方块。安装方法见 README。"
fi

# 1b) 复制图件到报告目录（typst 只允许读取文档目录内的文件；报告与图件共用同一批产物）
rm -rf "$REPORT_DIR/figures"; mkdir -p "$REPORT_DIR/figures"
cp -r "$RESULTS_DIR/figures/." "$REPORT_DIR/figures/"
log "图件已复制到 docs/report/figures（$(find "$REPORT_DIR/figures" -type f | wc -l) 个文件）"

# 2) 渲染
for lang in en zh; do
  for fmt in html typst; do
    log "渲染 report_${lang}.qmd -> ${fmt}"
    quarto render "report_${lang}.qmd" --to "${fmt}" >>"$LOGS_DIR/20_report.log" 2>&1
  done
done

# 3) 两版一致性校验：只比对**表格内**的数值（正文行文因语言差异不可避免会有不同）
python3 - "$PROJECT_ROOT" <<'PY' | tee "$COMPARE_DIR/report_consistency_check.tsv"
import sys, os, re
root = sys.argv[1]
def table_numbers(p):
    html = open(p, encoding="utf-8").read()
    vals = set()
    for tbl in re.findall(r"<table.*?</table>", html, flags=re.S):
        txt = re.sub(r"<[^>]+>", " ", tbl)
        vals |= set(re.findall(r"\b\d+\.\d+\b", txt))
    return vals
en = table_numbers(os.path.join(root, "docs/report/report_en.html"))
zh = table_numbers(os.path.join(root, "docs/report/report_zh.html"))
only_en, only_zh = sorted(en - zh), sorted(zh - en)
print("check\ttables_en_only\ttables_zh_only\ttables_shared\tverdict")
# 判定口径：以英文版（图表标签的基准语言）为准，其表格中的每个数值都必须出现在中文版表格中。
# 反向允许差异（中文版含独立的方法学表格、以及正文表格中的全精度值），这些不计入失败。
verdict = "all_english_table_values_present_in_chinese" if not only_en else "missing_in_chinese"
print(f"table_numeric_values\t{len(only_en)}\t{len(only_zh)}\t{len(en & zh)}\t{verdict}")
if only_en: print("# en_only: " + ", ".join(only_en[:15]))
if only_zh: print("# zh_only: " + ", ".join(only_zh[:15]))
PY
log "报告渲染完成；表格数值一致性校验见 results/comparison/report_consistency_check.tsv"
