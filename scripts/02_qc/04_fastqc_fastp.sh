#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 04_fastqc_fastp.sh   （PLAN 5.3）
# 目的：判定测序质量编码 -> 轻度修剪 -> 修剪前后质量对比报告
# 输入：data/raw/fastq/*.fastq.gz
# 输出：results/qc/fastqc/{raw,trimmed}/、results/qc/fastp/、results/qc/multiqc_report.html
#       data/processed/fastq/*.fastp.fastq.gz
#       docs/decisions/phred_encoding.md
# 说明：2014 年 GA II 数据可能是 Phred64（Illumina 1.3–1.7 管线），必须先判定再修剪，
#       误判会造成系统性"假低质量"与过度修剪。
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"

QC_RAW="$QC_DIR/fastqc/raw"; QC_TRIM="$QC_DIR/fastqc/trimmed"
FASTP_DIR="$QC_DIR/fastp"
SUB="$REFS_CUSTOM/adjudication/subsample"     # 复用裁决实验的 1 M reads 子采样做编码判定与修剪后 QC
mkdir -p "$QC_RAW" "$QC_TRIM" "$FASTP_DIR" "$PROC_FASTQ" "$DECISIONS_DIR"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/04_qc.log"; }
log "=== 读段 QC 开始 ==="

# --- 1. 用子采样 reads 做 FastQC，判定 Phred 编码 ---
log "FastQC（原始，子采样）"
for f in "$SUB"/*.1M.fq.gz; do
  [ -s "$QC_RAW/$(basename "${f%.fq.gz}")_fastqc.html" ] && continue
  fastqc -q -o "$QC_RAW" "$f" 2>>"$LOGS_DIR/04_qc.log"
done

# 从 FastQC 的 zip 中直接读取 Encoding（不依赖是否已解包；grep 失败不致脚本退出）
ENC="$(python3 - "$QC_RAW" <<'PYEOF'
import glob, os, sys, zipfile
vals = []
for z in glob.glob(os.path.join(sys.argv[1], "*_fastqc.zip")):
    with zipfile.ZipFile(z) as zf:
        for n in zf.namelist():
            if n.endswith("fastqc_data.txt"):
                for line in zf.read(n).decode("utf-8", "replace").splitlines():
                    if line.startswith("Encoding"):
                        vals.append(line.split("\t")[1].strip() if "\t" in line else line.split(None, 1)[1].strip())
from collections import Counter
print(Counter(vals).most_common(1)[0][0] if vals else "")
PYEOF
)"
ENC="${ENC:-Sanger / Illumina 1.9}"
log "FastQC 报告的全部编码取值：$(python3 -c "
import glob,os,zipfile
from collections import Counter
v=[]
for z in glob.glob('$QC_RAW/*_fastqc.zip'):
    with zipfile.ZipFile(z) as zf:
        for n in zf.namelist():
            if n.endswith('fastqc_data.txt'):
                for l in zf.read(n).decode('utf-8','replace').splitlines():
                    if l.startswith('Encoding'): v.append(l.split(chr(9))[1].strip())
print(dict(Counter(v)))
")"
log "判定质量编码：$ENC"
if printf '%s' "$ENC" | grep -q "1.5\|Illumina 1.3\|Illumina 1.4"; then
  PHRED_FLAG="--phred64"; PHRED_NOTE="Phred64（Illumina 1.3–1.7 管线）"
else
  PHRED_FLAG=""; PHRED_NOTE="Phred33（Sanger / Illumina 1.8+，或 fastqc 报告为 ASCII 33 起）"
fi

cat > "$DECISIONS_DIR/phred_encoding.md" <<EOF
# 决策记录：测序质量编码判定

**日期**：$(date '+%F %T')
**判定依据**：FastQC 对每样本 1 M reads 子采样的 \`Encoding\` 字段（见 \`results/qc/fastqc/raw/*/\`）

| 项 | 值 |
| --- | --- |
| FastQC 报告的 Encoding | \`$ENC\` |
| 结论 | $PHRED_NOTE |
| fastp 参数 | \`$PHRED_FLAG\`（空表示使用默认 Phred33） |

判定命令：\`grep '^Encoding' results/qc/fastqc/raw/*/fastqc_data.txt\`
EOF
log "编码判定写入 docs/decisions/phred_encoding.md"

# --- 2. fastp 轻度修剪（全长数据）---
log "fastp 修剪（单端，轻度：length_required=36；SE 接头检测为 fastp 默认行为）$PHRED_FLAG"
while IFS=$'\t' read -r sample_id gsm run condition replicate layout platform; do
  src="$RAW_FASTQ/$run.fastq.gz"
  out="$PROC_FASTQ/${run}_${sample_id}_${condition}.fastp.fastq.gz"
  [ -s "$src" ] || { log "!! 缺少 $src"; continue; }
  if [ -s "$out" ] && [ -s "$FASTP_DIR/${run}.json" ]; then log "已存在，跳过：$(basename "$out")"; continue; fi
  fastp -i "$src" -o "$out" \
    $PHRED_FLAG \
    --length_required 36 --qualified_quality_phred 20 --unqualified_percent_limit 40 \
    --cut_tail --cut_tail_mean_quality 20 \
    --thread 8 --json "$FASTP_DIR/${run}.json" --html "$FASTP_DIR/${run}.html" \
    >>"$LOGS_DIR/04_qc.log" 2>&1
  log "修剪完成 $run ($sample_id/$condition)"
done < <(tail -n +2 "$SAMPLESHEET")

# --- 3. 修剪后 FastQC（子采样 1 M reads，控制耗时）---
log "FastQC（修剪后，1 M 子采样）"
for f in "$PROC_FASTQ"/*.fastp.fastq.gz; do
  base="$(basename "${f%.fastp.fastq.gz}")"
  tmp="$QC_DIR/.tmp_${base}.1M.fq.gz"
  [ -s "$tmp" ] || seqtk sample -s100 "$f" 1000000 | gzip -c > "$tmp"
  [ -s "$QC_TRIM/${base}_fastqc.html" ] && continue
  fastqc -q -o "$QC_TRIM" "$tmp" 2>>"$LOGS_DIR/04_qc.log"
done
rm -f "$QC_DIR"/.tmp_*.1M.fq.gz

# --- 4. MultiQC 汇总 ---
log "MultiQC 汇总"
multiqc -q -f -o "$QC_DIR" -n multiqc_report "$QC_RAW" "$QC_TRIM" "$FASTP_DIR" \
  >>"$LOGS_DIR/04_qc.log" 2>&1
log "MultiQC 报告：results/qc/multiqc_report.html"

# --- 5. 修剪前后摘要表 ---
{
  printf 'run\tsample\tcondition\treads_before\treads_after\tpct_passed\tpct_low_quality\tpct_too_short\tadapter_pct\tq30_before\tq30_after\n'
  python3 - "$FASTP_DIR" "$SAMPLESHEET" <<'PY'
import json, os, sys
fastp_dir, sheet = sys.argv[1], sys.argv[2]
for line in open(sheet).read().splitlines()[1:]:
    f = line.split("\t")
    if len(f) < 4:
        continue
    sample, run, cond = f[0], f[2], f[3]
    p = os.path.join(fastp_dir, f"{run}.json")
    if not os.path.exists(p):
        continue
    d = json.load(open(p))
    s, fl = d.get("summary", {}), d.get("filtering_result", {})
    b, a = s.get("before_filtering", {}), s.get("after_filtering", {})
    tb = b.get("total_reads", 0) or 1
    print("\t".join(str(x) for x in [
        run, sample, cond, b.get("total_reads"), a.get("total_reads"),
        f"{100*a.get('total_reads',0)/tb:.2f}", fl.get("low_quality_reads", 0),
        fl.get("too_short_reads", 0), fl.get("passed_filter_reads", 0),
        f"{b.get('q30_rate',0):.4f}", f"{a.get('q30_rate',0):.4f}"]))
PY
} > "$QC_DIR/qc_summary.tsv"
log "修剪摘要：results/qc/qc_summary.tsv"
log "=== 读段 QC 结束 ==="
