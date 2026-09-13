#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 10_salmon_quant.sh   （PLAN 5.4 + 5.5 的参数敏感性）
# 目的：对 6 个样本做 salmon 定量（主用 k=31；随后跑 k=25 作 50 bp 短读的参数敏感性对照）
# 输入：data/processed/fastq/*.fastp.fastq.gz、refs/index/salmon/k{31,25}_decoy_host
# 输出：results/quantification/salmon{,_k25}/<run>_<sample>_<condition>/{quant.sf,lib_format_counts.json,logs}
#       results/qc/salmon_libtype.tsv（逐样本 libType 推断与比对率）
# 参数：--libType A（自动推断，单端）、--validateMappings、--gcBias、--seqBias、--dumpEq
#       并发 2（每个 8 线程，控制内存峰值）
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"
MAXJOBS="${MAXJOBS:-2}"
mkdir -p "$QUANT_DIR" "$LOGS_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/10_quant.log"; }

run_set() {           # $1 = kmer (31|25)
  local K="$1" IDX="$SALMON_INDEX/k${K}_decoy_host" OUT="$QUANT_DIR/salmon"
  [ "$K" = "25" ] && OUT="$QUANT_DIR/salmon_k25"
  mkdir -p "$OUT"
  local n=0
  while IFS=$'\t' read -r sample_id gsm run condition replicate layout platform; do
    local fq="$PROC_FASTQ/${run}_${sample_id}_${condition}.fastp.fastq.gz"
    local out="$OUT/${run}_${sample_id}_${condition}"
    [ -s "$fq" ] || { log "!! 缺少 $fq"; continue; }
    if [ -s "$out/quant.sf" ] && [ ! -f "$out/.stale" ]; then log "[k$K] 已存在，跳过：$(basename "$out")"; continue; fi
    # 单端数据的片段长度 = 读长；必须给出 FLD 先验，否则 salmon 2.7 会退回默认 250±25，
    # 导致 EffectiveLength 失真（实测 2472 bp 转录本报 40563，TPM 不可用）。
    FLD_MEAN="$(python3 - "$FASTP_DIR/${run}.json" <<'PYEOF'
import json, sys
d = json.load(open(sys.argv[1]))
print(int(round(d["summary"]["after_filtering"].get("read1_mean_length", 46))))
PYEOF
)"
    ( salmon quant -i "$IDX" -l A -r "$fq" -p 8 \
        --validateMappings --gcBias --seqBias --dumpEq \
        --fldMean "$FLD_MEAN" --fldSD 3 \
        -o "$out" > "$out.log" 2>&1 ) &
    n=$((n + 1))
    [ $((n % MAXJOBS)) -eq 0 ] && wait
  done < <(tail -n +2 "$SAMPLESHEET")
  wait
  log "[k$K] 本组定量完成"
}

for K in 31 25; do
  log "=== salmon 定量 k=$K 开始 ==="
  run_set "$K"
  log "=== salmon 定量 k=$K 结束 ==="
done

# --- 汇总 libType 与比对率（salmon 2.7 rust 版的 JSON 键名与 1.x 不同，故以日志为准）---
{
  printf 'kmer\tsample\tlibType\tobserved\tmapped\tpct_mapped\tmean_frag_len\tn_ec\tmultimap_or_orphan\n'
  for K in 31 25; do
    OUT="$QUANT_DIR/salmon"; [ "$K" = "25" ] && OUT="$QUANT_DIR/salmon_k25"
    for d in "$OUT"/*/; do
      [ -s "$d/quant.sf" ] || continue
      python3 - "$d" "$K" "$(basename "$d")" <<'PYEOF'
import json, os, re, sys
d, k, sample = sys.argv[1], sys.argv[2], sys.argv[3]
fmt = "?"
j = os.path.join(d, "lib_format_counts.json")
if os.path.exists(j):
    fmt = json.load(open(j)).get("expected_format", "?")
obs = mapped = rate = fld = ec = orph = ""
lg = os.path.join(d, "logs", "salmon_quant.log")
if os.path.exists(lg):
    t = open(lg).read()
    def g(pat):
        m = re.search(pat, t)
        return m.group(1) if m else ""
    obs = g(r"observed fragments:\s*(\d+)")
    mapped = g(r"mapped fragments:\s*(\d+)")
    rate = g(r"mapping rate:\s*([\d.]+)%")
    fld = g(r"fragment length mean \(sd\):\s*([\d.]+)")
    ec = g(r"number of equivalence classes:\s*(\d+)")
    orph = g(r"fragments mapped to multiple transcripts[^\d]*(\d+)")
print("\t".join([k, sample, fmt, obs, mapped, rate, fld, ec, orph]))
PYEOF
    done
  done
} > "$QC_DIR/salmon_libtype.tsv"
log "libType/比对率汇总：results/qc/salmon_libtype.tsv"
log "=== 定量全部结束 ==="
