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
    if [ -s "$out/quant.sf" ]; then log "[k$K] 已存在，跳过：$(basename "$out")"; continue; fi
    ( salmon quant -i "$IDX" -l A -r "$fq" -p 8 \
        --validateMappings --gcBias --seqBias --dumpEq \
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

# --- 汇总 libType 与比对率 ---
{
  printf 'kmer\tsample\tlibType\tnum_processed\tnum_mapped\tpct_mapped\tnum_decoy\tpct_decoy\tnum_ec\n'
  for K in 31 25; do
    OUT="$QUANT_DIR/salmon"; [ "$K" = "25" ] && OUT="$QUANT_DIR/salmon_k25"
    for d in "$OUT"/*/; do
      j="$d/lib_format_counts.json"
      [ -s "$j" ] || continue
      python3 - "$j" "$K" "$(basename "$d")" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
p = d.get("num_processed", 0) or 1
print("\t".join(str(x) for x in [
    sys.argv[2], sys.argv[3], d.get("expected_format", "?"), d.get("num_processed"),
    d.get("num_mapped"), f"{100*d.get('num_mapped',0)/p:.2f}",
    d.get("num_decoy_fragments", 0), f"{100*d.get('num_decoy_fragments',0)/p:.2f}",
    d.get("num_ec_fragments", 0)]))
PY
    done
  done
} > "$QC_DIR/salmon_libtype.tsv"
log "libType/比对率汇总：results/qc/salmon_libtype.tsv"
log "=== 定量全部结束 ==="
