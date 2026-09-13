#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 07_transgene_adjudication.sh   （PLAN 3.4.3 第 1 步：序列裁决实验）
# 目的：用每样本 1 M reads 判定样本中实际存在的外源序列版本
#       —— 天然 gh1-1 vs Figure S8 的 gh1-1a；标签为 eGFP 还是 sfGFP。
# 方法：候选集（含阴性对照）建索引 -> bowtie2 -k 5 --very-sensitive-local
#       -> 每条 read 取 NM 最小的候选归属 -> 汇总为判别矩阵。
# 输出：refs/custom/adjudication/{candidates.fa,candidates.tsv,*.bam,*.bam.bai,subsample/}
#       refs/custom/sequence_adjudication.tsv
#       results/qc/transgene_assembly/adjudication_detail.tsv
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"

ADJ="$REFS_CUSTOM/adjudication"
SUB="$ADJ/subsample"
mkdir -p "$ADJ" "$SUB" "$LOGS_DIR"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/07_adjudication.log"; }
log "=== 外源序列裁决开始 ==="

# --- 1. 候选序列集（含 Figure S8 阴性对照）---
python3 "$PROJECT_ROOT/scripts/08_original_study/07a_build_adjudication_candidates.py" \
  2>&1 | tee -a "$LOGS_DIR/07_adjudication.log"
md5sum "$ADJ/candidates.fa" | tee -a "$LOGS_DIR/07_adjudication.log"

# --- 2. 每样本 1 M reads 子采样（并发 3）---
log "子采样 1 M reads/样本"
n=0
while IFS=$'\t' read -r sample_id gsm run condition replicate layout platform; do
  out="$SUB/${run}_${sample_id}_${condition}.1M.fq.gz"
  if [ -s "$out" ]; then log "已存在，跳过：$(basename "$out")"; continue; fi
  [ -s "$RAW_FASTQ/$run.fastq.gz" ] || { log "!! 缺少 $run.fastq.gz"; continue; }
  ( seqtk sample -s100 "$RAW_FASTQ/$run.fastq.gz" 1000000 | gzip -c > "$out" ) &
  n=$((n + 1))
  [ $((n % 3)) -eq 0 ] && wait
done < <(tail -n +2 "$SAMPLESHEET")
wait
ls -la "$SUB" | tail -7 | tee -a "$LOGS_DIR/07_adjudication.log"

# --- 3. 索引与比对 ---
log "构建裁决索引"
bowtie2-build -q --threads "$THREADS" "$ADJ/candidates.fa" "$ADJ/cand" >>"$LOGS_DIR/07_adjudication.log" 2>&1

log "比对（-k 5 --very-sensitive-local）"
for f in "$SUB"/*.fq.gz; do
  sample="$(basename "$f" .1M.fq.gz)"
  [ -s "$ADJ/$sample.bam" ] && { log "已存在，跳过：$sample.bam"; continue; }
  bowtie2 -x "$ADJ/cand" -U "$f" -p "$THREADS" -k 5 --very-sensitive-local --no-unal \
    2> "$ADJ/$sample.bowtie2.log" | samtools sort -@ 4 -o "$ADJ/$sample.bam" -
  samtools index "$ADJ/$sample.bam"
  m=$(grep -oE "[0-9.]+% overall alignment rate" "$ADJ/$sample.bowtie2.log" || true)
  log "$sample 比对完成（$m）"
done

# --- 4. 汇总判别矩阵 ---
log "汇总判别矩阵"
python3 "$PROJECT_ROOT/scripts/08_original_study/07b_summarize_adjudication.py" \
  2>&1 | tee -a "$LOGS_DIR/07_adjudication.log"
log "=== 裁决结束 ==="
