#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 05_postmap_diagnostics.sh   （PLAN 5.6：比对诊断与 read fate）
# 目的：在没有全量 BAM 的前提下，为三类问题提供实测证据
#       ① 被排除的质粒骨架 reads 去了哪里；② rRNA/重复区占比；③ 线粒体转录本检出情况
# 方法：每样本 1 M reads 子采样（复用裁决实验的同一批子采样，保证不同分析口径可比）
#       一级：bowtie2 -> 宿主基因组 + 外源转录本（-k 2，统计多映射）
#       二级：一级未比对的 reads -> 诊断参考（pRS426 骨架 / 天然 2μ / TruSeq adapter）
# 输出：results/qc/postmap/{read_fate,chr_distribution,mt_gene_detection,rrdna_ty_hotspots}.tsv
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"
mkdir -p "$QC_DIR/postmap" "$REFS_INDEX/bowtie2" "$LOGS_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/05_postmap.log"; }
log "=== read fate 诊断开始 ==="

# --- 1. 一级参考：宿主基因组 + 外源转录本 ---
REF1="$REFS_INDEX/bowtie2/host_plus_transgenes.fa"
if [ ! -s "$REF1" ]; then
  log "构建一级参考（宿主基因组 + 外源转录本）"
  zcat "$REFS_HOST/${REFSEQ_ASM_DIR}_genomic.fna.gz" > "$REF1"
  cat "$REFS_CUSTOM/transgenes.fa" >> "$REF1"
fi
IDX1="$REFS_INDEX/bowtie2/host_plus_transgenes"
[ -s "${IDX1}.1.bt2" ] || bowtie2-build -q --threads "$THREADS" "$REF1" "$IDX1" >>"$LOGS_DIR/05_postmap.log" 2>&1

# --- 2. 二级参考：诊断序列（非基因组来源）---
DIAG="$REFS_DIAG/diagnostic_ref.fa"
if [ ! -s "$DIAG" ]; then
  log "构建二级诊断参考"
  {
    echo ">plasmid_pRS426_U03451"
    grep -v "^>" "$REFS_DIAG/plasmid_prs426_U03451.fa" | tr -d '\n'; echo
    echo ">plasmid_native_2micron_NC001398"
    grep -v "^>" "$REFS_DIAG/plasmid_2micron_NC001398.fa" | tr -d '\n'; echo
    echo ">adapter_TruSeq_read1"
    echo "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    echo ">adapter_TruSeq_read2"
    echo "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  } > "$DIAG"
fi
[ -s "${DIAG}.1.bt2" ] || bowtie2-build -q --threads "$THREADS" "$DIAG" "$DIAG" >>"$LOGS_DIR/05_postmap.log" 2>&1

# --- 3. 一级比对（复用 1 M 子采样）---
SUB="$REFS_CUSTOM/adjudication/subsample"
mkdir -p "$QC_DIR/postmap/bam" "$QC_DIR/postmap/unmapped"
while IFS=$'\t' read -r sample_id gsm run condition replicate layout platform; do
  fq="$SUB/${run}_${sample_id}_${condition}.1M.fq.gz"
  bam="$QC_DIR/postmap/bam/${run}_${sample_id}_${condition}.bam"
  [ -s "$fq" ] || { log "!! 缺少子采样 $fq"; continue; }
  if [ -s "$bam" ]; then log "已存在，跳过：$(basename "$bam")"; continue; fi
  bg="${run}_${sample_id}_${condition}"
  bowtie2 -x "$IDX1" -U "$fq" -p "$THREADS" -k 2 --sensitive --no-unal \
    --un-gz "$QC_DIR/postmap/unmapped/${bg}.fq.gz" \
    2> "$QC_DIR/postmap/bam/${run}.bowtie2.log" | samtools sort -@ 4 -o "$bam" -
  samtools index "$bam"
  log "$(basename "$bam") 比对完成：$(grep -oE '[0-9.]+% overall alignment rate' "$QC_DIR/postmap/bam/${run}.bowtie2.log" || true)"
done < <(tail -n +2 "$SAMPLESHEET")

# --- 4. 汇总指标 ---
log "汇总 read fate 指标"
python3 "$PROJECT_ROOT/scripts/02_qc/05a_read_fate_metrics.py" "$QC_DIR"/postmap/bam/*.bam \
  2>&1 | tee -a "$LOGS_DIR/05_postmap.log"
log "=== 诊断结束 ==="
