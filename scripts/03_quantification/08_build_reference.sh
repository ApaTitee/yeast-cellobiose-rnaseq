#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 08_build_reference.sh   （PLAN 5.2 / 3.5）
# 目的：构建自定义参考体系（宿主转录本 + 外源转录本 + decoy）与 ID 映射损失表
# 说明：宿主转录本用 gffread 从 GFF 抽取（RefSeq 的 rna_from_genomic 不含线粒体蛋白编码 mRNA，
#       直接使用会使 5.6 的 mtDNA 检出判定失效）。
# 输出：refs/host/host_transcripts_gffread.fa
#       refs/custom/{transcripts.fa,transgenes.fa,tx2gene.tsv,decoys.txt,gentrome.fa}
#       refs/checks/id_map_loss.tsv
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"
mkdir -p "$REFS_CUSTOM" "$REFS_CHECKS" "$LOGS_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/08_reference.log"; }
log "=== 构建自定义参考 ==="

# --- 1. gffread 抽取宿主转录本（含线粒体 mRNA）---
HOST_TR="$REFS_HOST/host_transcripts_gffread.fa"
if [ ! -s "$HOST_TR" ]; then
  log "gffread 抽取宿主转录本"
  TMP_FNA="$REFS_HOST/.tmp_genomic.fna"; TMP_GFF="$REFS_HOST/.tmp_genomic.gff"
  zcat "$REFS_HOST/${REFSEQ_ASM_DIR}_genomic.fna.gz" > "$TMP_FNA"
  zcat "$REFS_HOST/${REFSEQ_ASM_DIR}_genomic.gff.gz" > "$TMP_GFF"
  gffread "$TMP_GFF" -g "$TMP_FNA" -w "$HOST_TR" >>"$LOGS_DIR/08_reference.log" 2>&1
  rm -f "$TMP_FNA" "$TMP_GFF" "$TMP_FNA.fai"
fi
log "宿主转录本: $(grep -c '^>' "$HOST_TR") 条"

# --- 2. 组装参考 ---
python3 "$PROJECT_ROOT/scripts/03_quantification/08a_build_custom_reference.py" 2>&1 | tee -a "$LOGS_DIR/08_reference.log"

# --- 3. 自检与校验值 ---
{
  echo "--- 产物自检 ---"
  python3 - "$REFS_CUSTOM/transcripts.fa" <<'PY'
import sys
seqs, cur = {}, None
for line in open(sys.argv[1], encoding="utf-8"):
    if line.startswith(">"):
        cur = line[1:].split()[0]; seqs[cur] = []
    elif cur:
        seqs[cur].append(line.strip())
for k in ("cdt-1EGFP", "gh1-1"):
    s = "".join(seqs.get(k, []))
    print(f"  {k}: {len(s)} bp, 读框 {'OK' if len(s) % 3 == 0 else '错误'}")
print(f"  转录本总条目: {len(seqs)}")
PY
  md5sum "$REFS_CUSTOM/transcripts.fa" "$REFS_CUSTOM/tx2gene.tsv" "$REFS_CUSTOM/gentrome.fa"
} 2>&1 | tee -a "$LOGS_DIR/08_reference.log"
log "=== 构建完成 ==="
