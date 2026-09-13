#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 03_fetch_metadata.sh
# 目的：把附录 A 引用的外部事实以原始记录形式落库，使元数据结论离线可复现。
# 输出：data/metadata/geo_soft/{GSE54825,GSM*.txt}
#       data/metadata/original_study/geo/GSE54825_Cellobiose_versus_Glucose.xlsx.gz
#       checksums/metadata.md5
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"

mkdir -p "$META_GEOSOFT" "$META_ORIG_GEO" "$LOGS_DIR" "$CHECKSUM_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/03_metadata.log"; }

# --- 1. GEO SOFT：系列 + 6 个样本 ---
curl -fsSL --retry 3 -o "$META_GEOSOFT/${GEO_SERIES}.txt" \
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GEO_SERIES}&targ=self&form=text&view=brief"
log "获取 ${GEO_SERIES}.txt"

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_id gsm run condition replicate layout platform; do
  out="$META_GEOSOFT/${gsm}.txt"
  [ -s "$out" ] || curl -fsSL --retry 3 -o "$out" \
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}&targ=self&form=text&view=brief"
  log "获取 ${gsm}.txt ($sample_id, $condition)"
done

# --- 2. GEO 补充文件（逐基因 FC / FDR）---
xlsx="GSE54825_Cellobiose_versus_Glucose.xlsx.gz"
[ -s "$META_ORIG_GEO/$xlsx" ] && [ -s "$META_ORIG_GEO/$xlsx" ] || \
  curl -fsSL --retry 3 -o "$META_ORIG_GEO/$xlsx" "$GEO_SUPPL_BASE/$xlsx"
log "获取 $xlsx"

# --- 3. 校验值 ---
find "$META_GEOSOFT" "$META_ORIG_GEO" -type f -print0 | sort -z | xargs -0 md5sum \
  | sed "s|$PROJECT_ROOT/||" > "$CHECKSUM_DIR/metadata.md5"
log "md5 写入 checksums/metadata.md5（$(wc -l < "$CHECKSUM_DIR/metadata.md5") 个文件）"
