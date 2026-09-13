#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 09_salmon_index.sh   （PLAN 5.4）
# 目的：构建 decoy-aware 的 Salmon 索引（k=31 主用，k=25 供 50 bp 短读敏感性分析）
# 输入：refs/custom/gentrome.fa（转录本 + 宿主基因组）、refs/custom/decoys.txt
# 输出：refs/index/salmon/{k31_decoy_host, k25_decoy_host}/
# 记录：索引参数与 salmon 版本写入 refs/VERSIONS.md 与日志
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"
mkdir -p "$SALMON_INDEX" "$LOGS_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/09_index.log"; }
log "=== 构建 Salmon 索引 ==="
log "salmon 版本: $(salmon --version 2>&1 | head -1)"
for K in 31 25; do
  IDX="$SALMON_INDEX/k${K}_decoy_host"
  if [ -s "$IDX/info.json" ]; then log "已存在，跳过：k${K}"; continue; fi
  log "构建 k=${K} 索引 -> $IDX"
  salmon index -t "$REFS_CUSTOM/gentrome.fa" -d "$REFS_CUSTOM/decoys.txt" \
    -i "$IDX" -k "$K" -p "$THREADS" >>"$LOGS_DIR/09_index.log" 2>&1
  log "k=${K} 完成（索引大小 $(du -sh "$IDX" | cut -f1)）"
done
{
  echo
  echo "## Salmon 索引"
  echo
  echo "- salmon 版本：$(salmon --version 2>&1 | head -1)"
  echo "- 输入：\`refs/custom/gentrome.fa\`（md5 $(md5sum "$REFS_CUSTOM/gentrome.fa" | cut -d' ' -f1)）+ \`refs/custom/decoys.txt\`（$(wc -l < "$REFS_CUSTOM/decoys.txt") 条宿主基因组序列作 decoy）"
  echo "- 索引：k=31（主用，对应 ~50 bp reads）与 k=25（敏感性对照）"
  echo "- decoy 中不含质粒骨架与 URA3（见 PLAN 3.5 的槽位说明）"
} >> "$REFS_DIR/VERSIONS.md"
log "=== 索引构建完成 ==="
