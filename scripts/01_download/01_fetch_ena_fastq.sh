#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_fetch_ena_fastq.sh
# 目的：从 ENA 下载 samplesheet 中全部 run 的 FASTQ，并以 ENA 公布的 md5 校验。
# 输入：config/samplesheet.tsv
# 输出：data/raw/fastq/*.fastq.gz
#       data/metadata/ena/ena_filereport.tsv
#       checksums/raw.md5                （自算值）
#       checksums/ena_published.md5      （ENA 公布值）
#       checksums/md5_verification.tsv   （逐文件比对结论）
# 幂等：自算 md5 == ENA md5 时跳过；部分文件用 curl -C - 续传。
# 并发：MAXJOBS（默认 3）——ENA 单连接限速，多文件并发可提升总吞吐。
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"

MAXJOBS="${MAXJOBS:-3}"
FIELDS="run_accession,experiment_accession,sample_accession,sample_title,library_layout,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
FILEREPORT="$META_ENA/ena_filereport.tsv"
VDIR="$CHECKSUM_DIR/.verify"

mkdir -p "$META_ENA" "$RAW_FASTQ" "$LOGS_DIR" "$CHECKSUM_DIR" "$VDIR"
rm -f "$VDIR"/*.tsv   # 清理上一轮的陈旧校验记录，避免汇总表混入过期结果

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/01_download.log"; }

# 互斥锁：防止两个实例并发写同一个 FASTQ（并发写会损坏文件）
exec 9>"$LOGS_DIR/.01_fetch_ena_fastq.lock"
if ! flock -n 9; then log "另一个下载实例正在运行，退出"; exit 1; fi

# --- 1. ENA filereport（含官方 md5、read_count、字节数）---
if [ ! -s "$FILEREPORT" ]; then
  log "获取 ENA filereport: $ENA_PROJECT"
  curl -fsSL --retry 3 --retry-delay 5 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_PROJECT}&result=read_run&fields=${FIELDS}&format=tsv" \
    -o "$FILEREPORT"
fi

# --- 2. 单个 run 的下载 + 校验 ---
fetch_one() {
  local sample_id="$1" run="$2" cond="$3"
  local row ftp md5_ena n_reads bytes_ena url fname dest md5_self bytes status
  row="$(awk -F'\t' -v r="$run" '$1==r {print; exit}' "$FILEREPORT")"
  if [ -z "$row" ]; then log "!! $run 不在 ENA filereport 中，跳过"; return 0; fi
  ftp="$(printf '%s' "$row" | cut -f9 | cut -d';' -f1)"
  md5_ena="$(printf '%s' "$row" | cut -f10 | cut -d';' -f1)"
  n_reads="$(printf '%s' "$row" | cut -f7)"
  bytes_ena="$(printf '%s' "$row" | cut -f11 | cut -d';' -f1)"
  if [ -z "$ftp" ]; then log "!! $run 无 FASTQ 直链，跳过"; return 0; fi

  url="https://${ftp}"
  fname="$(basename "$ftp")"
  dest="$RAW_FASTQ/$fname"

  if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5_ena" ]; then
    log "校验通过，跳过下载：$fname"
  else
    log "下载 $run ($cond) -> $fname  reads=$n_reads  ENA_bytes=${bytes_ena:-NA}"
    # --speed-limit/--speed-time：连接速率低于 100 KB/s 持续 60 s 即中断，交给 --retry 重连
    # （应对 ENA 偶发的单连接限速；实测曾出现 54 KB/s 的停滞连接）
    curl -fL --retry 8 --retry-delay 10 --retry-all-errors \
         --speed-limit 102400 --speed-time 60 \
         -C - -sS --no-progress-meter -o "$dest" "$url" \
      2>>"$LOGS_DIR/01_download.err.log" || true
  fi

  md5_self="$(md5sum "$dest" | cut -d' ' -f1)"
  bytes="$(stat -c%s "$dest")"
  status="FAIL"; [ "$md5_self" = "$md5_ena" ] && status="PASS"
  # 额外一致性检查：字节数应与 ENA 公布的 fastq_bytes 相等
  if [ "$status" = "PASS" ] && [ -n "$bytes_ena" ] && [ "$bytes" != "$bytes_ena" ]; then
    status="FAIL_SIZE"; log "!! $fname 字节数 $bytes != ENA 公布 $bytes_ena"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$run" "$sample_id" "$cond" "$fname" "$bytes" "$md5_self" "$md5_ena" "$status" \
    > "$VDIR/$run.tsv"
  log "$status $fname ($(numfmt --to=iec "$bytes" 2>/dev/null || echo "$bytes"))"
}
export -f fetch_one
export FILEREPORT RAW_FASTQ LOGS_DIR CHECKSUM_DIR VDIR

# --- 3. 并发调度（每批 MAXJOBS 个）---
log "开始下载（并发 $MAXJOBS）"
batch=0
while IFS=$'\t' read -r sample_id gsm run condition replicate layout platform; do
  fetch_one "$sample_id" "$run" "$condition" &
  batch=$((batch + 1))
  if [ $((batch % MAXJOBS)) -eq 0 ]; then wait; fi
done < <(tail -n +2 "$SAMPLESHEET")
wait

# --- 4. 汇总与校验值 ---
{
  printf 'run\tsample_id\tcondition\tfile\tbytes\tmd5_self\tmd5_ena\tstatus\n'
  for f in "$VDIR"/*.tsv; do [ -s "$f" ] && cat "$f"; done
} > "$CHECKSUM_DIR/md5_verification.tsv"

: > "$CHECKSUM_DIR/raw.md5"; : > "$CHECKSUM_DIR/ena_published.md5"
while IFS=$'\t' read -r run sample_id cond fname bytes md5_self md5_ena status; do
  [ "$run" = "run" ] && continue
  printf '%s  %s\n' "$md5_self" "$fname" >> "$CHECKSUM_DIR/raw.md5"
  printf '%s  %s\n' "$md5_ena"  "$fname" >> "$CHECKSUM_DIR/ena_published.md5"
done < "$CHECKSUM_DIR/md5_verification.tsv"

n_all=$(( $(wc -l < "$CHECKSUM_DIR/md5_verification.tsv") - 1 ))
n_fail=$(awk -F'\t' 'NR>1 && $8!="PASS"' "$CHECKSUM_DIR/md5_verification.tsv" | wc -l)
log "校验汇总：$((n_all - n_fail))/$n_all 通过"
[ "$n_fail" -eq 0 ] || { log "!! 存在校验失败文件，见 checksums/md5_verification.tsv（重跑本脚本会自动重下或续传）"; exit 1; }
log "完成。自算值 checksums/raw.md5，ENA 公布值 checksums/ena_published.md5"
