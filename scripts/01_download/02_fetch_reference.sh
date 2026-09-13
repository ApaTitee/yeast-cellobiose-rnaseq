#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02_fetch_reference.sh
# 目的：冻结宿主参考、功能注释与外源候选/诊断序列，并生成 refs/VERSIONS.md。
# 原则：URL 写死版本，禁止 "latest" 类动态地址；全部文件记录 md5；NCBI 文件与官方
#       md5checksums.txt 比对，其余记录自算 md5。
# 输出：
#   refs/host/           RefSeq GCF_000146045.2_R64（genomic/rna/cds/gff/assembly_report）
#   refs/host/annotation sgd.gaf.gz、goslim_yeast.obo、go_slim_mapping.tab
#   refs/diagnostic/     外源候选序列 + 质粒骨架 + 天然 2μ（供 3.4 裁决与 5.6 read fate）
#   refs/VERSIONS.md     版本与校验值记录（附录 B 的凭证）
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"

ANNO_DIR="$REFS_HOST/annotation"
mkdir -p "$REFS_HOST" "$ANNO_DIR" "$REFS_DIAG" "$LOGS_DIR" "$CHECKSUM_DIR"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOGS_DIR/02_reference.log"; }

# ===========================================================================
# 1. 宿主参考（RefSeq，固定组装版本）
# ===========================================================================
REFSEQ_FILES=(
  "${REFSEQ_ASM_DIR}_genomic.fna.gz"
  "${REFSEQ_ASM_DIR}_rna_from_genomic.fna.gz"
  "${REFSEQ_ASM_DIR}_cds_from_genomic.fna.gz"
  "${REFSEQ_ASM_DIR}_genomic.gff.gz"
  "${REFSEQ_ASM_DIR}_assembly_report.txt"
  "md5checksums.txt"
)
for f in "${REFSEQ_FILES[@]}"; do
  [ -s "$REFS_HOST/$f" ] || { log "下载 $f"; curl -fsSL --retry 3 -o "$REFS_HOST/$f" "$NCBI_GENOMES_BASE/$f"; }
done

# 与 NCBI 官方 md5 比对（md5checksums.txt 内为相对文件名）
log "校验 RefSeq 文件（对照 NCBI md5checksums.txt）"
(cd "$REFS_HOST" && grep -E "$(printf '%s|' "${REFSEQ_FILES[@]}" | sed 's/|$//')" md5checksums.txt > .md5sub && md5sum -c .md5sub) \
  | tee -a "$LOGS_DIR/02_reference.log"
rm -f "$REFS_HOST/.md5sub"

# 注释发布信息（写死记录，供 VERSIONS.md）
ANNO_JSON="$REFS_HOST/annotation_release.json"
[ -s "$ANNO_JSON" ] || curl -fsSL --retry 3 -o "$ANNO_JSON" \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${REFSEQ_ASM}/dataset_report"

# ===========================================================================
# 2. 功能注释（GO / SGD 官方文件）
# ===========================================================================
log "下载功能注释文件"
[ -s "$ANNO_DIR/sgd.gaf.gz" ] || curl -fsSL --retry 3 -o "$ANNO_DIR/sgd.gaf.gz" \
  "https://current.geneontology.org/annotations/sgd.gaf.gz"
[ -s "$ANNO_DIR/goslim_yeast.obo" ] || curl -fsSL --retry 3 -o "$ANNO_DIR/goslim_yeast.obo" \
  "https://current.geneontology.org/ontology/subsets/goslim_yeast.obo"
[ -s "$ANNO_DIR/go_slim_mapping.tab" ] || curl -fsSL --retry 3 -o "$ANNO_DIR/go_slim_mapping.tab" \
  "https://downloads.yeastgenome.org/curation/literature/go_slim_mapping.tab"

# ===========================================================================
# 3. 外源候选序列与诊断序列
# ===========================================================================
efetch_nt() {  # $1=accession  $2=rettype  $3=输出文件
  [ -s "$3" ] || curl -fsSL --retry 3 -o "$3" \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$1&rettype=$2&retmode=text"
}
log "下载外源候选与诊断序列"
efetch_nt XM_958708.2     fasta_cds_na "$REFS_DIAG/ncrassa_cdt-1_NCU00801_candidate_CDS.fa"
efetch_nt XM_011395456.1  fasta_cds_na "$REFS_DIAG/ncrassa_gh1-1_NCU00130_candidate_CDS.fa"
efetch_nt PX636966.1      fasta_cds_na "$REFS_DIAG/sfgfp_candidate_PX636966.fa"
efetch_nt MW132720.1      fasta_cds_na "$REFS_DIAG/sfgfp_candidate_MW132720.fa"
efetch_nt JQ341914.1      fasta_cds_na "$REFS_DIAG/sfgfp_candidate_JQ341914.fa"
efetch_nt PX636966.1      fasta        "$REFS_DIAG/sfgfp_record_PX636966.fa"
efetch_nt U55762.1        fasta_cds_na "$REFS_DIAG/egfp_candidates_U55762.fa"
efetch_nt U55762.1        gb           "$REFS_DIAG/pegfpN1_record_U55762.gb"
efetch_nt U03451.1        fasta        "$REFS_DIAG/plasmid_prs426_U03451.fa"
efetch_nt U03451.1        gb           "$REFS_DIAG/plasmid_prs426_U03451.gb"
efetch_nt NC_001398.1     fasta        "$REFS_DIAG/plasmid_2micron_NC001398.fa"

# 候选序列裁决前的完整性自检：长度 + 翻译起点
python3 - "$REFS_DIAG" <<'PY'
import os, re, sys
d = sys.argv[1]
codon = {}
b = "TCAG"; aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
k = 0
for x in b:
    for y in b:
        for z in b:
            codon[x + y + z] = aas[k]; k += 1
tr = lambda s: "".join(codon.get(s[i:i+3], "X") for i in range(0, len(s) - 2, 3))
expect = {
    "ncrassa_cdt-1_NCU00801_candidate_CDS.fa": (1740, None),
    "ncrassa_gh1-1_NCU00130_candidate_CDS.fa": (1431, None),
    "sfgfp_candidate_PX636966.fa":             (717, None),   # 该记录 N 端为 MRKGEELFTGVV
    "sfgfp_candidate_MW132720.fa":             (717, None),   # 该记录 N 端为 MSKGEELFTGVV
    "sfgfp_candidate_JQ341914.fa":             (717, None),
    "egfp_candidates_U55762.fa":               (None, "MVSKGEELFTGV"),
    "plasmid_prs426_U03451.fa":                (5726, None),
    "plasmid_2micron_NC001398.fa":             (6318, None),
}
print("  --- 序列自检 ---")
for fn, (ln, aa) in expect.items():
    p = os.path.join(d, fn)
    if not os.path.exists(p):
        print(f"  {fn}: 缺失"); continue
    recs = [r for r in open(p).read().split(">") if r.strip()]
    for r in recs:
        head, *body = r.split("\n")
        s = "".join(body).upper()
        ok_len = "" if ln is None else ("len=OK" if len(s) == ln else f"len={len(s)}(期望{ln})")
        ok_aa = "" if aa is None else ("start=OK" if tr(s).startswith(aa) else f"start={tr(s)[:12]}")
        print(f"  {fn}: {len(s)} bp {ok_len} {ok_aa} | {head[:70]}")
PY

# ===========================================================================
# 4. VERSIONS.md（附录 B 凭证）
# ===========================================================================
ANNO_NAME="$(python3 -c "
import json
d=json.load(open('$ANNO_JSON'))
a=d['reports'][0].get('annotation_info',{})
print(f\"{a.get('name','NA')} | release_date={a.get('release_date','NA')} | provider={a.get('provider','NA')}\")" 2>/dev/null || echo "NA")"

{
  echo "# 参考体系版本与校验值"
  echo
  echo "生成时间：$(date '+%F %T %Z')"
  echo
  echo "## 宿主参考"
  echo
  echo "- Assembly: \`${REFSEQ_ASM}\` (S288C R64-1-1；16 染色体 + 线粒体 \`NC_001224.1\`，不含天然 2μ)"
  echo "- RefSeq 注释：${ANNO_NAME}"
  echo "- 来源：\`${NCBI_GENOMES_BASE}/\`"
  echo
  echo "| 文件 | 字节 | md5（NCBI 官方） |"
  echo "| --- | --- | --- |"
  (cd "$REFS_HOST" && md5sum "${REFSEQ_FILES[@]}" 2>/dev/null) | while read -r m f; do
    printf '| `%s` | %s | `%s` |\n' "$f" "$(stat -c%s "$REFS_HOST/$f")" "$m"
  done
  echo
  echo "## 功能注释"
  echo
  echo "| 文件 | 来源 | 字节 | md5 |"
  echo "| --- | --- | --- | --- |"
  (cd "$ANNO_DIR" && md5sum *) | while read -r m f; do
    case "$f" in
      sgd.gaf.gz)          src="https://current.geneontology.org/annotations/sgd.gaf.gz" ;;
      goslim_yeast.obo)    src="https://current.geneontology.org/ontology/subsets/goslim_yeast.obo" ;;
      go_slim_mapping.tab) src="https://downloads.yeastgenome.org/curation/literature/go_slim_mapping.tab" ;;
      *)                   src="(local)" ;;
    esac
    printf '| `%s` | %s | %s | `%s` |\n' "$f" "$src" "$(stat -c%s "$ANNO_DIR/$f")" "$m"
  done
  echo
  echo "## 外源候选与诊断序列"
  echo
  echo "来源：NCBI E-utilities \`efetch\`（nuccore）"
  echo
  echo "| 文件 | 记录 | 用途 | 字节 | md5 |"
  echo "| --- | --- | --- | --- | --- |"
  (cd "$REFS_DIAG" && md5sum *) | while read -r m f; do
    case "$f" in
      ncrassa_cdt-1*) rec="XM_958708.2 (NCU00801)" ;;
      ncrassa_gh1-1*) rec="XM_011395456.1 (NCU00130)" ;;
      sfgfp*)         rec="PX636966.1" ;;
      egfp*)          rec="U55762.1 (pEGFP-N1)" ;;
      plasmid_prs426*)rec="U03451.1 (pRS426, URA3 marker)" ;;
      plasmid_2micron*)rec="NC_001398.1 (2μ circle)" ;;
      *)              rec="-" ;;
    esac
    case "$f" in
      *candidate*) use="3.4 序列裁决候选" ;;
      *egfp_candidates*) use="3.4 标签裁决候选（按翻译起点筛选 EGFP）" ;;
      *record*|*.gb) use="接头/边界核对" ;;
      plasmid*) use="5.6 read fate 诊断参考" ;;
      *) use="-" ;;
    esac
    printf '| `%s` | %s | %s | %s | `%s` |\n' "$f" "$rec" "$use" "$(stat -c%s "$REFS_DIAG/$f")" "$m"
  done
  echo
  echo "## 说明"
  echo
  echo "- 所有 URL 均写死版本，不使用 \`latest\` 类动态地址。"
  echo "- 脚本：\`scripts/01_download/02_fetch_reference.sh\`。"
  echo "- 完整 md5 清单：\`checksums/reference.md5\`。"
} > "$REFS_DIR/VERSIONS.md"
log "已生成 refs/VERSIONS.md"

find "$REFS_HOST" "$ANNO_DIR" "$REFS_DIAG" -type f ! -name "*.md5" -print0 | sort -z | xargs -0 md5sum \
  | sed "s|$PROJECT_ROOT/||" > "$CHECKSUM_DIR/reference.md5"
log "md5 清单写入 checksums/reference.md5（$(wc -l < "$CHECKSUM_DIR/reference.md5") 个文件）"
log "完成"
