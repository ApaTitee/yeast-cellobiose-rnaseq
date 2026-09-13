#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_all.sh — 一键顺序重跑全流程（幂等）
# 用法：
#   bash run_all.sh              # 跑完所有已实现的步骤，缺失步骤给出警告
#   bash run_all.sh --strict     # 任一缺失步骤即失败
#   bash run_all.sh 02 03        # 只跑指定编号前缀的步骤
# 每个步骤的输出日志见 logs/<编号>_*.log
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config/paths.sh"

STRICT=0
SELECT=()
for a in "$@"; do
  case "$a" in
    --strict) STRICT=1 ;;
    *) SELECT+=("$a") ;;
  esac
done

# 有序步骤表：编号 -> 脚本
STEP_SCRIPTS=(
  "00|scripts/00_setup/00_init_project.sh"
  "01|scripts/01_download/01_fetch_ena_fastq.sh"
  "02|scripts/01_download/02_fetch_reference.sh"
  "03|scripts/01_download/03_fetch_metadata.sh"
  "04|scripts/02_qc/04_fastqc_fastp.sh"
  "05|scripts/02_qc/05_postmap_diagnostics.sh"
  "06|scripts/08_original_study/06_parse_original_study.sh"
  "07|scripts/08_original_study/07_transgene_adjudication.sh"
  "08|scripts/03_quantification/08_build_reference.sh"
  "09|scripts/03_quantification/09_salmon_index.sh"
  "10|scripts/03_quantification/10_salmon_quant.sh"
  "11|scripts/03_quantification/11_tximport_gene_matrix.R"
  "12|scripts/03_quantification/12_param_sensitivity.sh"
  "13|scripts/04_differential_expression/13_sample_level.R"
  "14|scripts/04_differential_expression/14_deseq2_de.R"
  "15|scripts/04_differential_expression/15_model_diagnostics.R"
  "16|scripts/05_enrichment/16_enrichment.R"
  "17|scripts/06_visualization/17_figures.R"
  "18|scripts/07_comparison/18_consistency_with_original.R"
  "19|scripts/07_comparison/19_hypothesis_scorecard.R"
  "20|scripts/07_comparison/20_render_report.sh"
  "21|scripts/00_setup/00_sanitize_paths.sh"
)

run_one() {
  local num="$1" rel="$2" script="$PROJECT_ROOT/$rel"
  if [ ! -f "$script" ]; then
    printf '[skip] %s 尚未实现：%s\n' "$num" "$rel"
    [ "$STRICT" -eq 1 ] && return 1 || return 0
  fi
  printf '\n===== [%s] %s =====\n' "$num" "$rel"
  case "$script" in
    *.R) Rscript "$script" 2>&1 | tee -a "$LOGS_DIR/run_all.log" ;;
    *)   bash "$script"     2>&1 | tee -a "$LOGS_DIR/run_all.log" ;;
  esac
}

mkdir -p "$LOGS_DIR"
for entry in "${STEP_SCRIPTS[@]}"; do
  num="${entry%%|*}"; rel="${entry##*|}"
  if [ "${#SELECT[@]}" -gt 0 ]; then
    keep=0; for s in "${SELECT[@]}"; do [ "$num" = "$s" ] && keep=1; done
    [ "$keep" -eq 1 ] || continue
  fi
  run_one "$num" "$rel"
done

printf '\n全部选定步骤执行完毕。产物见 results/ ，日志见 logs/run_all.log\n'
