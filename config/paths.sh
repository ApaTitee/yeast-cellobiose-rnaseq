#!/usr/bin/env bash
# 全流程统一路径变量。所有脚本通过 source 本文件获取路径，禁止在脚本内硬编码绝对路径。
# 用法： source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PROJECT_ROOT

# --- 配置与元数据 ---
export CONFIG_DIR="$PROJECT_ROOT/config"
export SAMPLESHEET="$CONFIG_DIR/samplesheet.tsv"
export DOCS_DIR="$PROJECT_ROOT/docs"
export LITERATURE_DIR="$DOCS_DIR/literature/original_study"
export DECISIONS_DIR="$DOCS_DIR/decisions"
export REPORT_DIR="$DOCS_DIR/report"

# --- 数据 ---
export DATA_DIR="$PROJECT_ROOT/data"
export RAW_FASTQ="$DATA_DIR/raw/fastq"
export PROC_FASTQ="$DATA_DIR/processed/fastq"
export META_DIR="$DATA_DIR/metadata"
export META_ENA="$META_DIR/ena"
export META_GEOSOFT="$META_DIR/geo_soft"
export META_ORIG="$META_DIR/original_study"
export META_ORIG_DERIVED="$META_ORIG/derived"
export META_ORIG_GEO="$META_ORIG/geo"

# --- 参考体系 ---
export REFS_DIR="$PROJECT_ROOT/refs"
export REFS_HOST="$REFS_DIR/host"
export REFS_CUSTOM="$REFS_DIR/custom"
export REFS_DIAG="$REFS_DIR/diagnostic"
export REFS_INDEX="$REFS_DIR/index"
export REFS_CHECKS="$REFS_DIR/checks"
export SALMON_INDEX="$REFS_INDEX/salmon"

# --- 结果 / 日志 / 校验 ---
export RESULTS_DIR="$PROJECT_ROOT/results"
export QC_DIR="$RESULTS_DIR/qc"
export FASTP_DIR="$QC_DIR/fastp"        # fastp 的 json/html 输出（读长等参数由此读取）
export QUANT_DIR="$RESULTS_DIR/quantification"
export DE_DIR="$RESULTS_DIR/differential_expression"
export ENRICH_DIR="$RESULTS_DIR/enrichment"
export COMPARE_DIR="$RESULTS_DIR/comparison"
export INTERP_DIR="$RESULTS_DIR/interpretation"
export FIG_DIR="$RESULTS_DIR/figures"
export LOGS_DIR="$PROJECT_ROOT/logs"
export CHECKSUM_DIR="$PROJECT_ROOT/checksums"
export ENV_DIR="$PROJECT_ROOT/env"

# --- 外部资源（写死版本，禁止 "latest" 类动态 URL）---
export ENA_PROJECT="PRJNA237759"
export REFSEQ_ASM="GCF_000146045.2"
export REFSEQ_ASM_DIR="GCF_000146045.2_R64"
export NCBI_GENOMES_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/${REFSEQ_ASM_DIR}"
export GEO_SERIES="GSE54825"
export GEO_SUPPL_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/${GEO_SERIES}/suppl"

# --- 计算资源 ---
export THREADS="${THREADS:-8}"
