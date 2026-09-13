#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# record_tool_versions.sh
# 目的：记录 CLI 工具版本与 R 包版本快照，作为可复现凭证（附录 B）。
# 输出：env/tool_versions.txt
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config/paths.sh"

OUT="$ENV_DIR/tool_versions.txt"
{
  echo "# Tool versions — generated $(date '+%F %T %Z')"
  echo "# host: $(uname -srm) | conda env: ${CONDA_DEFAULT_ENV:-none}"
  echo
  echo "## OS / scheduler"
  echo "kernel: $(uname -r)"
  nproc --all | sed 's/^/cores: /'
  free -g | awk 'NR==2{print "mem_total_GB: "$2}'
  echo
  echo "## CLI tools"
  for c in "fastqc --version" "fastp --version" "multiqc --version" "salmon --version" \
           "bowtie2 --version" "samtools --version" "seqtk" "gffread --version" \
           "datasets --version" "seqtk 2>/dev/null" "quarto --version" "typst --version" \
           "python3 --version" "md5sum --version"; do
    name="${c%% *}"
    if command -v "$name" >/dev/null 2>&1; then
      printf '%-10s %s\n' "$name" "$($c 2>&1 | head -1)"
    else
      printf '%-10s %s\n' "$name" "NOT INSTALLED"
    fi
  done
  echo
  echo "## R / Bioconductor"
  if command -v Rscript >/dev/null 2>&1; then
    Rscript -e 'cat("R: ", R.version.string, "\n", sep="")' 2>/dev/null
    Rscript -e 'if (requireNamespace("BiocManager", quietly=TRUE)) cat("Bioconductor: ", as.character(BiocManager::version()), "\n", sep="")' 2>/dev/null
    Rscript -e 'cat("renv: ", as.character(utils::packageVersion("renv")), "\n", sep="")' 2>/dev/null
    echo
    echo "## R packages"
    Rscript -e 'ip <- as.data.frame(installed.packages()[, c("Package","Version")]); ip <- ip[order(ip$Package), ]; write.table(ip, row.names=FALSE, sep="\t", quote=FALSE)' 2>/dev/null
  else
    echo "R: NOT INSTALLED"
  fi
  echo
  echo "## sessionInfo()"
  Rscript -e 'print(sessionInfo())' 2>/dev/null || echo "N/A"
} > "$OUT"
echo "已写入 $OUT"
