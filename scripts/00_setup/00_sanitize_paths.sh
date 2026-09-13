#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 00_sanitize_paths.sh —— 提交前把机器相关信息（绝对路径、用户名、主目录）替换为占位符
# 目的：仓库会被公开（GitHub），任何含本机绝对路径/用户名的文本都不应入库。
# 处理范围：受限的文本文件类型（跳过二进制、数据文件、报告渲染产物）
# 幂等：可重复运行；只替换已知的机器前缀。
# 用法： bash scripts/00_setup/00_sanitize_paths.sh [--check]
#   --check  仅报告存在机器相关字符串的文件，不做修改（返回码 1 表示存在）
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../config/paths.sh"

CHECK=0
[ "${1:-}" = "--check" ] && CHECK=1

# 需要替换的机器相关前缀（按长度从长到短替换，避免部分覆盖）
PROJ="$PROJECT_ROOT"
HOME_DIR="$(dirname "$(dirname "$PROJ")")"     # 一般为 /home/<user> 或 /Users/<user>
USER_NAME="$(basename "$HOME_DIR")"

declare -a PAIRS=(
  "$PROJ|<PROJECT_ROOT>"
  "$HOME_DIR|<HOME>"
)

# 只处理文本类文件（按扩展名白名单 + 无扩展名的脚本）
EXTS=("md" "qmd" "R" "sh" "py" "tsv" "csv" "json" "log" "txt" "yml" "yaml" "bib" "toml" "html" "tex")
INCLUDE_ARGS=()
for e in "${EXTS[@]}"; do INCLUDE_ARGS+=(--include="*.$e"); done

# 只检查**纳入版本控制**的文件（未跟踪/被忽略的产物如 logs/ 不在提交范围内）
if git -C "$PROJECT_ROOT" rev-parse --git-dir >/dev/null 2>&1; then
  mapfile -t FILES < <(cd "$PROJECT_ROOT" && git ls-files -z | xargs -0 grep -lE "$PROJ|$HOME_DIR" 2>/dev/null || true)
else
  mapfile -t FILES < <(grep -rIl -E "$PROJ|$HOME_DIR" "$PROJECT_ROOT" \
    --exclude-dir=.git --exclude-dir=.quarto --exclude-dir=figures "${INCLUDE_ARGS[@]}" 2>/dev/null || true)
fi

if [ "${#FILES[@]}" -eq 0 ]; then
  echo "  [ok] 未发现机器相关路径"
  exit 0
fi

echo "  发现 ${#FILES[@]} 个文件含机器相关路径："
for f in "${FILES[@]}"; do echo "    ${f#"$PROJECT_ROOT"/}"; done

if [ "$CHECK" -eq 1 ]; then
  echo "  [--check] 存在机器相关路径，未修改"
  exit 1
fi

for f in "${FILES[@]}"; do
  for pair in "${PAIRS[@]}"; do
    from="${pair%%|*}"; to="${pair##*|}"
    sed -i "s|$from|$to|g" "$f"
  done
  # 兜底：常见用户名残留（如路径片段 <user>/projects）
  sed -i "s|/${USER_NAME}/|<HOME>/|g" "$f"
done

echo "  [done] 已替换为 <PROJECT_ROOT> / <HOME>（${#FILES[@]} 个文件）"
echo "  提示：这些文件重新生成后（如 salmon 日志）可能再次带出绝对路径，提交前请重跑本脚本。"
