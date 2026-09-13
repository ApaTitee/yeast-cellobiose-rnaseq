# 决策记录：R 环境锁定方式（conda 精确锁定，不使用 renv 项目激活）

**日期**：2026-09-13
**状态**：已采纳
**影响文件**：`environment.yml`、`env/conda-explicit.txt`、`env/r_packages.tsv`、`docs/decisions/r_environment.md`

## 背景

`PLAN.md` 原定 R 环境由 `renv.lock` 锁定（`renv::restore()` 重建）。实测在本机环境下该方案不可用：

1. `renv::init()` 会在项目根写入 `.Rprofile`（`source("renv/activate.R")`）并创建空的 `renv/library`。
2. R 启动时 renv 激活，`.libPaths()` 指向项目自有库，**遮蔽 conda 环境下的 `lib/R/library`**。
3. 实测结果：`library(DESeq2)` 直接报错、进程 `Execution halted`；`renv.lock` 因项目内尚无 R 代码而只记录 6 个包（含 renv 自身），既不可用也不完整。
4. 若要让 renv 真正可用，须在 renv 项目库内重新从源码编译安装全部 Bioconductor 依赖（DESeq2/clusterProfiler 及其 170+ 依赖），在 15 GB RAM 的单机 WSL 上耗时且易因系统库缺失失败。

## 决策

**R 环境与 CLI 工具统一由 conda 精确锁定，不使用 renv 项目激活。**

| 凭证 | 内容 |
| --- | --- |
| `environment.yml` | 可读的依赖声明（CLI + R 4.4.3 + Bioconductor 包） |
| `env/conda-explicit.txt` | **精确锁定**：`conda list --explicit` 输出的全部包 URL 与版本（376 行），可逐位重建 |
| `env/r_packages.tsv` | 177 个 R 包的名称/版本/库路径快照 |
| `env/tool_versions.txt` | CLI 工具版本 + `BiocManager::version()` + `sessionInfo()` |

`r-renv` 仍保留在 `environment.yml` 中（不影响上述流程），以备其他项目使用；本项目的 R 脚本一律直接运行，不激活 renv。

## 重建方式

```bash
conda env create -f environment.yml          # 从声明重建
# 或按精确锁定逐位重建（最严格）：
conda create -n yeast-cellobiose-rnaseq --file env/conda-explicit.txt
```

## 对验收标准的影响

`PLAN.md` §1.3 的"环境可锁定"要求仍满足，且锁定强度不低于 renv：
`environment.yml`（声明）+ `env/conda-explicit.txt`（精确）+ `env/r_packages.tsv`（R 包快照）
+ 每次运行的 `sessionInfo()`。已同步更新 `PLAN.md` 的 §1.3、§5.0、§6 与附录 B。

## 复现验证

```bash
conda run -n yeast-cellobiose-rnaseq Rscript -e \
  'suppressMessages(library(DESeq2)); suppressMessages(library(clusterProfiler));
   cat(as.character(packageVersion("DESeq2")), as.character(BiocManager::version()), "\n")'
# -> 1.46.0 3.20   （2026-09-13 实测通过）
```
