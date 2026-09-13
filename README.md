# yeast-cellobiose-rnaseq

Reproducible re-analysis of **GEO GSE54825** — RNA-seq of engineered *Saccharomyces cerevisiae*
grown anaerobically on **cellobiose vs glucose** (n = 3 vs 3) — together with a **quantitative
comparison against the original 2014 study** (CLC Genomics Workbench + a private reference),
including an explicit account of what is *not* assessable and why.

Reference dataset: Lin Y, *et al.* *Biotechnol Biofuels.* 2014;7:126. PMID 25435910
([doi:10.1186/s13068-014-0126-6](https://doi.org/10.1186/s13068-014-0126-6)).

---

## What this repository demonstrates

| | |
| --- | --- |
| **Reproducibility engineering** | pinned environment, per-step scripts, one-command rerun (`run_all.sh`), md5 checksums for every downloaded artifact, provenance records for every non-public sequence |
| **Reference construction** | decoy-aware Salmon index on a frozen RefSeq assembly plus **manually reconstructed transgene sequences** (host-vector-insert junctions adjudicated from the reads themselves) |
| **Honest statistics** | unpaired primary model + paired sensitivity model, three threshold tiers plus an original-threshold-matched tier, explicit statement of which LFC (shrunken vs MLE) is used for thresholding |
| **Cross-study comparison** | per-gene log2FC correlation, DEG overlap (Jaccard / hypergeometric / threshold sweep), direction concordance, enrichment-layer and count-level reconciliation against the published datasets |
| **Failure accounting** | read-fate diagnostics, mitochondrial-transcript detectability, plasmid copy-number proxy, crosstalk checks, and a `not_assessable` list with reasons |

## Results at a glance

| | |
| --- | --- |
| Differential expression | **644 DEGs** (primary: BH-adjusted p < 0.05 and \|log2FC\| > 1); **632** at thresholds matched to the original study |
| Reproducibility of the original study | per-gene log2FC **Spearman 0.973**; **496 of the 519** original DEGs recovered (Jaccard 0.757, hypergeometric p ≈ 0); **100%** directional concordance |
| Functional layer | respiratory chain / oxidative phosphorylation / ATP synthesis up; amino-acid and **thiamine** biosynthesis down — the same processes reported in 2014 |
| Transgene work | cdt-1/gh1-1 reference sequences reconstructed and verified from the reads; the published codon-optimised `gh1-1a` is **not** the allele in these samples |
| Explicitly not assessable | mitochondrially encoded genes (poly-A selection), absolute transgene abundance, CLC-specific mapping statistics |

**Reports** (bilingual; all figure and table labels in English):

- English: [`docs/report/report_en.html`](docs/report/report_en.html) · [`report_en.pdf`](docs/report/report_en.pdf)
- 中文：[`docs/report/report_zh.html`](docs/report/report_zh.html) · [`report_zh.pdf`](docs/report/report_zh.pdf)

Stage acceptance records are in [`docs/acceptance/`](docs/acceptance/); parameter and scope decisions
(with their evidence) are in [`docs/decisions/`](docs/decisions/).

## Dataset

| | |
| --- | --- |
| GEO series | [GSE54825](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54825) / SRA `SRP037533` / BioProject `PRJNA237759` |
| Runs | `SRR1166442`–`SRR1166447` (ENA, public md5 used for verification) |
| Platform | Illumina Genome Analyzer II, single-end 50 bp, TruSeq barcoded multiplexed libraries |
| Design | anaerobic, exponential phase, single carbon source, 80 g/L; 2 conditions × 3 biological replicates |
| Strain | *S. cerevisiae* BY4742 + plasmid pRS426-BT (heterologous cellobiose-utilization pathway) |

## Reproducibility — how to rerun

```bash
# 1. environment (CLI tools + R/Bioconductor)
conda env create -f environment.yml
conda activate yeast-cellobiose-rnaseq

# 2. full pipeline (idempotent; skips completed steps)
bash run_all.sh

# 3. single step, e.g. download + checksum verification only
bash run_all.sh 01 03

# 4. render the bilingual reports (HTML + PDF)
bash scripts/07_comparison/20_render_report.sh
```

For the Chinese PDF, a CJK font must be available (Typst reads system fonts directly):

```bash
mkdir -p ~/.local/share/fonts
curl -fsSL -o ~/.local/share/fonts/NotoSansSC-Regular.otf \
  https://github.com/notofonts/noto-cjk/raw/main/Sans/SubsetOTF/SC/NotoSansSC-Regular.otf
fc-cache -f
```

Requirements: Linux (developed on WSL2 Ubuntu), ~15 GB free disk for raw data + indices,
internet access to NCBI/ENA.

Verification artifacts produced by the pipeline:

- `checksums/raw.md5`, `checksums/ena_published.md5`, `checksums/md5_verification.tsv`
- `refs/VERSIONS.md` (assembly accession, annotation release, download date, md5)
- `refs/custom/PROVENANCE.md` (origin of every non-public sequence)
- `env/tool_versions.txt` (tool versions + R package snapshot + `sessionInfo()`)
- `docs/decisions/` (parameter and scope decisions with reasons)

## Repository layout

```text
config/     samplesheet + shared path definitions
docs/       literature (read-only inputs + checksums), report sources, decision records
data/       raw/processed sequencing data (not tracked) and metadata provenance
refs/       host reference, reconstructed transgene reference, diagnostic sequences, indices
scripts/    numbered, runnable analysis steps (00_setup … 08_original_study)
results/    tables and figures produced by the pipeline
env/        environment lock and tool-version snapshots
logs/       run logs (not tracked)
checksums/  md5 records
```

The full execution plan, acceptance criteria and limitations are in [`PLAN.md`](PLAN.md).

## Data attribution

Sequencing data are public (SRA/ENA, BioProject `PRJNA237759`); the original study's supplementary
files are redistributed here under their original licenses for comparison purposes and are recorded
verbatim with md5 checksums (`docs/literature/CHECKSUMS.md`).

---

# 中文说明

**项目定位**：以**可复现流程**为主线的 RNA-seq 重分析项目，并给出与原研究（CLC + 私有参考体系）的
**量化一致性评估**。原研究：Lin et al. 2014, *Biotechnol Biofuels* 7:126。

**本仓库体现的能力**：可锁定环境与一键重跑、逐步骤可追溯脚本、下载物与序列的 md5/provenance 记录、
自定义参考体系（含**从测序数据裁决重建**的外源转基因序列）、分层统计口径与阈值敏感性分析、
以及对"不可评估项"的显式标注。

**重跑方式**：

```bash
conda env create -f environment.yml && conda activate yeast-cellobiose-rnaseq
bash run_all.sh          # 一键顺序重跑（幂等）
```

**执行依据**：全部范围、参数、阈值与验收标准见 [`PLAN.md`](PLAN.md)；参数决策记录见 `docs/decisions/`。
