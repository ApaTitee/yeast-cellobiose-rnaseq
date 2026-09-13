# Project 01 — 工程化酵母 RNA-seq 重分析：执行计划

**定位**：以**可复现流程**为主线的生信分析项目。全部工作在同一 WSL2 Ubuntu 环境内完成，本文件是该环境内的唯一执行依据。

**数据集**：GEO `GSE54825` — 工程化 *Saccharomyces cerevisiae*，Glucose vs Cellobiose，n = 3 vs 3。

---

## 0. 项目状态与执行前提

### 0.1 工作目录

- **项目根目录 = `<PROJECT_ROOT>`**（WSL 原生 ext4，可用 934 GB；**不得**迁到 `/mnt/c`）。
- 该目录当前不是 git 仓库 → P0 首件事 `git init`，并配置 GitHub 远端（见 5.0）。
- 资源：16 核 / 15 GB RAM（可用 ~11 GB）/ 无 GPU。

### 0.2 环境现状

| 项 | 现状 | P0 处理 |
| --- | --- | --- |
| conda | `26.1.1`（`<HOME>/miniconda3`） | 直接 `conda env create`；**不安装 mamba** |
| mamba / R / samtools / seqtk / salmon | 均未安装 | 全部由 `environment.yml` 提供 |
| Snakemake | 未安装 | P5（可选） |
| 网络 | 可用（NCBI / ENA / EBI 直连已验证） | 下载脚本一律用 `https://`，不用 `ftp://` |

### 0.3 本地已有材料（只读输入）

| 文件 | 内容 | md5 |
| --- | --- | --- |
| `s13068-014-0126-6.pdf` | Lin et al. 2014 正文 | `82a64a8127fecb9df6f3a4777001fe75` |
| `13068_2014_126_MOESM1_ESM.xlsx` | **Dataset S1** 全基因组 normalized RPKM（6351 条，含 2 条外源） | `98a937e749414d5dbd0145914569b149` |
| `13068_2014_126_MOESM2_ESM.xlsx` | **Dataset S2** 519 DEG + 19 TF + 7 GRN | `dd09b0b5c6a0b1d8e4c6ee6b04c24a82` |
| `13068_2014_126_MOESM3_ESM.xlsx` | **Dataset S3** FunSpec GO 结果（UP 244 / DOWN 256） | `98ec1472114123a85096cb2dcd87c613` |
| `13068_2014_126_MOESM4_ESM.pdf` | Supporting information：Fig S1–S8 + Table S1/S2 | `2ce0d449ba7caf012d065cfcf02775e4` |
| `science.1192838.pdf` | Galazka et al. 2010 *Science* 正文 | `b072b556036459ccbce66a5ce6cb0646` |
| `galazka-som.pdf` | Galazka et al. 2010 SOM（含全部克隆引物） | `5d3c1197ff3a118632376c84d5cf3c6d` |

P0 将这 7 个文件 `git mv` 至 `docs/literature/original_study/`，md5 写入 `docs/literature/CHECKSUMS.md`。

### 0.4 已确定的执行参数（固化进脚本，并在 `docs/decisions/` 留档）

| 项 | 取值 |
| --- | --- |
| 远端仓库 | GitHub 私有/公开待定，仓库名 `<GITHUB_REPO>`（账号信息由项目负责人提供） |
| 外源序列路线 | **路线 A（证据驱动重建）**；路线 B 仅作背景备选，不作为关键路径 |
| 报告 | **中英双语两版**（同一数据源自动生成），图表标签与图注**统一英文**，排版参照正规期刊 |
| DEG 主阈值 | `padj < 0.05` 且 `|log2FC| > 1`（阈值判定一律用**未收缩 LFC**） |
| 放宽 / 收紧阈值 | `padj < 0.1` / `padj < 0.01` 且 `|log2FC| > 1.5` |
| 原文对齐阈值 | `|FC| ≥ 2` 且 `FDR ≤ 0.001`（`results(alpha = 0.001)`，未收缩 LFC） |
| 富集 universe | 在 ≥3 个样本中 `count ≥ 10` 的基因 |
| 定量路线 B（STAR） | **关闭** |
| 参数敏感性 | `--kmerLen` 31 vs 25；修剪 vs 未修剪；decoy 是否含质粒骨架 |

---

## 1. 研究问题与目标

### 1.1 研究问题

工程化 *S. cerevisiae* 在 **cellobiose** 与 **glucose** 两种单一碳源、厌氧指数期条件下，全局转录差异是什么？差异基因涉及哪些生物学过程与代谢通路？其中哪些线索指向可工程化的候选基因与通路？

### 1.2 目标

1. 从原始测序数据出发，建立一条可重复执行的完整分析流程（下载 → QC → 定量 → 表达分析 → 差异表达 → 功能分析 → 解读）。
2. 鉴定差异表达基因，并分析其生物学过程、分子功能与代谢通路。
3. 解析工程化酵母在 cellobiose 条件下的转录响应。
4. 结合代谢工程背景整理潜在的工程化候选基因与通路。
5. 与本研究的原始发表结果进行**量化一致性评估**。

### 1.3 硬性验收标准

- 在干净环境中，按 `README.md` 指令可从头重跑，重新生成全部表与图。
- 环境可锁定：`environment.yml` + `env/conda-explicit.txt`（精确锁定，含 R/Bioconductor）+ `env/r_packages.tsv` + 基因组/注释/软件版本全部记录。
- 每个中间产物可追溯（脚本 + 参数 + 输入文件），关键文件带 `md5`。
- 报告中出现的每一个数字、阈值、版本，都能在产物文件中找到来源。
- 与原研究给出一致性的**量化结论**，并显式标注哪些结论**不可评估**及其原因。
- 结论与限制同等清晰。

### 1.4 非目标

- 湿实验验证、菌株构建、发酵实验。
- 外源通路的工程改造设计。
- 变异检测与位点级分析（唯一例外：3.4 的外源序列裁决）。
- 从零开始的软件方法学开发。
- 长期 portfolio 阶梯（Transcriptomics → Protein engineering → Structure-guided design → ML-assisted engineering）仅作路线图，不属本项目验收范围。

---

## 2. 数据资源

### 2.1 数据集事实

| 项目 | 内容 |
| --- | --- |
| 数据集 | GEO **GSE54825**（提交 2014-02-10，Public 2015-04-10） |
| 平台 | Illumina **Genome Analyzer II**（GPL9377），**单端 50 bp**，TruSeq barcoded 多重文库 |
| 实验设计 | 厌氧、指数期、单一碳源；2 组 × 3 生物学重复；80 g/L 起始糖，收获于残余糖 ≈ 50 g/L |
| 菌株 | **BY4742**（MATα, his3Δ1, leu2Δ0, lys2Δ0, ura3Δ0）+ 质粒 **pRS426-BT** |
| 原始数据 | SRA **SRP037533** / BioProject **PRJNA237759**；run `SRR1166442–SRR1166447`；每样本 34.5–45.5 M reads，合计 **11.93 Gb** |
| 原文参考体系 | S288C **R64-1-1**（16 染色体 + 线粒体）+ 手工注释的 *N. crassa* **gh1-1** 与 **eGFP 标签化 cdt-1**，合计 12.17 Mb。该参考文件未公开 |
| 原文分析软件 | **CLC Genomics Workbench 6.5**（比对与定量均在其内完成，未使用外部比对器） |
| 原文注释来源 | **SGD gene association file** |
| 原文定量与统计 | RPKM（"By totals" 再归一化）；unpaired two-group comparison + **Baggerley's test**；阈值 \|FC\| ≥ 2 且 FDR p ≤ 0.001 |
| 原文富集分析 | GO biological process 富集，**FunSpec（P 0.01 + Bonferroni）** |
| 原文映射统计 | 约 3,800 万 reads/文库（≈156× 覆盖）；**92.3% 导入 / 83.4% 可比对 / 76.3% 唯一比对** |
| 原文主要结果 | **519 / 6351 基因（8.2%）显著差异**；cellobiose 上调富集于**线粒体相关过程**（ATP 合成、电子与质子传递、TCA 循环）；下调富集于**氨基酸合成**（Met/Cys/Arg/His）与**硫胺素（VB1）合成** |
| 原文延伸结论 | 19 个转录因子被扰动；**SUT1 过表达** 与 **HAP4 缺失** 稳定提升纤维二糖发酵；外源 cdt-1/gh1-1 mRNA 在 cellobiose 下反而升高（尽管由 P_PGK1 驱动） |
| 论文 | Lin Y, Chomvong K, Acosta-Sampson L, Estrela R, Galazka JM, Kim SR, Jin YS, Cate JH. *Biotechnol Biofuels.* 2014;7:126. PMID 25435910；PMC PMC4243952；doi 10.1186/s13068-014-0126-6 |
| 关键序列标识 | cdt-1 = *N. crassa* **NCU00801**（XM_958708.2，CDS 1740 bp / 580 aa）；gh1-1 = **NCU00130**（XM_011395456.1，CDS 1431 bp / 477 aa） |

### 2.2 已核验的补充事实与证据

**(a) 样本表已双源验证**：ENA filereport（PRJNA237759）的 `sample_title` = `JCYL001B…JCYL003D`，与 GEO SOFT（GSM1324496–GSM1324501）的 `carbon source: Glucose / Cellobiose` 完全对应，与 2.3 表的 run↔GSM↔条件一致；全部 `SINGLE` / `Illumina Genome Analyzer II`，每条 read 恰为 50 bp。

**(b) 原文事实逐条命中**（正文核对）：92.3%/83.4%/76.3%、519/6351（8.2%）、\|FC\| ≥ 2.0 且 FDR ≤ 0.001、FunSpec P 0.01 + Bonferroni、CLC GW 6.5、R64-1-1 + 线粒体、SGD GAF、"12.17 Mb"、SUT1 过表达 / HAP4 缺失、P_PGK1 驱动的 cdt-1/gh1-1 在 cellobiose 下升高而内源 PGK1 无显著升高。

**(c) "6351" 的构成**：Dataset S1 数据行为 6351 条 = **6349 条宿主条目 + 2 条外源条目**（`cdt-1EGFP`、`gh1-1`）。故"519/6351"中的 519 **包含 1 条外源条目**（`gh1-1`，FC = 2.16，FDR p = 0）。报告统一表述为"518 个宿主基因 + 1 个外源条目"。

**(d) 原文 Dataset S2/S3 的符号约定**：

- `Fold Change` 列 = **带符号倍数**：上调为 `+(C8/G8)`，下调为 `-(G8/C8)`（如 THI4 记 −64.17）。
- `log` 列 = **`log2(C8/G8)`**（THI4 = −6.0038）。
- Dataset S3：`244 UP genes` + `256 DOWN genes`；UP 侧未映射 15 个基因（`ARG5,6`、`MFALPHA292`、Ty 元件等），DOWN 侧未映射 3 个（`COX26`、`MOS1`、`SHH3`）。
- **对账已闭合**（2026-09-13，由 `scripts/08_original_study/06_parse_original_study.py` 复算，见 `data/metadata/original_study/derived/reconciliation.tsv`）：244 + 256 + 15 + 3 = **518**，与 519 相差 **1** 条，即外源条目 `gh1-1`（无 GO 注释故未进入 FunSpec 分析）。此结论仍须在报告 5.12 中以数字形式复述。

**(e) 原文锚定值（Dataset S1/S2 的逐样本 RPKM）**：

| 条目 | C8（001D/002D/003D） | C8 均值 | G8（001B/002B/003B） | G8 均值 | FC(C8/G8) | log2 |
| --- | --- | --- | --- | --- | --- | --- |
| `cdt-1EGFP` | 1604.81 / 1839.95 / 1477.35 | 1640.70 | 970.54 / 950.36 / 709.29 | 876.73 | 1.871 | 0.904 |
| `gh1-1` | 1995.68 / 2186.65 / 1856.90 | 2013.08 | 997.33 / 958.49 / 839.41 | 931.74 | 2.160 | 1.111 |
| `PGK1` | 4334.94 / 4675.23 / 4335.91 | 4448.69 | 6123.11 / 5712.44 / 6081.18 | 5972.24 | 0.745 | −0.425 |
| `CYC1` | 17.66 / 15.63 / 28.22 | 20.50 | 7.78 / 7.83 / 8.66 | 8.09 | 2.534 | 1.341 |
| `URA3` | 646.84 / 581.86 / 687.18 | 638.63 | 419.38 / 433.69 / 562.01 | 471.69 | 1.354 | 0.437 |
| `HAP4` | 138.81 / 139.71 / 152.82 | 143.78 | 22.65 / 25.34 / 22.45 | 23.48 | 6.124 | 2.614 |
| `SUT1` | 164.55 / 164.21 / 154.86 | 161.21 | 46.08 / 47.88 / 48.34 | 47.43 | 3.399 | 1.765 |

由此得到两条**预注册判据**（写入 5.11/5.12，在查看本研究结果之前确定）：

- **URA3 拷贝数代理**：质粒来源 URA3 在 cellobiose 中为 glucose 的 **1.354 倍**。若外源基因升高完全来自质粒拷贝数，则 cdt-1/gh1-1 的 FC 应 ≈ 1.35；原文实测为 1.87 / 2.16 → 拷贝数最多解释一部分（URA3 归一化后：cdt-1 = 1.382，gh1-1 = 1.596）。
- **同启动子反证**：外源基因与内源 `PGK1` 同为 P_PGK1 驱动但方向相反（外源升、PGK1 降 0.745）→ 启动子活性不能解释外源升高。`CYC1`（终止子来源基因）本身在 cellobiose 中升高 2.53 倍（属 H1 的 OXPHOS 基因），**因此 CYC1 的串扰无法从生物学变化中分离**，须在 Limitations 声明。

**(f) 参考总量核算**：宿主 R64-1-1 = 12,071,326 bp（BK006934–BK006949）+ 线粒体 NC_001224.1 = 85,779 bp → 12,157,105 bp；加两条外源转录本 ≈ 3.9 kb → ≈ 12.161 Mb，与原文"12.17 Mb"在 ±10 kb 内一致；原文参考未纳入质粒骨架。

### 2.3 样本表与原始数据

`config/samplesheet.tsv`：

| sample_id | gsm | run | condition | replicate | layout | platform | ENA read_count | ENA md5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JCYL001B | GSM1324496 | SRR1166442 | glucose | 1 | SE 50bp | GAII | 40,585,197 | `68e674b8a20e88f962492481671f619e` |
| JCYL002B | GSM1324497 | SRR1166443 | glucose | 2 | SE 50bp | GAII | 44,937,781 | `dcc91f00361255629d2ec0cc8ef6083d` |
| JCYL003B | GSM1324498 | SRR1166444 | glucose | 3 | SE 50bp | GAII | 34,778,661 | `f58e37e44319f4f45ba871d760c8ab24` |
| JCYL001D | GSM1324499 | SRR1166445 | cellobiose | 1 | SE 50bp | GAII | 45,543,917 | `a201447b2ddad46b08f85e1ad7311e34` |
| JCYL002D | GSM1324500 | SRR1166446 | cellobiose | 2 | SE 50bp | GAII | 38,223,389 | `028b8c4edc3bbdbd8db69ca9f2edf758` |
| JCYL003D | GSM1324501 | SRR1166447 | cellobiose | 3 | SE 50bp | GAII | 34,542,268 | `cbb987dbfaa3fda2da823ff583f9d675` |

- `condition` 是唯一比较因子，`glucose` 设为 reference level；`replicate` 仅用于敏感性分析（4.2）。
- ENA 公布的 md5 为校验基准（自算 md5 只用于比对，不作为"通过"依据）。
- 下载用 `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/00X/SRRxxxxxxx/…`（由 ENA filereport 的 `fastq_ftp` 字段换算），不使用 SRA Toolkit。预计 FASTQ.gz 合计 ≈ 6–9 GB。

### 2.4 原文结果资源

| 资源 | 本地状态 | 落库位置 | 用途 |
| --- | --- | --- | --- |
| Dataset S1（RPKM 全表） | 已有 | `docs/literature/original_study/` → `data/metadata/original_study/derived/dataset_s1_rpkm.tsv` | 定量层对照 + 外源锚定值 |
| Dataset S2（519 DEG / 19 TF / 7 GRN） | 已有 | 同上 → `dataset_s2_deg.tsv` 等 | DEG 重叠、方向一致率、TF 集合 |
| Dataset S3（FunSpec GO） | 已有 | 同上 → `dataset_s3_goslim_funspec.tsv` | 富集层对照（含 244/256 与 "not found" 清单） |
| Supporting information（Fig S1–S8、Table S1/S2） | 已有 | 同上 + `table_s1_plasmids.tsv`、`table_s2_primers.tsv` | 菌株/质粒结构、引物、Fig S8 序列（3.4 的阴性对照） |
| Galazka 2010 正文 + SOM | 已有 | 同上 + `galazka2010_cloning_primers.tsv` | 外源序列重建依据 |
| GEO 补充文件 | 缺 | `data/metadata/original_study/geo/` | 逐基因 FC/FDR 对照 |

GEO 补充文件直链（P0 下载并记 md5）：
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54825/suppl/GSE54825_Cellobiose_versus_Glucose.xlsx.gz`

---

## 3. 参考体系

### 3.1 宿主参考（冻结）

- **RefSeq `GCF_000146045.2`（R64）**：17 条序列 = 16 染色体（BK006934–BK006949 / NC_001133–NC_001148）+ 线粒体（**NC_001224.1，85,779 bp**）。该组装**不含天然 2μ 质粒**，此事实为 3.5 与 5.6 的前提。
- 冻结文件：`GCF_000146045.2_R64_genomic.fna.gz`、`_rna_from_genomic.fna.gz`、`_cds_from_genomic.fna.gz`、`_genomic.gff.gz`。转录本骨架用 RefSeq，功能注释用 SGD GAF + GO slim。
- **冻结规范**：
  1. `refs/VERSIONS.md` 记录 **assembly accession + RefSeq annotation release + 下载日期 + 每个文件 md5**；
  2. 脚本内写死文件名与 md5，**禁止** "latest" 类动态 URL；
  3. 跨步骤 join key 统一为 **systematic ORF name**，并输出 **ID 映射损失表**（5.2）；
  4. `refs/host/` 下的 GFF/GTF 与 ID 映射表入库（体积小），基因组 FASTA 与索引不入库。

### 3.2 偏差声明（写入报告的固定文本）

> 本研究分析的是 **BY4742** 衍生工程菌（MATα, his3Δ1, leu2Δ0, lys2Δ0, ura3Δ0），而参考基因组为 **S288C R64-1-1**（`GCF_000146045.2`）。两者高度接近但非完全相同：存在少量株系间 SNP/indel，且 4 个营养缺陷标记位点在 BY4742 中为缺失等位，其序列与参考不符。因此：(i) 基因水平表达与差异分析受影响极小；(ii) `URA3`/`LEU2`/`HIS3`/`LYS2` 四个位点的计数不反映宿主基因表达，解读时予以排除；(iii) 本研究不进行位点级或变异层面分析（唯一例外是 3.4 的外源序列裁决，且不进主结论）。

**URA3 代理指标的边界条件**：BY4742 的 `ura3Δ0` 为完全缺失等位（工作假设：自 −223 至 +880，相对 ATG；P0 以 SGD 等位页核实并留档）。代理指标成立的**充分操作条件**是"**宿主侧不产生 URA3 转录本**"，即缺失区间覆盖 URA3 的整个转录单元（5'UTR–CDS–3'UTR）。若核实结果为仅删除 ORF 而保留完整 UTR，则在 Limitations 追加"残留 UTR 可能贡献极少量映到 URA3 的 reads"，并将代理降级为半定量。参考中保留野生型 `YEL021W` 序列（不额外加入质粒 URA3 条目，见 3.5）。

### 3.3 关键位点核查清单

逐项检查并填入 `refs/checks/locus_check.tsv`：

| 检查项 | 预期 | 说明 |
| --- | --- | --- |
| 线粒体基因组 | 参考含 `NC_001224.1` 且注释含 mt 编码基因 | H1 的核/线粒体限定依赖此项 |
| 天然 2μ 质粒 | 参考中不存在（17 条序列） | 决定"2μ 来源 reads 只能进 read-fate 桶"（5.6） |
| 营养缺陷标记 | `URA3` YEL021W、`LEU2` YCL018W、`HIS3` YOR202W、`LYS2` YBR115C | 系统性名以 SGD 核对后写入 `VERSIONS.md` |
| `URA3` 转录单元边界 | 记录 YEL021W 的 5'UTR/CDS/3'UTR 坐标 | 配合 3.2 的边界条件 |
| `ura3Δ0` 缺失区间 | −223 … +880（待 SGD 核实） | 操作判据见 3.2 |
| **URA3 的特殊性** | pRS426 以 URA3 为选择标记，宿主该位点缺失 | 质粒来源 URA3 reads 唯一落到参考 URA3；计数作**质粒拷贝数代理**；预注册值见 2.2(e) |
| Ty 元件 / rDNA（RDN1） | 多重映射热点 | 在 5.6 统计多映射比例，解释比对率缺口 |
| 标准名/系统名映射 | 与 SGD 一致 | 输出 6351 → RefSeq 的映射损失表（5.2） |
| `PGK1` / `CYC1` 基因体与侧翼 | 记录坐标与 3' 端 | 3.5 串扰检查义务的对象 |

### 3.4 外源序列重建（P1 前置任务）

#### 3.4.1 已确证的事实

1. **RNA-seq 样本的质粒是 `pRS426-BT`**（Methods 原文："Anaerobic cultures of strain BY4742 expressing the cellobiose-utilizing pathway from plasmid pRS426-BT…"）。Table S1：`pRS426-BT` = *N. crassa* **gh1-1 与 cdt-1**，P_PGK1 驱动 + **T_CYC1** 终止，骨架 pRS426（2μ, URA3），"This study"。
2. **外源基因来自 cDNA，不是密码子优化版本**。Galazka 2010 SOM 原文："**All *N. crassa* genes were amplified by PCR from cDNA synthesized from mRNA isolated from *N. crassa* (FGSC 2489)**"。→ 参考序列应为**去内含子的天然 CDS**。
3. **必须区分两个版本的 gh1-1**：
   - 样本中的 `gh1-1`（pRS426-BT）= 天然 cDNA CDS（1431 bp）+ C 端 6×His；
   - 补充材料 **Figure S8 的 `gh1-1a`（1449 bp）**是本研究后期优化实验合成、用于 `pRS316-BaT` / `pRS315` 的版本，**不在 RNA-seq 样本中**。实测 `gh1-1a` 与天然 NCU00130 CDS 的**核酸同一性仅 74.5%**（蛋白 100% 一致 + His6），**不得作为定量参考**，仅作 3.4.3 判别实验的阴性对照。
4. **两个基因的天然 CDS 与 2010 年克隆引物逐位吻合**：gh1-1 正向引物 3' 端命中天然 CDS 第 5 bp 起，反向引物命中 3' 末端；cdt-1 的 GFP 版反向引物反补后为天然 CDS 倒数第 20 bp 起的 `GTCCGAAGCTATCGTTGCT`。→ 骨架序列可从公共库精确取得。
5. **标签/接头信息可由 Galazka SOM 引物直接导出**：
   - `cdt-1`（NCU00801）：存在 **Myc 版**（BamHI–EcoRI，C 端 Myc）与 **GFP 融合版**（BamHI–EcoRI 插入已含 `ClaI–GSGS–sfGFP–SalI` 的 P_PGK1-pRS426 骨架）；`sfGFP` 为 **superfolder GFP**，N 端带 **Gly-Ser-Gly-Ser** 接头，反向引物含终止密码子。
   - `gh1-1`（NCU00130）：Kozak 优化 + **C 端 6×His**（引物 `…TTA ATG×6 GTCCTTCTTGATCAAAGAGTCA AAG`），蛋白末端为 `…KPLFDSLIKKD-HHHHHH-*`（与天然 CDS 末端 `…AAGAAGGAC TAA` 一致）。
   - 终止子为 **CYC1 终止子**（XhoI–KpnI）；启动子为酵母基因组 **PGK1 启动子**（可由宿主参考直接取得）。
6. **由数据裁决的部分（已于 2026-09-13 完成，见 `refs/custom/PROVENANCE.md`）**：
   - pRS426-BT 中 cdt-1 为天然 CDS + 标签融合体；标签经 reads 裁决为 **superfolder GFP（sfGFP）**，其密码子接近人源化 eGFP
     （reads 共识与 pEGFP-N1 的 EGFP CDS 仅 7 nt 差异；翻译后为 eGFP 骨架 + sfGFP 全部 6 个 superfolder 替换）；
     因此参考中使用 **reads 共识序列**，不使用任何 GenBank 的 sfGFP 记录（密码子不同，会丢失比对）。
   - `gh1-1` 为 **天然 CDS（0 处差异）+ C 端 6×His**；Figure S8 的 `gh1-1a` 在样本中不存在（阴性对照，仅 2 reads）。
   - **仍待完成**：cdt-1–标签接头的精确碱基、两条转录本的 5'/3' 边界（由软剪切 reads 与覆盖断点界定）。

#### 3.4.2 重建方案

**路线 A（主路线）：证据驱动的确定性重建**

| 元件 | 来源 | 确定性 |
| --- | --- | --- |
| cdt-1 CDS（1740 bp） | NCBI RefSeq `XM_958708.2` 的去内含子 CDS | 引物已逐位验证，高 |
| gh1-1 CDS（1431 bp） | NCBI RefSeq `XM_011395456.1` 的去内含子 CDS | 引物已逐位验证，高 |
| 标签 | 由 Galazka SOM 引物翻译（6×His；Myc 或 GSGS-sfGFP） | 中高（变体待裁决） |
| P_PGK1 / T_CYC1 | 冻结的宿主参考基因组（记录坐标） | 高 |
| 接头与 5'/3' 边界 | 由 reads 裁决（3.4.3） | 需实测 |

**路线 B（背景备选，不进关键路径）**：向通讯作者（J. H. D. Cate / Y.-S. Jin）索取 pRS426-BT 图谱或测序文件；并在 Addgene/NCBI 检索该质粒的衍生登记。无论是否回复，均按路线 A 推进。

**路线 C（降级，仅当路线 A 裁决失败时启用）**：

- **F1**：仅用宿主参考跑全流程；报告声明"外源基因表达无法定量"，一致性评估中外源相关结论标 `not_assessable`；
- **F2**：以 *N. crassa* 转录本作独立旁证定量，结论标注"探索性、非与原研究可比"；
- **F3**：外源基因完全移出主分析，仅作 Limitations 条目。

#### 3.4.3 序列裁决与接头组装

1. **判别实验（先做）**：每样本取 1 M reads，`bowtie2 --very-sensitive-local`（或 `minimap2 -ax sr`）比对到候选序列，输出 `refs/custom/sequence_adjudication.tsv`，字段 = 候选（`gh1-1_native` / `gh1-1a_FigS8` / `cdt-1_Myc` / `cdt-1_GSGS-sfGFP` / `cdt-1_GSGS-eGFP` / `eGFP` / `sfGFP` / `yEGFP`）× 指标（比对率、覆盖度均匀性、错配数、错配位点）。
   - `gh1-1_native` 应接近 100% 覆盖且零系统性错配；`gh1-1a_FigS8` 应在约 25% 位置出现错配（阴性对照）。
   - 标签裁决：sfGFP 与 eGFP 仅差数个碱基（F64L/S65T、F99S/M153T/V163A 等），在 ~150× 覆盖下可无歧义区分；输出差异位点的 per-base 计数作为证据。
2. **接头与边界组装**：对未比对/软剪切 reads 用 `rnaSPAdes` 或 MEGAHIT 组装；以 `eGFP/sfGFP` 起始与终止、`cdt-1` 终止、`CYC1` 序列为锚，拼出跨接头 contig。**验收**：接头两侧各 ≥ 100 bp 支持，跨接头 reads ≥ 20 条；5'/3' 端以覆盖断点定义，若无法确定则转录本条目只取 CDS 主体并在 PROVENANCE 标注"边界为观测估计"。
3. **残差核验**：以组装序列为参考做 `bcftools mpileup`，要求 **CDS 区域错配率 < 0.1% 且无系统性错配**；所有非参考位点逐一判读留档（`results/qc/transgene_assembly/`）。
4. **蛋白层核验**：cdt-1 ORF 翻译 = NCU00801 蛋白（580 aa）+ 融合标签，帧内无提前终止；gh1-1 ORF 翻译 = 477 aa + His6。
5. **循环性声明（写入报告）**：外源参考部分来源于样本自身 reads → **外源的绝对 TPM 不可独立验证**；但 C8/G8 **比值**对参考偏差不敏感（同一参考、同一偏差用于两条件）。因此 **H7 只报比值、配对统计与置信区间，不报绝对丰度**。

#### 3.4.4 记录与验收

- `refs/custom/PROVENANCE.md`：每元件一行，字段 = 元件 / 来源类型（公共库 | 文献引物推导 | 数据组装）/ accession 或文献位置 / 长度 / md5 / 裁决依据；
- `refs/custom/sequence_adjudication.tsv`：3.4.3(1) 的判别矩阵；
- `results/qc/transgene_assembly/`：接头 contig、跨接头 read 计数、错配清单；
- **P1 出口条件**：上述三项齐备，或明确转入路线 C 并在报告中相应标注。

### 3.5 自定义参考的构建规范

**量化单元 = 外源转录本条目**：`cdt-1EGFP`（CDS + 标签 + 接头，按 3.4.3 定界）与 `gh1-1`（CDS + His6），各附 `tx2gene` 映射。

**"量化转录本"与"decoy"是两个不同槽位**：

| 序列 | 作为量化转录本 | 作为 decoy | 理由 |
| --- | --- | --- | --- |
| 宿主 16 染色体 + 线粒体 | 用其转录本（`rna_from_genomic`） | **是**（gentrome 的 decoy） | 吸收未注释/基因间区 reads，降低假阳性 |
| cdt-1 / gh1-1 外源转录本 | **是** | 否 | 主分析对象 |
| 质粒骨架：pBluescript 部分 | 否 | 可选（建议加入） | 大肠杆菌来源，作为 decoy 把这批 reads 从"未比对"变为"decoy"，对任何基因计数零影响，仅改善比对率口径 |
| 质粒骨架：2μ 部分 | 否 | 可选（建议加入） | 参考中不存在天然 2μ 序列，不构成重复目标；加入 decoy 只是让这批 reads 有归属 |
| **URA3** | 否（野生型位点已在宿主中） | **禁止加入** | 一旦加入，质粒来源的 URA3 reads 会被 decoy 吸走，破坏拷贝数代理指标（3.3 / 2.2(e)） |
| P_PGK1 / T_CYC1（质粒拷贝） | 否 | 可选 | 序列与基因组相同；启动子不进入 mRNA，终止子会进入 mRNA，属已知微量串扰，由 5.12 检查义务处理 |

- **主分析索引**：宿主转录本 + 两条外源转录本，decoy = 宿主基因组。
- **对照索引**：decoy 额外加入质粒骨架（剔除 URA3 区段），作为比对率口径的敏感性对照（预期计数完全一致；若不一致必须查明）。两次索引的参数与结果记入 `refs/VERSIONS.md` 与 `results/comparison/decoy_sensitivity.tsv`。
- **诊断参考（供 5.6 使用）**：`refs/diagnostic/` 备 ① 质粒骨架（pRS426：pBluescript + 2μ，来自公共库，记录来源与 md5）② 天然 2μ 质粒序列 ③ 胞质 rRNA（RDN1 的 35S/5S）④ 线粒体 rRNA（21S/15S）⑤ adapter 序列（TruSeq）。这些序列只用于 read fate 分类，不进入主索引。

---

## 4. 实验设计与统计模型

### 4.1 设计

- 6 个样本，2 组 × 3 生物学重复，单一批次，厌氧指数期（80 g/L 起始糖，取样于残余糖 ≈ 50 g/L）。
- 样本命名中 `JCYL00xB / 00xD` 共享编号，是否存在配对/批次关系在公开元数据中无法确证；原文在 Methods 中明确使用非配对两组比较。

### 4.2 主模型与敏感性模型

| 角色 | 设计公式 | 残差自由度 | 说明 |
| --- | --- | --- | --- |
| **主分析** | `~ condition` | 4 | 与原文一致，保守 |
| **敏感性分析** | `~ replicate + condition` | 3 | 3×2 平衡设计，数学上恒合法 |

两个模型都必须运行并对比：一致 → 结论对设计假设稳健；不一致 → 报告必须讨论差异并说明取舍。配对并非生物学确证，仅作敏感性，理由记入 `docs/decisions/`。

### 4.3 设计与方法学注意事项

- **2μ 高拷贝质粒**：外源基因计数受拷贝数影响，**不等于表达强度**；质粒还带来额外的转录负担。以 URA3 代理（2.2(e)）作定量约束。
- **无 ERCC spike-in**：无法评估绝对定量与动态范围。
- **单批次**：无独立批次效应评估空间；n = 3 vs 3 下 SVA/RUV 不可用（其自由度会吞掉条件效应），报告中须写明"已考虑但不适用"。
- **50 bp 单端短读**：映射率与异构体分辨低于现代 100 bp 双端数据。
- **poly-A 建库与线粒体转录本**：原文仅写"TruSeq RNA Sample Prep Kit"，未写"poly-A"；GEO 的 `library_selection` 记为 `cDNA`。TruSeq RNA Sample Prep Kit 的标准流程含 oligo(dT) 富集，故"poly-A 富集"属**基于试剂盒的推定**，报告中必须标注为推定。
  线粒体方面：**芽殖酵母线粒体 mRNA 缺少可供 oligo(dT) 高效捕获的长 poly(A) 尾；其成熟 3' 端由保守十二聚体元件 AAUAA(U/C)AUUCUU 定义，且该细胞器内未鉴定到 PAP 活性**（Chang & Tong 2012 原文："Budding yeast (*S. cerevisiae*) mitochondrial mRNAs are not polyadenylated, and no PAP activity has been identified in this organelle in yeast. Instead, they carry a conserved dodecamer sequence…"）。文献中存在方向相反或细化的报道（1982 *MCB* 的 oligoadenylate；2024 *RNA* 31:208 的 3' 端加工），故报告统一表述为"捕获效率低"，不使用绝对化的"不 poly(A) 化"。
  据此：**H1 的检验主体限定为核编码的线粒体功能基因**；mtDNA 编码基因的检出情况作为显式检查项（5.6），若接近 0 则判定 `not_assessable`。
- **测序质量编码**：2014 年 GA II 数据可能为 Phred64（Illumina 1.3–1.7 管线）。必须先由 FastQC 判定编码，fastp 相应加 `--phred64`，判定依据留档（`docs/decisions/phred_encoding.md`）。误判会造成系统性"假低质量"与过度修剪。
- **文库链特异性**：TruSeq RNA Sample Prep Kit v1 与后续 stranded 版本行为不同。以 `salmon --libType A` 的逐样本推断结果为准；若样本间推断不一致，必须查明原因后再继续（5.4）。

---

## 5. 分析流程

每一步给出：目的 / 输入 / 工具与参数 / 输出 / 验收。

### 5.0 环境初始化与仓库

- **目的**：建立可重建的分析环境、项目骨架与版本控制。
- **工具**：WSL2 Ubuntu；conda 26.1.1（不装 mamba）；R 4.4.3 + Bioconductor 3.20（由 conda 提供；不使用 renv 项目激活，理由见 `docs/decisions/r_environment.md`）；Git + GitHub。
- **步骤**：
  1. `git init`，配置 `user.name` / `user.email`；
  2. 添加 GitHub 远端（`git remote add origin <GITHUB_REPO>`），首次推送；凭据方式（SSH key 或 HTTPS + PAT）由项目负责人提供，凭据本身不入库；
  3. 建立目录骨架（§6），`git mv` 归档 7 个文献文件并生成 `docs/literature/CHECKSUMS.md`；
  4. 编写 `environment.yml` / `renv.lock` / `.gitignore` / `README.md` / `run_all.sh` / `config/paths.sh`。
- **输出**：`environment.yml`、`env/conda-explicit.txt`、`env/r_packages.tsv`、`.gitignore`、`README.md`、`run_all.sh`、`config/paths.sh`、`env/tool_versions.txt`、`.git` + 远端。
- **验收**：新 shell 中 `conda env create -f environment.yml` 成功；`conda create --file env/conda-explicit.txt` 可逐位重建；R 关键包（DESeq2/tximport/apeglm/clusterProfiler/fgsea/org.Sc.sgd.db/GO.db）可加载（脚本断言）；`git log` 有首个 commit 且远端可见；`env/tool_versions.txt` 含 Bioconductor release。
- **注意**：`data/raw/`、`refs/index/`、`*.bam`、`logs/` 不入库；凭据、token 一律不写入仓库。

### 5.1 数据获取与校验

- **目的**：取得全部原始测序数据与元数据并确保完整性。
- **输入**：`config/samplesheet.tsv`
- **工具**：`curl`（https）+ `md5sum`；ENA Portal API filereport；GEO SOFT。
- **输出**：`data/raw/fastq/*.fastq.gz`、`checksums/raw.md5`、`data/metadata/ena/ena_filereport.tsv`、`data/metadata/geo_soft/GSM*.txt`、`data/metadata/original_study/geo/GSE54825_Cellobiose_versus_Glucose.xlsx.gz`、`logs/download.log`
- **验收**：6 个文件齐全；**自算 md5 与 ENA 公布 md5 逐一相等**（值见 2.3）；read 数与 ENA filereport 一致；GEO/ENA 元数据文件落库。

### 5.2 参考体系构建与核查

- **目的**：构建冻结的宿主参考、含外源序列的自定义参考，以及 5.6 的诊断参考。
- **输入**：`GCF_000146045.2`（genomic fna + rna_from_genomic + cds_from_genomic + gff）、SGD GAF + GO slim、3.4 的外源序列、pRS426/2μ/rRNA 序列（诊断用）。
- **工具**：`datasets` CLI（锁定版本）/ `gffread` / `samtools faidx` / 序列追加脚本 / `salmon index`。
- **输出**：`refs/host/`、`refs/custom/`（`transcripts.fa` + `tx2gene.tsv` + `PROVENANCE.md`）、`refs/diagnostic/`、`refs/index/salmon/{k31_decoy_host, k25_decoy_host, k31_decoy_host_plus_backbone}/`、`refs/VERSIONS.md`、`refs/checks/locus_check.tsv`、`refs/checks/id_map_loss.tsv`。
- **ID 映射损失表**：输出 `6351（Dataset S1）→ RefSeq 基因 ID`、`519 → RefSeq` 的映射成功率与未映射清单。已知 Dataset S3 中 FunSpec 有 18 个 "not found"（含 `ARG5,6`、`MFALPHA292`、多个 Ty 元件），本项目须给出自己的数字并与之对账。
- **验收**：3.3 清单逐项通过并留档；外源转录本可被 `tx2gene` 正确映射；索引参数（k-mer、decoy、salmon 版本）被记录；ID 映射损失表产出。

### 5.3 读段 QC 与修剪

- **目的**：评估测序质量并做轻度修剪。
- **输入**：`data/raw/fastq/`
- **工具与参数**：`FastQC`（先判定 Phred 编码）→ `fastp`（单端，`--detect_adapter_for_se`，如判定 Phred64 则加 `--phred64`；轻度修剪，`--length_required 36`）；`MultiQC`。
- **输出**：`results/qc/fastqc/`、`results/qc/fastp/`（含 `*_fastp.json` 的 adapter 比例）、`results/qc/multiqc_report.html`、`data/processed/fastq/`、`docs/decisions/phred_encoding.md`。
- **验收**：per-base 质量、adapter 含量、重复率、长度分布齐备；修剪前后可对比；Phred 编码判定有留档；adapter 残留率 < 1%。
- **附加**：保留一份**未修剪**数据，供 5.5 使用。

### 5.4 定量（Salmon，主路线）

- **目的**：产出基因水平表达矩阵。
- **输入**：修剪后 FASTQ（SE 50 bp）+ `refs/index/salmon/k31_decoy_host/`
- **工具与参数**（按锁定的 Salmon 版本核对有效 flag 后写死在脚本中；`--validateMappings` / `--gcBias` / `--seqBias` 的默认值在 0.14→1.10 之间发生过变化，不得照抄教程）：
  - `salmon quant --libType A --gcBias --seqBias --threads 16`（若该版本已废弃 `--validateMappings` 则删除该 flag 并记录）；
  - `tximport(..., type = "salmon", countsFromAbundance = "lengthScaledTPM")` → `gene_counts.tsv`（README 中写明该矩阵为 length-scaled TPM 计数，非原始 read 计数）；`gene_tpm.tsv` 由对应 abundance 直接给出；
  - `DESeqDataSetFromTximport()` 承接下游；`countsFromAbundance` 的选择记入 `docs/decisions/`。
- **输出**：`results/quantification/salmon/<sample>/quant.sf` + `lib_format_counts.json`、`gene_counts.tsv`、`gene_tpm.tsv`、`results/qc/salmon_libtype.tsv`。
- **验收**：
  - 各样本映射率、libType 推断结果、检出基因数入表；样本间 libType 不一致则必须停下查明；
  - `lib_format_counts.json` 的多映射/decoy 比例作为**一级比对诊断**（5.6 为二级归因）；
  - TPM 与原文 RPKM 的相关性作为 sanity check（数值进 5.12）；
  - 外源转录本检出情况单独记录，并与 2.2(e) 的锚定值并列。

### 5.5 定量参数敏感性

- **目的**：确认定量结果不依赖于关键参数选择。
- **比较项**：`--kmerLen`（31 vs 25）；修剪 vs 未修剪；decoy 是否含质粒骨架（3.5）。
- **输出**：`results/comparison/param_sensitivity.tsv`（各设置下 log2FC 相关性、DEG 重叠率、检出基因数）、`results/comparison/decoy_sensitivity.tsv`。
- **验收**：明确回答"结论是否依赖参数选择"，给出数字而非断言。

### 5.6 比对诊断与 read fate

- **目的**：在没有全量 BAM 的前提下，为三类问题提供实测证据——被排除的质粒骨架去往何处、rRNA 占比、线粒体转录本检出情况。
- **方法**：每样本用 `seqtk` 子采样 2–5 M reads，两级比对：
  1. 比对到宿主基因组 + 外源转录本（`bowtie2` 或 `minimap2`）→ 比对率、线粒体序列 reads 占比与覆盖度、mtDNA 编码基因检出情况、重复率、多重映射热点（Ty / rDNA）、天然 2μ 来源 reads；
  2. 取第 1 级未比对 reads，比对到诊断参考（质粒骨架 + 天然 2μ + 胞质 rRNA + 线粒体 rRNA + adapter，见 3.5）→ reads 去向分类。
- **输出**：`results/qc/postmap/`（read fate 表、chrM 统计、mt 基因检出表、多映射热点统计）。
- **验收与判定规则（已按实测修正，2026-09-14）**：
  - read fate 采用三级阶梯归因：端到端比对 → 局部模式回收 → 诊断参考（质粒骨架/天然 2μ/adapter）。
    实测（每样本 ~1.14 M reads）：端到端比对 **94.5–95.5%**；未比对 4.5–5.5%，其中
    **19.9–22.8% 可由局部模式回收**（源于端到端模式偏严 + GA II 读段错误）、
    **35.6–40.4% 命中质粒骨架/天然 2μ/adapter**；**残余 2.2–3.0%（占全库）未完全归因**，
    候选原因在 Limitations 列明（读段错误、E. coli 等污染、R64 参考未覆盖的亚端粒/Ty1-Y′ 区、
    以及启动子–CDS、CDS–终止子等参考中不存在的嵌合 reads）。
    **因此原定"各桶合计 ≥ 95%"不适用于未比对桶**；改为：全库口径下已归因比例 ≥ 97%，
    且残余部分必须给出量化占比与候选原因。
  - 若 mtDNA 编码基因计数接近 0 → 判定 `not_assessable`，H1 仅在核编码范围检验并写明原因；
  - 子采样诊断用于定性归因，不与原文 CLC 映射率作严格跨工具对标；原文 92.3%/83.4%/76.3% 仅作量级背景参考。

### 5.7 样本水平分析

- **目的**：评估样本整体结构与重复一致性。
- **输入**：`gene_counts.tsv`
- **方法与参数**：归一化（VST/rlog）；样本相关性（Spearman/Pearson）；PCA；层次聚类。
- **输出**：归一化矩阵；`results/figures/pca.png`、`sample_correlation.png`、`clustering.png`。
- **验收**：两条件是否分离有明确结论；同组内重复一致性有明确结论（给出相关系数区间）。n = 3 vs 3，PCA 仅作定性描述。

### 5.8 差异表达

- **目的**：鉴定 cellobiose 相对 glucose 的差异表达基因。
- **输入**：`gene_counts.tsv`
- **主分析**：`DESeq2`，设计 `~ condition`，独立过滤开启。
- **LFC 口径（写死在脚本中）**：
  - `results()` 的 `padj` 基于 **MLE（未收缩）LFC**；`lfcShrink(type = "apeglm")` 只替换 `log2FoldChange` 列；
  - **阈值判定一律使用未收缩 LFC**；收缩 LFC 仅用于排序与 MA/volcano/GSEA 展示；
  - 该口径记入 `docs/decisions/` 并在脚本注释中标注。
- **主阈值**：`padj < 0.05` 且 `|log2FC| > 1`；同时输出放宽（`padj < 0.1`）与收紧（`padj < 0.01` 且 `|log2FC| > 1.5`）两档。
- **敏感性分析**：`~ replicate + condition`（配对模型）；**阈值对齐原文**版本（`|FC| ≥ 2` 且 `FDR ≤ 0.001`，`results(alpha = 0.001)`，未收缩 LFC）。原文阈值作用于 RPKM 倍数，本研究作用于归一化计数倍数，属口径差异，须在 5.12 注明。
- **输出**：`DEG_primary.tsv`、`DEG_relaxed.tsv`、`DEG_stringent.tsv`、`DEG_paired.tsv`、`DEG_threshold_matched.tsv`（每表含 baseMean / LFC(MLE) / lfcSE / stat / pvalue / padj / LFC(apeglm)）。
- **验收**：上/下调基因数完整（分宿主/外源两类计数）；与原文 519 的重叠已计算；主分析与配对模型的结论差异有明确说明。

### 5.9 模型诊断

- **内容**：p 值直方图、dispersion 拟合图、MA plot（shrunken LFC）、Cook's distance 与离群样本检查。
- **输出**：`results/figures/model_diagnostics/`。
- **验收**：p 值分布形态合理（无系统性反保守）；无样本被 Cook's cutoff 剔除而未被讨论。

### 5.10 功能分析

- **目的**：解析差异基因涉及的生物学过程与通路。
- **工具与参数**：
  - ORA：`clusterProfiler`；**GO slim 需自建映射**——用 SGD `go_slim_mapping.tab`（或等价文件）构建 `TERM2GENE` 后走 `enricher()`；`org.Sc.sgd.db` 仅作便利来源（注释有滞后），主路径以 3.1 冻结的 SGD GAF + GO slim 为准；
  - 通路：KEGG（`sce`），**原始结果缓存入库**（`results/enrichment/kegg_raw/`），保证离线可复现；
  - **GSEA**：以收缩 LFC 排序，`fgsea` / `clusterProfiler::GSEA`。
- **背景集（universe）**：定义为"**在 ≥3 个样本中 count ≥ 10 的基因**"；该集合的基因数、与 6351 的差异写入报告；外源条目是否入 universe 单独声明。
- **与原文的方法学差异**：原文 **FunSpec + Bonferroni（P 0.01）**（且含 18 个未能映射的基因）；本研究 clusterProfiler + BH。该差异须在报告与 5.12 中注明。
- **输出**：`results/enrichment/`（GO slim / KEGG / GSEA 表）、`results/figures/enrichment/`。
- **验收**：背景集定义、校正方法、富集阈值均可见；结果可独立重跑。

### 5.11 生物学与工程解读

在查看结果**之前**固定以下假设清单，之后逐条判定为「复现 / 未复现 / 新发现 / 无法评估」。**每条判定必须附集合水平统计量**（该基因集 log2FC 的 Wilcoxon 检验 vs 全基因背景，或 fgsea NES/p）。

| ID | 假设（cellobiose 相对 glucose） | 代表基因 / 集合 |
| --- | --- | --- |
| H0 | 两条件整体转录谱可区分 | 全基因组（PCA / 聚类，见 5.7） |
| H1 | 线粒体功能激活：TCA 循环、电子传递链、OXPHOS、ATP 合成上调 | CIT1/2/3, ACO1, IDH1/2, KGD1/2, LSC1/2, SDH1-4, FUM1, MDH1, COX4-9, QCR*, CYT1, RIP1, ATP1-7（核编码基因为检验主体） |
| H2 | 氨基酸与硫胺素合成下调 | MET*, CYS*, ARG*, HIS*, THI* |
| H3 | 葡萄糖感知与信号通路仅部分激活 | PKA（CYR1, GPA2, RAS1/2, TPK1-3, BCY1, IRA1/2）；Snf3-Rgt2-Rgt1（SNF3, RGT2, RGT1, MTH1, STD1）；Snf1-Mig1（SNF1, MIG1, MIG2, HXK2, GRR1） |
| H4 | 己糖转运蛋白家族重排（缺少胞外葡萄糖） | HXT1-7, GAL2, SNF3, MPH2/3 |
| H5 | 厌氧 / 固醇 / 缺氧相关基因变化 | SUT1, ERG*, DAN/TIR/PAU |
| H6 | 储存碳与糖异生相关基因变化 | TPS1/2, TSL1, GSY1/2, GLC3, PCK1, FBP1 |
| H7 | 外源通路基因（cdt-1 / gh1-1）行为可被定量 | 预注册判据见下 |
| H8 | 由 H3 / H5 可导出工程候选靶点，且覆盖原文线索 | SUT1 上调（原文 FC 3.40）；**HAP4 转录上调（原文 FC 6.12），而"HAP4 缺失改善发酵"是表型层结论，方向非单调**；DAL80（原文 FC 2.17） |
| H9 | （增补，见 5.15）cellobiose 与 glucose 的差异部分由生长速率 / 葡萄糖阻遏解释 | 核糖体蛋白基因（RPL*/RPS*）、rRNA 加工基因作为生长速率代理 |

**H7 的预注册判据**：

- 检验对象：`cdt-1EGFP` 与 `gh1-1` 的 **C8/G8 比值**（不报绝对 TPM，理由见 3.4.3(5)）；
- 原文锚定值：cdt-1EGFP = 1.871（log2 0.904）、gh1-1 = 2.160（log2 1.111）；
- 方向一致性判据：本研究比值 > 1 且 log2FC 与锚定值之差在 n = 3 的置信区间内；
- 混杂控制：同时报告 (i) 质粒拷贝数代理 URA3 的比值（原文 1.354）；(ii) **URA3 归一化后的外源比值**（原文推导值：cdt-1 = 1.382、gh1-1 = 1.596）；(iii) 内源 `PGK1` 的比值（原文 0.745）作为同启动子反证；
- 判定为「复现」的条件：方向一致且 URA3 归一化后仍 > 1；若归一化后 ≈ 1，则判定为「由质粒拷贝数差异解释，非转录调控」，属新发现而非复现。

- **输出**：`results/interpretation/hypothesis_scorecard.tsv`（含每条假设的统计量）+ 报告章节。
- **验收**：每条假设都有明确判定、数字依据与未复现项说明。

### 5.12 与原研究的一致性评估

- **目的**：给出量化的可比性结论，而非定性印象。

| 指标 | 说明与要求 |
| --- | --- |
| per-gene log2FC 相关性 | 共同基因上的 Spearman / Pearson（本研究 vs Dataset S1 与 GEO xlsx）；在统一符号约定后计算（原文 `log` 列 = log2(C8/G8)，见 2.2(d)） |
| DEG 集合重叠 | 阈值对齐版 vs Dataset S2 的 519 条：Venn / 混淆矩阵 / Jaccard / 超几何检验 p；**外源条目（gh1-1）单独标记**，不计入宿主基因分子分母 |
| 跨阈值比较 | 输出 precision/recall 随阈值变化曲线（单点重叠在 n = 3 下噪声过大） |
| 方向一致率 | 双方皆显著基因中方向一致的比例（按 2.2(d) 换算） |
| 数字对账 | 复现 `519 = 244 UP + 256 DOWN + 19`？逐项对账并解释差额（含 FunSpec 的 18 个 not found） |
| 定量层一致性 | 本研究 TPM/count 与原文 RPKM 的线性关系与离群基因 |
| 富集层对照 | 本研究 GO 结果 vs Dataset S3（FunSpec + Bonferroni 0.01），UP/DOWN 两侧分别对照 |
| 结论层对照 | 519 基因、线粒体上调、氨基酸与硫胺素下调是否可重现 |
| 映射层背景对照 | 原文 92.3% / 83.4% / 76.3% 仅作量级参考；跨工具严格对标不可行，需在报告中说明 |
| 可复现性边界 | CLC 的多重比对处理为软件内部算法（原文仅 76.3% 唯一比对），因此**绝对 RPKM/TPM 不会与原文完全一致**，可复现的是基因水平趋势 |
| 质粒拷贝数代理 | URA3 的 C8/G8 比值；URA3 归一化后的外源比值；URA3 vs cdt-1/gh1-1 的逐样本散点图 |
| 串扰检查义务 | `PGK1` 与 `CYC1` 的逐样本计数与条件偏移；`CYC1` 因本身是 H1 基因（原文 FC 2.53）**无法分离**，须在 Limitations 明写 |
| ID 映射损失 | 6351 / 519 → RefSeq 的映射率与未映射清单（对应 5.2 的 `id_map_loss.tsv`） |

- **输出**：`results/comparison/`、报告对应章节。
- **验收**：给出一致性的数值与图，并显式列出 `not_assessable` 项及原因。

### 5.13 定量路线 B（可选，本项目关闭）

- **启用条件**：需要 BAM 级诊断（IGV/覆盖度检查）或需要方法学对比章节，且磁盘与时间允许。本机 RAM 15 GB，**本项目不启用**。
- **工具与参数**（若后续启用）：`STAR`（`--outSAMtype BAM SortedByCoordinate`，可选 `--twopassMode Basic`）→ `samtools index` → `featureCounts`（`-t exon -g gene_id`，单端不加 `-p`；`-s` 由实测链特异性决定）。
- **替代方案**：诊断能力由 5.6 的子采样轻量比对提供。
- **若启用**，比较指标：比对率、多映射比例、检出基因数、与 Salmon 的 log2FC 相关性、DEG 重叠率、wall time、峰值内存、磁盘占用。

### 5.14 报告产出

- **形式**：**中英双语两版**（`docs/report/report_zh.qmd`、`docs/report/report_en.qmd`），由同一数据源自动生成，两版的数字、表、图完全一致；渲染为 HTML + PDF。
- **图表**：所有图的坐标轴标签、图例、图注、表头、表注**统一使用英文**（两版共用同一批图文件）。
- **排版**：参照正规期刊（编号章节、Vancouver 编号引用、LaTeX 排版，图注由 Quarto 自动编号）。参考文献条目由 BibTeX 统一维护；PDF 由 **xelatex** 生成（中文版通过 xeCJK 指定 Noto Sans SC 静态字体）。
- **验收**：`quarto render` 一次命令生成两版；两版数字一致性由脚本断言校验（同一 TSV 输入）。

### 5.15 增补分析（P2.5，主结果产出后执行）

> 本组事项不影响流程能否跑通，但与结论质量直接相关。默认在 P2 完成、主结果落盘之后执行；执行前在 `docs/decisions/` 记录决策。

| ID | 事项 | 预期产出 |
| --- | --- | --- |
| E1 | **生长速率 / 葡萄糖阻遏混杂分析（H9）** | 以核糖体蛋白基因与 rRNA 加工基因作生长速率代理，与 H1/H2 基因集做重叠与效应量对照；给出"差异中有多少可归因于生长速率/葡萄糖阻遏"的量化讨论 |
| E2 | **LFC 感知检验（`lfcThreshold = 1` 的 Wald 检验）** | 作为统计上站得住的推断集，与双重过滤集并列报告（后者降级为可比性用途） |
| E3 | **URA3 内参化的外源定量重分析** | 逐样本 URA3 归一化后的外源比值 + 置信区间 |
| E4 | **与原文的偏离表** | 单表列出参考体系、比对器、定量器、归一化、统计、富集六个层面本研究与原文做法的差异及理由 |
| E5 | **仓库归档（Zenodo DOI）** | 可引用的版本化发布 |

---

## 6. 目录结构

```text
<PROJECT_ROOT>/          # 项目根
├── README.md                     # 如何从零重跑
├── PLAN.md                       # 本文件
├── environment.yml
├── run_all.sh                    # 顺序执行 00→08，幂等（存在即跳过 + 校验）
├── .gitignore                    # 排除 data/raw、refs/index、*.bam、logs/*
│
├── config/
│   ├── samplesheet.tsv
│   └── paths.sh                  # 全流程统一路径变量
│
├── docs/
│   ├── literature/
│   │   ├── original_study/       # 7 个只读输入文件
│   │   └── CHECKSUMS.md
│   ├── acceptance/               # 阶段出口验收记录（P0.md …）
│   ├── report/                   # report_zh.qmd / report_en.qmd / refs.bib / 渲染产物
│   └── decisions/                # phred 编码、countsFromAbundance、LFC 口径、5.13 启停、E1–E5
│
├── data/
│   ├── raw/fastq/                # 不入库
│   ├── processed/fastq/
│   └── metadata/
│       ├── ena/ena_filereport.tsv
│       ├── geo_soft/GSM*.txt
│       └── original_study/
│           ├── geo/GSE54825_Cellobiose_versus_Glucose.xlsx.gz
│           └── derived/          # dataset_s1_rpkm.tsv、dataset_s2_deg.tsv、dataset_s3_goslim.tsv、
│                                 # table_s1_plasmids.tsv、table_s2_primers.tsv、galazka2010_cloning_primers.tsv
│
├── refs/
│   ├── host/                     # GCF_000146045.2（fasta 不入库，gff/gtf 与 md5 入库）
│   ├── custom/                   # transcripts.fa + tx2gene.tsv + PROVENANCE.md + sequence_adjudication.tsv
│   ├── diagnostic/               # 质粒骨架、天然 2μ、rRNA、adapter（仅用于 5.6）
│   ├── index/                    # salmon 索引（不入库）
│   ├── checks/                   # locus_check.tsv、id_map_loss.tsv
│   └── VERSIONS.md
│
├── env/tool_versions.txt         # CLI + R + Bioconductor 版本快照
├── env/conda-explicit.txt        # 精确锁定（conda list --explicit）
├── env/r_packages.tsv            # R 包版本清单
│
├── scripts/
│   ├── 00_setup/  01_download/  02_qc/  03_quantification/
│   ├── 04_differential_expression/  05_enrichment/
│   ├── 06_visualization/  07_comparison/
│   └── 08_original_study/        # 原文补充材料解析、外源序列裁决与重建
│
├── notebooks/                    # 探索性分析；不承载最终结论
│
├── results/
│   ├── qc/{fastqc,fastp,postmap,salmon_libtype.tsv,transgene_assembly/}
│   ├── quantification/   differential_expression/   enrichment/
│   ├── comparison/       interpretation/            figures/
│
├── logs/
└── checksums/
```

---

## 7. 项目产物

### 7.1 数据与分析文件

samplesheet；参考与注释版本记录（含 annotation release 与 md5）；位点核查表；ID 映射损失表；外源序列 provenance + 序列裁决表；QC 报告（含 read fate、libType、线粒体检出）；基因计数与 TPM 矩阵；归一化矩阵；DEG 表（主/放宽/收紧/配对/阈值对齐）；富集表（GO slim/KEGG/GSEA）；一致性表（含超几何 p、Jaccard、跨阈值曲线、URA3 归一化结果）；参数与 decoy 敏感性表；假设判定表（含集合水平统计量）。

### 7.2 图表

QC summary（含 read fate、libType、线粒体检出）；PCA、样本相关性热图、层次聚类；MA plot（shrunken）、volcano plot、DEG 热图；模型诊断图组；GO slim / KEGG / GSEA 富集图；锚定图（外源逐样本 RPKM/TPM 对比、URA3 vs 外源散点、HAP4/SUT1 的 per-sample 表达、葡萄糖感知通路热图）；一致性图（log2FC 散点相关、DEG 重叠 Venn、跨阈值 PR 曲线）。所有图的标签与图注为英文。

### 7.3 报告结构（中英双语，同一结构）

```text
1.  Background
2.  Research Question
3.  Dataset
4.  Experimental Design
5.  Methods（版本、参数、阈值、软件、LFC 口径）
6.  Quality Control（读段 QC + libType + read fate + 线粒体检出 + 比对率背景对照）
7.  Quantification（Salmon 主路线 + 参数/decoy 敏感性）
8.  Sample-level Analysis
9.  Differential Expression（含模型诊断与配对敏感性）
10. Functional Enrichment（注明与 FunSpec 的方法学差异）
11. Biological Interpretation（H0–H9 判定表 + 集合水平统计量）
12. Comparison with Original Study（量化一致性 + 偏离表 + 可复现性边界）
13. Limitations（参考偏差、外源可评估性与循环性、2μ 拷贝数、短读长、无 ERCC、单批次、poly-A 推定与线粒体、PGK1/CYC1 串扰、CYC1 不可分离）
14. Engineering Implications
15. Reproducibility（环境、命令、校验、如何重跑）
16. Conclusion
```

### 7.4 代码仓库

源码、流程脚本、图、报告（双语）、环境文件、README，以及全部可复现凭证（CHECKSUMS、VERSIONS、PROVENANCE、SOFT/ENA 元数据、ID 映射表、决策记录）。大文件（FASTQ、索引、BAM）不入库，仅保留下载与构建脚本及校验值。远端仓库用于对外展示项目结构与流程完整度。

---

## 8. 风险与限制

| 风险 / 限制 | 影响 | 应对 |
| --- | --- | --- |
| 外源标签变体不确定（sfGFP vs eGFP、Myc vs GFP） | 外源定量偏差 | 3.4.3 判别实验以 reads 裁决；候选并行评估 |
| 外源参考部分来自样本自身 reads | 绝对 TPM 不可独立验证 | 只报比值与配对统计；Limitations 声明循环性 |
| 外源参考版本误用（`gh1-1a` / Figure S8） | 比对率骤降、结论错误 | 3.4.2 明确禁用；3.4.3 判别实验含阴性对照 |
| 50 bp 单端短读 | 映射率与异构体分辨下降 | 5.5 参数敏感性；报告中说明 |
| Phred 编码误判（GA II 时代风险） | 过度修剪、假低质量 | 5.3 显式判定 + 留档 |
| Salmon flag/默认值随版本漂移 | 参数无效或不可复现 | 锁定版本 + 按版本核对有效 flag |
| n = 3，残差自由度 4 | 统计功效有限，阈值敏感 | 主模型 + 配对敏感性 + 三档阈值 + GSEA + 跨阈值 PR |
| 配对关系无法确证 | 模型选择不确定 | 两个模型都跑并对比；非配对为主分析 |
| 生长速率 / 葡萄糖阻遏混杂 | H1/H2 的"复现"可能只是慢生长效应 | E1（5.15）：生长速率代理基因集对照，量化归因 |
| 2μ 质粒拷贝数波动 | 外源基因定量不可靠 | 不将外源计数当表达强度；URA3 代理 + 预注册判据 |
| 质粒酵母终止子串扰 | 对 PGK1 / CYC1 产生微量计数偏移 | 3.5 的槽位设计 + 5.12 检查义务 + Limitations |
| `CYC1` 串扰不可分离 | 无法纯化该项 | 明写 Limitations |
| 参考中无天然 2μ | 2μ 来源 reads 成为未比对 | 5.6 read-fate 显式统计 |
| poly-A 推定对线粒体转录本的盲区 | mtDNA 编码基因计数接近 0，易被误读为"无变化" | 5.6 显式检查；H1 限核编码；表述按 4.3 |
| 株系与参考差异 / `ura3Δ0` 边界未核实 | 少数位点计数异常；URA3 代理强度不确定 | 3.2 边界条件 + 3.3 核查 |
| 原文统计与富集框架不同 | 直接比较不公平 | 阈值对齐 + 偏离表 + 可复现性边界 |
| KEGG 在线依赖 | 富集结果不可复现 | 缓存原始结果入库 |
| ID 映射损失（Ty 元件、合并名 `ARG5,6` 等） | 重叠率被低估 | 5.2 映射损失表 + 5.12 对账 |
| 范围蔓延 | 项目无法收尾 | 1.4 非目标 + 9.1 交付范围 + E1–E5 延后 |

---

## 9. 执行顺序与验收门

### 9.1 交付范围

- **本次交付 = P0–P3**（含 3.4 的外源重建；若转入路线 C 则按 F1 交付并标注）。
- **MVP 完成线**：原始数据校验通过 → 参考与 ID 映射留档 → 计数/TPM 矩阵 → DEG 五表 → 富集三表 → 假设记分卡（带统计量）→ 一致性量化结论（含 `not_assessable` 清单）→ 中英双语报告可渲染。
- **P4 / P5 / 5.15(E1–E5) 不属本次交付**；E1 在 P2 之后立即执行。

| 阶段 | 内容 | 出口验收 |
| --- | --- | --- |
| **P0** 准备 | git init + GitHub 远端、目录、7 个文献文件归档 + md5、环境、samplesheet、原始数据与元数据下载校验 | 环境可从零重建；6 个 FASTQ 与 ENA md5 一致；SOFT/ENA/filereport/GEO xlsx 落库；远端可见首个 commit |
| **P1** 参考与定量 | 3.3 位点核查 + ID 映射损失表、3.4 外源序列裁决与重建、读段 QC（含 Phred 判定）、Salmon 定量、read fate 与线粒体检出 | 计数与 TPM 矩阵；`PROVENANCE.md` + `sequence_adjudication.tsv` + 接头证据齐备；read fate 解释 ≥ 95% 未比对 reads |
| **P2** 主分析 | 参数/decoy 敏感性、样本水平分析、差异表达（主 + 敏感性 + 阈值对齐）、模型诊断、功能分析 | DEG 表、诊断图、富集表齐备且可重跑；背景集与 LFC 口径写入决策记录 |
| **P2.5** 增补 | E1（H9 生长速率混杂）、E2（lfcThreshold 检验） | 见 5.15 |
| **P3** 对比与报告 | 一致性评估（超几何/Jaccard/PR 曲线/URA3 归一化）、假设判定、中英双语报告 | 报告含量化一致性结论、偏离表、可复现性边界与 Limitations |
| **P4** 扩展 | GSE69404（4 种糖 × 好氧/厌氧） | 独立目录，复用同一流程骨架 |
| **P5** 工程化封装 | Snakemake 薄封装 | 封装输出与编号脚本输出一致（校验通过） |

**执行要点**：3.4 是 P1 的前置任务，其判据是"**先用 1 M reads 裁决序列，再建索引**"。

---

## 10. 扩展

### 10.1 GSE69404（多因子扩展）

- **内容**：同一实验室的工程菌（cdt-1 + 密码子优化 gh1-1 + XYL1/XYL2/XKS1 + cdt-2 + gh43-2 + gh43-7），glucose / cellobiose / xylose / xylodextrins × 好氧 / 厌氧，生物学三重复。
- **价值**：把单因子设计扩展为双因子设计（碳源 × 氧），引入"protein biosynthesis 下调"这一新维度，并可部分缓解 H9 的生长速率混杂。
- **方式**：独立 samplesheet 与独立报告，不与 GSE54825 混合计算。

### 10.2 Snakemake 封装

策略：编号脚本作为可读的 ground truth；项目末尾加一层薄封装，并验证其输出与脚本输出一致。

### 10.3 长期路线图（非本项目验收内容）

```text
Transcriptomics → Protein engineering → Structure-guided design → ML-assisted engineering
```

---

## 附录 A：外部事实来源

**一手记录**

- **GEO Series GSE54825**（SOFT）：系列摘要、总体设计、样本列表、平台 GPL9377、补充文件、SRA 关系。
- **GEO Samples GSM1324496–GSM1324501**（SOFT，`form=text&view=brief`）：`Sample_title` 与 `carbon source` 字段（条件与 run 的对应关系）。
- **ENA Portal API filereport（PRJNA237759）**：6 个 run 的 accession、`library_layout=SINGLE`、`instrument_model=Illumina Genome Analyzer II`、read/base count、FASTQ https 直链、官方 md5（值见 2.3）。
- **NCBI RefSeq `GCF_000146045.2` 组装报告**：17 条序列（BK006934–BK006949 + NC_001224.1），宿主 12,071,326 bp，线粒体 85,779 bp；不含 2μ 质粒。

**文献**

- **Lin Y, et al.** *Biotechnol Biofuels.* 2014;7:126. PMID 25435910；PMC PMC4243952；doi 10.1186/s13068-014-0126-6（正文、Methods、Dataset S1–S3、Fig S1–S8、Table S1–S2）。
- **Galazka JM, Tian C, Beeson WT, Martinez B, Glass NL, Cate JH.** *Science.* 2010;330(6000):84-86. doi 10.1126/science.1192838 及其 SOM：外源基因"amplified by PCR from cDNA"；pRS426 + P_PGK1 构建；cdt-1 的 Myc 版与 GSGS–superfolder GFP 融合版引物；gh1-1 的 C 端 6×His 引物；所有构建体用 CYC1 终止子。
- **Sikorski RS, Hieter P.** *Genetics.* 1989;122:19-27（pRS 系列质粒来源）。
- **Tian C, et al.** *PNAS.* 2009;106:22157-22162（*N. crassa* 植物细胞壁降解系统分析；cdt-1/gh1-1 的基因鉴定来源）。
- **NCBI RefSeq**：`XM_958708.2`（cdt-1 / NCU00801，CDS 1740 bp）、`XM_011395456.1`（gh1-1 / NCU00130，CDS 1431 bp）。
- **SGD**（yeastgenome.org）：基因系统性名与标准名、gene association file、GO slim、`ura3Δ0` 等位描述。
- **酵母线粒体 RNA 3' 端与 poly(A)**：Chang DD, Tong AH. 2012（PMC3307840）："Budding yeast (*S. cerevisiae*) mitochondrial mRNAs are not polyadenylated, and no PAP activity has been identified in this organelle in yeast. Instead, they carry a conserved dodecamer sequence, AAUAA(U/C)AUUCUU, at their 3′-ends."；另见 1982 *Mol Cell Biol* 2:450（oligoadenylate）与 2024 *RNA* 31:208（芽殖酵母 3' 端加工的物种特异性元件）。
- **TruSeq RNA Sample Prep Kit**（Illumina）标准流程含 oligo(dT) 富集 → 本项目记为"poly-A 富集（基于试剂盒推定）"，原文未明写。

## 附录 B：版本与可复现性记录规范

| 文件 | 记录内容 |
| --- | --- |
| `refs/VERSIONS.md` | 参考 assembly accession、RefSeq annotation release、下载日期、每个文件 md5；两次 decoy 索引的参数与结果 |
| `refs/custom/PROVENANCE.md` | 每个外源元件的来源类型、accession/文献位置、长度、md5、裁决依据 |
| `refs/custom/sequence_adjudication.tsv` | 3.4.3(1) 的候选×指标判别矩阵（含 `gh1-1a` 阴性对照） |
| `refs/checks/locus_check.tsv` | 3.3 清单逐项检查结果（含 `ura3Δ0` 区间核实结论） |
| `refs/checks/id_map_loss.tsv` | 6351 / 519 → RefSeq 的映射率与未映射清单 |
| `docs/literature/CHECKSUMS.md` | 7 个只读输入文件的 md5 |
| `data/metadata/{ena,geo_soft,original_study}/` | 使附录 A 的每条事实离线可复现 |
| `env/tool_versions.txt` | fastp / FastQC / MultiQC / Salmon / bowtie2 / samtools / seqtk / R 及包版本 + Bioconductor release |
| `env/conda-explicit.txt` | 精确环境锁定（conda list --explicit）；与 `environment.yml` 共同构成环境凭证 |
| `env/r_packages.tsv` | R 包名称/版本/库路径快照（177 个包） |
| `checksums/raw.md5` | 原始数据自算 md5（与 ENA 公布值比对结论） |
| `logs/` | 各步骤运行日志 |
| 运行结束的 `sessionInfo()` | R 包版本快照 |
| `docs/decisions/` | phred 编码、`countsFromAbundance`、LFC 口径、5.13 启停、E1–E5 决策与阈值选择依据 |
