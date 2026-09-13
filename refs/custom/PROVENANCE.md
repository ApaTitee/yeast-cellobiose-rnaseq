# 外源序列 PROVENANCE（证据驱动的重建记录）

**生成时间**：2026-09-13
**方法**：`PLAN.md` §3.4.2 路线 A —— 以公共序列为候选骨架，用**样本自身 reads** 裁决版本与标签，
并对接头的存在性做覆盖度检验。裁决实验：`scripts/08_original_study/07_transgene_adjudication.sh`
（每样本 1 M reads，bowtie2 `-k 5 --very-sensitive-local`，按最小编辑距离归属）。

---

## 1. 裁决结论（一句话）

| 元件 | 结论 | 关键证据 |
| --- | --- | --- |
| `gh1-1` | **天然 NCU00130 cDNA CDS（1431 bp，去内含子、非密码子优化）+ C 端 6×His** | 与公共 CDS **逐位一致（0 处差异）**；His6 特有区（1432–1449）覆盖 **31.9×** |
| `gh1-1a`（本论文 Figure S8） | **样本中不存在**（阴性对照） | 仅 2/6,000,000 reads 归属，且均为错配；与天然 CDS 核酸同一性 74.5% |
| `cdt-1` | **天然 NCU00801 cDNA CDS（1740 bp，去内含子）** | 与公共 CDS **逐位一致（0 处差异）**，平均覆盖 69.5× |
| 标签 | **superfolder GFP（sfGFP）**，密码子接近人源化 eGFP，**不是** GenBank 中的 sfGFP 记录，也**不是**裸 eGFP | reads 共识与 pEGFP-N1 的 EGFP CDS 仅 7 nt 差异；翻译后为 eGFP 骨架 **+ sfGFP 的全部 6 个 superfolder 替换**（见 §3） |
| 质粒骨架 / 天然 2μ | 样本中存在（预期内），仅用于 5.6 read fate 归类 | 见 §4 |

> 该结论**取代**了基于文献的推断。`PLAN.md` §3.4.1 原先记录的两项待裁决问题（标签变体、gh1-1 版本）至此闭合。

## 2. 逐位一致性（reads 共识 vs 公共候选）

共识由 `07c_consensus_from_reads.py` 生成（mpileup 多数碱基，要求 ≥60% 支持，否则记 N）。

| 元件 | 候选来源 | 长度 | 平均覆盖 | **与候选的差异位点** | 未定(N) |
| --- | --- | --- | --- | --- | --- |
| `cdt-1` | XM_958708.2 CDS（NCU00801） | 1740 bp | 69.5×（cellobiose）/ 42.0×（glucose） | **0** | 4 |
| `gh1-1` | XM_011395456.1 CDS（NCU00130） | 1431 bp | 41.1× / 21.6× | **0** | 4–5 |
| 标签 | U55762.1（pEGFP-N1）egfp CDS | 720 bp | 69.8× / 37.0× | **7**（6 处非同义） | 4–5 |

两个条件、两个样本的结论完全一致（标签的 7 处差异在同一位置重复出现），排除单样本假象。

## 3. 标签身份的判定过程

| 候选 | 归属 reads（6 样本 × 1 M） | NM=0 | 判定 |
| --- | --- | --- | --- |
| pEGFP-N1 EGFP CDS（U55762.1） | **4,964** | 2,771 | 骨架最接近 |
| `sfGFP` JQ341914.1 | 8 | 5 | 密码子差异过大 |
| `sfGFP` MW132720.1 | 0 | 0 | 同上 |
| `sfGFP` PX636966.1 | 0 | 0 | 同上（且 N 端为 MRKGEELFTGVV） |

共识翻译与 EGFP 的 6 处氨基酸差异：

| 位置（U55762 EGFP 帧） | EGFP → 共识 | 密码子 | **sfGFP 编号** | sfGFP 已知突变 |
| --- | --- | --- | --- | --- |
| 31 | S → R | TCC → CGC | 30 | **S30R** ✓ |
| 40 | Y → N | TAC → AAC | 39 | **Y39N** ✓ |
| 106 | N → T | AAC → ACC | 105 | **N105T** ✓ |
| 146 | Y → F | TAC → TTC | 145 | **Y145F** ✓ |
| 172 | I → V | ATC → GTC | 171 | **I171V** ✓ |
| 207 | A → V | GCC → GTC | 206 | **A206V** ✓ |

**6/6 命中 sfGFP 的 superfolder 替换集合**（编号偏移 1 由 N 端 Met-Val 的差异解释），同时保留 EGFP 的 F64L/S65T。
这与 Galazka et al. 2010 SOM 的描述（"superfolder GFP with an N-terminal linker of Gly-Ser-Gly-Ser"）一致；
2014 年论文中的 "eGFP" 属宽松表述。

**结论**：参考中应使用 **reads 共识序列**作为标签条目，而不是任何 GenBank 记录
（GenBank 的 sfGFP 记录使用不同密码子，会系统性丢失比对）。

## 4. 未做/待做（进入 `refs/custom/transcripts.fa` 前必须完成）

| 项 | 状态 | 说明 |
| --- | --- | --- |
| `gh1-1` C 端 6×His | **已确认存在** | His6 特有区覆盖 31.9×，且有 31 条比对越过天然 CDS 终止子 |
| 标签与 `cdt-1` 的**接头精确碱基** | **已确定** | reads 的软剪切片段直接给出：`…GTTGCT` + **`ATCGAT`（ClaI）** + **`GGTAGTGGTAGT`（Gly-Ser-Gly-Ser）** + `GTGAGC…`；与 Galazka 2010 SOM 的克隆设计（ClaI 位点 + GSGS 接头）完全一致；证据片段见 §4.1 |
| 两条转录本的 **5'/3' 边界**（P_PGK1 → CDS；CDS/tag → T_CYC1） | **待界定** | 以覆盖断点定义；无法确定时只取 CDS 主体并在此处标注 |
| 质粒骨架中 URA3 区段 | 不进参考 | 保持 URA3 的"唯一映射"状态以维持拷贝数代理（3.5） |

### 4.1 接头证据（软剪切 reads 片段）

在标签上带 5' 软剪切的片段（参考方向）中直接观察到：

| 片段 | 解读 |
| --- | --- |
| `GTTGCTATCGATGGTAGTGGTAGT` | cdt-1 末 6 bp（`GTTGCT`）+ ClaI（`ATCGAT`）+ GSGS（`GGTAGTGGTAGT`） |
| `CGATGGTAGTGGTAGT` | 同上，起点前移 2 bp |
| `TGGTAGT` | GSGS 的尾部 |

另有两条相互独立的证据支持该结构：
1. cdt-1 的**终止密码子不存在于任何 read**（含终止子的 20-mer 命中 0 次），而"去终止子"版的 20-mer 命中 20 次 → 融合确实去掉了 stop；
2. 标签的比对起始位置集中在公共 eGFP 参考的**第 4 位**（12 reads），无任何 read 起始于第 1 位 → 构建体中标签不带起始 ATG，与"载体提供 Met/接头"的融合结构一致。

### 已确认的关键序列与校验值

| 序列 | 文件 | 长度 | md5 |
| --- | --- | --- | --- |
| 标签（reads 共识，sfGFP） | `refs/custom/adjudication/consensus/tag_eGFP_U55762.SRR1166445_JCYL001D_cellobiose.fa` | 720 bp | `c4957a75897c70cf86f134d65597ecfb` |
| `gh1-1a`（Figure S8，阴性对照） | `refs/diagnostic/gh1-1a_FigureS8.fa` | 1449 bp | `6477425dc9c09fad97ef54d5d050dc16` |
| `cdt-1` 候选 CDS | `refs/diagnostic/ncrassa_cdt-1_NCU00801_candidate_CDS.fa` | 1740 bp | `9c48707606ad99406c470f4e3b137861` |
| `gh1-1` 候选 CDS | `refs/diagnostic/ncrassa_gh1-1_NCU00130_candidate_CDS.fa` | 1431 bp | `606919bd13e4e43f85a2cfa6b300452c` |
| 候选集（10 条） | `refs/custom/adjudication/candidates.fa` | — | `12713f7168a26d01b51f158a7185da79` |

## 5. 裁决实验的副产物（对 5.11/5.12 直接有用）

每样本 1 M reads 的归属计数（`results/qc/transgene_assembly/adjudication_detail.tsv`）：

| 指标 | glucose 均值 | cellobiose 均值 | 比值 | 原文锚定值（2.2(e)） |
| --- | --- | --- | --- | --- |
| `cdt-1` reads / 1 M | 1,392.7 | 2,625.3 | **1.89** | 1.871 |
| `gh1-1` reads / 1 M | 1,205.0 | 2,679.7 | **2.22** | 2.160 |

→ 用完全不同的工具链（bowtie2 + 1 M 子采样读数计数）**独立复现了原文的外源基因升高幅度**
（1.89 vs 1.871；2.22 vs 2.160）。这是 H7 的第一条独立证据，先于主定量流程产生。

另注：工程质粒骨架来源的 reads 在 cellobiose 中也升高（≈1.65×），而天然 2μ 来源 reads 基本不变
（≈1.0×）—— 提示升高的是**工程质粒的拷贝数**而非天然 2μ 质粒，与 2.2(e) 中 URA3 代理的
预注册值（1.354×）方向一致。该观察将在 5.12 用主定量结果（URA3 计数）复核后写入报告。

## 6. 循环性声明（写入报告）

标签序列与 C 端接头部分来源于样本自身的 reads。因此：

- 外源转录本的**绝对** TPM 不可独立验证；
- 但 C8/G8 **比值**对参考偏差不敏感（同一参考、同一偏差用于两条件），故 H7 只报比值与配对统计；
- 本文件的每一项结论都同时给出了**可独立获得的反证**（公共 CDS 逐位一致、Figure S8 阴性对照、
  覆盖度证据、两条件一致），以降低循环论证的风险。
