# 决策记录：测序质量编码判定

**日期**：2026-09-13 23:36:18
**判定依据**：FastQC 对每样本 1 M reads 子采样的 `Encoding` 字段（见 `results/qc/fastqc/raw/*/`）

| 项 | 值 |
| --- | --- |
| FastQC 报告的 Encoding | `Sanger / Illumina 1.9` |
| 结论 | Phred33（Sanger / Illumina 1.8+，或 fastqc 报告为 ASCII 33 起） |
| fastp 参数 | ``（空表示使用默认 Phred33） |

判定命令：`grep '^Encoding' results/qc/fastqc/raw/*/fastqc_data.txt`
