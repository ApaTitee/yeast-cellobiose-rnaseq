# 决策记录：salmon 2.7（rust 版）EffectiveLength 异常及其影响评估

**日期**：2026-09-14
**状态**：已评估，结论为"不影响本项目的相对定量结论"，继续使用，并在报告中声明
**影响文件**：`scripts/03_quantification/10_salmon_quant.sh`、`docs/decisions/`

## 现象

1. **单端数据默认 FLD 错误**：salmon 2.7（rust 版）在未给先验时对 SE 50 bp 数据使用
   `fragment length mean (sd) = 250.00 (25.00)`，导致 `EffectiveLength` 完全失真
   （实测 2472 bp 转录本报 40563）。**已修复**：按样本实测读长传入 `--fldMean 49 --fldSD 3`。
2. **修复后 EffectiveLength 仍与转录本长度不成合理关系**：
   - 宿主转录本：中位 `efflen/len ≈ 0.465`（50 bp 单端理论上应 ≈ 1 − 50/L，即 >0.9）；
   - 外源转录本：`efflen` **大于** 长度 10 倍以上（cdt-1EGFP：2472 bp → 26488）。

## 影响评估（三项独立检验）

| 检验 | 结果 | 结论 |
| --- | --- | --- |
| 有/无 decoy 对照（同一样本、同一参数） | 共同 6295 条目；`efflen/len` 中位 **0.465 vs 0.465**；NumReads 比值中位 0.992；TPM 的 log2 差异中位 **0.027**，5960 条中仅 17 条 >1 | decoy 设计安全，且异常与 decoy 无关 |
| 与原文独立定量（RPKM）的相关性 | 5796 个可比基因 **Spearman = 0.880** | 相对定量整体可信 |
| 关键基因点对点核对 | PGK1 6220 vs 5972；THI4 7979 vs 6772；MET17 5402 vs 4538；**URA3 432 vs 472**；ATP1 187 vs 171 | 无系统性偏移 |

## 决策

1. **继续使用 salmon 2.7 的 TPM/计数**，因为它与独立来源（CLC RPKM）高度一致，且对 decoy 选择不敏感。
2. 报告中使用"**TPM 表观语义异常**"的措辞，并把上述三项检验写入 Methods/Limitations。
3. **外源基因（cdt-1EGFP / gh1-1）只报 C8/G8 比值**（H7 预注册判据），并用三套相互独立的估计交叉验证：
   - 本记录 TPM 比值；
   - 1 M reads 子采样的读数比值（裁决实验，1.89 / 2.22）；
   - 原文锚定值（1.871 / 2.160）。
4. 若后续需要绝对 TPM 的严格解释，可改用 `--noLengthCorrection` 或换用 salmon 1.x 复算；
   该项列入 P2.5 增补分析（见 PLAN 5.15）。

## 修复后的可复现命令

```bash
salmon quant -i refs/index/salmon/k31_decoy_host -l A \
  -r data/processed/fastq/<sample>.fastp.fastq.gz -p 8 \
  --validateMappings --gcBias --seqBias --dumpEq \
  --fldMean 49 --fldSD 3 -o results/quantification/salmon/<sample>
```
