# 只读输入文件校验值

原始文献与补充材料，位于 `docs/literature/original_study/`。这些文件是本研究"与原研究对照"的唯一依据，
内容不得修改；如需变更，必须同时更新本文件的 md5 并在 `docs/decisions/` 记录原因。

| 文件 | 内容 | 大小 (bytes) | md5 |
| --- | --- | --- | --- |
| `s13068-014-0126-6.pdf` | Lin et al. 2014, *Biotechnol Biofuels* 7:126（正文） | 2152116 | `82a64a8127fecb9df6f3a4777001fe75` |
| `13068_2014_126_MOESM1_ESM.xlsx` | Additional file 1 = Dataset S1（全基因组 normalized RPKM，6351 条） | 631122 | `98a937e749414d5dbd0145914569b149` |
| `13068_2014_126_MOESM2_ESM.xlsx` | Additional file 2 = Dataset S2（519 DEG + 19 TF + 7 GRN） | 65584 | `dd09b0b5c6a0b1d8e4c6ee6b04c24a82` |
| `13068_2014_126_MOESM3_ESM.xlsx` | Additional file 3 = Dataset S3（FunSpec GO，UP 244 / DOWN 256） | 12954 | `98ec1472114123a85096cb2dcd87c613` |
| `13068_2014_126_MOESM4_ESM.pdf` | Additional file 4 = Supporting information（Fig S1–S8、Table S1–S2） | 14239883 | `2ce0d449ba7caf012d065cfcf02775e4` |
| `science.1192838.pdf` | Galazka et al. 2010, *Science* 330:84-86（正文） | 308948 | `b072b556036459ccbce66a5ce6cb0646` |
| `galazka-som.pdf` | Galazka et al. 2010 Supporting Online Material（含全部克隆引物） | 887516 | `5d3c1197ff3a118632376c84d5cf3c6d` |

校验命令：

```bash
cd docs/literature/original_study && md5sum -c <(sed -n 's/^| `\([^`]*\)`.*| `\([0-9a-f]\{32\}\)` |$/\2  \1/p' ../CHECKSUMS.md)
```
