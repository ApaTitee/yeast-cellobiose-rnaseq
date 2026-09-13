# 参考体系版本与校验值

生成时间：2026-09-13 22:48:14 HKT

## 宿主参考

- Assembly: `GCF_000146045.2` (S288C R64-1-1；16 染色体 + 线粒体 `NC_001224.1`，不含天然 2μ)
- RefSeq 注释：SGD R64-5-1 | release_date=2026-07-10 | provider=SGD
- 来源：`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/`

| 文件 | 字节 | md5（NCBI 官方） |
| --- | --- | --- |
| `GCF_000146045.2_R64_genomic.fna.gz` | 3843460 | `88c38b957b721dfc50e6c4df03b6242e` |
| `GCF_000146045.2_R64_rna_from_genomic.fna.gz` | 3054486 | `605c359a9417cd65358a1d65c1ffdff3` |
| `GCF_000146045.2_R64_cds_from_genomic.fna.gz` | 3026445 | `cdd083d7b301e312804ea773f12953d4` |
| `GCF_000146045.2_R64_genomic.gff.gz` | 2203839 | `1fddbd976c4ce61e8a36a3908d46c25d` |
| `GCF_000146045.2_R64_assembly_report.txt` | 2603 | `dc747ca264896255e63aa1f05d56dcee` |
| `md5checksums.txt` | 6046 | `8a1d8c22193c28803703017d10aa4746` |

## 功能注释

| 文件 | 来源 | 字节 | md5 |
| --- | --- | --- | --- |
| `go_slim_mapping.tab` | https://downloads.yeastgenome.org/curation/literature/go_slim_mapping.tab | 4197139 | `cc1a6091823a465f5ce06c7d120ec7ce` |
| `goslim_yeast.obo` | https://current.geneontology.org/ontology/subsets/goslim_yeast.obo | 152859 | `588ae838289c141334f6aae78e192f77` |
| `sgd.gaf.gz` | https://current.geneontology.org/annotations/sgd.gaf.gz | 2602501 | `bd1438fe970c0fcd735f3f17df83a281` |

## 外源候选与诊断序列

来源：NCBI E-utilities `efetch`（nuccore）

| 文件 | 记录 | 用途 | 字节 | md5 |
| --- | --- | --- | --- | --- |
| `egfp_candidates_U55762.fa` | U55762.1 (pEGFP-N1) | 3.4 序列裁决候选 | 1811 | `7ffaefdb6be8f77be12ecaed8be7bbf4` |
| `ncrassa_cdt-1_NCU00801_candidate_CDS.fa` | XM_958708.2 (NCU00801) | 3.4 序列裁决候选 | 1950 | `9c48707606ad99406c470f4e3b137861` |
| `ncrassa_gh1-1_NCU00130_candidate_CDS.fa` | XM_011395456.1 (NCU00130) | 3.4 序列裁决候选 | 1641 | `606919bd13e4e43f85a2cfa6b300452c` |
| `pegfpN1_record_U55762.gb` | - | 接头/边界核对 | 9370 | `cd9645eb3ba9a72fb664adc31f44a7e4` |
| `plasmid_2micron_NC001398.fa` | NC_001398.1 (2μ circle) | 5.6 read fate 诊断参考 | 6497 | `538bd526349a400a36aa763fce6607d1` |
| `plasmid_prs426_U03451.fa` | U03451.1 (pRS426, URA3 marker) | 5.6 read fate 诊断参考 | 5884 | `d1ef803fbf7de917c3747a291682efb7` |
| `plasmid_prs426_U03451.gb` | U03451.1 (pRS426, URA3 marker) | 接头/边界核对 | 8646 | `448f5c031a96bf2001e005268f06c7c4` |
| `sfgfp_candidate_JQ341914.fa` | PX636966.1 | 3.4 序列裁决候选 | 842 | `cbc54b3d7dd1619c24cc67ccc4bc63c4` |
| `sfgfp_candidate_MW132720.fa` | PX636966.1 | 3.4 序列裁决候选 | 854 | `8baba229c5bade5fb3f637d9c4c4fd44` |
| `sfgfp_candidate_PX636966.fa` | PX636966.1 | 3.4 序列裁决候选 | 842 | `602103f50e6fb605f1aef402e0a03e7f` |
| `sfgfp_record_PX636966.fa` | PX636966.1 | 接头/边界核对 | 796 | `a26f96363717cf6c9d4046e9f0c459e3` |

## 说明

- 所有 URL 均写死版本，不使用 `latest` 类动态地址。
- 脚本：`scripts/01_download/02_fetch_reference.sh`。
- 完整 md5 清单：`checksums/reference.md5`。
