#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05a_read_fate_metrics.py — 汇总 read fate 诊断的两级比对结果（PLAN 5.6）

一级比对：宿主基因组 + 外源转录本（bowtie2，-k 2 以便统计多映射）
二级比对：一级未比对的 reads -> 诊断参考（质粒骨架 pRS426 / 天然 2μ / adapter / rRNA 等）

输出（results/qc/postmap/）：
  read_fate.tsv        逐样本读数去向分类
  chr_distribution.tsv 逐染色体（含线粒体）reads 占比
  mt_gene_detection.tsv mtDNA 编码基因的检出判定（决定 H1 的可评估范围）
  rrdna_ty_hotspots.tsv rDNA / Ty 等重复区域的 reads 占比（解释比对率缺口与多映射）

用法： python3 05a_read_fate_metrics.py <bam> <bam> ...
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "qc", "postmap")
HOST_GFF = os.path.join(ROOT, "refs", "host", "GCF_000146045.2_R64_genomic.gff.gz")
DIAG_FA = os.path.join(ROOT, "refs", "diagnostic", "diagnostic_ref.fa")
UNMAPPED_DIR = os.path.join(ROOT, "results", "qc", "postmap", "unmapped")
IDX1 = os.path.join(ROOT, "refs", "index", "bowtie2", "host_plus_transgenes")


def sh(cmd: list[str]) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def parse_gff_features():
    """从 GFF 取：线粒体基因坐标、rDNA 与 Ty 区域坐标（用于热点解释）。"""
    import gzip
    mt, hotspots = [], []
    with gzip.open(HOST_GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, ftype, start, end, attrs = f[0], f[2], int(f[3]), int(f[4]), f[8]
            if chrom == "NC_001224.1" and ftype in ("gene", "mRNA"):
                m = re.search(r"(?:gene|locus_tag)=([^;]+)", attrs)
                if m:
                    mt.append((m.group(1), start, end, ftype))
            if ftype == "gene":
                m = re.search(r"(?:gene|Name)=([^;]+)", attrs)
                name = m.group(1) if m else ""
                if name in ("RDN1", "RDN5-1", "RDN5-2", "RDN5-3", "RDN5-4", "RDN5-5",
                            "RDN18-1", "RDN25-1", "RDN58-1") or name.startswith("TY"):
                    hotspots.append((chrom, name, start, end))
    # 合并同一基因的多个 feature（取最长的 mRNA）
    mt_genes = {}
    for name, s, e, t in mt:
        if t == "mRNA":
            mt_genes[name] = (s, e)
    return mt_genes, hotspots


def depth_regions(bam: str, regions: list[tuple]) -> list[float]:
    """返回每个区域的平均深度（samtools depth）。区域可为 (chrom, start, end) 或 (chrom, name, start, end)。"""
    out = []
    for reg in regions:
        chrom, s, e = reg[0], reg[-2], reg[-1]
        try:
            txt = sh(["samtools", "depth", "-a", "-r", f"{chrom}:{s}-{e}", bam])
        except subprocess.CalledProcessError:
            out.append(0.0)
            continue
        vals = [int(l.split("\t")[2]) for l in txt.splitlines() if l.strip()]
        out.append(sum(vals) / len(vals) if vals else 0.0)
    return out


def main() -> int:
    bams = sys.argv[1:]
    if not bams:
        raise SystemExit(__doc__)
    os.makedirs(OUT, exist_ok=True)
    mt_genes, hotspots = parse_gff_features()

    fate_rows, chr_rows, mt_rows, hot_rows = [], [], [], []
    for bam in bams:
        sample = os.path.basename(bam)[:-4]
        mapped_in_bam = int(sh(["samtools", "view", "-c", bam]).strip())
        # 唯一/多重：-k 2 时 XS 标签表示存在次优比对
        n_multi = int(sh(["bash", "-lc",
                          f"samtools view {bam} | grep -c 'XS:i:' || true"]).strip() or 0)
        idx = sh(["samtools", "idxstats", bam])
        mapped = 0
        unmapped = 0
        for line in idx.splitlines():
            f = line.split("\t")
            if f[0] == "*":
                unmapped = int(f[3])
                continue
            n = int(f[2])
            mapped += n
            chr_rows.append([sample, f[0], n, ""])
        # 线粒体
        mt_reads = int(sh(["bash", "-lc", f"samtools idxstats {bam} | awk '$1==\"NC_001224.1\"{{print $3}}'"]).strip() or 0)
        # 一级补：把一级未比对的 reads 用**局部模式**重比对到同一参考
        # （端到端模式对 GA II 的老数据偏严；局部模式可量化"因比对模式/测序错误而未比对"的部分）
        import gzip
        n_unmapped_used = 0
        n_local_recovered = 0
        ufq = os.path.join(UNMAPPED_DIR, f"{sample}.fq.gz")
        residual2 = os.path.join(OUT, f".{sample}.residual2.fq.gz")
        if os.path.exists(ufq):
            with gzip.open(ufq, "rt") as fh:
                n_unmapped_used = sum(1 for _ in fh) // 4
            if n_unmapped_used > 0 and os.path.exists(IDX1 + ".1.bt2"):
                sam1b = os.path.join(OUT, f".{sample}.local1b.sam")
                sh(["bash", "-lc",
                    f"bowtie2 -x {IDX1} -U {ufq} --very-sensitive-local -k 1 -p 8 "
                    f"--un-gz {residual2} -S {sam1b} 2>/dev/null"])
                n_local_recovered = sum(1 for l in open(sam1b)
                                        if not l.startswith("@") and l.split("\t")[2] != "*")
                os.remove(sam1b)

        # 二级：残余未比对 reads 的去向（质粒骨架 / 天然 2μ / adapter）
        diag_counts: Counter = Counter()
        ufq = residual2 if os.path.exists(residual2) else ufq
        n_residual = 0
        if os.path.exists(ufq):
            with gzip.open(ufq, "rt") as fh:
                n_residual = sum(1 for _ in fh) // 4
        if n_residual > 0 and os.path.exists(DIAG_FA + ".1.bt2"):
            out_sam = os.path.join(OUT, f".{sample}.diag.sam")
            sh(["bash", "-lc",
                f"bowtie2 -x {DIAG_FA} -U {ufq} --very-sensitive-local -k 1 -p 8 --no-unal -S {out_sam} 2>/dev/null"])
            for line in open(out_sam):
                if line.startswith("@"):
                    continue
                f = line.split("\t")
                if len(f) > 2:
                    diag_counts[f[2]] += 1
            os.remove(out_sam)

        total = mapped_in_bam + n_unmapped_used
        mapped = mapped_in_bam
        unmapped = n_unmapped_used
        fate_rows.append([sample, total, mapped, unmapped, n_local_recovered, n_multi,
                          round(100.0 * mapped / total, 2) if total else 0,
                          round(100.0 * n_multi / total, 2) if total else 0,
                          mt_reads, round(100.0 * mt_reads / total, 2) if total else 0,
                          n_residual,
                          sum(diag_counts.values()),
                          round(100.0 * sum(diag_counts.values()) / n_residual, 2) if n_residual else 0,
                          ";".join(f"{k}={v}" for k, v in diag_counts.most_common())])

        # 线粒体基因检出
        if mt_genes:
            names = list(mt_genes)
            depths = depth_regions(bam, [("NC_001224.1", *mt_genes[n]) for n in names])
            for n, d in zip(names, depths):
                mt_rows.append([sample, n, mt_genes[n][0], mt_genes[n][1], round(d, 2),
                                "detected" if d >= 1.0 else ("trace" if d > 0 else "not_detected")])

        # 重复热点（rDNA / Ty）
        if hotspots:
            hd = depth_regions(bam, hotspots)
            for (chrom, name, s, e), d in zip(hotspots, hd):
                if d > 0:
                    hot_rows.append([sample, chrom, name, s, e, round(d, 2)])

    def write(path, header, rows):
        with open(os.path.join(OUT, path), "w", encoding="utf-8") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        print(f"  {path}: {len(rows)} 行")

    write("read_fate.tsv",
          ["sample", "reads_total", "mapped_genome_e2e", "unmapped_e2e", "recovered_local", "multi_XS",
           "pct_mapped", "pct_multi", "chrM_reads", "pct_chrM", "residual_after_local",
           "diag_hits", "pct_of_residual_diag", "diag_breakdown"], fate_rows)
    write("chr_distribution.tsv", ["sample", "seqname", "reads", "note"], chr_rows)
    write("mt_gene_detection.tsv",
          ["sample", "mt_gene", "start", "end", "mean_depth", "status"], mt_rows)
    write("rrdna_ty_hotspots.tsv",
          ["sample", "chrom", "feature", "start", "end", "mean_depth"], hot_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
