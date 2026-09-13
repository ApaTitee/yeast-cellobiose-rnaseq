#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07b_summarize_adjudication.py
---------------------------------------------------------------------------
目的：把裁决比对的 BAM 汇总为"候选 × 样本"判别矩阵（PLAN 3.4.3 第 1 步的可读产物）。
方法：bowtie2 以 -k N 报出每条 read 的多个候选比对；对每条 read 取 **NM 最小**（并列时 AS 最大）
      的候选作为该 read 的归属；再按候选汇总。真版本应表现为 NM=0 的 read 占绝对多数，
      假版本则出现系统性错配。
输入：refs/custom/adjudication/{cand}.bam（每样本一个）
输出：refs/custom/sequence_adjudication.tsv
      results/qc/transgene_assembly/adjudication_detail.tsv
---------------------------------------------------------------------------
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ADJ = os.path.join(ROOT, "refs", "custom", "adjudication")
DETAIL_DIR = os.path.join(ROOT, "results", "qc", "transgene_assembly")
OUT = os.path.join(ROOT, "refs", "custom", "sequence_adjudication.tsv")
NM = re.compile(r"NM:i:(\d+)")
AS = re.compile(r"AS:i:(-?\d+)")
CIGAR_SOFT = re.compile(r"^(\d+)S")


def parse_bam(bam: str):
    """返回 (每样本统计, read->最优候选 计数)。"""
    best: dict[str, tuple[int, int, str]] = {}   # read -> (nm, -as, ref)
    per_ref = defaultdict(lambda: {"hits": 0, "nm0": 0})
    total = 0
    proc = subprocess.Popen(["samtools", "view", "-F", "4", bam], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:                      # type: ignore[union-attr]
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        total += 1
        read, ref = f[0], f[2]
        nm = int(NM.search(line).group(1)) if NM.search(line) else 999
        asc = int(AS.search(line).group(1)) if AS.search(line) else -9999
        per_ref[ref]["hits"] += 1
        if nm == 0:
            per_ref[ref]["nm0"] += 1
        key = (nm, -asc, ref)
        if read not in best or key < best[read]:
            best[read] = key
    proc.wait()                                   # type: ignore[union-attr]
    winner = defaultdict(int)
    winner_nm0 = defaultdict(int)
    for _, (nm, _, ref) in best.items():
        winner[ref] += 1
        if nm == 0:
            winner_nm0[ref] += 1
    return total, per_ref, winner, winner_nm0, len(best)


def main() -> int:
    os.makedirs(DETAIL_DIR, exist_ok=True)
    bams = sorted(f for f in os.listdir(ADJ) if f.endswith(".bam"))
    if not bams:
        raise SystemExit("!! 未找到 BAM，请先运行 07_transgene_adjudication.sh")
    rows, detail = [], []
    for b in bams:
        sample = b[:-4]
        total, per_ref, winner, winner_nm0, n_reads = parse_bam(os.path.join(ADJ, b))
        for ref in sorted(set(list(per_ref) + list(winner))):
            rows.append([sample, ref, winner.get(ref, 0), winner_nm0.get(ref, 0),
                         per_ref.get(ref, {}).get("hits", 0),
                         per_ref.get(ref, {}).get("nm0", 0)])
        detail.append((sample, total, n_reads))
        print(f"  {sample}: 命中行 {total}, 可归属 read {n_reads}")

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("sample\tcandidate\treads_best_NM\ttotal_reads\treads_best_NM0\n")
        agg: dict[str, list[int]] = defaultdict(lambda: [0, 0])
        for r in rows:
            agg[r[1]][0] += r[2]
            agg[r[1]][1] += r[3]
        for cand in sorted(agg):
            fh.write(f"ALL\t{cand}\t{agg[cand][0]}\t\t{agg[cand][1]}\n")

    with open(os.path.join(DETAIL_DIR, "adjudication_detail.tsv"), "w", encoding="utf-8") as fh:
        fh.write("sample\tcandidate\treads_best_NM\treads_best_NM0\tall_hits\thits_NM0\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"\n== 全样本汇总（reads_best_NM = 按最小编辑距离归属到该候选的 read 数）==")
    for cand in sorted(agg):
        print(f"  {cand:26s} {agg[cand][0]:10d}  (其中 NM=0: {agg[cand][1]})")
    print(f"\n已写入:\n  {os.path.relpath(OUT, ROOT)}\n  {os.path.relpath(os.path.join(DETAIL_DIR, 'adjudication_detail.tsv'), ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
