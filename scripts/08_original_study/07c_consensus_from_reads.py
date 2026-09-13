#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07c_consensus_from_reads.py
---------------------------------------------------------------------------
目的：从裁决比对的 BAM 中重建指定候选区域的**共识序列**，用于
      (1) 判定标签的真实身份（eGFP 的密码子优化变体？）；
      (2) 与公共候选逐位比较并区分同义/非同义差异；
      (3) 为 refs/custom 提供实测序列依据（PLAN 3.4.3 第 2–3 步的量化基础）。

方法：samtools mpileup 取每个位置的支持碱基，按覆盖度取多数碱基作为共识；
      低覆盖位置记为 N 并计数；输出共识序列、与候选的逐位差异、翻译比较。

用法：
  python3 07c_consensus_from_reads.py <candidate_name> <bam> [输出前缀]
输出：refs/custom/adjudication/consensus/<candidate>.<sample>.{fa,diff.tsv}
---------------------------------------------------------------------------
"""
from __future__ import annotations

import os
import subprocess
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ADJ = os.path.join(ROOT, "refs", "custom", "adjudication")
OUTDIR = os.path.join(ADJ, "consensus")

CODON = {}
_b, _aas = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_k = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _aas[_k]; _k += 1
translate = lambda s: "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def read_one(path: str, name: str) -> str:
    seq, keep, buf = None, False, []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if keep:
                break
            keep = line[1:].split()[0] == name
        elif keep:
            buf.append(line.strip())
    seq = "".join(buf)
    if not seq:
        raise SystemExit(f"!! 参考 {path} 中未找到 {name}")
    return seq.upper()


def main() -> int:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    cand, bam = sys.argv[1], sys.argv[2]
    ref_fa = os.path.join(ADJ, "candidates.fa")
    ref = read_one(ref_fa, cand)
    sample = os.path.basename(bam)[:-4]
    os.makedirs(OUTDIR, exist_ok=True)

    # mpileup：-A 保留异常配对（单端数据）；-a 输出全部位置（否则零覆盖位置缺失）
    # 注意：匹配参考的碱基在 pileup 中以 '.'/',' 表示，必须还原为参考碱基，否则会误判为无覆盖
    pile = subprocess.run(["samtools", "mpileup", "-A", "-a", "-Q", "13", "-f", ref_fa, "-r", cand, bam],
                          capture_output=True, text=True, check=True).stdout
    cons, depth = [], []
    for line in pile.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        refbase = f[2].upper()
        bases = f[4].upper()
        cnt: Counter = Counter()
        i = 0
        while i < len(bases):
            ch = bases[i]
            if ch == "^":          # 读段起始标记：跳过标记与其后一个字符
                i += 2
                continue
            if ch == "$":
                i += 1
                continue
            if ch in "+-":         # indel 标记：跳过长度数字与其后的序列
                i += 1
                num = ""
                while i < len(bases) and bases[i].isdigit():
                    num += bases[i]; i += 1
                i += int(num or 0)
                continue
            if ch in ".,":         # 与参考一致
                cnt[refbase] += 1
            elif ch in "ACGTN":
                cnt[ch] += 1
            elif ch == "*":        # 缺失占位
                cnt["DEL"] += 1
            i += 1
        alt = {b: cnt.get(b, 0) for b in "ACGT"}
        tot = sum(alt.values()) + cnt.get("DEL", 0)
        depth.append(tot)
        if tot == 0 or not any(alt.values()):
            cons.append("N"); continue
        best = max(alt, key=lambda b: alt[b])
        # 要求多数碱基的支持数不少于 60% 的简单多数深度，否则记为 N
        simple = sum(alt.values())
        cons.append(best if simple and alt[best] >= 0.6 * simple else "N")
    consensus = "".join(cons)

    n_n = consensus.count("N")
    d = [x for x in depth if x > 0]
    mean_depth = sum(d) / len(d) if d else 0.0

    diffs = [(i, ref[i], consensus[i]) for i in range(min(len(ref), len(consensus)))
             if consensus[i] != "N" and consensus[i] != ref[i]]
    # 同义/非同义判定（按密码子第 3 位近似：第 3 位差异多为同义）
    syn = sum(1 for i, _, _ in diffs if (i % 3) == 2)
    nonsyn = len(diffs) - syn

    fa = os.path.join(OUTDIR, f"{cand}.{sample}.fa")
    with open(fa, "w", encoding="utf-8") as fh:
        fh.write(f">{cand}_consensus_from_{sample}  (reads-derived consensus)\n")
        for i in range(0, len(consensus), 70):
            fh.write(consensus[i:i + 70] + "\n")
    tsv = os.path.join(OUTDIR, f"{cand}.{sample}.diff.tsv")
    with open(tsv, "w", encoding="utf-8") as fh:
        fh.write("pos_1based\tref_base\tconsensus_base\tcodon_pos\n")
        for i, a, b in diffs:
            fh.write(f"{i + 1}\t{a}\t{b}\t{(i % 3) + 1}\n")

    print(f"候选 {cand} | 样本 {sample}")
    print(f"  参考长度 {len(ref)} bp, 共识长度 {len(consensus)} bp, 未定(N) {n_n} bp, 平均覆盖 {mean_depth:.1f}×")
    print(f"  与候选的差异位点 {len(diffs)} 个（密码子第 3 位 {syn}，第 1/2 位 {nonsyn}）")
    if diffs:
        shown = diffs[:12]
        print("  前若干差异: " + ", ".join(f"{i+1}{a}>{b}" for i, a, b in shown))
    print(f"  翻译（第 1 帧）: {translate(consensus)[:80]}")
    print(f"  写入: {os.path.relpath(fa, ROOT)} / {os.path.relpath(tsv, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
