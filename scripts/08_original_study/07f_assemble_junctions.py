#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07f_assemble_junctions.py — 用软剪切 reads 组装 cdt-1 与标签之间的接头序列

原理：跨越接头的 read 无法完整比对到任一候选，bowtie2（local）会把它软剪切。
      收集 (a) 比对到 cdt-1 且 3' 端软剪切的片段、(b) 比对到标签且 5' 端软剪切的片段，
      两者分别提供"接头左侧"与"接头右侧"的证据；再用 overlap-extension 从 cdt-1 3' 端
      向标签 5' 端延伸，直到抵达标签起点，并以跨接头 reads 数验证。

用法：
  python3 07f_assemble_junctions.py <bam> [<bam> ...]
输出：
  refs/custom/junctions.tsv                   （接头序列 + 支持证据）
  refs/custom/adjudication/junction_cdt1_tag.fa
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ADJ = os.path.join(ROOT, "refs", "custom", "adjudication")

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(COMP)[::-1]


def read_one(path: str, name: str) -> str:
    buf, keep = [], False
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if keep:
                break
            keep = line[1:].split()[0] == name
        elif keep:
            buf.append(line.strip())
    s = "".join(buf).upper()
    if not s:
        raise SystemExit(f"!! 参考中未找到 {name}")
    return s


def clips_of(bam: str, ref: str, min_mapq: int = 5):
    """返回该参考上的 (左软剪切序列列表, 右软剪切序列列表)，方向已统一为参考方向。"""
    left, right = [], []
    p = subprocess.Popen(["samtools", "view", "-F", "0x904", bam, ref],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:                       # type: ignore[union-attr]
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        try:
            if int(f[4]) < min_mapq:
                continue
        except ValueError:
            continue
        seq = f[9].upper()
        if int(f[1]) & 16:
            seq = rc(seq)
        ops = [(int(n), op) for n, op in CIGAR.findall(f[5])]
        i = 0
        if ops and ops[0][1] == "S":
            left.append(seq[:ops[0][0]])
            i = ops[0][0]
        if len(ops) > 1 and ops[-1][1] == "S":
            right.append(seq[len(seq) - ops[-1][0]:])
    p.wait()                                    # type: ignore[union-attr]
    return left, right


def assemble(anchor: str, tails: list[str], target_start: str, min_overlap: int = 18,
             min_support: int = 3, max_steps: int = 12):
    """从 anchor 的末端向 target 方向做 overlap-extension。"""
    contig = anchor[-60:]
    used, steps = set(), []
    for _ in range(max_steps):
        cands = defaultdict(list)
        for idx, t in enumerate(tails):
            if idx in used:
                continue
            for k in range(min(30, len(t) - 1), min_overlap - 1, -1):
                if contig.endswith(t[:k]):
                    cands[(k, t[k:])].append(idx)
                    break
        if not cands:
            break
        # 先取支持数 >= min_support 的候选中重叠最长、延伸最长者；否则取重叠最长者
        good = {c: v for c, v in cands.items() if len(v) >= min_support and c[1]}
        pool = good or {c: v for c, v in cands.items() if c[1]}
        if not pool:
            break
        (k, ext), idxs = max(pool.items(), key=lambda kv: (kv[0][0], len(kv[0][1])))
        contig += ext
        used.update(idxs)
        steps.append((k, ext, len(idxs)))
        if contig.endswith(target_start[:min(25, len(target_start))]):
            break
    reached = contig.endswith(target_start[:min(25, len(target_start))])
    return contig, steps, reached


def main() -> int:
    bams = sys.argv[1:]
    if not bams:
        raise SystemExit(__doc__)
    ref_fa = os.path.join(ADJ, "candidates.fa")
    cdt = read_one(ref_fa, "cdt-1_native_CDS")
    tag = read_one(ref_fa, "tag_eGFP_U55762")
    gh1 = read_one(ref_fa, "gh1-1_native_CDS")

    # --- cdt-1 -> tag 接头 ---
    left_tails, right_tails, ok_span = [], [], 0
    for b in bams:
        # cdt-1 上的 3' 软剪切（接头右侧证据）；限制在 cdt-1 末端附近，避免其他来源的剪切片段
        p = subprocess.Popen(["samtools", "view", "-F", "0x904", b, "cdt-1_native_CDS"],
                             stdout=subprocess.PIPE, text=True)
        for line in p.stdout:                   # type: ignore[union-attr]
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            try:
                if int(f[4]) < 5:
                    continue
            except ValueError:
                continue
            ops = [(int(n), op) for n, op in CIGAR.findall(f[5])]
            ref_consumed = sum(n for n, op in ops if op in "MDN=X")
            seq = f[9].upper()
            if int(f[1]) & 16:
                seq = rc(seq)
            if ops and ops[0][1] == "S":
                # 左剪切：其右侧序列应与 cdt-1 末端相连 -> 属于 cdt-1 之后的序列
                if int(f[3]) + ref_consumed >= 1700:
                    right_tails.append(seq[:ops[0][0]] + seq[ops[0][0]:])
            if len(ops) > 1 and ops[-1][1] == "S" and int(f[3]) + ref_consumed >= 1650:
                right_tails.append(seq[-(ops[-1][0] + 0):])
        p.wait()                                # type: ignore[union-attr]
        # 标签上的 5' 软剪切（接头左侧证据）
        l2, _ = clips_of(b, "tag_eGFP_U55762")
        left_tails.extend(l2)

    contig, steps, reached = assemble(cdt, right_tails + left_tails, tag)
    junction = contig[len(cdt[-60:]):] if contig.startswith(cdt[-60:]) else ""
    # 跨接头 reads 计数：在原始 BAM 中找同时包含 cdt-1 末端 20 bp 与标签起始 20 bp 的 reads
    k1, k2 = cdt[-20:], tag[:20]
    span = 0
    for b in bams:
        p = subprocess.Popen(["samtools", "view", "-F", "0x904", b], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:                   # type: ignore[union-attr]
            s = line.split("\t")[9].upper()
            if k1 in s or k1 in rc(s):
                if k2 in s or k2 in rc(s):
                    span += 1
        p.wait()                                # type: ignore[union-attr]

    rows = [["cdt-1->tag", contig, junction, str(reached), str(span),
             ";".join(f"k{k}+{len(e)}(n={n})" for k, e, n in steps)]]
    out_tsv = os.path.join(ROOT, "refs", "custom", "junctions.tsv")
    with open(out_tsv, "w", encoding="utf-8") as fh:
        fh.write("junction\tassembled_contig\tjunction_sequence\treached_target\tspanning_reads\tassembly_steps\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    fa = os.path.join(ADJ, "junction_cdt1_tag.fa")
    with open(fa, "w", encoding="utf-8") as fh:
        fh.write(">junction_cdt1_tag  (cdt-1 3'end + linker; assembled from soft-clipped reads)\n")
        fh.write(contig + "\n")

    print("=== cdt-1 -> tag 接头组装 ===")
    print(f"  右侧证据片段（cdt-1 端）: {len(right_tails)}  左侧证据片段（tag 端）: {len(left_tails)}")
    print(f"  组装步数: {len(steps)}  抵达标签起点: {reached}")
    print(f"  组装序列 ({len(contig)} bp): {contig}")
    print(f"  推定接头 (cdt-1 之后、标签之前, {len(junction)} bp): {junction}")
    print(f"  跨接头 reads（同时含 cdt-1 末 20 bp 与标签首 20 bp）: {span}")
    print(f"  写出: refs/custom/junctions.tsv, {os.path.relpath(fa, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
