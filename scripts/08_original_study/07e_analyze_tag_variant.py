#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07e_analyze_tag_variant.py — 把 reads 共识与公共 eGFP 的差异翻译到氨基酸层面，
并与常见 GFP 变体的突变集合对照，判定标签的真实身份。

用法： python3 07e_analyze_tag_variant.py <consensus_of_tag.fa> <egfp_candidate.fa>
输出： refs/custom/adjudication/consensus/tag_variant_analysis.tsv
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

CODON = {}
_b, _aas = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_k = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _aas[_k]; _k += 1
tr = lambda s: "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))

# 常见 GFP 变体相对 avGFP/EGFP 的替换（用于身份比对；序列来源见附录 A 文献）
KNOWN = {
    "EGFP (GFPmut1 + humanized)": {"F64L", "S65T"},
    "GFPmut2": {"F64L", "S65T", "V68L", "S72A"},
    "cycle3": {"F99S", "M153T", "V163A"},
    "GFPmut3": {"S65G", "S72A", "T203Y"},
    "sfGFP (superfolder)": {"F64L", "S65T", "F99S", "M153T", "V163A", "S30R", "Y39N", "N105T", "Y145F", "I171V", "A206V"},
}


def read_one(path: str) -> str:
    buf, keep = [], False
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            keep = True
        elif keep:
            buf.append(line.strip())
    return "".join(buf).upper()


def find_frame(nt_query: str, nt_ref: str) -> int:
    """在 ref 的 3 个翻译帧中寻找与 query 蛋白最相似者，返回帧号。"""
    q = tr(nt_query)
    best, best_frame = -1, 0
    for f in range(3):
        r = tr(nt_ref[f:])
        n = min(len(q), len(r))
        same = sum(1 for i in range(n) if q[i] == r[i])
        if same > best:
            best, best_frame = same, f
    return best_frame


def main() -> int:
    cons_path, ref_path = sys.argv[1], sys.argv[2]
    cons, ref = read_one(cons_path), read_one(ref_path)
    frame = find_frame(cons, ref)
    r = ref[frame:]
    q = cons
    p_ref, p_con = tr(r), tr(q)

    # 逐密码子比较（忽略 N）
    changes = []
    for i in range(min(len(r), len(q)) // 3):
        rc, qc = r[i * 3:i * 3 + 3], q[i * 3:i * 3 + 3]
        if "N" in qc or rc == qc:
            continue
        a, b = p_ref[i], p_con[i]
        if a != b:
            changes.append((i + 1, a, b, rc, qc))

    subs = {f"{a}{pos}{b}" for pos, a, b, _, _ in changes}
    out = os.path.join(os.path.dirname(cons_path), "tag_variant_analysis.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("codon_index\teGFP_aa\tconsensus_aa\teGFP_codon\tconsensus_codon\n")
        for row in changes:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"共识 {os.path.basename(cons_path)} vs 参考 {os.path.basename(ref_path)}（参考帧偏移 {frame}）")
    print(f"  核酸差异（不含 N）: {sum(1 for a, b in zip(r, q) if a != b and b != 'N')}")
    print(f"  氨基酸替换（{len(changes)} 处）:")
    for pos, a, b, rc, qc in changes:
        print(f"    {a}{pos}{b}   ({rc} -> {qc})")
    print("\n  与已知变体的替换集合对照:")
    for name, known in KNOWN.items():
        inter = subs & known
        print(f"    {name:28s} 匹配 {len(inter):2d}/{len(subs):2d} -> {sorted(inter)}")
    print(f"\n写入 {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
