#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07d_compare_consensus.py — 把 reads 共识序列与公共候选逐一比较（蛋白/核酸同一性）

用法： python3 07d_compare_consensus.py <consensus.fa> <candidates.fa> [输出tsv]
输出： refs/custom/adjudication/consensus/tag_identity.tsv
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
translate = lambda s: "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def read_fasta(path: str) -> dict[str, str]:
    out, cur, buf = {}, None, []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if cur:
                out[cur] = "".join(buf).upper()
            cur, buf = line[1:].split()[0], []
        elif cur:
            buf.append(line.strip())
    if cur:
        out[cur] = "".join(buf).upper()
    return out


def main() -> int:
    cons_path, cand_path = sys.argv[1], sys.argv[2]
    out_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(
        os.path.dirname(cons_path), "tag_identity.tsv")
    cons = read_fasta(cons_path)
    cands = read_fasta(cand_path)
    cname = next(iter(cons))
    cseq = cons[cname]
    pcons = translate(cseq)

    rows = []
    for nm, s in cands.items():
        p = translate(s)
        n = min(len(p), len(pcons))
        prot_id = sum(1 for i in range(n) if p[i] == pcons[i]) / n if n else 0.0
        m = min(len(s), len(cseq))
        nt_id = sum(1 for i in range(m) if s[i] == cseq[i]) / m if m else 0.0
        rows.append((nm, len(s), len(p), f"{prot_id:.4f}", f"{nt_id:.4f}",
                     sum(1 for i in range(m) if s[i] != cseq[i])))

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("candidate\tnt_len\taa_len\tprotein_identity\tnt_identity\tnt_diff\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"共识序列 {cname} ({len(cseq)} bp, 翻译 {len(pcons)} aa)")
    print(f"{'候选':26s} {'蛋白同一性':>10s} {'核酸同一性':>10s} {'核酸差异':>8s}")
    for nm, ln, lp, pid, nid, nd in rows:
        print(f"{nm:26s} {pid:>10s} {nid:>10s} {nd:>8d}")
    print(f"\n写入 {os.path.relpath(out_path, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
