#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07a_build_adjudication_candidates.py
---------------------------------------------------------------------------
目的：构建"外源序列裁决"的候选序列集（PLAN 3.4.3 第 1 步），并记录每个候选的来源与 md5。
      同时候选集中包含阴性对照（Figure S8 的 gh1-1a、质粒骨架、天然 2μ），
      使"样本中实际是哪个版本"成为**被测量的量**，而不是引用文献的推断。

输入：refs/diagnostic/*、docs/literature/original_study/13068_2014_126_MOESM4_ESM.pdf
输出：refs/custom/adjudication/candidates.fa
      refs/custom/adjudication/candidates.tsv     （候选名/长度/来源/md5/预期）
      refs/diagnostic/gh1-1a_FigureS8.fa          （从 PDF 提取，附长度与 md5）
---------------------------------------------------------------------------
"""
from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DIAG = os.path.join(ROOT, "refs", "diagnostic")
OUTDIR = os.path.join(ROOT, "refs", "custom", "adjudication")
PDF = os.path.join(ROOT, "docs", "literature", "original_study", "13068_2014_126_MOESM4_ESM.pdf")

CODON = {}
_b, _aas = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_k = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _aas[_k]; _k += 1
translate = lambda s: "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
md5 = lambda b: hashlib.md5(b).hexdigest()


def read_fasta(path: str) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def extract_figS8() -> str:
    """从补充材料 PDF 中提取 Figure S8 的 gh1-1a 序列。
    Figure S8 的序列块在版式上位于图题之前，故先取 'Figure S8' 之前的文本块再收集纯 ACGT 行。"""
    txt = subprocess.run(["pdftotext", "-layout", PDF, "-"], capture_output=True, text=True, check=True).stdout
    i = txt.find("Figure S8")
    if i < 0:
        raise SystemExit("!! 未在补充材料中找到 Figure S8")
    before = txt[max(0, i - 6000):i]
    seq = "".join(l.strip() for l in before.split("\n") if re.fullmatch(r"[ACGT]{20,}", l.strip()))
    if not seq:
        raise SystemExit("!! Figure S8 序列提取失败（版式可能变化）")
    if len(seq) % 3 != 0:
        raise SystemExit(f"!! Figure S8 序列长度 {len(seq)} 不是 3 的倍数")
    return seq


def main() -> int:
    os.makedirs(OUTDIR, exist_ok=True)
    cands: list[tuple[str, str, str, str]] = []  # (name, seq, source, expectation)

    # --- 1. Figure S8 的 gh1-1a（阴性对照：样本中不应为此版本）---
    s8 = extract_figS8()
    p_s8 = os.path.join(DIAG, "gh1-1a_FigureS8.fa")
    with open(p_s8, "w", encoding="utf-8") as fh:
        fh.write(">gh1-1a_FigureS8  codon-optimized gh1-1a with C-terminal His6 (DNA2.0)\n")
        for i in range(0, len(s8), 70):
            fh.write(s8[i:i + 70] + "\n")
    print(f"  Figure S8 提取: {len(s8)} bp -> {os.path.relpath(p_s8, ROOT)}")

    # --- 2. gh1-1 候选 ---
    gh_native = read_fasta(os.path.join(DIAG, "ncrassa_gh1-1_NCU00130_candidate_CDS.fa"))[0][1].upper()
    gh_orf = gh_native[:-3] if gh_native.endswith("TAA") else gh_native   # 去终止密码子
    cands += [
        ("gh1-1_native_CDS", gh_native, "XM_011395456.1 (NCU00130) CDS", "样本预期为真"),
        ("gh1-1_native_His6", gh_orf + "CAT" * 6 + "TAA",
         "XM_011395456.1 CDS + C 端 6xHis（Galazka 2010 SOM 引物推导）", "样本预期为真"),
        ("gh1-1a_FigureS8", s8, "本论文 Figure S8（DNA2.0 密码子优化）", "阴性对照：不应有 reads"),
    ]

    # --- 3. cdt-1 候选 ---
    cdt = read_fasta(os.path.join(DIAG, "ncrassa_cdt-1_NCU00801_candidate_CDS.fa"))[0][1].upper()
    cands.append(("cdt-1_native_CDS", cdt, "XM_958708.2 (NCU00801) CDS", "样本预期为真（eGFP 融合体的一部分）"))

    # --- 4. 标签候选：eGFP 与 3 个 sfGFP 记录 ---
    egfp_recs = read_fasta(os.path.join(DIAG, "egfp_candidates_U55762.fa"))
    egfp = next((s.upper() for _, s in egfp_recs if translate(s.upper()).startswith("MVSKGEELFTGV")), None)
    if egfp is None:
        raise SystemExit("!! 未能在 U55762 中找到 EGFP CDS（翻译起点匹配失败）")
    cands.append(("tag_eGFP_U55762", egfp, "U55762.1 (pEGFP-N1) egfp CDS", "标签候选"))
    for rec in ("PX636966", "MW132720", "JQ341914"):
        p = os.path.join(DIAG, f"sfgfp_candidate_{rec}.fa")
        s = read_fasta(p)[0][1].upper()
        cands.append((f"tag_sfGFP_{rec}", s, f"{rec} superfolder GFP CDS", "标签候选"))

    # --- 5. 阴性/背景对照 ---
    for fn, nm, src in [("plasmid_prs426_U03451.fa", "CTRL_pRS426_backbone", "U03451.1 (pRS426)"),
                        ("plasmid_2micron_NC001398.fa", "CTRL_native_2micron", "NC_001398.1 (2μ circle)")]:
        s = read_fasta(os.path.join(DIAG, fn))[0][1].upper()
        cands.append((nm, s, src, "诊断对照：用于 read fate 归类"))

    # --- 输出 ---
    fa = os.path.join(OUTDIR, "candidates.fa")
    with open(fa, "w", encoding="utf-8") as fh:
        for name, seq, _, _ in cands:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    tsv = os.path.join(OUTDIR, "candidates.tsv")
    with open(tsv, "w", encoding="utf-8") as fh:
        fh.write("candidate\tlength_bp\tsource\texpectation\tmd5\n")
        for name, seq, src, exp in cands:
            fh.write(f"{name}\t{len(seq)}\t{src}\t{exp}\t{md5(seq.encode())}\n")

    print(f"  候选集: {len(cands)} 条 -> {os.path.relpath(fa, ROOT)}")
    for name, seq, _, exp in cands:
        print(f"    {name:24s} {len(seq):6d} bp   {exp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
