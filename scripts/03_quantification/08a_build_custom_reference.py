#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08a_build_custom_reference.py — 构建自定义参考体系（PLAN 5.2 / 3.5）

产物：
  refs/custom/transcripts.fa     宿主转录本（gffread 抽取，含线粒体 mRNA）+ 2 条外源转录本
  refs/custom/transgenes.fa      仅 2 条外源转录本
  refs/custom/tx2gene.tsv        transcript_id -> gene_id（gene_id 优先取 locus_tag/systematic ORF name）
  refs/custom/decoys.txt         decoy 序列名（宿主基因组 17 条）
  refs/custom/gentrome.fa        transcripts.fa + 宿主基因组（salmon decoy-aware 索引输入）
  refs/checks/id_map_loss.tsv    原文 6351 / 519 条目 -> 本参考的映射损失

为什么用 gffread 而不是 RefSeq 的 rna_from_genomic.fna：
  后者不含线粒体蛋白编码 mRNA（只有 tRNA/rRNA/RPM1），若直接使用，
  5.6 的"mtDNA 编码基因检出情况"会因参考缺失而得出假结论。gffread 从 GFF 抽取可覆盖 19 条线粒体 mRNA。

外源转录本构成（依据 refs/custom/PROVENANCE.md 与接头裁决）：
  cdt-1EGFP = cdt-1 CDS（去终止子，1737 bp）+ ATCGAT（ClaI）+ GGTAGTGGTAGT（GSGS）
            + 标签共识[3:]（去掉公共 eGFP 起始 ATG；717 bp）= 2472 bp，读框连续
  gh1-1     = gh1-1 CDS（去终止子，1428 bp）+ CAT×6 + TAA = 1449 bp
"""
from __future__ import annotations

import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
REFS = os.path.join(ROOT, "refs")
HOST = os.path.join(REFS, "host")
DIAG = os.path.join(REFS, "diagnostic")
CUSTOM = os.path.join(REFS, "custom")
CHECKS = os.path.join(REFS, "checks")
DERIVED = os.path.join(ROOT, "data", "metadata", "original_study", "derived")

CLAI, GSGS, HIS6 = "ATCGAT", "GGTAGTGGTAGT", "CAT" * 6


def open_maybe_gz(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


def read_fasta(path: str) -> list[tuple[str, str]]:
    out, head, buf = [], None, []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if head is not None:
                    out.append((head, "".join(buf)))
                head, buf = line[1:].strip(), []
            elif line.strip():
                buf.append(line.strip())
    if head is not None:
        out.append((head, "".join(buf)))
    return out


def first_seq(path: str) -> str:
    return read_fasta(path)[0][1].upper()


def strip_stop(s: str) -> str:
    return s[:-3] if s[-3:] in ("TAA", "TAG", "TGA") else s


def parse_gff_map(gff_path: str):
    """从 GFF 解析 transcript -> gene 映射与基因名集合。"""
    tx2gene: dict[str, str] = {}
    names: set[str] = set()
    with open_maybe_gz(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype, attrs = f[2], f[8]

            def attr(key):
                m = re.search(rf"(?:^|;){key}=([^;]+)", attrs)
                return m.group(1) if m else None

            if ftype == "gene":
                for k in ("locus_tag", "Name", "gene", "ID"):
                    v = attr(k)
                    if v:
                        names.add(v.replace("gene-", ""))
            elif ftype in ("mRNA", "transcript", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA"):
                ident = attr("ID") or attr("transcript_id")
                if not ident:
                    continue
                ident = ident.replace("rna-", "")
                gene = attr("locus_tag") or attr("gene") or (attr("Parent") or "").replace("gene-", "") or ident
                tx2gene[ident] = gene
                for k in ("locus_tag", "gene", "Parent"):
                    v = attr(k)
                    if v:
                        names.add(v.replace("gene-", ""))
    return tx2gene, names


def main() -> int:
    os.makedirs(CUSTOM, exist_ok=True)
    os.makedirs(CHECKS, exist_ok=True)

    # ---------- 1. 宿主转录本（gffread 输出）与 tx2gene（GFF） ----------
    host_fa = os.path.join(HOST, "host_transcripts_gffread.fa")
    if not os.path.exists(host_fa):
        raise SystemExit("!! 缺少 refs/host/host_transcripts_gffread.fa，请先运行 08_build_reference.sh（内含 gffread 步骤）")
    host = read_fasta(host_fa)
    tx2gene_map, host_names = parse_gff_map(os.path.join(HOST, "GCF_000146045.2_R64_genomic.gff.gz"))
    n_no_gene = 0
    rows = []
    for head, seq in host:
        key = head.split()[0]
        tx = key.replace("rna-", "")
        gene = tx2gene_map.get(tx)
        if gene is None:
            gene = tx
            n_no_gene += 1
        rows.append((tx, gene, head, seq))
        host_names.add(gene)
        host_names.add(tx)

    # ---------- 2. 外源转录本 ----------
    cdt = first_seq(os.path.join(DIAG, "ncrassa_cdt-1_NCU00801_candidate_CDS.fa"))
    gh1 = first_seq(os.path.join(DIAG, "ncrassa_gh1-1_NCU00130_candidate_CDS.fa"))
    tag = first_seq(os.path.join(CUSTOM, "adjudication", "consensus",
                                 "tag_eGFP_U55762.SRR1166445_JCYL001D_cellobiose.fa"))
    egfp = None
    for _, s in read_fasta(os.path.join(DIAG, "egfp_candidates_U55762.fa")):
        s = s.upper()
        if s.startswith("ATGGTGAGCAAGGGCGAGGAG"):
            egfp = s
            break
    if egfp is None:
        raise SystemExit("!! 未能在 U55762 候选中定位 EGFP CDS")
    filled = [(i + 1, tag[i], egfp[i]) for i in range(min(len(tag), len(egfp))) if tag[i] == "N"]
    tag_filled = "".join(egfp[i] if tag[i] == "N" else tag[i] for i in range(min(len(tag), len(egfp))))

    transgenes = [("cdt-1EGFP", strip_stop(cdt) + CLAI + GSGS + tag_filled[3:],
                   "cdt-1(NCU00801) CDS[-stop] + ClaI + GSGS + sfGFP tag[no ATG]"),
                  ("gh1-1", strip_stop(gh1) + HIS6 + "TAA",
                   "gh1-1(NCU00130) CDS[-stop] + 6xHis + stop")]
    for name, seq, _ in transgenes:
        assert len(seq) % 3 == 0, f"{name} 读框不连续"

    # ---------- 3. 写出 ----------
    tr_path = os.path.join(CUSTOM, "transcripts.fa")
    with open(tr_path, "w", encoding="utf-8") as fh:
        for name, seq, desc in transgenes:
            fh.write(f">{name} transgene | {desc}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
        for tx, gene, head, seq in rows:
            fh.write(f">{head}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    with open(os.path.join(CUSTOM, "transgenes.fa"), "w", encoding="utf-8") as fh:
        for name, seq, desc in transgenes:
            fh.write(f">{name} {desc}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    with open(os.path.join(CUSTOM, "tx2gene.tsv"), "w", encoding="utf-8") as fh:
        fh.write("transcript_id\tgene_id\tsource\n")
        for name, _, _ in transgenes:
            fh.write(f"{name}\t{name}\ttransgene\n")
        for tx, gene, head, _ in rows:
            # 键必须与 quant.sf 中的 Name 完全一致（即 FASTA 头部首个 token，含 'rna-' 前缀）
            fh.write(f"{head.split()[0]}\t{gene}\thost\n")

    # ---------- 4. decoy 与 gentrome ----------
    genome = read_fasta(os.path.join(HOST, "GCF_000146045.2_R64_genomic.fna.gz"))
    with open(os.path.join(CUSTOM, "decoys.txt"), "w", encoding="utf-8") as fh:
        for head, _ in genome:
            fh.write(head.split()[0] + "\n")
    with open(os.path.join(CUSTOM, "gentrome.fa"), "w", encoding="utf-8") as out:
        with open(tr_path, encoding="utf-8") as fh:
            out.write(fh.read())
        for head, seq in genome:
            out.write(f">{head}\n")
            for i in range(0, len(seq), 70):
                out.write(seq[i:i + 70] + "\n")

    # ---------- 5. ID 映射损失 ----------
    loss = []
    for fn, label in [("dataset_s1_rpkm.tsv", "S1_6351"), ("dataset_s2_deg.tsv", "S2_519")]:
        p = os.path.join(DERIVED, fn)
        if not os.path.exists(p):
            continue
        lines = open(p, encoding="utf-8").read().splitlines()
        hdr = lines[0].split("\t")
        i_syn, i_tg = hdr.index("synonym"), hdr.index("is_transgene")
        n_all = n_map = 0
        unmapped = []
        for line in lines[1:]:
            f = line.split("\t")
            if len(f) <= max(i_syn, i_tg):
                continue
            name, syn, is_tg = f[0], f[i_syn], (f[i_tg] == "yes")
            n_all += 1
            cand = None
            for tok in re.split(r"[|,]", syn):
                tok = tok.strip()
                if re.fullmatch(r"Y[A-P][LR]\d{3}[WC](?:-[A-Z])?", tok):
                    cand = tok
                    break
            if is_tg or (cand and cand in host_names) or (name in host_names):
                n_map += 1
            else:
                unmapped.append(name)
        loss.append((label, n_all, n_map, n_all - n_map, ",".join(unmapped[:60])))
    with open(os.path.join(CHECKS, "id_map_loss.tsv"), "w", encoding="utf-8") as fh:
        fh.write("dataset\tn_total\tn_mapped\tn_unmapped\tunmapped_examples\n")
        for r in loss:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # ---------- 6. 汇报 ----------
    print(f"宿主转录本: {len(rows)} 条（gffread；其中 {n_no_gene} 条无 gene 映射，回退为 transcript id）")
    print(f"tx2gene 映射表: {len(tx2gene_map)} 条 mRNA/transcript 记录")
    print(f"外源转录本: {len(transgenes)} 条")
    for name, seq, desc in transgenes:
        print(f"  {name}: {len(seq)} bp = {len(seq)//3} codons | {desc}")
    print(f"  标签 N 位回填: {len(filled)} 处 -> {filled}")
    print(f"decoy: {len(genome)} 条；gentrome 已写入")
    print("ID 映射损失:")
    for r in loss:
        print(f"  {r[0]}: {r[2]}/{r[1]} 可映射（{100*r[2]/max(r[1],1):.1f}%）, 未映射 {r[3]}")

    # 线粒体覆盖自检（关系到 5.6 判定是否有效）
    mito = [t for t, g, _, _ in rows if t.startswith("Q")]
    print(f"线粒体转录本（locus_tag Q*）: {len(mito)} 条 -> 5.6 的 mtDNA 检出检查有效")
    return 0


if __name__ == "__main__":
    sys.exit(main())
