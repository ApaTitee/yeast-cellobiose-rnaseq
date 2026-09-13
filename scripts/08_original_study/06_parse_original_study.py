#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_parse_original_study.py
---------------------------------------------------------------------------
目的：把原研究补充材料（Dataset S1–S3）解析为结构化 TSV，作为一致性评估的对照基准。
      仅用 Python 标准库（zipfile + xml.etree），无需第三方依赖，可在任何环境下重跑。

输入：docs/literature/original_study/13068_2014_126_MOESM{1,2,3}_ESM.xlsx
输出：data/metadata/original_study/derived/
        dataset_s1_rpkm.tsv          6351 条（6349 宿主 + 2 外源）逐样本 RPKM + 复算比值
        dataset_s2_deg.tsv           519 条 DEG（含原文 reported 值与复算值）
        dataset_s2_tf.tsv            19 个转录因子
        dataset_s2_grn.tsv           7 个 GRN 调控基因
        dataset_s3_goslim.tsv        FunSpec GO 结果（UP 244 / DOWN 256）
        dataset_s3_notfound.tsv      FunSpec 未能映射的基因
        reconciliation.tsv           数字对账（519 = 244 UP + 256 DOWN + ? 等）

口径说明（关键）：
  * 原文 Dataset S2 的 `Fold Change` 列是**带符号倍数**：上调 = +(C8/G8)，下调 = -(G8/C8)；
    `log` 列是 log2(C8/G8)。本脚本保留原值，并在 *_recomputed 列给出复算结果以便交叉校验。
  * 基因标识：原文 Feature ID 为标准名（standard name）；Synonym 列含 systematic ORF name。
    本脚本不做名称映射（映射属 5.2 的 id_map_loss.tsv），只原样保留两列。
---------------------------------------------------------------------------
"""
from __future__ import annotations

import math
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
LIT = os.path.join(ROOT, "docs", "literature", "original_study")
OUT = os.path.join(ROOT, "data", "metadata", "original_study", "derived")

TRANSGENES = {"cdt-1EGFP", "gh1-1"}


# --------------------------------------------------------------------- xlsx
def load_sheet(path: str, sheet_index: int) -> list[dict[str, str]]:
    """把 xlsx 的某个 worksheet 读成 [{列字母: 值}]，值统一为字符串。"""
    with zipfile.ZipFile(path) as z:
        try:
            shared = [t.text or "" for t in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(NS + "t")]
        except KeyError:
            shared = []
        root = ET.fromstring(z.read(f"xl/worksheets/sheet{sheet_index}.xml"))
        rows = []
        for row in root.iter(NS + "row"):
            cells: dict[str, str] = {}
            for c in row.iter(NS + "c"):
                ref = c.get("r") or ""
                col = "".join(ch for ch in ref if ch.isalpha())
                v = c.find(NS + "v")
                if v is None:
                    cells[col] = ""
                elif c.get("t") == "s":
                    cells[col] = shared[int(v.text)]
                else:
                    cells[col] = v.text or ""
            rows.append(cells)
    return rows


def fnum(x: str):
    """安全转 float；空/非数值返回 None。"""
    try:
        if x is None or x == "":
            return None
        return float(x)
    except ValueError:
        return None


def fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def write_tsv(path: str, header: list[str], rows: list[list[str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    print(f"  [write] {os.path.relpath(path, ROOT)}  ({len(rows)} 行)")


# ------------------------------------------------------------------ S1 解析
def parse_s1() -> list[list[str]]:
    """Dataset S1：A=Feature ID, B=Database object name, C=Synonym,
    D/E/F=Cellobiose 三重复, G=C8 mean, H/I/J=Glucose 三重复, K=G8 mean。
    数据自第 5 行开始（第 1 行是作者遗留标签 'cdt-1'）。"""
    rows = load_sheet(os.path.join(LIT, "13068_2014_126_MOESM1_ESM.xlsx"), 1)
    out = []
    for i, r in enumerate(rows, start=1):
        name = (r.get("A") or "").strip()
        if not name or name in ("Feature ID",):
            continue
        c8 = [fnum(r.get(c, "")) for c in "DEF"]
        g8 = [fnum(r.get(c, "")) for c in "HIJ"]
        c8m, g8m = fnum(r.get("G", "")), fnum(r.get("K", ""))
        if not any(v is not None for v in c8 + g8):
            continue  # 表头行 / 作者遗留标签行（无任何数值）
        # 复算均值（原文 G/K 列为空或为四舍五入值时的交叉校验）
        c8m_rc = sum(c8) / 3 if all(v is not None for v in c8) else None
        g8m_rc = sum(g8) / 3 if all(v is not None for v in g8) else None
        ratio = (c8m / g8m) if (c8m and g8m and g8m > 0) else None
        log2fc = math.log2(ratio) if ratio and ratio > 0 else None
        out.append([
            name,
            (r.get("B") or "").strip(),
            (r.get("C") or "").strip(),
            fmt(c8[0]), fmt(c8[1]), fmt(c8[2]), fmt(c8m), fmt(c8m_rc),
            fmt(g8[0]), fmt(g8[1]), fmt(g8[2]), fmt(g8m), fmt(g8m_rc),
            fmt(ratio), fmt(log2fc),
            "yes" if name in TRANSGENES else "no",
        ])
    return out


# ------------------------------------------------------------------ S2 解析
def parse_s2_sheet(sheet_index: int, has_pathway: bool) -> tuple[list[list[str]], list[str]]:
    rows = load_sheet(os.path.join(LIT, "13068_2014_126_MOESM2_ESM.xlsx"), sheet_index)
    out = []
    for i, r in enumerate(rows, start=1):
        name = (r.get("A") or "").strip()
        if not name or name in ("Feature ID",):
            continue
        if has_pathway:  # A=ID, B=Pathway, C=FC, D=log, E=FDRp, F=C8mean, G=G8mean, H=desc, I=synonym
            fc, logv = fnum(r.get("C", "")), fnum(r.get("D", ""))
            p = r.get("E", "")
            c8m, g8m = fnum(r.get("F", "")), fnum(r.get("G", ""))
            desc, syn, pathway = r.get("H", ""), r.get("I", ""), r.get("B", "")
        else:  # A=ID, B=FC, C=log, D=FDRp, E=C8mean, F=G8mean, G=desc, H=synonym
            fc, logv = fnum(r.get("B", "")), fnum(r.get("C", ""))
            p = r.get("D", "")
            c8m, g8m = fnum(r.get("E", "")), fnum(r.get("F", ""))
            desc, syn, pathway = r.get("G", ""), r.get("H", ""), ""
        if fc is None and logv is None and c8m is None and g8m is None:
            continue  # 表头行 / 遗留标签行
        # 复算：原文 FC 为带符号倍数，log 为 log2(C8/G8)
        ratio = (c8m / g8m) if (c8m and g8m and g8m > 0) else None
        log2_rc = math.log2(ratio) if ratio and ratio > 0 else None
        fc_rc = None if ratio is None else (ratio if ratio >= 1 else -1.0 / ratio)
        out.append([
            name, pathway,
            fmt(fc), fmt(fc_rc), fmt(logv), fmt(log2_rc),
            p, fmt(c8m), fmt(g8m),
            "yes" if name in TRANSGENES else "no",
            desc.strip(), syn.strip(),
        ])
    return out, [r[0] for r in out]


# ------------------------------------------------------------------ S3 解析
def parse_s3():
    """Dataset S3：两个区块（UP / DOWN）。结构为：区块标题行（含该侧基因总数）+ GO 结果行
    + 'N genes were not found:' 计数行与其基因清单行。清单行不一定紧跟计数行
    （DOWN 侧的清单排在 GO 表之后），因此采用"等待消费"状态机。
    A=GO process, B=p(Bonferroni), C=k, D=f, E=In Category from Cluster。"""
    rows = load_sheet(os.path.join(LIT, "13068_2014_126_MOESM3_ESM.xlsx"), 1)
    go_rows: list[list[str]] = []
    notfound: list[list[str]] = []
    section_sizes: dict[str, int] = {}
    section = None
    pending: tuple[str, int] | None = None  # (方向, 计数)
    for r in rows:
        a = (r.get("A") or "").strip()
        if not a:
            continue
        m_sec = re.match(r"^GO process for\s+(\d+)\s+(UP|DOWN)\s+genes", a)
        if m_sec:
            section = m_sec.group(2)
            section_sizes[section] = int(m_sec.group(1))
            continue
        m_nf = re.match(r"^(\d+)\s+genes were not found", a)
        if m_nf:
            pending = (section or "NA", int(m_nf.group(1)))
            continue
        if "[GO:" in a:
            term = a[: a.rfind("[GO:")].strip()
            go_id = a[a.rfind("[GO:") + 1:].rstrip("]").strip()
            genes = (r.get("E") or "").strip()
            go_rows.append([
                section or "NA", term, go_id,
                (r.get("B") or "").strip(), (r.get("C") or "").strip(),
                (r.get("D") or "").strip(), genes,
                str(len(genes.split())),
            ])
            continue
        if pending is not None:  # 未映射基因清单行
            glist = [g.strip() for g in a.split(",") if g.strip()]
            notfound.append([pending[0], a, str(len(glist)), str(pending[1])])
            pending = None
    return go_rows, notfound, section_sizes


# --------------------------------------------------------------------- main
def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    print("解析原研究补充材料 ->", os.path.relpath(OUT, ROOT))

    s1 = parse_s1()
    write_tsv(os.path.join(OUT, "dataset_s1_rpkm.tsv"),
              ["feature_id", "sgd_description", "synonym",
               "C8_JCYL001D", "C8_JCYL002D", "C8_JCYL003D", "C8_mean", "C8_mean_recomputed",
               "G8_JCYL001B", "G8_JCYL002B", "G8_JCYL003B", "G8_mean", "G8_mean_recomputed",
               "ratio_C8_over_G8_recomputed", "log2FC_recomputed", "is_transgene"], s1)

    header_s2 = ["feature_id", "pathway", "FC_reported", "FC_recomputed",
                 "log2FC_reported", "log2FC_recomputed", "FDR_p_reported",
                 "C8_mean", "G8_mean", "is_transgene", "sgd_description", "synonym"]
    s2_deg, deg_ids = parse_s2_sheet(1, has_pathway=False)
    write_tsv(os.path.join(OUT, "dataset_s2_deg.tsv"), header_s2, s2_deg)
    s2_tf, tf_ids = parse_s2_sheet(2, has_pathway=False)
    write_tsv(os.path.join(OUT, "dataset_s2_tf.tsv"), header_s2, s2_tf)
    s2_grn, grn_ids = parse_s2_sheet(3, has_pathway=True)
    write_tsv(os.path.join(OUT, "dataset_s2_grn.tsv"),
              ["feature_id", "pathway", "FC_reported", "FC_recomputed",
               "log2FC_reported", "log2FC_recomputed", "FDR_p_reported",
               "C8_mean", "G8_mean", "is_transgene", "sgd_description", "synonym"], s2_grn)

    s3_go, s3_nf, s3_sizes = parse_s3()
    write_tsv(os.path.join(OUT, "dataset_s3_goslim.tsv"),
              ["direction", "go_term", "go_id", "p_bonferroni", "k_in_category",
               "f_category_size", "genes", "n_genes"], s3_go)
    write_tsv(os.path.join(OUT, "dataset_s3_notfound.tsv"),
              ["direction", "genes", "n_genes", "n_reported"], s3_nf)

    # ---- 数字对账 ----
    n_up_terms = sum(1 for r in s3_go if r[0] == "UP")
    n_dn_terms = sum(1 for r in s3_go if r[0] == "DOWN")
    nf_up = sum(int(r[3]) for r in s3_nf if r[0] == "UP")
    nf_dn = sum(int(r[3]) for r in s3_nf if r[0] == "DOWN")
    n_up = s3_sizes.get("UP", 0)
    n_dn = s3_sizes.get("DOWN", 0)
    n_trans_in_deg = sum(1 for r in s2_deg if r[9] == "yes")
    recon = [
        ["dataset_s1_rows_total", str(len(s1)), "Dataset S1 数据行数（含 2 条外源）"],
        ["dataset_s1_host_genes", str(len(s1) - sum(1 for r in s1 if r[15] == "yes")), "宿主基因条目数"],
        ["dataset_s1_transgenes", str(sum(1 for r in s1 if r[15] == "yes")), "外源条目数（cdt-1EGFP, gh1-1）"],
        ["dataset_s2_deg_rows", str(len(s2_deg)), "原文报告的 519 个 DEG"],
        ["dataset_s2_deg_transgenes", str(n_trans_in_deg), "519 中的外源条目数（gh1-1）"],
        ["dataset_s2_tf_rows", str(len(s2_tf)), "19 个转录因子"],
        ["dataset_s2_grn_rows", str(len(s2_grn)), "7 个 GRN 调控基因"],
        ["dataset_s3_go_terms_up", str(n_up_terms), "UP 侧 GO 条目数"],
        ["dataset_s3_go_terms_down", str(n_dn_terms), "DOWN 侧 GO 条目数"],
        ["dataset_s3_up_genes", str(n_up), "Dataset S3 区块标题声明的 UP 基因数"],
        ["dataset_s3_down_genes", str(n_dn), "Dataset S3 区块标题声明的 DOWN 基因数"],
        ["dataset_s3_notfound_up", str(nf_up), "FunSpec UP 未映射基因数"],
        ["dataset_s3_notfound_down", str(nf_dn), "FunSpec DOWN 未映射基因数"],
        ["check_sum_of_parts", str(n_up + n_dn + nf_up + nf_dn), "对账：UP + DOWN + 两侧 not found"],
        ["check_gap_vs_519", str(519 - (n_up + n_dn + nf_up + nf_dn)),
         "与 519 的差额；需在报告 5.12 解释（候选：外源条目 gh1-1 无 GO 注释）"],
    ]
    write_tsv(os.path.join(OUT, "reconciliation.tsv"), ["item", "value", "note"], recon)
    return 0


if __name__ == "__main__":
    sys.exit(main())
