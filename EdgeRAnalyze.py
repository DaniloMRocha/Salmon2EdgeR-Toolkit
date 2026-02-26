#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
edgeR Analyzer (GUI) — A vs B (padronizado)

Padronização
-----------
- edgeR logFC esperado: log2(A/B)
- Direção:
    logFC > 0  => A-up
    logFC < 0  => B-up
- Para gráficos divergentes por categoria:
    plotFC = -logFC  (A-up fica negativo; B-up fica positivo)

Requisitos
----------
  pip install pandas numpy matplotlib

Rodar
-----
  python3 edger_analyzer_gui_AB.py
"""

from __future__ import annotations

import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk


# ----------------------------
# Helpers: leitura e colunas
# ----------------------------

def _read_tabular(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", dtype=str, encoding="utf-8")
    except UnicodeDecodeError:
        return pd.read_csv(path, sep="\t", dtype=str, encoding="latin-1")


def _detect_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def _as_numeric_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def sanitize_basename(name: str) -> str:
    name = name.strip()
    if not name:
        return ""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


# ----------------------------
# Filtro de descrições genéricas (opcional)
# ----------------------------

GENERIC_DESC_PATTERNS = [
    r"\buncharacterized\b",
    r"\bhypothetical\b",
    r"\bpredicted\b",
    r"\bunknown function\b",
    r"\bprotein of unknown function\b",
    r"\bunnamed protein product\b",
    r"\bputative uncharacterized\b",
    r"\bconserved unknown\b",
    r"\bno functional annotation\b",
]

def is_generic_description(desc: str) -> bool:
    if desc is None:
        return True
    d = str(desc).strip().lower()
    if d == "" or d == "nan":
        return True
    for p in GENERIC_DESC_PATTERNS:
        if re.search(p, d):
            return True
    if d.startswith("uncharacterized"):
        return True
    return False


def filter_informative_only(df: pd.DataFrame) -> pd.DataFrame:
    """Remove genes cuja descrição é genérica."""
    if "Description" not in df.columns:
        return df.copy()
    mask = ~df["Description"].apply(is_generic_description)
    return df.loc[mask].copy()


# ----------------------------
# Categorias (regex)
# ----------------------------

CATEGORY_RULES = [
    ("HSP_HSF_Heat", [
        r"\bHSP\b", r"\bHSF\b", r"heat", r"shock", r"chaperon", r"Hsc70", r"HSP70", r"HSP90", r"DNAJ", r"BiP"
    ]),
    ("UPR_ER_Stress", [
        r"\bUPR\b", r"unfolded protein", r"endoplasmic reticulum", r"\bER\b",
        r"IRE1", r"bZIP60", r"calnexin", r"calreticulin", r"protein disulfide isomerase", r"\bPDI\b"
    ]),
    ("ROS_Antioxidant", [
        r"peroxidase", r"catalase", r"superoxide dismutase", r"\bSOD\b",
        r"glutathione", r"ascorbate", r"thioredoxin", r"peroxiredoxin", r"oxidative", r"\bAPX\b", r"ascorbate peroxidase"
    ]),
    ("Metabolism", [
        r"glycolysis", r"gluconeogenesis", r"tricarboxylic", r"\bTCA\b", r"citric acid cycle",
        r"respiration", r"mitochondrial", r"oxidative phosphorylation", r"\bATP synthase\b",
        r"malate dehydrogenase", r"isocitrate dehydrogenase", r"succinate dehydrogenase",
        r"pyruvate kinase", r"hexokinase", r"phosphofructokinase",
        r"shikimate", r"phenylpropanoid", r"flavonoid", r"terpen", r"carotenoid",
        r"fatty acid", r"lipid", r"beta-oxidation"
    ]),
    ("Ubiquitin_Proteasome", [
        r"ubiquitin", r"\bE3\b", r"F-box", r"RING", r"BTB", r"proteasome", r"SCF", r"ubiquitin ligase"
    ]),
    ("Cytoskeleton_Polarity_Trafficking", [
        r"microtubule", r"tubulin", r"\bTUB\b", r"kinesin", r"dynein", r"myosin", r"MAP65", r"katanin",
        r"\bactin\b", r"\bACTIN\b",
        r"polarity", r"\bpolar\b", r"\bROP\b", r"\bRAC\b", r"\bCDC42\b", r"SOSEKI", r"DUF966",
        r"\bRAB\b", r"\bARF\b", r"SNARE", r"exocyst", r"clathrin", r"endocyt", r"vesicle", r"traffick"
    ]),
    ("CellWall_Callose_Pectin", [
        r"callose", r"\bGSL\b", r"glucan", r"beta-1,3", r"cell wall", r"pectin",
        r"expansin", r"xyloglucan", r"cellulose", r"lignin", r"extensin"
    ]),
    ("Hormone_Development", [
        r"auxin", r"\bIAA\b", r"gibberellin", r"\bGA\b", r"cytokinin",
        r"abscisic", r"\bABA\b", r"ethylene", r"\bACC\b", r"jasmon", r"salicylic",
        r"brassin", r"floral", r"pollen", r"anther", r"tapetum", r"microspore"
    ]),
    ("TF_Chromatin", [
        r"\bMYB\b", r"\bbHLH\b", r"\bNAC\b", r"\bWRKY\b", r"\bMADS\b", r"\bbZIP\b",
        r"chromatin", r"histone", r"SWI/SNF", r"methyltransferase", r"acetyltransferase", r"deacetylase"
    ]),
]

def categorize_row(symbol: str, desc: str) -> list[str]:
    s = ((symbol or "") + " " + (desc or "")).lower()
    cats = []
    for cat, patterns in CATEGORY_RULES:
        for p in patterns:
            if re.search(p.lower(), s):
                cats.append(cat)
                break
    if not cats:
        cats = ["Other"]
    return cats


# ----------------------------
# DETECÇÃO ROBUSTA DAS AMOSTRAS (3, 4 ou 5; A e B com mesmo n)
# ----------------------------

def detect_and_normalize_sample_cols(counts_df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Mantém exatamente o comportamento para 3 replicatas e expande para 4/5 quando existir.
    Regras:
      - aceita n = 3, 4 ou 5
      - A e B DEVEM ter o mesmo n
      - aceita colunas A1..An,B1..Bn OU Counts_A1..Counts_Bn (case-insensitive)
    Retorna:
      - df com colunas renomeadas para A1..An,B1..Bn
      - sample_cols na ordem A1..An,B1..Bn
    """
    cols = list(counts_df.columns)

    def wanted_for(n: int) -> list[str]:
        return [f"A{i}" for i in range(1, n + 1)] + [f"B{i}" for i in range(1, n + 1)]

    # 1) Caso já esteja normalizado (A1.. / B1..) — prioriza 5, depois 4, depois 3
    for n in (5, 4, 3):
        wanted = wanted_for(n)
        if all(c in cols for c in wanted):
            return counts_df.copy(), wanted

    # 2) Caso esteja com prefixo Counts_ (Counts_A1..Counts_Bn)
    for n in (5, 4, 3):
        wanted = wanted_for(n)
        mapping = {}
        for w in wanted:
            for pref in ("Counts_", "counts_"):
                key = f"{pref}{w}"
                if key in cols:
                    mapping[key] = w
        if len(mapping) == len(wanted):
            df = counts_df.copy().rename(columns=mapping)
            return df, wanted

    # 3) Caso misto / variações de caixa: regex ^(?:Counts_)?([AB][1-5])$
    regex = re.compile(r"^(?:Counts_)?([AB][1-5])$", flags=re.IGNORECASE)
    found = {}
    for c in cols:
        m = regex.match(c)
        if m:
            found[c] = m.group(1).upper()

    # escolher o maior n em {5,4,3} que tenha conjunto completo e pareado
    for n in (5, 4, 3):
        wanted = set(wanted_for(n))
        if set(found.values()) == wanted:
            df = counts_df.copy().rename(columns=found)
            return df, wanted_for(n)

    raise ValueError(
        "Counts: não consegui identificar colunas para A vs B com n=3,4 ou 5 replicatas.\n"
        f"Colunas encontradas: {cols}\n\n"
        "Garanta colunas exatamente:\n"
        "  - 'A1 A2 A3 B1 B2 B3' (ou até A5/B5),\n"
        "OU\n"
        "  - 'Counts_A1 ... Counts_B3' (ou até Counts_A5/Counts_B5).\n"
        "E lembre: A e B devem ter o mesmo número de replicatas (3, 4 ou 5)."
    )


# ----------------------------
# CPM / logCPM / z-score
# ----------------------------

def compute_log2cpm_from_counts(counts_df: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    X = counts_df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(float)
    libsize = X.sum(axis=0).replace(0, np.nan)
    cpm = (X / libsize) * 1e6
    log2cpm = np.log2(cpm + 1.0)
    out = counts_df[["GeneID"]].copy()
    out[sample_cols] = log2cpm
    return out


def zscore_rows(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    X = df[cols].astype(float).values
    mu = np.nanmean(X, axis=1, keepdims=True)
    sd = np.nanstd(X, axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    out = df.copy()
    out[cols] = Z
    return out


# ----------------------------
# Labels: SEMPRE descrição (fallback símbolo)
# ----------------------------

def _build_desc_labels(desc_series: pd.Series, symbol_series: pd.Series | None = None,
                      max_chars: int = 140) -> pd.Series:
    dsc = desc_series.fillna("").astype(str).replace("nan", "").str.strip()
    sym = symbol_series.fillna("").astype(str).replace("nan", "").str.strip() if symbol_series is not None else pd.Series([""] * len(dsc))
    labels = dsc.copy()
    empty = labels.eq("") | labels.str.lower().eq("nan")
    labels.loc[empty] = sym.loc[empty]
    labels = labels.replace("", "NA")
    labels = labels.apply(lambda x: x if len(x) <= max_chars else x[: max_chars - 3] + "...")
    return labels


# ----------------------------
# PLOTS — HEATMAPS
# ----------------------------

def save_heatmap(
    matrix_df: pd.DataFrame,
    out_png: Path,
    title: str,
    sample_cols: list[str],
    row_labels: pd.Series,
    max_rows: int = 300,
    also_pdf: bool = True,
):
    # NÃO MEXER: usa SOMENTE sample_cols
    X = matrix_df[sample_cols].astype(float).values

    if X.shape[0] > max_rows:
        var = np.nanvar(X, axis=1)
        idx = np.argsort(var)[::-1][:max_rows]
        X = X[idx, :]
        row_labels = row_labels.iloc[idx]

    n_rows = X.shape[0]
    n_cols = len(sample_cols)

    fig_w = max(12, 1.3 * n_cols + 11)
    fig_h = max(10, 0.25 * n_rows)

    plt.figure(figsize=(fig_w, fig_h))
    plt.imshow(X, aspect="auto")
    plt.colorbar(label="value")
    plt.xticks(np.arange(n_cols), sample_cols, rotation=45, ha="right")
    plt.yticks(np.arange(n_rows), row_labels.astype(str).values)

    if n_rows <= 40:
        fs = 10
    elif n_rows <= 80:
        fs = 8
    elif n_rows <= 140:
        fs = 6
    else:
        fs = 5

    plt.gca().tick_params(axis="y", labelsize=fs)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    if also_pdf:
        plt.savefig(out_png.with_suffix(".pdf"))
    plt.close()


def save_corr_heatmap(corr: pd.DataFrame, out_png: Path, title: str):
    X = corr.values.astype(float)
    labels = corr.columns.tolist()

    plt.figure(figsize=(9, 8))
    plt.imshow(X, aspect="auto", vmin=-1, vmax=1)
    plt.colorbar(label="correlation")
    plt.xticks(np.arange(len(labels)), labels, rotation=45, ha="right")
    plt.yticks(np.arange(len(labels)), labels)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_png.with_suffix(".pdf"))
    plt.close()


# ----------------------------
# BARPLOT DIVERGENTE POR CATEGORIA (UM ÚNICO)
# ----------------------------

def save_category_diverging_logfc_plot(
    df_cat: pd.DataFrame,
    out_png: Path,
    title: str,
    max_label_chars: int = 140,
):
    """
    Um único gráfico por categoria:
    - valor plotado = plotFC = -logFC (A-up fica negativo, B-up fica positivo)
    - y = descrição (fallback símbolo)
    - ordenado pelo plotFC
    """
    if df_cat.empty:
        return

    d = df_cat.copy()
    d["logFC"] = pd.to_numeric(d["logFC"], errors="coerce")
    d = d.dropna(subset=["logFC"]).copy()
    if d.empty:
        return

    d["plotFC"] = -d["logFC"]  # <<< ESSENCIAL: A-up vira esquerda, B-up vira direita

    labels = _build_desc_labels(d["Description"], d.get("Symbol"), max_chars=max_label_chars)
    d["_label"] = labels

    d = d.sort_values("plotFC", ascending=True)

    n = len(d)
    fig_h = max(8, 0.42 * n)
    fig_w = 24
    plt.figure(figsize=(fig_w, fig_h))

    y = np.arange(n)
    plt.barh(y, d["plotFC"].astype(float).values)

    plt.yticks(y, d["_label"].astype(str).values, fontsize=8)
    plt.axvline(0, linewidth=1)

    plt.xlabel("plotFC = -log2FC(A/B)  (A-up à esquerda; B-up à direita)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_png.with_suffix(".pdf"))
    plt.close()


def save_topN_up_plots_separate(
    filt_sig: pd.DataFrame,
    out_dir: Path,
    base: str,
    top_n: int,
):
    """
    TopN A-up e TopN B-up em gráficos separados.
    """
    df = filt_sig.copy()
    df["logFC"] = pd.to_numeric(df["logFC"], errors="coerce")
    df = df.dropna(subset=["logFC"]).copy()

    # A-up
    dA = df[df["logFC"] > 0].sort_values("logFC", ascending=False).head(top_n).copy()
    if not dA.empty:
        labA = _build_desc_labels(dA["Description"], dA.get("Symbol"), max_chars=140)
        dA["_label"] = labA
        dA = dA.sort_values("logFC", ascending=True)
        nA = len(dA)
        plt.figure(figsize=(24, max(7, 0.55*nA)))
        y = np.arange(nA)
        plt.barh(y, dA["logFC"].astype(float).values)
        plt.yticks(y, dA["_label"].values, fontsize=8)
        plt.xlabel("log2FC (A/B)")
        plt.title(f"Top {top_n} A-up (SIG)")
        plt.tight_layout()
        outA = out_dir / f"{base}_Top{top_n}_Aup.png"
        plt.savefig(outA, dpi=300)
        plt.savefig(outA.with_suffix(".pdf"))
        plt.close()

    # B-up
    dB = df[df["logFC"] < 0].sort_values("logFC", ascending=True).head(top_n).copy()
    if not dB.empty:
        labB = _build_desc_labels(dB["Description"], dB.get("Symbol"), max_chars=140)
        dB["_label"] = labB
        dB["mag"] = (-dB["logFC"]).astype(float)
        dB = dB.sort_values("mag", ascending=True)
        nB = len(dB)
        plt.figure(figsize=(24, max(7, 0.55*nB)))
        y = np.arange(nB)
        plt.barh(y, dB["mag"].values)
        plt.yticks(y, dB["_label"].values, fontsize=8)
        plt.xlabel("magnitude = -log2FC (A/B)  (B-up)")
        plt.title(f"Top {top_n} B-up (SIG)")
        plt.tight_layout()
        outB = out_dir / f"{base}_Top{top_n}_Bup.png"
        plt.savefig(outB, dpi=300)
        plt.savefig(outB.with_suffix(".pdf"))
        plt.close()


# ----------------------------
# Pipeline principal
# ----------------------------

def run_all_analyses(
    edger_path: Path,
    counts_path: Path,
    out_dir: Path,
    base_name: str,
    fdr_cutoff: float = 0.05,
    logfc_abs_cutoff: float = 1.0,
    top_n_bar: int = 20,
    max_rows_heatmap: int = 300,
    drop_uncharacterized: bool = True,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    base = sanitize_basename(base_name)
    if not base:
        raise ValueError("Nome base inválido.")

    # ---- Ler edgeR
    ed = _read_tabular(edger_path)
    geneid_col = _detect_column(ed, ["NCBI.GeneID", "GeneID", "geneid", "GeneId", "Gene ID"])
    sym_col    = _detect_column(ed, ["Symbol", "GeneName", "gene", "Name"])
    desc_col   = _detect_column(ed, ["Description", "descrição", "Descricao", "Desc"])
    logfc_col  = _detect_column(ed, ["logFC", "log2FC", "log2FoldChange"])
    fdr_col    = _detect_column(ed, ["FDR", "adj.P.Val", "padj"])
    logcpm_col = _detect_column(ed, ["logCPM", "AveLogCPM", "logCPM."])

    missing = [("GeneID", geneid_col), ("Symbol", sym_col), ("Description", desc_col),
               ("logFC", logfc_col), ("FDR", fdr_col)]
    missing = [k for k, v in missing if v is None]
    if missing:
        raise ValueError(
            f"edgeR: faltando colunas obrigatórias: {missing}.\n"
            f"Colunas encontradas: {list(ed.columns)}"
        )

    ed = ed.rename(columns={
        geneid_col: "GeneID",
        sym_col: "Symbol",
        desc_col: "Description",
        logfc_col: "logFC",
        fdr_col: "FDR",
        **({logcpm_col: "logCPM"} if logcpm_col else {}),
    })

    ed["GeneID"] = ed["GeneID"].astype(str)
    ed["logFC"] = _as_numeric_series(ed["logFC"])
    ed["FDR"] = _as_numeric_series(ed["FDR"])
    if "logCPM" in ed.columns:
        ed["logCPM"] = _as_numeric_series(ed["logCPM"])

    # ---- Ler counts
    ct = _read_tabular(counts_path)
    ct_geneid = _detect_column(ct, ["GeneID", "NCBI.GeneID", "geneid", "Gene ID"])
    if ct_geneid is None:
        ct_geneid = ct.columns[0]
    ct = ct.rename(columns={ct_geneid: "GeneID"})
    ct["GeneID"] = ct["GeneID"].astype(str)

    # Detectar e normalizar A1..B3/4/5
    ct, sample_cols = detect_and_normalize_sample_cols(ct)

    # manter só GeneID + amostras
    ct = ct[["GeneID"] + sample_cols].copy()
    ct[sample_cols] = ct[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # ---- Merge edgeR + counts
    m = ed.merge(ct, on="GeneID", how="left")
    for c in sample_cols:
        if c not in m.columns:
            m[c] = 0.0
        m[c] = pd.to_numeric(m[c], errors="coerce").fillna(0.0)

    # ---- logCPM (se edgeR não tiver)
    if "logCPM" not in m.columns:
        log2cpm = compute_log2cpm_from_counts(m[["GeneID"] + sample_cols], sample_cols)
        tmp = m.merge(log2cpm, on="GeneID", how="left", suffixes=("", "_log2cpm"))
        m["logCPM"] = tmp[sample_cols].mean(axis=1)

    # ---- Direção
    m["Direction"] = np.where(m["logFC"] > 0, "A-up", np.where(m["logFC"] < 0, "B-up", "0"))

    # ---- Categorias
    m["CategoryList"] = [
        categorize_row(sym, desc)
        for sym, desc in zip(m["Symbol"].fillna("").astype(str), m["Description"].fillna("").astype(str))
    ]
    m["Category"] = m["CategoryList"].apply(lambda x: x[0] if isinstance(x, list) and x else "Other")

    # ---- Filtro estrito SIG
    filt_sig = m[(m["FDR"] < fdr_cutoff) & (m["logFC"].abs() >= logfc_abs_cutoff)].copy()

    # ---- Excluir descrições genéricas
    if drop_uncharacterized:
        filt_sig = filter_informative_only(filt_sig)

    # ---- Saídas tabulares (SIG)
    out_all = out_dir / f"{base}_SIG_FDR{fdr_cutoff:g}_logFC{logfc_abs_cutoff:g}_ALL.tabular"
    out_A   = out_dir / f"{base}_SIG_FDR{fdr_cutoff:g}_logFC{logfc_abs_cutoff:g}_Aup.tabular"
    out_B   = out_dir / f"{base}_SIG_FDR{fdr_cutoff:g}_logFC{logfc_abs_cutoff:g}_Bup.tabular"

    head_cols = ["GeneID", "Symbol", "Description", "logFC", "FDR", "logCPM", "Direction", "Category"]
    rest_cols = [c for c in filt_sig.columns if c not in head_cols]
    filt_sig = filt_sig[head_cols + rest_cols]

    filt_sig.to_csv(out_all, sep="\t", index=False)
    filt_sig[filt_sig["Direction"] == "A-up"].to_csv(out_A, sep="\t", index=False)
    filt_sig[filt_sig["Direction"] == "B-up"].to_csv(out_B, sep="\t", index=False)

    summary = (
        filt_sig.groupby(["Category", "Direction"], as_index=False)
               .size()
               .pivot(index="Category", columns="Direction", values="size")
               .fillna(0)
               .astype(int)
               .reset_index()
    )
    out_sum = out_dir / f"{base}_summary_category_by_direction.tabular"
    summary.to_csv(out_sum, sep="\t", index=False)

    # ---- TOP A-up e TOP B-up (separados)
    save_topN_up_plots_separate(
        filt_sig=filt_sig,
        out_dir=out_dir,
        base=base,
        top_n=top_n_bar,
    )

    # ---- Gráficos por categoria (UM gráfico divergente por categoria)
    cats = sorted(filt_sig["Category"].dropna().unique().tolist())
    for cat in cats:
        sub = filt_sig[filt_sig["Category"] == cat].copy()
        if len(sub) < 5:
            continue
        out_png = out_dir / f"{base}_Category_{cat}_logFC_diverging.png"
        save_category_diverging_logfc_plot(
            df_cat=sub,
            out_png=out_png,
            title=f"{cat} (SIG) — plotFC=-logFC (A-up esquerda; B-up direita) — n={len(sub)}",
            max_label_chars=160,
        )

    # ---- Heatmaps 
    log2cpm_mat = compute_log2cpm_from_counts(m[["GeneID"] + sample_cols], sample_cols)
    hm = filt_sig[["GeneID", "Symbol", "Description", "Category", "Direction"]].merge(
        log2cpm_mat, on="GeneID", how="left"
    ).fillna(0.0)

    hm_z = zscore_rows(hm, sample_cols)

    # Heatmap SIG geral
    rowlab_sig = _build_desc_labels(hm_z["Description"], hm_z["Symbol"], max_chars=140)
    save_heatmap(
        hm_z,
        out_png=out_dir / f"{base}_heatmap_SIG_zscore.png",
        title=f"Heatmap (z-score por gene) — SIG: FDR<{fdr_cutoff:g} & |logFC|>={logfc_abs_cutoff:g}",
        sample_cols=sample_cols,
        row_labels=rowlab_sig,
        max_rows=max_rows_heatmap,
    )

    # Heatmaps por categoria (SIG)
    for cat in sorted(hm_z["Category"].unique()):
        sub = hm_z[hm_z["Category"] == cat].copy()
        if len(sub) < 10:
            continue
        rowlab_cat = _build_desc_labels(sub["Description"], sub["Symbol"], max_chars=140)
        save_heatmap(
            sub,
            out_png=out_dir / f"{base}_heatmap_{cat}_zscore.png",
            title=f"Heatmap (z-score) — {cat} (n={len(sub)})",
            sample_cols=sample_cols,
            row_labels=rowlab_cat,
            max_rows=max_rows_heatmap,
        )

    # Forçar heatmap ROS + Metabolism
    for forced_cat in ["ROS_Antioxidant", "Metabolism"]:
        sub = hm_z[hm_z["Category"] == forced_cat].copy()
        if len(sub) >= 2:
            rowlab_forced = _build_desc_labels(sub["Description"], sub["Symbol"], max_chars=140)
            save_heatmap(
                sub,
                out_png=out_dir / f"{base}_heatmap_{forced_cat}_zscore.png",
                title=f"Heatmap (z-score) — {forced_cat} (n={len(sub)})",
                sample_cols=sample_cols,
                row_labels=rowlab_forced,
                max_rows=max_rows_heatmap,
            )

    # ---- Perfis por categoria (SIG) e correlação
    profiles = []
    for cat in sorted(hm_z["Category"].unique()):
        sub = hm_z[hm_z["Category"] == cat]
        if len(sub) < 5:
            continue
        mean_profile = sub[sample_cols].astype(float).mean(axis=0)
        row = {"Category": cat, **{k: float(v) for k, v in mean_profile.items()}, "n_genes": int(len(sub))}
        profiles.append(row)

    prof_df = pd.DataFrame(profiles)
    out_prof = out_dir / f"{base}_category_profiles_zscore_SIG.tabular"
    if not prof_df.empty:
        prof_df.to_csv(out_prof, sep="\t", index=False)

        corr = prof_df.set_index("Category")[sample_cols].T.corr()
        out_corr = out_dir / f"{base}_category_profile_correlation.tabular"
        corr.to_csv(out_corr, sep="\t")
        save_corr_heatmap(
            corr,
            out_png=out_dir / f"{base}_category_profile_correlation.png",
            title="Correlação entre categorias (perfil médio z-score, genes SIG)",
        )

        plt.figure(figsize=(12, 6))
        for _, r in prof_df.iterrows():
            plt.plot(sample_cols, [r[c] for c in sample_cols], marker="o",
                     label=f"{r['Category']} (n={r['n_genes']})")
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("mean z-score")
        plt.title("Perfis médios por categoria (genes SIG)")
        plt.tight_layout()
        plt.legend(fontsize=8, loc="best")
        plt.savefig(out_dir / f"{base}_category_profiles_lines.png", dpi=300)
        plt.savefig((out_dir / f"{base}_category_profiles_lines.png").with_suffix(".pdf"))
        plt.close()

    return {
        "n_total_edger": int(len(m)),
        "n_sig_after_filters": int(len(filt_sig)),
        "drop_uncharacterized": bool(drop_uncharacterized),
        "out_dir": str(out_dir),
        "base": base,
        "sample_cols": sample_cols,
    }


# ----------------------------
# GUI
# ----------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("edgeR Analyzer — A vs B (GUI)")

        self.edger_var = tk.StringVar()
        self.counts_var = tk.StringVar()
        self.outdir_var = tk.StringVar()
        self.base_var = tk.StringVar(value="A_vs_B")

        self.fdr_var = tk.StringVar(value="0.05")
        self.logfc_var = tk.StringVar(value="1.0")
        self.topn_var = tk.StringVar(value="20")
        self.maxrows_var = tk.StringVar(value="300")

        self.drop_unchar_var = tk.BooleanVar(value=True)

        pad = {"padx": 8, "pady": 6}

        frm = ttk.Frame(self)
        frm.pack(fill="both", expand=True)

        ttk.Label(frm, text="edgeR result (.tabular):").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(frm, textvariable=self.edger_var, width=78).grid(row=0, column=1, **pad)
        ttk.Button(frm, text="Selecionar…", command=self.pick_edger).grid(row=0, column=2, **pad)

        ttk.Label(frm, text="Counts (GeneID + A1..B3/4/5) ou (GeneID + ... + Counts_A1..Counts_B3/4/5):").grid(row=1, column=0, sticky="w", **pad)
        ttk.Entry(frm, textvariable=self.counts_var, width=78).grid(row=1, column=1, **pad)
        ttk.Button(frm, text="Selecionar…", command=self.pick_counts).grid(row=1, column=2, **pad)

        ttk.Label(frm, text="Pasta de saída:").grid(row=2, column=0, sticky="w", **pad)
        ttk.Entry(frm, textvariable=self.outdir_var, width=78).grid(row=2, column=1, **pad)
        ttk.Button(frm, text="Selecionar…", command=self.pick_outdir).grid(row=2, column=2, **pad)

        ttk.Label(frm, text="Nome base do experimento:").grid(row=3, column=0, sticky="w", **pad)
        ttk.Entry(frm, textvariable=self.base_var, width=78).grid(row=3, column=1, **pad)

        pfrm = ttk.LabelFrame(frm, text="Parâmetros (padrão: estrito)")
        pfrm.grid(row=4, column=0, columnspan=3, sticky="we", **pad)

        ttk.Label(pfrm, text="FDR cutoff:").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(pfrm, textvariable=self.fdr_var, width=10).grid(row=0, column=1, sticky="w", **pad)

        ttk.Label(pfrm, text="|logFC| cutoff:").grid(row=0, column=2, sticky="w", **pad)
        ttk.Entry(pfrm, textvariable=self.logfc_var, width=10).grid(row=0, column=3, sticky="w", **pad)

        ttk.Label(pfrm, text="Top N (Top A-up e Top B-up):").grid(row=1, column=0, sticky="w", **pad)
        ttk.Entry(pfrm, textvariable=self.topn_var, width=10).grid(row=1, column=1, sticky="w", **pad)

        ttk.Label(pfrm, text="Max genes no heatmap:").grid(row=1, column=2, sticky="w", **pad)
        ttk.Entry(pfrm, textvariable=self.maxrows_var, width=10).grid(row=1, column=3, sticky="w", **pad)

        ttk.Checkbutton(
            pfrm,
            text="Excluir descrições genéricas (uncharacterized/hypothetical/predicted/unknown...)",
            variable=self.drop_unchar_var
        ).grid(row=2, column=0, columnspan=4, sticky="w", **pad)

        bfrm = ttk.Frame(frm)
        bfrm.grid(row=5, column=0, columnspan=3, sticky="e", **pad)

        ttk.Button(bfrm, text="Rodar análises", command=self.run).grid(row=0, column=0, **pad)
        ttk.Button(bfrm, text="Clear", command=self.clear).grid(row=0, column=1, **pad)

        self.status = tk.StringVar(value="Pronto.")
        ttk.Label(frm, textvariable=self.status).grid(row=6, column=0, columnspan=3, sticky="w", **pad)

    def pick_edger(self):
        p = filedialog.askopenfilename(
            title="Selecione o arquivo edgeR (.tabular)",
            filetypes=[("Tabular/TSV", "*.tabular *.tsv *.txt"), ("All files", "*.*")]
        )
        if p:
            self.edger_var.set(p)
            if not self.base_var.get().strip():
                self.base_var.set(Path(p).stem)

    def pick_counts(self):
        p = filedialog.askopenfilename(
            title="Selecione o arquivo counts (edgeReady)",
            filetypes=[("Tabular/TSV", "*.tabular *.tsv *.txt"), ("All files", "*.*")]
        )
        if p:
            self.counts_var.set(p)

    def pick_outdir(self):
        p = filedialog.askdirectory(title="Selecione a pasta de saída")
        if p:
            self.outdir_var.set(p)

    def clear(self):
        self.edger_var.set("")
        self.counts_var.set("")
        self.outdir_var.set("")
        self.base_var.set("A_vs_B")
        self.fdr_var.set("0.05")
        self.logfc_var.set("1.0")
        self.topn_var.set("20")
        self.maxrows_var.set("300")
        self.drop_unchar_var.set(True)
        self.status.set("Pronto.")

    def run(self):
        try:
            edger = Path(self.edger_var.get().strip()).expanduser()
            counts = Path(self.counts_var.get().strip()).expanduser()
            outdir_str = self.outdir_var.get().strip()
            base = self.base_var.get().strip()

            if not edger.exists():
                raise FileNotFoundError("edgeR: arquivo não encontrado.")
            if not counts.exists():
                raise FileNotFoundError("Counts: arquivo não encontrado.")
            if not outdir_str:
                raise ValueError("Selecione a pasta de saída.")
            outdir = Path(outdir_str).expanduser()

            fdr = float(self.fdr_var.get().strip().replace(",", "."))
            logfc = float(self.logfc_var.get().strip().replace(",", "."))
            topn = int(float(self.topn_var.get().strip().replace(",", ".")))
            maxrows = int(float(self.maxrows_var.get().strip().replace(",", ".")))

            if fdr <= 0 or fdr >= 1:
                raise ValueError("FDR cutoff deve ser entre 0 e 1.")
            if logfc < 0:
                raise ValueError("|logFC| cutoff deve ser >= 0.")
            if topn < 1:
                raise ValueError("Top N deve ser >= 1.")
            if maxrows < 10:
                raise ValueError("Max genes no heatmap deve ser >= 10.")

            self.status.set("Rodando…")
            self.update_idletasks()

            res = run_all_analyses(
                edger_path=edger,
                counts_path=counts,
                out_dir=outdir,
                base_name=base,
                fdr_cutoff=fdr,
                logfc_abs_cutoff=logfc,
                top_n_bar=topn,
                max_rows_heatmap=maxrows,
                drop_uncharacterized=self.drop_unchar_var.get(),
            )

            self.status.set("Concluído.")
            msg_extra = "SIM" if res["drop_uncharacterized"] else "NÃO"
            messagebox.showinfo(
                "OK",
                "Concluído!\n\n"
                f"Genes edgeR: {res['n_total_edger']}\n"
                f"Genes SIG (após filtros): {res['n_sig_after_filters']}\n"
                f"Excluir descrições genéricas: {msg_extra}\n"
                f"Amostras detectadas: {', '.join(res['sample_cols'])}\n\n"
                f"Pasta:\n{res['out_dir']}\n\n"
                "Saídas principais:\n"
                f"- {res['base']}_Top{topn}_Aup.png\n"
                f"- {res['base']}_Top{topn}_Bup.png\n"
                f"- {res['base']}_Category_<Categoria>_logFC_diverging.png\n"
                f"- {res['base']}_heatmap_SIG_zscore.png\n"
                f"- {res['base']}_heatmap_<Categoria>_zscore.png\n"
            )

        except Exception as e:
            self.status.set("Erro.")
            messagebox.showerror("Erro", str(e))


if __name__ == "__main__":
    App().mainloop()