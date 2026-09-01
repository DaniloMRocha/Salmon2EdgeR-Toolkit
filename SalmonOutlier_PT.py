#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
salmon_qc_gui_cli.py

Modo GUI (para alunos):
  python salmon_qc_gui_cli.py

Modo CLI (reprodutível / pipeline):
  python salmon_qc_gui_cli.py run \
    --inputs A.tabular B.tabular C.tabular D.tabular \
    --names A B C D \
    --fasta reference.fa \
    --gene_table genes.tsv \
    --outdir qc_out \
    --top_outlier_genes 50 \
    --min_tpm 1

Reprodutibilidade:
  - Salva run_manifest.json e reproduce_command.txt no outdir

Entradas:
  - Até 5 arquivos Salmon (TSV/tabular) com colunas: Name, TPM
  - FASTA com headers contendo [GeneID=12345] (mesma lógica do seu salmongenemapper)
  - Tabela de genes (opcional): colunas "NCBI GeneID", "Symbol", "Description"

Saídas:
  - correlation_matrix_gene.csv
  - pca_scores_gene.csv
  - pca_PC1_PC2_gene.png
  - corr_heatmap_gene.png
  - outlier_report_gene.txt
  - top_outlier_genes_annotated.csv
  - run_manifest.json
  - reproduce_command.txt
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# GUI
import tkinter as tk
from tkinter import filedialog, messagebox


# -----------------------------
# Parsing/mapping (fasta -> geneid)
# -----------------------------
HEADER_RE = re.compile(r"^>(?P<tx>\S+)\s+(?P<gene>\S+).*?\[GeneID=(?P<geneid>\d+)\]")
REQUIRED_SALMON_COLS = {"Name", "TPM"}


def build_map_from_fasta(fasta_path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Retorna: tx_to_gene, tx_to_geneid"""
    tx_to_gene: Dict[str, str] = {}
    tx_to_geneid: Dict[str, str] = {}

    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m = HEADER_RE.match(line.strip())
            if not m:
                continue
            tx = m.group("tx")
            tx_to_gene[tx] = m.group("gene")
            tx_to_geneid[tx] = m.group("geneid")

    if not tx_to_geneid:
        raise ValueError(
            "Não consegui extrair GeneID dos headers do FASTA.\n"
            "Precisa conter algo como: >TRANSCRIPTO GENE ... [GeneID=12345]"
        )
    return tx_to_gene, tx_to_geneid


def load_gene_table(gene_table_path: Path) -> pd.DataFrame:
    """
    Espera colunas:
      - NCBI GeneID
      - Symbol
      - Description
    """
    try:
        gt = pd.read_csv(gene_table_path, sep="\t", dtype=str, encoding="utf-8")
    except UnicodeDecodeError:
        gt = pd.read_csv(gene_table_path, sep="\t", dtype=str, encoding="latin-1")

    required = {"NCBI GeneID", "Symbol", "Description"}
    missing = required - set(gt.columns)
    if missing:
        raise ValueError(f"Tabela de genes: faltando colunas {sorted(missing)}")

    gt = gt[["NCBI GeneID", "Symbol", "Description"]].copy()
    gt.rename(columns={"NCBI GeneID": "GeneID"}, inplace=True)
    gt["GeneID"] = gt["GeneID"].astype(str)
    return gt


def read_salmon(path: Path) -> pd.DataFrame:
    """Lê TSV/tabular do Salmon. Espera Name e TPM."""
    df = pd.read_csv(path, sep="\t")

    # alguns exports mudam o nome da primeira coluna
    if "Name" not in df.columns:
        df.rename(columns={df.columns[0]: "Name"}, inplace=True)

    missing = REQUIRED_SALMON_COLS - set(df.columns)
    if missing:
        raise ValueError(
            f"{path.name}: faltando colunas {sorted(missing)}. "
            f"Colunas encontradas: {list(df.columns)}"
        )

    df = df[["Name", "TPM"]].copy()
    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce").fillna(0.0)
    return df


def salmon_to_gene_tpm(quant_path: Path, sample_label: str,
                       tx_to_gene: Dict[str, str],
                       tx_to_geneid: Dict[str, str]) -> pd.DataFrame:
    """
    Converte 1 quant para gene-level (TPM somado).
    Retorna: GeneID, GeneName, TPM_<label>
    """
    df = read_salmon(quant_path)

    df["GeneName"] = df["Name"].map(tx_to_gene)
    df["GeneID"] = df["Name"].map(tx_to_geneid)

    dfg = df[df["GeneID"].notna()].copy()
    dfg["GeneID"] = dfg["GeneID"].astype(str)

    gene_df = (
        dfg.groupby(["GeneID", "GeneName"], as_index=False)
           .agg({"TPM": "sum"})
    )

    gene_df = gene_df.rename(columns={"TPM": f"TPM_{sample_label}"})
    return gene_df


def merge_gene_tables(per_sample: List[pd.DataFrame]) -> pd.DataFrame:
    """Merge por GeneID e resolve GeneName."""
    wide = per_sample[0]
    for nxt in per_sample[1:]:
        wide = wide.merge(nxt, on=["GeneID"], how="outer", suffixes=("", "_dup"))
        if "GeneName_dup" in wide.columns:
            wide["GeneName"] = wide["GeneName"].fillna(wide["GeneName_dup"])
            wide.drop(columns=["GeneName_dup"], inplace=True)
    return wide


def add_annotation(wide: pd.DataFrame, gene_table_path: Optional[Path]) -> pd.DataFrame:
    if gene_table_path is None:
        wide["simbolo"] = wide["GeneName"]
        wide["descricao"] = pd.NA
        return wide

    gt = load_gene_table(gene_table_path)
    wide["GeneID"] = wide["GeneID"].astype(str)
    wide = wide.merge(gt, on="GeneID", how="left")
    wide["simbolo"] = wide["Symbol"].fillna(wide["GeneName"])
    wide["descricao"] = wide["Description"]
    wide.drop(columns=[c for c in ["Symbol", "Description"] if c in wide.columns], inplace=True)
    return wide


# -----------------------------
# PCA / Outlier / Genes discrepantes
# -----------------------------
def detect_outlier_pca(scores: pd.DataFrame) -> Tuple[str, pd.Series]:
    """Outlier = maior distância média no espaço PC1-PC2."""
    X = scores[["PC1", "PC2"]].to_numpy()
    n = X.shape[0]
    d = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            d[i, j] = np.linalg.norm(X[i] - X[j])
    mean_dist = (d.sum(axis=1) - np.diag(d)) / (n - 1)
    s = pd.Series(mean_dist, index=scores["sample"].tolist()).sort_values(ascending=False)
    return s.index[0], s


def compute_outlier_genes(expr_log: pd.DataFrame, outlier: str, top_n: int) -> pd.DataFrame:
    """
    expr_log: index=GeneID, colunas=samples com log2(TPM+1)
    z = (x_outlier - mean_others)/sd_others
    """
    x = expr_log[outlier]
    others = expr_log.drop(columns=[outlier])

    mu = others.mean(axis=1)
    sd = others.std(axis=1).replace(0, np.nan)

    z = (x - mu) / sd
    z = z.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    top = z.abs().sort_values(ascending=False).head(top_n)
    out = pd.DataFrame({
        "GeneID": top.index,
        "z_abs": top.values,
        "z": z.loc[top.index].values,
        "outlier_log2TPM1": x.loc[top.index].values,
        "mean_others_log2TPM1": mu.loc[top.index].values,
        "sd_others_log2TPM1": sd.loc[top.index].values,
    })
    return out


def corr_heatmap(corr: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr.values, aspect="auto")
    ax.set_xticks(range(corr.shape[1]))
    ax.set_yticks(range(corr.shape[0]))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticklabels(corr.index)
    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            ax.text(j, i, f"{corr.values[i, j]:.3f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("Correlation (Pearson) on log2(TPM+1) [gene-level]")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ✅ ÚNICA MUDANÇA: outlier em vermelho, sem texto extra
def pca_plot(scores: pd.DataFrame, explained: np.ndarray, out_png: Path, outlier: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))

    # Pontos: outlier vermelho, demais azul
    for _, r in scores.iterrows():
        if r["sample"] == outlier:
            ax.scatter([r["PC1"]], [r["PC2"]], s=110, color="red")
        else:
            ax.scatter([r["PC1"]], [r["PC2"]], s=80)

        # Nome aparece 1x (incluindo o outlier) — sem "OUTLIER:" em cima
        ax.text(r["PC1"], r["PC2"], f" {r['sample']}", va="center")

    ax.set_xlabel(f"PC1 ({explained[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1]*100:.1f}%)")
    ax.set_title("PCA (gene-level log2(TPM+1), standardized)")
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# -----------------------------
# Core runner (shared by GUI/CLI)
# -----------------------------
@dataclass
class RunConfig:
    inputs: List[str]
    names: List[str]
    fasta: str
    gene_table: Optional[str]
    outdir: str
    top_outlier_genes: int
    min_tpm: float


def ensure_outdir(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)


def save_manifest(cfg: RunConfig, outdir: Path, script_path: str) -> None:
    manifest = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "script": str(Path(script_path).resolve()),
        "inputs": [str(Path(p).resolve()) for p in cfg.inputs],
        "names": cfg.names,
        "fasta": str(Path(cfg.fasta).resolve()),
        "gene_table": str(Path(cfg.gene_table).resolve()) if cfg.gene_table else None,
        "outdir": str(outdir.resolve()),
        "parameters": {
            "top_outlier_genes": cfg.top_outlier_genes,
            "min_tpm": cfg.min_tpm,
        },
        "versions": {
            "python": sys.version,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        }
    }
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    cmd = [
        "python", str(Path(script_path).name), "run",
        "--inputs", *manifest["inputs"],
        "--names", *cfg.names,
        "--fasta", manifest["fasta"],
        "--outdir", manifest["outdir"],
        "--top_outlier_genes", str(cfg.top_outlier_genes),
        "--min_tpm", str(cfg.min_tpm),
    ]
    if cfg.gene_table:
        cmd.extend(["--gene_table", manifest["gene_table"]])

    (outdir / "reproduce_command.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")


def run_analysis(cfg: RunConfig, script_path: str) -> Dict[str, str]:
    if not (2 <= len(cfg.inputs) <= 5):
        raise ValueError("Selecione entre 2 e 5 arquivos Salmon.")

    if len(cfg.names) != len(cfg.inputs):
        raise ValueError("Número de nomes deve ser igual ao número de arquivos.")

    outdir = Path(cfg.outdir)
    ensure_outdir(outdir)

    fasta = Path(cfg.fasta)
    if not fasta.exists():
        raise FileNotFoundError(f"FASTA não encontrado: {fasta}")

    gene_table = Path(cfg.gene_table) if cfg.gene_table else None
    if gene_table and not gene_table.exists():
        raise FileNotFoundError(f"Tabela de genes não encontrada: {gene_table}")

    save_manifest(cfg, outdir, script_path)

    tx_to_gene, tx_to_geneid = build_map_from_fasta(fasta)

    per_sample = []
    for qp, lab in zip([Path(p) for p in cfg.inputs], cfg.names):
        per_sample.append(salmon_to_gene_tpm(qp, lab, tx_to_gene, tx_to_geneid))

    wide = merge_gene_tables(per_sample).fillna(0.0)
    wide = add_annotation(wide, gene_table)

    tpm_cols = [f"TPM_{lab}" for lab in cfg.names]
    for c in tpm_cols:
        if c not in wide.columns:
            wide[c] = 0.0
        wide[c] = pd.to_numeric(wide[c], errors="coerce").fillna(0.0)

    expr = wide.set_index("GeneID")[tpm_cols].copy()
    expr.columns = cfg.names

    if cfg.min_tpm and cfg.min_tpm > 0:
        keep = (expr >= cfg.min_tpm).sum(axis=1) >= 2
        expr = expr.loc[keep]

    expr_log = np.log2(expr + 1.0)

    corr = expr_log.corr(method="pearson")
    corr_path = outdir / "correlation_matrix_gene.csv"
    corr.to_csv(corr_path)

    X = expr_log.T
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    pca = PCA(n_components=2, random_state=0).fit(Xs)
    pcs = pca.transform(Xs)

    scores = pd.DataFrame({"sample": X.index, "PC1": pcs[:, 0], "PC2": pcs[:, 1]})
    scores_path = outdir / "pca_scores_gene.csv"
    scores.to_csv(scores_path, index=False)

    outlier, mean_dist = detect_outlier_pca(scores)

    top_genes = compute_outlier_genes(expr_log, outlier, cfg.top_outlier_genes)
    annot = wide[["GeneID", "GeneName", "simbolo", "descricao"]].drop_duplicates("GeneID")
    top_genes = top_genes.merge(annot, on="GeneID", how="left")

    top_path = outdir / "top_outlier_genes_annotated.csv"
    top_genes.to_csv(top_path, index=False)

    pca_png = outdir / "pca_PC1_PC2_gene.png"
    corr_png = outdir / "corr_heatmap_gene.png"
    pca_plot(scores, pca.explained_variance_ratio_, pca_png, outlier=outlier)
    corr_heatmap(corr, corr_png)

    mean_corr = corr.mean(axis=1).sort_values()
    worst_corr = mean_corr.index[0]

    report = []
    report.append("=== PCA outlier (gene-level) ===\n")
    report.append(f"Samples: {', '.join(cfg.names)}\n")
    report.append(f"Genes (após filtro): {expr.shape[0]}\n")
    report.append(f"Outlier (PCA mean distance PC1-PC2): {outlier}\n\n")
    report.append("Mean distance (PC1-PC2):\n")
    for s, v in mean_dist.items():
        report.append(f"  {s}\t{v:.4f}\n")
    report.append("\nWorst mean correlation:\n")
    report.append(f"  {worst_corr}\n")
    report.append("Mean correlations:\n")
    for s, v in mean_corr.items():
        report.append(f"  {s}\t{v:.4f}\n")

    report_path = outdir / "outlier_report_gene.txt"
    report_path.write_text("".join(report), encoding="utf-8")

    return {
        "outdir": str(outdir.resolve()),
        "outlier": outlier,
        "pca_png": str(pca_png.resolve()),
        "corr_png": str(corr_png.resolve()),
        "top_genes_csv": str(top_path.resolve()),
        "report": str(report_path.resolve()),
    }


# -----------------------------
# GUI
# -----------------------------
class App(tk.Tk):
    def __init__(self, script_path: str):
        super().__init__()
        self.title("Salmon QC (PCA + Outlier + Genes) - GUI")
        self.geometry("820x560")
        self.script_path = script_path

        self.salmon_files: List[str] = []
        self.sample_names: List[str] = []
        self.fasta_path: Optional[str] = None
        self.gene_table_path: Optional[str] = None
        self.outdir: Optional[str] = None

        self._build()

    def _build(self):
        pad = {"padx": 10, "pady": 6}

        frm1 = tk.LabelFrame(self, text="1) Arquivos Salmon (2 a 5) [Name + TPM]")
        frm1.pack(fill="x", **pad)

        btn_sel = tk.Button(frm1, text="Selecionar arquivos Salmon...", command=self.pick_salmon_files)
        btn_sel.pack(anchor="w", padx=10, pady=6)

        self.files_box = tk.Text(frm1, height=6, width=100)
        self.files_box.pack(padx=10, pady=6)

        frm2 = tk.LabelFrame(self, text="2) Referência (FASTA) e Tabela de genes (opcional)")
        frm2.pack(fill="x", **pad)

        btn_fa = tk.Button(frm2, text="Selecionar FASTA (obrigatório)", command=self.pick_fasta)
        btn_fa.grid(row=0, column=0, padx=10, pady=6, sticky="w")

        self.fasta_lbl = tk.Label(frm2, text="(nenhum)")
        self.fasta_lbl.grid(row=0, column=1, padx=10, pady=6, sticky="w")

        btn_gt = tk.Button(frm2, text="Selecionar gene_table TSV (opcional)", command=self.pick_gene_table)
        btn_gt.grid(row=1, column=0, padx=10, pady=6, sticky="w")

        self.gt_lbl = tk.Label(frm2, text="(nenhum)")
        self.gt_lbl.grid(row=1, column=1, padx=10, pady=6, sticky="w")

        frm3 = tk.LabelFrame(self, text="3) Saída e parâmetros")
        frm3.pack(fill="x", **pad)

        btn_out = tk.Button(frm3, text="Selecionar pasta de saída", command=self.pick_outdir)
        btn_out.grid(row=0, column=0, padx=10, pady=6, sticky="w")

        self.out_lbl = tk.Label(frm3, text="(nenhuma)")
        self.out_lbl.grid(row=0, column=1, padx=10, pady=6, sticky="w")

        tk.Label(frm3, text="Top genes discrepantes:").grid(row=1, column=0, padx=10, pady=6, sticky="w")
        self.topn_var = tk.StringVar(value="50")
        tk.Entry(frm3, textvariable=self.topn_var, width=10).grid(row=1, column=1, padx=10, pady=6, sticky="w")

        tk.Label(frm3, text="Filtro min_TPM (>= em pelo menos 2 amostras; 0 = desliga):").grid(row=2, column=0, padx=10, pady=6, sticky="w")
        self.mintpm_var = tk.StringVar(value="1")
        tk.Entry(frm3, textvariable=self.mintpm_var, width=10).grid(row=2, column=1, padx=10, pady=6, sticky="w")

        frm4 = tk.Frame(self)
        frm4.pack(fill="x", **pad)

        self.run_btn = tk.Button(frm4, text="RODAR ANÁLISE", command=self.run_clicked, bg="#e6e6e6")
        self.run_btn.pack(side="left", padx=10, pady=10)

        self.status = tk.Label(frm4, text="Status: aguardando seleção de arquivos.")
        self.status.pack(side="left", padx=10)

        frm5 = tk.LabelFrame(self, text="O que sai (reprodutível)")
        frm5.pack(fill="both", expand=True, **pad)

        help_txt = (
            "O programa salva no outdir:\n"
            "- pca_PC1_PC2_gene.png (PCA; outlier em vermelho)\n"
            "- corr_heatmap_gene.png (correlação)\n"
            "- top_outlier_genes_annotated.csv (genes mais discrepantes no outlier)\n"
            "- outlier_report_gene.txt (quem é o outlier e métricas)\n"
            "- run_manifest.json (tudo para reproduzir)\n"
            "- reproduce_command.txt (comando pronto)\n"
        )
        self.help_box = tk.Text(frm5, height=9, width=100)
        self.help_box.insert("1.0", help_txt)
        self.help_box.configure(state="disabled")
        self.help_box.pack(padx=10, pady=10)

    def pick_salmon_files(self):
        files = filedialog.askopenfilenames(
            title="Selecione 2 a 5 arquivos Salmon (TSV/tabular)",
            filetypes=[("TSV/tabular", "*.tsv *.tabular *.txt *.*"), ("All files", "*.*")]
        )
        if not files:
            return
        if len(files) < 2 or len(files) > 5:
            messagebox.showerror("Erro", "Selecione entre 2 e 5 arquivos.")
            return

        self.salmon_files = list(files)
        self.sample_names = [Path(f).stem for f in self.salmon_files]

        self.files_box.configure(state="normal")
        self.files_box.delete("1.0", "end")
        for f, n in zip(self.salmon_files, self.sample_names):
            self.files_box.insert("end", f"{n}\t{f}\n")
        self.files_box.configure(state="disabled")
        self.status.configure(text="Status: Salmon selecionado. Agora selecione FASTA e pasta de saída.")

    def pick_fasta(self):
        f = filedialog.askopenfilename(
            title="Selecione o FASTA referência (headers com [GeneID=...])",
            filetypes=[("FASTA", "*.fa *.fasta *.fna *.*"), ("All files", "*.*")]
        )
        if not f:
            return
        self.fasta_path = f
        self.fasta_lbl.configure(text=f)

    def pick_gene_table(self):
        f = filedialog.askopenfilename(
            title="Selecione gene_table TSV (opcional)",
            filetypes=[("TSV", "*.tsv *.txt *.*"), ("All files", "*.*")]
        )
        if not f:
            return
        self.gene_table_path = f
        self.gt_lbl.configure(text=f)

    def pick_outdir(self):
        d = filedialog.askdirectory(title="Selecione a pasta de saída")
        if not d:
            return
        self.outdir = d
        self.out_lbl.configure(text=d)

    def run_clicked(self):
        try:
            if not self.salmon_files or len(self.salmon_files) < 2:
                messagebox.showerror("Erro", "Selecione 2 a 5 arquivos Salmon.")
                return
            if not self.fasta_path:
                messagebox.showerror("Erro", "Selecione o FASTA (obrigatório).")
                return
            if not self.outdir:
                messagebox.showerror("Erro", "Selecione a pasta de saída.")
                return

            topn = int(self.topn_var.get().strip())
            mintpm = float(self.mintpm_var.get().strip())

            cfg = RunConfig(
                inputs=self.salmon_files,
                names=self.sample_names,
                fasta=self.fasta_path,
                gene_table=self.gene_table_path,
                outdir=self.outdir,
                top_outlier_genes=topn,
                min_tpm=mintpm,
            )

            self.status.configure(text="Status: rodando... (pode levar alguns segundos)")
            self.update_idletasks()

            res = run_analysis(cfg, self.script_path)

            msg = (
                f"Concluído!\n\n"
                f"Outlier sugerido: {res['outlier']}\n"
                f"Saída em: {res['outdir']}\n\n"
                f"Arquivos principais:\n"
                f"- {Path(res['top_genes_csv']).name}\n"
                f"- {Path(res['pca_png']).name}\n"
                f"- {Path(res['corr_png']).name}\n"
                f"- {Path(res['report']).name}\n"
                f"\nReprodutibilidade:\n"
                f"- run_manifest.json\n"
                f"- reproduce_command.txt\n"
            )
            messagebox.showinfo("OK", msg)
            self.status.configure(text=f"Status: concluído. Outlier = {res['outlier']}")

        except Exception as e:
            messagebox.showerror("Erro", str(e))
            self.status.configure(text="Status: erro. Veja a mensagem.")


# -----------------------------
# CLI entry
# -----------------------------
def parse_args():
    ap = argparse.ArgumentParser(add_help=True)
    sub = ap.add_subparsers(dest="cmd")

    runp = sub.add_parser("run", help="Rodar análise via CLI (reprodutível)")
    runp.add_argument("--inputs", nargs="+", required=True, help="2 a 5 arquivos Salmon (Name, TPM)")
    runp.add_argument("--names", nargs="+", required=False, help="Nomes das amostras (se omitido usa basename)")
    runp.add_argument("--fasta", required=True, help="FASTA com headers contendo [GeneID=...]")
    runp.add_argument("--gene_table", required=False, help="TSV com NCBI GeneID / Symbol / Description (opcional)")
    runp.add_argument("--outdir", required=True, help="Pasta de saída")
    runp.add_argument("--top_outlier_genes", type=int, default=50, help="Top N genes discrepantes no outlier")
    runp.add_argument("--min_tpm", type=float, default=1.0, help="Filtro TPM (>= em >=2 amostras); 0 desliga")

    return ap.parse_args()


def main():
    args = parse_args()
    script_path = sys.argv[0]

    if args.cmd is None:
        app = App(script_path=script_path)
        app.mainloop()
        return

    if args.cmd == "run":
        inputs = args.inputs
        if not (2 <= len(inputs) <= 5):
            raise SystemExit("Forneça entre 2 e 5 arquivos em --inputs.")

        if args.names is None:
            names = [Path(p).stem for p in inputs]
        else:
            if len(args.names) != len(inputs):
                raise SystemExit("O número de --names deve ser igual ao número de --inputs.")
            names = args.names

        cfg = RunConfig(
            inputs=inputs,
            names=names,
            fasta=args.fasta,
            gene_table=args.gene_table,
            outdir=args.outdir,
            top_outlier_genes=args.top_outlier_genes,
            min_tpm=args.min_tpm,
        )

        res = run_analysis(cfg, script_path=script_path)
        print(f"[OK] Outdir: {res['outdir']}")
        print(f"[OK] Outlier (PCA mean distance): {res['outlier']}")
        print(f"[OK] Top genes: {res['top_genes_csv']}")
        print(f"[OK] Report: {res['report']}")
        return


if __name__ == "__main__":
    main()