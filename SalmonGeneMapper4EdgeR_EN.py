#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Salmon Gene Mapper (GUI) with 2 tabs:

Tab 1: (wide gene-level)
- Accepts 1 a 3 files do Salmon (tabular/quant.sf exportado em TSV)
- Uses FASTA (headers containing GeneID) to map transcript -> GeneID/GeneName
- (Optional) Uses gene table TSV (NCBI GeneID / Symbol / Description)
- Generates 1 wide file:
    GeneID / symbol / description / Counts_1 TPM_1 / Counts_2 TPM_2 / Counts_3 TPM_3

Tab 2: (2 conditions A vs B x replicates)
- Accepts 3 a 5 files (A1-A5 e B1-B5)
- Automatically generates 4 outputs from:
    (a) Output folder
    (b) Experiment base name
  Outputs:
    <base>_edgeR_counts.tabular        (with GeneID/symbol/description + Counts_A1..B5 as selected)
    <base>_edgeReady_counts.tabular    (only GeneID + A1..B5 as selected)  <- "pure edgeR"
    <base>_TPM.tabular                 (GeneID/symbol/description + TPM_A1..B5 as selected)
    <base>_TPM_actnorm.tabular         (TPM normalized by sum(TPM actinas) per sample)

Requirements:
  pip install pandas

Run:
  python3 salmon_gui_duas_abas.py
"""

import re
from pathlib import Path
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk


# ----------------------------
# Utilities (shared)
# ----------------------------

# Expected header example:
# >NM_001247092.2 PG2 [organism=...] [GeneID=544051]
HEADER_RE = re.compile(
    r"^>(?P<tx>\S+)\s+(?P<gene>\S+).*?\[GeneID=(?P<geneid>\d+)\]"
)


def sanitize_basename(name: str) -> str:
    """Makes the string safe for use as a filename."""
    name = name.strip()
    if not name:
        return ""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def build_map_from_fasta(fasta_path: Path):
    """Returns dicts: transcript_id -> GeneName, transcript_id -> GeneID (strings)."""
    tx_to_gene = {}
    tx_to_geneid = {}

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
        raise ValueError("Could not extract GeneID from the FASTA headers. Check the header format.")
    return tx_to_gene, tx_to_geneid


def load_gene_table(gene_table_path: Path):
    """
    Expected columns:
      - NCBI GeneID
      - Symbol
      - Description

    Returns a DataFrame with: GeneID (string), Symbol, Description
    """
    try:
        gt = pd.read_csv(gene_table_path, sep="\t", dtype=str, encoding="utf-8")
    except UnicodeDecodeError:
        gt = pd.read_csv(gene_table_path, sep="\t", dtype=str, encoding="latin-1")

    required = {"NCBI GeneID", "Symbol", "Description"}
    missing = required - set(gt.columns)
    if missing:
        raise ValueError(f"Gene table: missing columns {sorted(missing)}")

    gt = gt[["NCBI GeneID", "Symbol", "Description"]].copy()
    gt.rename(columns={"NCBI GeneID": "GeneID"}, inplace=True)
    gt["GeneID"] = gt["GeneID"].astype(str)
    return gt


def read_salmon_tabular(quant_path: Path) -> pd.DataFrame:
    """Reads Salmon TSV/tabular file. Expected columns Name/NumReads/TPM."""
    df = pd.read_csv(quant_path, sep="\t")

    if "Name" not in df.columns:
        df.rename(columns={df.columns[0]: "Name"}, inplace=True)

    for col in ("NumReads", "TPM"):
        if col not in df.columns:
            raise ValueError(f"Salmon file {quant_path.name}: column not found '{col}'.")
    return df


def salmon_to_gene_df(
    quant_path: Path,
    tx_to_gene: dict,
    tx_to_geneid: dict,
    sample_label: str,
) -> pd.DataFrame:
    """
    Converts one Salmon file to gene level by aggregating:
      Counts_{label} (NumReads somado)
      TPM_{label} (TPM somado)
    Returns columns: GeneID, GeneName, Counts_label, TPM_label
    """
    df = read_salmon_tabular(quant_path)

    df["GeneName"] = df["Name"].map(tx_to_gene)
    df["GeneID"] = df["Name"].map(tx_to_geneid)

    dfg = df[df["GeneID"].notna()].copy()
    dfg["GeneID"] = dfg["GeneID"].astype(str)

    dfg["NumReads"] = pd.to_numeric(dfg["NumReads"], errors="coerce")
    dfg["TPM"] = pd.to_numeric(dfg["TPM"], errors="coerce")

    gene_df = (
        dfg.groupby(["GeneID", "GeneName"], as_index=False)
           .agg({"NumReads": "sum", "TPM": "sum"})
    )

    gene_df = gene_df.rename(columns={
        "NumReads": f"Counts_{sample_label}",
        "TPM": f"TPM_{sample_label}",
    })
    return gene_df


def merge_gene_dfs(per_sample: list[pd.DataFrame]) -> pd.DataFrame:
    """Progressive outer merge by GeneID. Resolves duplicated GeneName columns."""
    wide = per_sample[0]
    for nxt in per_sample[1:]:
        wide = wide.merge(nxt, on=["GeneID"], how="outer", suffixes=("", "_dup"))
        if "GeneName_dup" in wide.columns:
            wide["GeneName"] = wide["GeneName"].fillna(wide["GeneName_dup"])
            wide.drop(columns=["GeneName_dup"], inplace=True)
    return wide


def add_symbol_description(wide: pd.DataFrame, gene_table_path: Path | None):
    """Adds symbol/description from the gene table when available; otherwise uses GeneName."""
    if gene_table_path is not None:
        gt = load_gene_table(gene_table_path)
        wide["GeneID"] = wide["GeneID"].astype(str)
        wide = wide.merge(gt, on="GeneID", how="left")
        wide["simbolo"] = wide["Symbol"].fillna(wide["GeneName"])
        wide["descrição"] = wide["Description"]
    else:
        wide["simbolo"] = wide["GeneName"]
        wide["descrição"] = pd.NA
    return wide


def ensure_numeric_fill0(df: pd.DataFrame, cols: list[str]):
    for c in cols:
        if c not in df.columns:
            df[c] = 0.0
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df


def actin_mask(symbol_series: pd.Series, desc_series: pd.Series) -> pd.Series:
    """
    Defines actin genes for normalization:
    - symbol contains 'ACT' OU
    - symbol/description contains 'actin' (case-insensitive)
    """
    sym = symbol_series.fillna("").astype(str)
    dsc = desc_series.fillna("").astype(str)
    m = (
        sym.str.contains("ACT", case=False, regex=False) |
        sym.str.contains("actin", case=False, regex=False) |
        dsc.str.contains("actin", case=False, regex=False)
    )
    return m


# ----------------------------
# TAB 1: Wide (up to 3)
# ----------------------------

def make_wide_table_tab1(
    quant_paths: list[Path],
    fasta_path: Path,
    gene_table_path: Path | None,
) -> pd.DataFrame:
    tx_to_gene, tx_to_geneid = build_map_from_fasta(fasta_path)

    per_sample = []
    for i, qp in enumerate(quant_paths, start=1):
        per_sample.append(salmon_to_gene_df(qp, tx_to_gene, tx_to_geneid, str(i)))

    wide = merge_gene_dfs(per_sample)
    wide = add_symbol_description(wide, gene_table_path)

    col_order = ["GeneID", "simbolo", "descrição"]
    for i in range(1, len(quant_paths) + 1):
        col_order += [f"Counts_{i}", f"TPM_{i}"]

    for c in col_order:
        if c not in wide.columns:
            wide[c] = pd.NA

    wide = wide[col_order]
    wide = ensure_numeric_fill0(wide, [c for c in wide.columns if c.startswith("Counts_") or c.startswith("TPM_")])

    sum_counts = wide[[c for c in wide.columns if c.startswith("Counts_")]].sum(axis=1)
    wide = wide.assign(_sumCounts=sum_counts).sort_values("_sumCounts", ascending=False).drop(columns=["_sumCounts"])
    return wide


class Tab1Frame(ttk.Frame):
    """Tab 1: 1 to 3 samples -> wide gene-level (output folder + base name)."""

    def __init__(self, parent):
        super().__init__(parent)

        self.quant1_var = tk.StringVar()
        self.quant2_var = tk.StringVar()
        self.quant3_var = tk.StringVar()
        self.fasta_var = tk.StringVar()
        self.gene_table_var = tk.StringVar()

        self.out_dir_var = tk.StringVar()
        self.base_name_var = tk.StringVar()

        pad = {"padx": 8, "pady": 6}

        ttk.Label(self, text="Salmon 1 (required):").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.quant1_var, width=72).grid(row=0, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_quant1).grid(row=0, column=2, **pad)

        ttk.Label(self, text="Salmon 2 (optional):").grid(row=1, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.quant2_var, width=72).grid(row=1, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_quant2).grid(row=1, column=2, **pad)

        ttk.Label(self, text="Salmon 3 (optional):").grid(row=2, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.quant3_var, width=72).grid(row=2, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_quant3).grid(row=2, column=2, **pad)

        ttk.Label(self, text="Reference FASTA (headers containing GeneID):").grid(row=3, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.fasta_var, width=72).grid(row=3, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_fasta).grid(row=3, column=2, **pad)

        ttk.Label(self, text="Gene table (TSV) (optional):").grid(row=4, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.gene_table_var, width=72).grid(row=4, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_gene_table).grid(row=4, column=2, **pad)

        ttk.Label(self, text="Output folder:").grid(row=5, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.out_dir_var, width=72).grid(row=5, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_out_dir).grid(row=5, column=2, **pad)

        ttk.Label(self, text="Experiment base name:").grid(row=6, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.base_name_var, width=72).grid(row=6, column=1, **pad)

        ttk.Button(self, text="Run", command=self.run).grid(row=7, column=1, sticky="e", **pad)
        ttk.Button(self, text="Clear", command=self.clear_all).grid(row=7, column=2, sticky="w", **pad)

    def _pick_quant(self, title: str):
        return filedialog.askopenfilename(
            title=title,
            filetypes=[("Tabular/TSV", "*.tabular *.tsv *.txt *.sf"), ("All files", "*.*")],
        )

    def pick_quant1(self):
        p = self._pick_quant("Select Salmon 1")
        if p:
            self.quant1_var.set(p)
            if not self.base_name_var.get():
                self.base_name_var.set(Path(p).stem)

    def pick_quant2(self):
        p = self._pick_quant("Select Salmon 2")
        if p:
            self.quant2_var.set(p)

    def pick_quant3(self):
        p = self._pick_quant("Select Salmon 3")
        if p:
            self.quant3_var.set(p)

    def pick_fasta(self):
        p = filedialog.askopenfilename(
            title="Select o FASTA",
            filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("All files", "*.*")],
        )
        if p:
            self.fasta_var.set(p)

    def pick_gene_table(self):
        p = filedialog.askopenfilename(
            title="Select a gene table (TSV)",
            filetypes=[("TSV", "*.tsv *.txt"), ("All files", "*.*")],
        )
        if p:
            self.gene_table_var.set(p)

    def pick_out_dir(self):
        p = filedialog.askdirectory(title="Select a output folder")
        if p:
            self.out_dir_var.set(p)

    def clear_all(self):
        self.quant1_var.set("")
        self.quant2_var.set("")
        self.quant3_var.set("")
        self.fasta_var.set("")
        self.gene_table_var.set("")
        self.out_dir_var.set("")
        self.base_name_var.set("")

    def run(self):
        try:
            q1 = self.quant1_var.get().strip()
            if not q1:
                raise ValueError("Select at least Salmon 1.")

            quant_paths = [Path(q1).expanduser()]
            for q in (self.quant2_var.get().strip(), self.quant3_var.get().strip()):
                if q:
                    quant_paths.append(Path(q).expanduser())

            fasta_str = self.fasta_var.get().strip()
            if not fasta_str:
                raise ValueError("Select o Reference FASTA.")
            fasta = Path(fasta_str).expanduser()

            gene_table_str = self.gene_table_var.get().strip()
            gene_table = Path(gene_table_str).expanduser() if gene_table_str else None

            out_dir_str = self.out_dir_var.get().strip()
            base = sanitize_basename(self.base_name_var.get())
            if not out_dir_str:
                raise ValueError("Select a output folder.")
            if not base:
                raise ValueError("Enter an experiment base name.")

            out_dir = Path(out_dir_str).expanduser()
            out_dir.mkdir(parents=True, exist_ok=True)

            out_path = out_dir / f"{base}_wide_gene_level.tabular"

            for qp in quant_paths:
                if not qp.exists():
                    raise FileNotFoundError(f"Salmon file not found: {qp}")
            if not fasta.exists():
                raise FileNotFoundError(f"FASTA not found: {fasta}")
            if gene_table is not None and not gene_table.exists():
                raise FileNotFoundError(f"Gene table not found: {gene_table}")

            wide = make_wide_table_tab1(quant_paths, fasta, gene_table)

            wide.to_csv(out_path, sep="\t", index=False)

            messagebox.showinfo(
                "OK",
                f"Completed!\n\nSamples: {len(quant_paths)}\nGenes: {len(wide)}\n\nOutput:\n{out_path}"
            )
        except Exception as e:
            messagebox.showerror("Error", str(e))


# ----------------------------
# TAB 2: 2 conditions x 3-5 replicates (with edgeReady)
# ----------------------------

def make_tables_tab2(
    A_paths: list[Path],
    B_paths: list[Path],
    fasta_path: Path,
    gene_table_path: Path | None,
    round_counts_for_edger: bool,
):
    """
    Accepts 3 (required) up to 4 ou 5 por condition.
    Returns:
      counts_df         (with symbol/description + Counts_A1..B5 as selected)
      tpm_df            (with symbol/description + TPM_A1..B5 as selected)
      tpm_actnorm_df    (TPM actin-normalized)
      edge_ready_df     (GeneID + A1..B5 as selected)  <- edgeR "puro"
      act_sums          (dict with actin TPM sums per sample)
    """
    tx_to_gene, tx_to_geneid = build_map_from_fasta(fasta_path)

    labels_A = [f"A{i}" for i in range(1, len(A_paths) + 1)]
    labels_B = [f"B{i}" for i in range(1, len(B_paths) + 1)]

    per_sample = []
    for p, lab in zip(A_paths, labels_A):
        per_sample.append(salmon_to_gene_df(p, tx_to_gene, tx_to_geneid, lab))
    for p, lab in zip(B_paths, labels_B):
        per_sample.append(salmon_to_gene_df(p, tx_to_gene, tx_to_geneid, lab))

    wide = merge_gene_dfs(per_sample)
    wide = add_symbol_description(wide, gene_table_path)

    count_cols = [f"Counts_{x}" for x in labels_A + labels_B]
    tpm_cols   = [f"TPM_{x}" for x in labels_A + labels_B]
    wide = ensure_numeric_fill0(wide, count_cols + tpm_cols)

    # Counts edgeR (with annotation)
    counts_df = wide[["GeneID", "simbolo", "descrição"] + count_cols].copy()

    # EdgeReady: GeneID + A1.. + B1.. (without prefix Counts_)
    edge_ready_df = wide[["GeneID"] + count_cols].copy()
    edge_ready_df = edge_ready_df.rename(columns={f"Counts_{x}": x for x in labels_A + labels_B})

    if round_counts_for_edger:
        for c in count_cols:
            counts_df[c] = counts_df[c].round().astype("Int64")
        for c in labels_A + labels_B:
            edge_ready_df[c] = edge_ready_df[c].round().astype("Int64")

    # TPM (with annotation)
    tpm_df = wide[["GeneID", "simbolo", "descrição"] + tpm_cols].copy()

    # TPM actin-normalized
    m_act = actin_mask(wide["simbolo"], wide["descrição"])
    act_sums = {c: float(wide.loc[m_act, c].sum()) for c in tpm_cols}

    tpm_actnorm_df = tpm_df.copy()
    for c in tpm_cols:
        denom = act_sums[c]
        if denom <= 0:
            tpm_actnorm_df[c] = 0.0
        else:
            tpm_actnorm_df[c] = pd.to_numeric(tpm_actnorm_df[c], errors="coerce").fillna(0.0) / denom

    return counts_df, tpm_df, tpm_actnorm_df, edge_ready_df, act_sums


class Tab2Frame(ttk.Frame):
    """Tab 2: 2 conditions (A, B) x 3-5 replicates -> edgeR/TPM/actin + edgeReady (output folder + base name)."""

    def __init__(self, parent):
        super().__init__(parent)

        self.A1 = tk.StringVar(); self.A2 = tk.StringVar(); self.A3 = tk.StringVar()
        self.A4 = tk.StringVar(); self.A5 = tk.StringVar()
        self.B1 = tk.StringVar(); self.B2 = tk.StringVar(); self.B3 = tk.StringVar()
        self.B4 = tk.StringVar(); self.B5 = tk.StringVar()

        self.fasta_var = tk.StringVar()
        self.gene_table_var = tk.StringVar()

        self.out_dir_var = tk.StringVar()
        self.base_name_var = tk.StringVar()

        self.round_counts_var = tk.BooleanVar(value=True)

        pad = {"padx": 8, "pady": 6}

        # Condition A
        ttk.Label(self, text="Condition A - Salmon 1:").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.A1, width=72).grid(row=0, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.A1, "Select A1")).grid(row=0, column=2, **pad)

        ttk.Label(self, text="Condition A - Salmon 2:").grid(row=1, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.A2, width=72).grid(row=1, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.A2, "Select A2")).grid(row=1, column=2, **pad)

        ttk.Label(self, text="Condition A - Salmon 3:").grid(row=2, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.A3, width=72).grid(row=2, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.A3, "Select A3")).grid(row=2, column=2, **pad)

        ttk.Label(self, text="Condition A - Salmon 4 (optional):").grid(row=3, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.A4, width=72).grid(row=3, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.A4, "Select A4")).grid(row=3, column=2, **pad)

        ttk.Label(self, text="Condition A - Salmon 5 (optional):").grid(row=4, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.A5, width=72).grid(row=4, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.A5, "Select A5")).grid(row=4, column=2, **pad)

        # Condition B
        ttk.Label(self, text="Condition B - Salmon 1:").grid(row=5, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.B1, width=72).grid(row=5, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.B1, "Select B1")).grid(row=5, column=2, **pad)

        ttk.Label(self, text="Condition B - Salmon 2:").grid(row=6, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.B2, width=72).grid(row=6, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.B2, "Select B2")).grid(row=6, column=2, **pad)

        ttk.Label(self, text="Condition B - Salmon 3:").grid(row=7, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.B3, width=72).grid(row=7, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.B3, "Select B3")).grid(row=7, column=2, **pad)

        ttk.Label(self, text="Condition B - Salmon 4 (optional):").grid(row=8, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.B4, width=72).grid(row=8, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.B4, "Select B4")).grid(row=8, column=2, **pad)

        ttk.Label(self, text="Condition B - Salmon 5 (optional):").grid(row=9, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.B5, width=72).grid(row=9, column=1, **pad)
        ttk.Button(self, text="Browse…", command=lambda: self.pick_quant(self.B5, "Select B5")).grid(row=9, column=2, **pad)

        # FASTA / Gene table
        ttk.Label(self, text="Reference FASTA (headers containing GeneID):").grid(row=10, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.fasta_var, width=72).grid(row=10, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_fasta).grid(row=10, column=2, **pad)

        ttk.Label(self, text="Gene table (TSV) (optional):").grid(row=11, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.gene_table_var, width=72).grid(row=11, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_gene_table).grid(row=11, column=2, **pad)

        # Output: folder + base name
        ttk.Label(self, text="Output folder:").grid(row=12, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.out_dir_var, width=72).grid(row=12, column=1, **pad)
        ttk.Button(self, text="Browse…", command=self.pick_out_dir).grid(row=12, column=2, **pad)

        ttk.Label(self, text="Experiment base name:").grid(row=13, column=0, sticky="w", **pad)
        ttk.Entry(self, textvariable=self.base_name_var, width=72).grid(row=13, column=1, **pad)

        ttk.Checkbutton(
            self,
            text="Round counts (recommended for edgeR)",
            variable=self.round_counts_var
        ).grid(row=14, column=0, columnspan=2, sticky="w", **pad)

        ttk.Button(self, text="Run", command=self.run).grid(row=15, column=1, sticky="e", **pad)
        ttk.Button(self, text="Clear", command=self.clear_all).grid(row=15, column=2, sticky="w", **pad)

    def pick_quant(self, var: tk.StringVar, title: str):
        p = filedialog.askopenfilename(
            title=title,
            filetypes=[("Tabular/TSV", "*.tabular *.tsv *.txt *.sf"), ("All files", "*.*")],
        )
        if p:
            var.set(p)
            # Suggest a base name from the first selected file
            if not self.base_name_var.get():
                self.base_name_var.set(Path(p).stem)

    def pick_fasta(self):
        p = filedialog.askopenfilename(
            title="Select o FASTA",
            filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("All files", "*.*")],
        )
        if p:
            self.fasta_var.set(p)

    def pick_gene_table(self):
        p = filedialog.askopenfilename(
            title="Select a gene table (TSV)",
            filetypes=[("TSV", "*.tsv *.txt"), ("All files", "*.*")],
        )
        if p:
            self.gene_table_var.set(p)

    def pick_out_dir(self):
        p = filedialog.askdirectory(title="Select a output folder")
        if p:
            self.out_dir_var.set(p)

    def clear_all(self):
        for v in (self.A1, self.A2, self.A3, self.A4, self.A5,
                  self.B1, self.B2, self.B3, self.B4, self.B5,
                  self.fasta_var, self.gene_table_var,
                  self.out_dir_var, self.base_name_var):
            v.set("")
        self.round_counts_var.set(True)

    def run(self):
        try:
            pathsA_all = [self.A1.get().strip(), self.A2.get().strip(), self.A3.get().strip(),
                          self.A4.get().strip(), self.A5.get().strip()]
            pathsB_all = [self.B1.get().strip(), self.B2.get().strip(), self.B3.get().strip(),
                          self.B4.get().strip(), self.B5.get().strip()]

            pathsA = [p for p in pathsA_all if p]
            pathsB = [p for p in pathsB_all if p]

            if len(pathsA) < 3 or len(pathsB) < 3:
                raise ValueError("You must select at least 3 files for condition A and 3 for condition B.")
            if len(pathsA) > 5 or len(pathsB) > 5:
                raise ValueError("Maximum of 5 files per condition (A and B).")

            A = [Path(p).expanduser() for p in pathsA]
            B = [Path(p).expanduser() for p in pathsB]

            fasta_str = self.fasta_var.get().strip()
            if not fasta_str:
                raise ValueError("Select o Reference FASTA.")
            fasta = Path(fasta_str).expanduser()

            gene_table_str = self.gene_table_var.get().strip()
            gene_table = Path(gene_table_str).expanduser() if gene_table_str else None

            out_dir_str = self.out_dir_var.get().strip()
            base = sanitize_basename(self.base_name_var.get())
            if not out_dir_str:
                raise ValueError("Select a output folder.")
            if not base:
                raise ValueError("Enter an experiment base name (e.g., TETRADxFREE).")

            out_dir = Path(out_dir_str).expanduser()
            out_dir.mkdir(parents=True, exist_ok=True)

            # Automatic names
            out_counts_p     = out_dir / f"{base}_edgeR_counts.tabular"
            out_edge_ready_p = out_dir / f"{base}_edgeReady_counts.tabular"
            out_tpm_p        = out_dir / f"{base}_TPM.tabular"
            out_tpm_act_p    = out_dir / f"{base}_TPM_actnorm.tabular"

            for p in A + B:
                if not p.exists():
                    raise FileNotFoundError(f"Salmon file not found: {p}")
            if not fasta.exists():
                raise FileNotFoundError(f"FASTA not found: {fasta}")
            if gene_table is not None and not gene_table.exists():
                raise FileNotFoundError(f"Gene table not found: {gene_table}")

            counts_df, tpm_df, tpm_act_df, edge_ready_df, act_sums = make_tables_tab2(
                A_paths=A,
                B_paths=B,
                fasta_path=fasta,
                gene_table_path=gene_table,
                round_counts_for_edger=self.round_counts_var.get(),
            )

            counts_df.to_csv(out_counts_p, sep="\t", index=False)
            edge_ready_df.to_csv(out_edge_ready_p, sep="\t", index=False)
            tpm_df.to_csv(out_tpm_p, sep="\t", index=False)
            tpm_act_df.to_csv(out_tpm_act_p, sep="\t", index=False)

            denom_msg = "\n".join([f"{k}: {v:.6g}" for k, v in act_sums.items()])

            messagebox.showinfo(
                "OK",
                "Completed!\n\n"
                f"Genes: {len(counts_df)}\n\n"
                f"Counts edgeR (with annotation):\n{out_counts_p}\n\n"
                f"edgeReady (edgeR puro):\n{out_edge_ready_p}\n\n"
                f"TPM:\n{out_tpm_p}\n\n"
                f"TPM / soma(actinas):\n{out_tpm_act_p}\n\n"
                f"Somas de TPM de actina (per sample):\n{denom_msg}"
            )

        except Exception as e:
            messagebox.showerror("Error", str(e))


# ----------------------------
# APP principal
# ----------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Salmon GUI - Wide + 2xTriplicatas")

        notebook = ttk.Notebook(self)
        notebook.pack(fill="both", expand=True)

        tab1 = Tab1Frame(notebook)
        tab2 = Tab2Frame(notebook)

        notebook.add(tab1, text="Tab 1: Wide (up to 3)")
        notebook.add(tab2, text="Tab 2: 2xTriplicatas (edgeR/TPM/Actina)")


if __name__ == "__main__":
    App().mainloop()