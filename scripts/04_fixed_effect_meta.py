#!/usr/bin/env python3
"""Fixed-effect (inverse-variance) meta-analysis of pPS -> T2D associations
across the MASC and SAZ cohorts, separately for the AFR and MA pPS panels.

Usage:
    python3 scripts/04_fixed_effect_meta.py meta_analysis/input meta_analysis/output

Inputs (in the input dir): {MASC,SAZ}_results_{AFR,MA}_pPS.xlsx with columns
Model, PRS, N, Cases, Controls, Beta, SE, OR, CI_lower, CI_upper, P.

Method
  w_i      = 1 / SE_i^2
  beta_FE  = sum(w_i * beta_i) / sum(w_i)
  SE_FE    = 1 / sqrt(sum(w_i))
  OR_FE    = exp(beta_FE), 95% CI = exp(beta_FE +/- 1.959964 * SE_FE)
  Q        = sum(w_i * (beta_i - beta_FE)^2), df = k - 1, P_Q ~ chi2(df)
  I^2      = max(0, (Q - df) / Q);  H^2 = Q / df
  tau^2_DL = max(0, (Q - df) / (sum(w) - sum(w^2)/sum(w)))   [sensitivity only]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Z975 = stats.norm.ppf(0.975)
COHORTS = ["MASC", "SAZ"]
PANELS = ["AFR", "MA"]
MODEL_ORDER = ["Base_age_sex_PCs", "BMI_adjusted", "Waist_adjusted"]


def load(indir: Path) -> pd.DataFrame:
    frames = []
    for c in COHORTS:
        for p in PANELS:
            d = pd.read_excel(indir / f"{c}_results_{p}_pPS.xlsx")
            # MASC labels the outcome column 'Outcome', SAZ labels it 'Analysis'
            d = d.rename(columns={"Analysis": "Outcome"})
            d["Cohort"], d["Panel"] = c, p
            d["Score"] = d["PRS"].str.replace(rf"_{p}$", "", regex=True)
            frames.append(d)
    return pd.concat(frames, ignore_index=True)


def qc(df: pd.DataFrame) -> list[str]:
    """Check internal consistency of each input row; returns messages."""
    msgs = []
    dup = df.duplicated(["Cohort", "Panel", "Model", "PRS"])
    msgs.append(f"duplicate rows: {int(dup.sum())}")
    or_err = np.abs(np.exp(df.Beta) - df.OR).max()
    lo_err = np.abs(np.exp(df.Beta - Z975 * df.SE) - df.CI_lower).max()
    hi_err = np.abs(np.exp(df.Beta + Z975 * df.SE) - df.CI_upper).max()
    p_err = np.abs(2 * stats.norm.sf(np.abs(df.Beta / df.SE)) - df.P).max()
    msgs.append(f"max |exp(Beta)-OR| = {or_err:.2e}")
    msgs.append(f"max |exp(Beta-1.96SE)-CI_lower| = {lo_err:.2e}")
    msgs.append(f"max |exp(Beta+1.96SE)-CI_upper| = {hi_err:.2e}")
    msgs.append(f"max |Wald P - reported P| = {p_err:.2e}")
    for p in PANELS:
        for m in MODEL_ORDER:
            sets = [set(df.query("Cohort==@c and Panel==@p and Model==@m").Score) for c in COHORTS]
            if sets[0] != sets[1]:
                msgs.append(f"MISMATCH {p}/{m}: only MASC={sets[0]-sets[1]}, only SAZ={sets[1]-sets[0]}")
    for c in COHORTS:
        msgs.append(f"{c} outcome label(s): {sorted(df.query('Cohort==@c').Outcome.unique())}")
    return msgs


def fe_meta(g: pd.DataFrame) -> pd.Series:
    b, se = g.Beta.to_numpy(), g.SE.to_numpy()
    w = 1 / se**2
    k = len(b)
    beta = np.sum(w * b) / np.sum(w)
    se_fe = 1 / np.sqrt(np.sum(w))
    z = beta / se_fe
    q = float(np.sum(w * (b - beta) ** 2))
    df_q = k - 1
    p_q = stats.chi2.sf(q, df_q) if df_q > 0 else np.nan
    i2 = max(0.0, (q - df_q) / q) if q > 0 else 0.0
    h2 = q / df_q if df_q > 0 else np.nan
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - df_q) / c) if df_q > 0 else 0.0
    w_re = 1 / (se**2 + tau2)
    beta_re = np.sum(w_re * b) / np.sum(w_re)
    se_re = 1 / np.sqrt(np.sum(w_re))

    by = {r.Cohort: r for r in g.itertuples()}
    out = {"k": k}
    for coh in COHORTS:
        r = by[coh]
        out.update({
            f"{coh}_N": r.N, f"{coh}_Cases": r.Cases, f"{coh}_Controls": r.Controls,
            f"{coh}_Beta": r.Beta, f"{coh}_SE": r.SE, f"{coh}_OR": r.OR,
            f"{coh}_CI_lower": r.CI_lower, f"{coh}_CI_upper": r.CI_upper, f"{coh}_P": r.P,
            f"{coh}_Weight_pct": 100 * (1 / r.SE**2) / np.sum(w),
        })
    out.update({
        "FE_Beta": beta, "FE_SE": se_fe,
        "FE_OR": np.exp(beta),
        "FE_CI_lower": np.exp(beta - Z975 * se_fe),
        "FE_CI_upper": np.exp(beta + Z975 * se_fe),
        "FE_Z": z, "FE_P": 2 * stats.norm.sf(abs(z)),
        "Q": q, "Q_df": df_q, "Q_P": p_q, "I2_pct": 100 * i2, "H2": h2, "tau2_DL": tau2,
        "Direction": "".join("+" if x > 0 else "-" for x in b),
        "RE_DL_OR": np.exp(beta_re),
        "RE_DL_CI_lower": np.exp(beta_re - Z975 * se_re),
        "RE_DL_CI_upper": np.exp(beta_re + Z975 * se_re),
        "RE_DL_P": 2 * stats.norm.sf(abs(beta_re / se_re)),
    })
    return pd.Series(out)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1)
    return out


def forest(res: pd.DataFrame, panel: str, model: str, path: Path) -> None:
    d = res.query("Panel==@panel and Model==@model").sort_values("FE_OR", ascending=False)
    rows = []  # (label, OR, lo, hi, kind, weight, annot)
    for r in d.itertuples():
        rows.append((r.Score.replace("_", " "), None, None, None, "header", None, ""))
        for coh in COHORTS:
            rows.append((f"   {coh}", getattr(r, f"{coh}_OR"), getattr(r, f"{coh}_CI_lower"),
                         getattr(r, f"{coh}_CI_upper"), coh, getattr(r, f"{coh}_Weight_pct"), ""))
        het = f"I²={r.I2_pct:.0f}%, P_Q={r.Q_P:.2f}"
        rows.append(("   Fixed effect", r.FE_OR, r.FE_CI_lower, r.FE_CI_upper, "fe", None, het))
        rows.append(("", None, None, None, "gap", None, ""))
    rows = rows[:-1]

    n = len(rows)
    fig, ax = plt.subplots(figsize=(11, 0.28 * n + 1.6))
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()
    colors = {"MASC": "#4C72B0", "SAZ": "#DD8452", "fe": "#222222"}
    all_lo = [r[2] for r in rows if r[2] is not None]
    all_hi = [r[3] for r in rows if r[3] is not None]
    xmin, xmax = min(all_lo) * 0.95, max(all_hi) * 1.05
    for y, (lab, o, lo, hi, kind, wt, het) in enumerate(rows):
        if kind == "header":
            ax.text(-0.01, y, lab, transform=ax.get_yaxis_transform(), ha="right",
                    va="center", fontsize=9, fontweight="bold")
            continue
        if kind == "gap":
            continue
        ax.text(-0.01, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8.5)
        txt = f"{o:.2f} ({lo:.2f}–{hi:.2f})"
        if kind == "fe":
            h = 0.32
            ax.fill([lo, o, hi, o], [y, y - h, y, y + h], color=colors["fe"], zorder=3)
            txt += f"   {het}"
            weight = "bold"
        else:
            ax.plot([lo, hi], [y, y], color=colors[kind], lw=1.4, zorder=2)
            ax.plot(o, y, "s", color=colors[kind], ms=3 + 6 * wt / 100, zorder=3)
            txt += f"   w={wt:.0f}%"
            weight = "normal"
        ax.text(1.01, y, txt, transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=8.5, fontweight=weight)
    ax.axvline(1, color="grey", ls="--", lw=0.8, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(xmin, xmax)
    ticks = [t for t in [0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 1.75] if xmin <= t <= xmax]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_yticks([])
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.set_xlabel("Odds ratio for T2D per 1-unit pPS (log scale)")
    ax.text(1.01, -0.9, "OR (95% CI)", transform=ax.get_yaxis_transform(), ha="left",
            va="center", fontsize=9, fontweight="bold")
    ax.set_title(f"{panel} pPS – {model} – inverse-variance fixed-effect meta-analysis (MASC + SAZ)",
                 fontsize=10.5, loc="left", pad=18)
    fig.subplots_adjust(left=0.2, right=0.62, top=1 - 0.9 / (0.28 * n + 1.6), bottom=0.6 / (0.28 * n + 1.6))
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(indir: str, outdir: str) -> None:
    indir, outdir = Path(indir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = load(indir)

    print("== QC ==")
    for m in qc(df):
        print("  ", m)

    res = (df.groupby(["Panel", "Model", "Score"])
             .apply(fe_meta, include_groups=False)
             .reset_index())
    res["k"] = res["k"].astype(int)
    # Multiple testing: BH-FDR within each panel x model family; Bonferroni likewise
    res["FE_P_FDR_BH"] = res.groupby(["Panel", "Model"])["FE_P"].transform(bh)
    res["FE_P_Bonferroni"] = res.groupby(["Panel", "Model"])["FE_P"].transform(
        lambda p: np.minimum(p * len(p), 1))
    res["Model"] = pd.Categorical(res["Model"], MODEL_ORDER, ordered=True)
    res = res.sort_values(["Panel", "Model", "FE_P"]).reset_index(drop=True)
    res.insert(3, "PRS", res.Score + "_" + res.Panel)

    csv = outdir / "fixed_effect_meta_MASC_SAZ_pPS.csv"
    res.to_csv(csv, index=False, float_format="%.6g")
    print(f"\nwrote {csv} ({len(res)} rows)")

    for p in PANELS:
        for m in MODEL_ORDER:
            f = outdir / f"forest_{p}_pPS_{m}.png"
            forest(res, p, m, f)
            print("wrote", f)

    cols = ["Panel", "Model", "PRS", "FE_OR", "FE_CI_lower", "FE_CI_upper", "FE_P",
            "FE_P_FDR_BH", "Q", "Q_P", "I2_pct", "Direction"]
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(res[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))


if __name__ == "__main__":
    main(*sys.argv[1:3])
