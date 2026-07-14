"""Figure 6: ChIP-seq peak classification and positional bias.

Builds the eight-panel figure:

  A-D  ROC curves separating NKX2.1 peaks from those of a negative control factor
       (GATA1, MYOD1, NKX2.5, RXRA), for the CORE, FLANK and ALL models and the PWM.
  E-H  Where each model places its highest-scoring 24 bp window inside a 500 bp
       NKX2.1 peak (ALL, CORE, FLANK, PWM).

Panels are lettered in the order the negative factors appear in ``NEG_TFS``:
A = GATA1, B = MYOD1, C = NKX2.5, D = RXRA.

``POS_TF`` and ``NEG_TFS`` keep the hyphenated gene symbols (``NKX2-1``, ``NKX2-5``)
because those name the ChIP-Atlas input files. The figure labels the proteins
(NKX2.1, NKX2.5), which is what the plotting code renders.

Inputs are the prediction tables written by the scoring step and the FIMO output:

    <root>/predictions/500/500/{all,core,flank}_mean/standardize/sequences_{TF}_VCNNBpnet.csv
    <root>/predictions/500/24/{all,core,flank}_mean/standardize/sequences_NKX2-1_VCNNBpnet.csv
    <root>/FIMO/500/{TF}/fimo.tsv
    <root>/FIMO/500/NKX2-1/best_site.narrowPeak

Usage::

    python scripts/figure6_chipseq_roc.py --root results --out results/plots_final
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_curve

MOTIF_ID = "MA1994.1"
POS_TF = "NKX2-1"
NEG_TFS = ["GATA1", "MYOD1", "NKX2-5", "RXRA"]
MODELS = ["all_mean", "core_mean", "flank_mean"]
MODEL_LABEL = {"all_mean": "ALL", "core_mean": "CORE", "flank_mean": "FLANK"}
COLORS = {"CORE": "tab:blue", "FLANK": "tab:orange", "ALL": "tab:green", "PWM": "blue"}
SEED = 42


def sturge_optimal_bins(data: np.ndarray) -> int:
    data = np.asarray(data)
    if data.size == 0:
        return 1
    width = 1.0 + np.log2(data.size)
    return max(1, math.ceil((data.max() - data.min()) / width))


class Paths:
    def __init__(self, root: Path):
        self.roc_pred = root / "predictions/500/500"
        self.hist_pred = root / "predictions/500/24"
        self.fimo = root / "FIMO/500"


def load_nn_max_scores(p: Paths, model: str, tf: str) -> pd.Series:
    path = p.roc_pred / model / "standardize" / f"sequences_{tf}_VCNNBpnet.csv"
    return pd.read_csv(path, index_col=0)["max_score"]


def load_pwm_scores(p: Paths, tf: str) -> pd.Series:
    df = pd.read_csv(p.fimo / tf / "fimo.tsv", sep="\t", header=0, comment="#")
    return df[df["motif_id"] == MOTIF_ID].groupby("sequence_name")["score"].max()


def load_best_positions_nn(p: Paths, model: str) -> np.ndarray:
    path = p.hist_pred / model / "standardize" / f"sequences_{POS_TF}_VCNNBpnet.csv"
    df = pd.read_csv(path)
    df["all_scores"] = df["all_scores"].apply(lambda x: [float(j) for j in x.strip("[]").split(",")])
    return df["all_scores"].apply(np.argmax).to_numpy()


def load_best_positions_pwm(p: Paths) -> np.ndarray:
    df = pd.read_csv(
        p.fimo / POS_TF / "best_site.narrowPeak", sep="\t", comment="#",
        names=["seq", "start", "stop", "motif", "color", "strand",
               "score", "p-value", "q-value", "peak"],
    )
    return df["start"].to_numpy()


def roc_of(pos: pd.Series, neg: pd.Series) -> tuple[np.ndarray, np.ndarray, float]:
    """Downsample the larger negative set so the classes are balanced, then score."""
    if len(neg) > len(pos):
        neg = neg.sample(n=len(pos), random_state=SEED)
    y_true = np.concatenate([np.ones_like(pos), np.zeros_like(neg)])
    y_score = np.concatenate([pos, neg])
    fpr, tpr, _ = roc_curve(y_true, y_score)
    return fpr, tpr, auc(fpr, tpr)


def panel_label(ax, letter: str) -> None:
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes,
            fontsize=18, fontweight="bold", va="top", ha="left")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True, help="directory holding predictions/ and FIMO/")
    ap.add_argument("--out", type=Path, default=Path("results/plots_final"))
    args = ap.parse_args()

    p = Paths(args.root)
    args.out.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
        "legend.fontsize": 9, "xtick.labelsize": 10, "ytick.labelsize": 10,
    })
    fig, axes = plt.subplots(4, 2, figsize=(12, 16))

    # ---- A-D: ROC ----
    pos_pwm = load_pwm_scores(p, POS_TF)
    pos_nn = {m: load_nn_max_scores(p, m, POS_TF) for m in MODELS}

    for letter, ax, neg_tf in zip("ABCD", [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]], NEG_TFS):
        for model in MODELS:
            label = MODEL_LABEL[model]
            fpr, tpr, a = roc_of(pos_nn[model], load_nn_max_scores(p, model, neg_tf))
            ax.plot(fpr, tpr, color=COLORS[label], label=f"{label} (AUC={a:.3f})", linewidth=1.6)

        fpr, tpr, a = roc_of(pos_pwm, load_pwm_scores(p, neg_tf))
        ax.plot(fpr, tpr, color=COLORS["PWM"], linestyle="--", label=f"PWM (AUC={a:.3f})", linewidth=1.6)

        ax.plot([0, 1], [0, 1], color="black", linestyle=":", linewidth=0.8)
        ax.set_title(f"ROC: NKX2.1 vs {neg_tf.replace('-', '.')}")
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.005)
        ax.grid(linestyle="--", alpha=0.5)
        ax.legend(loc="lower right", framealpha=0.9)
        panel_label(ax, letter)

    # ---- E-H: best-hit position histograms ----
    panels = [
        ("ALL", load_best_positions_nn(p, "all_mean")),
        ("CORE", load_best_positions_nn(p, "core_mean")),
        ("FLANK", load_best_positions_nn(p, "flank_mean")),
        ("PWM", load_best_positions_pwm(p)),
    ]
    for letter, ax, (name, pos) in zip("EFGH", [axes[2, 0], axes[2, 1], axes[3, 0], axes[3, 1]], panels):
        ax.hist(pos, alpha=0.75, density=True, bins=sturge_optimal_bins(pos) + 10,
                color=COLORS[name], edgecolor="white", linewidth=0.3)
        ax.set_title(f"{name}: best-hit position on NKX2.1 peaks")
        ax.set_xlabel("Position in 500 bp window")
        ax.set_ylabel("Density")
        ax.set_xlim(0, 500)
        ax.set_ylim(0, 0.008)
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        panel_label(ax, letter)

    fig.tight_layout(rect=(0, 0, 1, 0.985))
    for ext in ("png", "pdf"):
        path = args.out / f"figure6_chipseq_roc.{ext}"
        fig.savefig(path, dpi=200 if ext == "png" else None, bbox_inches="tight")
        print(f"Wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
