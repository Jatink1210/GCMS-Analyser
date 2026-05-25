"""Plot figures from the abundance / diversity tables."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

# Distinct, perceptually well-separated colours; cycled if there are
# more groups than colours.
_PALETTE = [
    "#2ca02c", "#98df8a", "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5", "#7f3fbf", "#1f77b4",
    "#aec7e8", "#ff7f0e", "#ffbb78", "#17becf",
    "#bcbd22", "#e377c2", "#8c564b", "#7f7f7f",
]


def palette_for(groups: list[str]) -> dict[str, str]:
    return {g: _PALETTE[i % len(_PALETTE)] for i, g in enumerate(groups)}


# ---------------------------------------------------------------------------
def fig_pca(mat: pd.DataFrame, save: Path) -> None:
    if mat.shape[1] < 2 or mat.shape[0] < 2:
        # Not enough groups or metabolites for ordination — skip silently
        return
    X = np.log10(1.0 + mat.T.values)
    Xz = StandardScaler().fit_transform(X)
    n_comp = min(2, *Xz.shape)
    if n_comp < 1:
        return
    pca = PCA(n_components=n_comp).fit(Xz)
    pc = pca.transform(Xz)
    var = pca.explained_variance_ratio_ * 100
    pal = palette_for(list(mat.columns))

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    for i, g in enumerate(mat.columns):
        x = pc[i, 0]
        y = pc[i, 1] if pc.shape[1] > 1 else 0.0
        ax.scatter(x, y, s=170, c=pal[g], edgecolor="black", linewidth=0.8,
                   label=g, zorder=3)
        ax.annotate(g, (x, y), xytext=(7, 7), textcoords="offset points",
                    fontsize=9)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    if len(var) > 1:
        ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    else:
        ax.set_ylabel("(only one PC available)")
    ax.set_title("PCA — GC-MS metabolome\n"
                 "(Intra + Extra merged, after re-identification)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
              fontsize=8, title="Group")
    fig.tight_layout()
    fig.savefig(save, dpi=300)
    plt.close(fig)


def fig_class_diversity_bar(class_rich: pd.DataFrame, save: Path) -> None:
    if class_rich.empty:
        return
    keep = class_rich[class_rich.max(axis=1) >= 2]
    if keep.empty:
        keep = class_rich
    keep = keep.loc[keep.sum(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(max(11, len(keep.columns) * 1.1), 7))
    bottom = np.zeros(len(keep.columns))
    cmap = plt.get_cmap("tab20", max(len(keep), 1))
    for i, cls in enumerate(keep.index):
        ax.bar(keep.columns, keep.loc[cls], bottom=bottom,
               label=cls, color=cmap(i), edgecolor="white", linewidth=0.4)
        bottom += keep.loc[cls].values
    ax.set_ylabel("Number of distinct metabolites")
    ax.set_title("Metabolite class diversity per group")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
              fontsize=8, title="Compound class")
    ax.tick_params(axis="x", rotation=30)
    for t in ax.get_xticklabels():
        t.set_horizontalalignment("right")
    fig.tight_layout()
    fig.savefig(save, dpi=300)
    plt.close(fig)


def fig_alpha_diversity(div: pd.DataFrame, save: Path) -> None:
    if div.empty:
        return
    pal = palette_for(list(div.index))
    fig, axes = plt.subplots(1, 3, figsize=(max(13, len(div) * 1.1), 4.5))
    cols = ["richness", "shannon", "simpson"]
    titles = ["Richness (S)", "Shannon (H')", "Simpson (1−D)"]
    for ax, c, t in zip(axes, cols, titles):
        bars = ax.bar(div.index, div[c],
                      color=[pal[g] for g in div.index],
                      edgecolor="black", linewidth=0.5)
        ax.set_title(t, fontsize=11)
        ax.tick_params(axis="x", rotation=35)
        for tl in ax.get_xticklabels():
            tl.set_horizontalalignment("right")
        for b, v in zip(bars, div[c]):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                    f"{v:g}", ha="center", va="bottom", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Alpha-diversity of GC-MS metabolome", fontsize=12)
    fig.tight_layout()
    fig.savefig(save, dpi=300)
    plt.close(fig)


def fig_heatmap_top(mat: pd.DataFrame, save: Path, top_n: int = 40) -> None:
    if mat.empty:
        return
    keep = mat.loc[mat.max(axis=1).sort_values(ascending=False).index[:top_n]]
    log = np.log10(keep + 0.01)
    fig, ax = plt.subplots(figsize=(max(8, len(mat.columns) * 0.8),
                                     max(7, len(keep) * 0.22)))
    sns.heatmap(log, ax=ax, cmap="viridis",
                cbar_kws={"label": "log10(area% + 0.01)"},
                linewidths=0.2, linecolor="white",
                yticklabels=[s if len(s) < 55 else s[:52] + "…"
                             for s in keep.index])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(f"Top-{top_n} metabolites (log10 area%)", fontsize=11)
    ax.tick_params(axis="x", rotation=30)
    for t in ax.get_xticklabels():
        t.set_horizontalalignment("right")
    fig.tight_layout()
    fig.savefig(save, dpi=300)
    plt.close(fig)


def fig_heatmap_classes(class_tab: pd.DataFrame, save: Path) -> None:
    if class_tab.empty:
        return
    log = np.log10(class_tab + 0.01)
    fig, ax = plt.subplots(figsize=(max(8, len(class_tab.columns) * 0.9), 7))
    sns.heatmap(log, ax=ax, cmap="rocket_r",
                cbar_kws={"label": "log10(sum area% + 0.01)"},
                linewidths=0.3, linecolor="white",
                annot=class_tab.round(1), fmt="",
                annot_kws={"fontsize": 7})
    ax.set_title("Compound-class abundance per group", fontsize=11)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=30)
    for t in ax.get_xticklabels():
        t.set_horizontalalignment("right")
    fig.tight_layout()
    fig.savefig(save, dpi=300)
    plt.close(fig)


def fig_unique_vs_core(buckets: dict[str, pd.DataFrame], conditions: list[str],
                        save: Path) -> None:
    counts = [len(buckets.get(f"unique_to_{c}", [])) for c in conditions]
    core_n = len(buckets.get("core_metabolites", []))
    pal = palette_for(conditions)
    colors = [pal[c] for c in conditions] + ["#7f7f7f"]
    fig, ax = plt.subplots(figsize=(max(7.5, len(conditions) * 1.4), 4.5))
    bars = ax.bar(conditions + [f"Core (all {len(conditions)})"],
                  counts + [core_n], color=colors,
                  edgecolor="black", linewidth=0.6)
    for b, v in zip(bars, counts + [core_n]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), str(v),
                ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("Number of distinct metabolites")
    ax.set_title("Unique-to-condition vs core metabolites")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(save, dpi=300)
    plt.close(fig)
