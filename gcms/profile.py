"""
Build the metabolite x group abundance matrix and the diversity tables.

All operations honour the user's intra+extra merging request: within a
(condition, day) group the abundance values from both fractions are
summed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _attach_group(reid: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Return reid with a guaranteed-fresh 'group' column.

    Some upstream callers already merge meta into reid; others don't.
    We always rebuild the group column from the metadata so the output is
    deterministic and there is no risk of duplicate / stale columns.
    """
    out = reid.drop(columns=[c for c in ("group", "condition", "day",
                                          "fraction") if c in reid.columns],
                     errors="ignore")
    out = out.merge(meta[["sample", "condition", "day", "fraction", "group"]],
                    on="sample", how="left")
    return out


def build_abundance_matrix(reid: pd.DataFrame, meta: pd.DataFrame,
                            group_order: list[str]) -> pd.DataFrame:
    """Return DataFrame indexed by metabolite, columns = experimental groups."""
    df = _attach_group(reid, meta)
    df = df[~df["final_class"].isin(["Unidentified"])]
    df = df[~df["is_contaminant_final"]]
    if df.empty:
        return pd.DataFrame(columns=group_order)
    mat = (df.groupby(["final_name", "group"], as_index=False)["area_pct"]
              .sum()
              .pivot(index="final_name", columns="group", values="area_pct")
              .fillna(0.0))
    mat.index.name = "metabolite"
    for g in group_order:
        if g not in mat.columns:
            mat[g] = 0.0
    return mat[group_order]


def diversity_indices(mat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g in mat.columns:
        x = mat[g].values
        present = x[x > 0]
        S = present.size
        if S == 0:
            rows.append(dict(group=g, richness=0, shannon=0, simpson=0,
                             evenness=0, total_signal=0))
            continue
        p = present / present.sum()
        H = float(-(p * np.log(p)).sum())
        D = float(1 - (p ** 2).sum())
        J = H / np.log(S) if S > 1 else 0.0
        rows.append(dict(group=g, richness=S, shannon=round(H, 3),
                         simpson=round(D, 3), evenness=round(J, 3),
                         total_signal=round(present.sum(), 1)))
    return pd.DataFrame(rows).set_index("group").reindex(mat.columns)


def class_abundance(reid: pd.DataFrame, meta: pd.DataFrame,
                     group_order: list[str]) -> pd.DataFrame:
    df = _attach_group(reid, meta)
    df = df[~df["final_class"].isin(["Unidentified"])]
    df = df[~df["is_contaminant_final"]]
    if df.empty:
        return pd.DataFrame(columns=group_order)
    return (df.groupby(["final_class", "group"], as_index=False)["area_pct"]
              .sum()
              .pivot(index="final_class", columns="group", values="area_pct")
              .fillna(0.0)
              .reindex(columns=group_order, fill_value=0.0))


def class_richness(reid: pd.DataFrame, meta: pd.DataFrame,
                    group_order: list[str]) -> pd.DataFrame:
    df = _attach_group(reid, meta)
    df = df[~df["final_class"].isin(["Unidentified"])]
    df = df[~df["is_contaminant_final"]]
    if df.empty:
        return pd.DataFrame(columns=group_order)
    return (df.groupby(["final_class", "group"])["final_name"]
              .nunique()
              .unstack(fill_value=0)
              .reindex(columns=group_order, fill_value=0))


def shared_unique(mat: pd.DataFrame) -> pd.DataFrame:
    presence = (mat > 0).astype(int)
    presence["#groups"] = presence.sum(axis=1)
    presence["groups_present"] = (mat > 0).apply(
        lambda r: ", ".join([c for c, v in r.items() if v]), axis=1)
    presence["mean_area_pct"] = mat.mean(axis=1).round(3)
    presence["max_area_pct"] = mat.max(axis=1).round(3)
    return presence.sort_values(["#groups", "max_area_pct"],
                                ascending=[True, False])


def core_vs_unique(mat: pd.DataFrame, conditions: list[str],
                    group_order: list[str]) -> dict[str, pd.DataFrame]:
    """Buckets keyed by 'core_metabolites', f'unique_to_{cond}' for cond in conditions."""
    cond_groups = {c: [g for g in group_order if g.split("_d")[0] == c]
                   for c in conditions}
    cond_pres = pd.DataFrame({
        c: (mat[gs].sum(axis=1) > 0).astype(int) for c, gs in cond_groups.items()
    })
    cond_pres["#conditions"] = cond_pres.sum(axis=1)
    out: dict[str, pd.DataFrame] = {}
    out["core_metabolites"] = cond_pres[
        cond_pres["#conditions"] == len(conditions)].copy()
    for c in conditions:
        others = [x for x in conditions if x != c]
        mask = (cond_pres[c] == 1) & (cond_pres[others].sum(axis=1) == 0)
        out[f"unique_to_{c}"] = cond_pres.loc[mask].copy()
    out["all_summary"] = cond_pres
    return out


def cond_unique_3way(mat: pd.DataFrame, conditions: list[str],
                      group_order: list[str]) -> dict[str, pd.DataFrame]:
    """For each subset combination, build a 'unique to this condition' table
    (excluding all conditions in `conditions` other than the target).

    Returned dict has key=cond, value=DataFrame indexed by metabolite with
    a single column 'abundance_area_pct'."""
    cond_groups = {c: [g for g in group_order if g.split("_d")[0] == c]
                   for c in conditions}
    cond_abund = pd.DataFrame({
        c: mat[gs].sum(axis=1) for c, gs in cond_groups.items()
    })
    cond_pres = (cond_abund > 0).astype(int)
    out: dict[str, pd.DataFrame] = {}
    for c in conditions:
        others = [x for x in conditions if x != c]
        mask = (cond_pres[c] == 1) & (cond_pres[others].sum(axis=1) == 0)
        sub = cond_abund.loc[mask, [c]].copy()
        sub.columns = ["abundance_area_pct"]
        sub["abundance_area_pct"] = sub["abundance_area_pct"].round(3)
        sub = sub.sort_values("abundance_area_pct", ascending=False)
        out[c] = sub
    return out


def top_abundant_overall(mat: pd.DataFrame, top_n: int = 40) -> pd.DataFrame:
    out = mat.copy()
    out["max_area_pct"] = out.max(axis=1).round(3)
    out["mean_area_pct"] = out.mean(axis=1).round(3)
    out["#groups_present"] = (mat > 0).sum(axis=1)
    out["argmax_group"] = mat.idxmax(axis=1)
    return (out[["max_area_pct", "mean_area_pct",
                 "#groups_present", "argmax_group"]]
            .sort_values("max_area_pct", ascending=False)
            .head(top_n))


def top_n_per_group(mat: pd.DataFrame, n: int = 5) -> pd.DataFrame:
    rows = []
    for g in mat.columns:
        s = mat[g].sort_values(ascending=False).head(n)
        for rank, (m, v) in enumerate(s.items(), start=1):
            rows.append(dict(group=g, rank=rank, metabolite=m,
                             area_pct=round(float(v), 3)))
    return pd.DataFrame(rows)
