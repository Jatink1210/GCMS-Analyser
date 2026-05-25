"""
Knowledge-based re-identification engine.

Takes the parsed peaks and hits DataFrames and returns a single
DataFrame with one row per peak, augmented with the corrected
identification, the reasoning code, the assigned compound class and the
contaminant / unidentified flags.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .config import (
    CLASS_RULES, CONTAMINANT_SUBSTRINGS, IMPLAUSIBLE_SUBSTRINGS,
    SYNONYM_MAP,
)


def is_implausible(name: str) -> bool:
    if not isinstance(name, str):
        return False
    n = name.lower()
    return any(s.lower() in n for s in IMPLAUSIBLE_SUBSTRINGS)


def is_contaminant(name: str) -> bool:
    if not isinstance(name, str):
        return False
    n = name.lower()
    return any(s.lower() in n for s in CONTAMINANT_SUBSTRINGS)


def classify(name: str) -> str:
    if not isinstance(name, str):
        return "Unclassified"
    for cls, pat in CLASS_RULES:
        if pat.search(name):
            return cls
    return "Other"


def strip_tms(name: str) -> str:
    """Remove `, TMS derivative` / `, 2TMS derivative` / `, TBDMS …` suffixes."""
    if not isinstance(name, str):
        return name
    out = re.sub(r",\s*\d*TMS derivative\s*$", "", name, flags=re.I)
    out = re.sub(r",\s*\d*TBDMS derivative\s*$", "", out, flags=re.I)
    out = re.sub(r",\s*\d*TBS derivative\s*$", "", out, flags=re.I)
    return out.strip()


def normalize_name(name: str) -> str:
    """Apply SYNONYM_MAP."""
    if not isinstance(name, str):
        return name
    n = name.strip()
    return SYNONYM_MAP.get(n, n)


def reidentify(peaks: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    """Walk candidate hits per peak and return augmented DataFrame."""
    out_rows: list[dict] = []
    for _, p in peaks.iterrows():
        cand = hits[(hits["sample"] == p["sample"]) &
                    (np.isclose(hits["rt"], p["rt"], atol=0.01))].sort_values(
                        ["hit", "si"], ascending=[True, False])

        sys_name = str(p["system_name"])
        sys_implausible = is_implausible(sys_name)
        sys_contaminant = is_contaminant(sys_name)

        chosen, chosen_reason = None, ""
        seen = set()
        for _, h in cand.iterrows():
            cname = str(h["compname"]).strip()
            if cname in seen:
                continue
            seen.add(cname)
            if is_implausible(cname):
                continue
            if is_contaminant(cname):
                if chosen is None:
                    chosen = h
                    chosen_reason = "contaminant_fallback"
                continue
            chosen = h
            chosen_reason = (
                "kept_system_pick" if h["hit"] == 1 and not sys_implausible
                else f"reassigned_from_hit{int(h['hit'])}"
            )
            break

        if chosen is None and not cand.empty:
            chosen = cand.iloc[0]
            chosen_reason = "no_plausible_alternative"

        # TMS / silyl artefact correction
        chosen_name = str(chosen["compname"]) if chosen is not None else sys_name
        if re.search(r"TMS derivative|TBDMS derivative|trimethylsilyl",
                     chosen_name, re.I):
            stripped = strip_tms(chosen_name)
            if stripped and stripped != chosen_name:
                chosen_name = stripped
                chosen_reason = (
                    chosen_reason + "+TMS_stripped"
                    if chosen_reason and chosen_reason != "kept_system_pick"
                    else "TMS_artefact_stripped"
                )

        chosen_name = normalize_name(chosen_name)

        out_rows.append(dict(
            sample=p["sample"],
            peak=p["peak"],
            rt=p["rt"],
            area=p["area"],
            area_pct=p["area_pct"],
            height=p["height"],
            height_pct=p["height_pct"],
            system_name=sys_name,
            system_si=int(cand.iloc[0]["si"]) if not cand.empty else np.nan,
            final_name=chosen_name if chosen is not None else sys_name,
            final_formula=str(chosen["formula"]) if chosen is not None else "",
            final_mw=int(chosen["molweight"]) if chosen is not None else np.nan,
            final_si=int(chosen["si"]) if chosen is not None else np.nan,
            final_retindex=int(chosen["retindex"]) if chosen is not None else np.nan,
            hit_chosen=int(chosen["hit"]) if chosen is not None else np.nan,
            reidentification=chosen_reason,
            is_system_implausible=sys_implausible,
            is_system_contaminant=sys_contaminant,
        ))

    df = pd.DataFrame(out_rows)
    df["final_class"] = df["final_name"].apply(classify)
    df["is_contaminant_final"] = df["final_name"].apply(is_contaminant)
    df["is_unidentified"] = df["final_name"].apply(is_implausible)
    df.loc[df["is_unidentified"], "final_class"] = "Unidentified"
    return df


def correction_summary(reid: pd.DataFrame) -> pd.DataFrame:
    """Aggregate reidentified peaks into (system_name -> final_name) pairs."""
    changed = reid[
        (reid["final_name"].astype(str) != reid["system_name"].astype(str))
        | reid["is_system_implausible"]
        | reid["is_system_contaminant"]
    ].copy()
    if changed.empty:
        return pd.DataFrame()
    agg = (changed
           .groupby(["system_name", "final_name", "final_class",
                     "is_system_implausible", "is_system_contaminant"],
                    dropna=False)
           .agg(occurrences=("peak", "count"),
                mean_area_pct=("area_pct", "mean"),
                samples=("sample", lambda s: ", ".join(sorted(set(s)))))
           .reset_index())
    agg["mean_area_pct"] = agg["mean_area_pct"].round(3)
    agg = agg.sort_values(["is_system_implausible", "occurrences"],
                          ascending=[False, False])
    return agg
